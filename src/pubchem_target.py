"""
PubChem BioAssay + ChEMBL Similarity 기반 타깃 보완 수집기
ChEMBL 데이터 없는 성분(paeoniflorin, pachymic acid 등)의 타깃 보완용

수집 전략:
  1순위: ChEMBL 실험 데이터 (collect.py 에서 기본 처리)
  2순위: PubChem BioAssay (인간 단백질 대상 실험 데이터)
  3순위: ChEMBL 유사도 검색 (Tanimoto ≥ 65%, 유사 분자 타깃 전이)
"""

import json
import time
import logging
import urllib.parse
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

log = logging.getLogger(__name__)


def _make_session():
    s = requests.Session()
    retry = Retry(total=3, backoff_factor=1.0, status_forcelist=(429, 500, 502, 503, 504))
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://",  HTTPAdapter(max_retries=retry))
    s.headers["User-Agent"] = "netpharm-research/2.0 (academic)"
    return s


SESSION = _make_session()


# ─────────────────────────────────────────────────────────────────────────────
# 1. UniProt ID → Gene symbol (인간 한정)
# ─────────────────────────────────────────────────────────────────────────────
_uniprot_cache: dict[str, str | None] = {}


def _uniprot_to_gene(uniprot_id: str, delay: float = 0.3) -> dict | None:
    """UniProt accession → {gene_symbol, uniprot_id, target_name} (Homo sapiens만)"""
    if uniprot_id in _uniprot_cache:
        return _uniprot_cache[uniprot_id]
    time.sleep(delay)
    try:
        r = SESSION.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=12,
        )
        if r.status_code != 200:
            _uniprot_cache[uniprot_id] = None
            return None
        d = r.json()
        org = d.get("organism", {}).get("scientificName", "")
        if "Homo sapiens" not in org:
            _uniprot_cache[uniprot_id] = None
            return None
        genes = d.get("genes", [])
        gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
        if not gene:
            _uniprot_cache[uniprot_id] = None
            return None
        prot_desc = (
            d.get("proteinDescription", {})
             .get("recommendedName", {})
             .get("fullName", {})
             .get("value", "")
        )
        result = {
            "gene_symbol":  gene.upper(),
            "uniprot_id":   uniprot_id,
            "target_name":  prot_desc,
            "action_type":  "",
            "pchembl_value": None,
        }
        _uniprot_cache[uniprot_id] = result
        return result
    except Exception as e:
        log.debug(f"[UniProt] {uniprot_id} 조회 실패: {e}")
        _uniprot_cache[uniprot_id] = None
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 2. PubChem BioAssay 타깃 수집
# ─────────────────────────────────────────────────────────────────────────────

def get_pubchem_targets(compound_name: str, cache: "CacheDB | None" = None) -> list[dict]:
    """
    화합물명 → PubChem BioAssay 기반 인간 타깃 리스트
    반환: [{gene_symbol, uniprot_id, target_name, action_type, source:'PubChem'}, ...]
    """
    cache_key = f"PUBCHEM:targets:{compound_name.lower()}"
    if cache:
        cached = cache.get(cache_key)
        if cached is not None:
            log.debug(f"[PubChem] 캐시: {compound_name}")
            return cached

    log.info(f"[PubChem] '{compound_name}' BioAssay 조회 시작")

    # 1) 이름 → CID
    try:
        enc_name = urllib.parse.quote(compound_name)
        r = SESSION.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc_name}/cids/JSON",
            timeout=15,
        )
        if r.status_code != 200:
            log.warning(f"[PubChem] CID 조회 실패: {compound_name}")
            return []
        cids = r.json().get("IdentifierList", {}).get("CID", [])
        if not cids:
            return []
        cid = cids[0]
    except Exception as e:
        log.error(f"[PubChem] CID 조회 오류: {e}")
        return []

    # 2) CID → BioAssay 요약 (UniProt 컬럼 수집)
    time.sleep(0.5)
    try:
        r2 = SESSION.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON",
            timeout=25,
        )
        if r2.status_code != 200:
            return []
        rows = r2.json().get("Table", {}).get("Row", [])
    except Exception as e:
        log.error(f"[PubChem] BioAssay 조회 오류: {e}")
        return []

    # UniProt ID 수집 (컬럼 5 = UniProt Accession)
    uniprot_ids = set()
    for row in rows:
        cells = row.get("Cell", [])
        if len(cells) > 5 and cells[5]:
            uid = str(cells[5]).strip()
            if 6 <= len(uid) <= 10 and uid.replace("_", "").isalnum():
                uniprot_ids.add(uid)

    if not uniprot_ids:
        log.info(f"[PubChem] '{compound_name}': UniProt ID 없음")
        return []

    # 3) UniProt → Gene symbol (인간 한정)
    targets = []
    for uid in list(uniprot_ids)[:30]:
        gene_info = _uniprot_to_gene(uid)
        if gene_info:
            targets.append({**gene_info, "source": "PubChem"})

    log.info(f"[PubChem] '{compound_name}': 인간 타깃 {len(targets)}개")

    if cache and targets is not None:
        cache.set(cache_key, targets)

    return targets


# ─────────────────────────────────────────────────────────────────────────────
# 3. ChEMBL 유사도 기반 타깃 전이 (Tanimoto ≥ cutoff)
# ─────────────────────────────────────────────────────────────────────────────

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get_similarity_targets(
    smiles: str,
    compound_name: str = "",
    cutoff: int = 65,
    max_similar: int = 5,
    cache: "CacheDB | None" = None,
) -> list[dict]:
    """
    SMILES → ChEMBL 유사도 검색 → 유사 분자의 인간 타깃 전이
    Tanimoto coefficient ≥ cutoff(%) 인 분자들의 타깃 수집

    반환: [{gene_symbol, uniprot_id, target_name, action_type, similarity, source:'ChEMBL_similarity'}, ...]
    """
    if not smiles:
        return []

    cache_key = f"CHEMBL_SIM:{cutoff}:{smiles[:80]}"
    if cache:
        cached = cache.get(cache_key)
        if cached is not None:
            log.debug(f"[Similarity] 캐시: {compound_name}")
            return cached

    log.info(f"[Similarity] '{compound_name}' 유사도 검색 (cutoff={cutoff}%)")

    enc_smiles = urllib.parse.quote(smiles)
    time.sleep(0.5)

    try:
        r = SESSION.get(
            f"{CHEMBL_BASE}/similarity/{enc_smiles}/{cutoff}",
            params={"format": "json", "limit": max_similar + 1},
            timeout=25,
        )
        if r.status_code != 200:
            return []
        mols = r.json().get("molecules", [])
    except Exception as e:
        log.error(f"[Similarity] 유사도 검색 오류: {e}")
        return []

    # 자기 자신(100% 유사도) 제외하고 타깃 수집
    all_targets: dict[str, dict] = {}
    for m in mols:
        sim = float(m.get("similarity", 0))
        if sim >= 99.9:  # 자기 자신 제외
            continue
        cid = m.get("molecule_chembl_id", "")
        if not cid:
            continue

        time.sleep(0.3)
        try:
            r2 = SESSION.get(
                f"{CHEMBL_BASE}/activity",
                params={
                    "molecule_chembl_id":     cid,
                    "target_organism":        "Homo sapiens",
                    "pchembl_value__isnull":  False,
                    "format":                 "json",
                    "limit":                  100,
                },
                timeout=20,
            )
            if r2.status_code != 200:
                continue
            activities = r2.json().get("activities", [])
        except Exception:
            continue

        for act in activities:
            tid = act.get("target_chembl_id", "")
            if not tid or tid in all_targets:
                continue
            # 타깃 상세
            time.sleep(0.2)
            try:
                rt = SESSION.get(f"{CHEMBL_BASE}/target/{tid}", params={"format": "json"}, timeout=15)
                if rt.status_code != 200:
                    continue
                tdata = rt.json()
                for comp in tdata.get("target_components", []):
                    uid = comp.get("accession", "")
                    gene = ""
                    for syn in comp.get("target_component_synonyms", []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            gene = syn.get("component_synonym", "")
                            break
                    if gene:
                        all_targets[tid] = {
                            "gene_symbol":   gene.upper(),
                            "uniprot_id":    uid,
                            "target_name":   tdata.get("pref_name", ""),
                            "action_type":   act.get("action_type", ""),
                            "pchembl_value": act.get("pchembl_value"),
                            "similarity":    sim,
                            "source":        "ChEMBL_similarity",
                        }
                        break
            except Exception:
                continue

    targets = list(all_targets.values())
    log.info(f"[Similarity] '{compound_name}': 유사도 기반 타깃 {len(targets)}개 (cutoff={cutoff}%)")

    if cache and targets is not None:
        cache.set(cache_key, targets)

    return targets


# ─────────────────────────────────────────────────────────────────────────────
# 4. 통합 보완 수집 함수 (collect.py 에서 호출)
# ─────────────────────────────────────────────────────────────────────────────

def get_curated_targets(chembl_id: str) -> list[dict]:
    """
    문헌 기반 curated 타깃 조회 (tcm_known_targets.json)
    ChEMBL/PubChem에 데이터가 없는 TCM 화합물용
    """
    import json
    curated_path = Path(__file__).parent / "tcm_known_targets.json"
    if not curated_path.exists():
        return []
    try:
        db = json.loads(curated_path.read_text(encoding="utf-8"))
        entry = db.get(chembl_id, {})
        results = []
        for t in entry.get("targets", []):
            results.append({
                "gene_symbol":   t.get("gene_symbol", ""),
                "uniprot_id":    t.get("uniprot_id", ""),
                "target_name":   t.get("target_name", ""),
                "action_type":   t.get("action_type", ""),
                "source":        "curated_literature",
                "source_pmid":   t.get("source_pmid", ""),
            })
        if results:
            log.info(f"[Curated] {chembl_id} → {len(results)}개 문헌 타깃 보완")
        return results
    except Exception as e:
        log.warning(f"[Curated] 로드 실패: {e}")
        return []


def supplement_targets(
    compound_name: str,
    smiles: str = "",
    chembl_id: str = "",
    herb_name: str = "",
    cache=None,
    similarity_cutoff: int = 65,
    use_pubchem: bool = True,
    use_similarity: bool = True,
    use_curated: bool = True,
    use_tcmsp: bool = True,
) -> list[dict]:
    """
    ChEMBL 타깃이 0개인 성분을 위한 보완 수집
    1단계: 문헌 curated DB (tcm_known_targets.json)
    2단계: PubChem BioAssay
    3단계: ChEMBL 유사도 기반 전이 (SMILES 있을 때만)
    4단계: TCMSP-E 스크래핑 (약재명 알 때)

    반환: 타깃 리스트 (source 필드로 출처 구분)
    """
    targets = []

    # 1단계: 문헌 curated DB
    if use_curated and chembl_id:
        curated = get_curated_targets(chembl_id)
        targets.extend(curated)

    # 2단계: PubChem BioAssay
    if not targets and use_pubchem:
        pc_targets = get_pubchem_targets(compound_name, cache=cache)
        targets.extend(pc_targets)

    # 3단계: PubChem에도 없고 SMILES가 있으면 유사도 검색
    if not targets and use_similarity and smiles:
        sim_targets = get_similarity_targets(
            smiles, compound_name=compound_name,
            cutoff=similarity_cutoff, cache=cache,
        )
        targets.extend(sim_targets)

    # 4단계: TCMSP-E 스크래핑 (약재명이 있을 때)
    if not targets and use_tcmsp and herb_name:
        try:
            from tcmsp_scraper import get_tcmsp_targets_for_compound
            tcmsp_targets = get_tcmsp_targets_for_compound(compound_name, herb_name)
            if tcmsp_targets:
                log.info(f"[TCMSP] {compound_name}: {len(tcmsp_targets)}개 타깃 보완")
            targets.extend(tcmsp_targets)
        except Exception as e:
            log.warning(f"[TCMSP] 스크래핑 실패: {e}")

    return targets
