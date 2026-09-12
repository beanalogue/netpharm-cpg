"""
Network Pharmacology Data Collector
범용 한약/양약-타겟 데이터 수집 모듈
SQLite 캐싱으로 반복 요청 최소화 (캐시 유효기간: 30일)

동작 확인된 데이터소스 (2025-08):
  - ChEMBL REST API  : 화합물-타겟 활성 데이터 (양약/한약 성분 공통)
  - STRING REST API  : 단백질-단백질 상호작용 (PPI)
  - Enrichr API      : GO/KEGG 농축분석 (enrichment.py 에서 사용)
  - UniProt API      : 유전자-단백질 정보 보완

수동 다운로드 필요:
  - TCMSP            : https://tcmsp-e.com (로그인 후 TSV 다운로드)
  - HERB             : https://herb.ac.cn (공개 API 미제공)

사용 예:
    collector = NetPharmCollector()
    herb_data = collector.collect_herb_compounds_chembl(["ginsenoside Rb1", "baicalein"])
    drug_data = collector.collect_drug_targets_chembl(["aspirin", "metformin"])
    ppi       = collector.collect_ppi(genes)
    collector.close()
"""

import sqlite3
import json
import time
import logging
from pathlib import Path
from datetime import datetime, timedelta

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_DIR = Path(__file__).parent.parent
DB_PATH  = BASE_DIR / "db" / "cache.db"
RAW_DIR  = BASE_DIR / "data" / "raw"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# 0. HTTP 세션 (지수 백오프 자동 재시도)
# ─────────────────────────────────────────────────────────────────────────────
def make_retry_session(
    retries:          int   = 3,
    backoff_factor:   float = 1.0,
    status_forcelist: tuple = (429, 500, 502, 503, 504),
) -> requests.Session:
    """
    HTTPAdapter + urllib3 Retry 조합
    - 429(Rate limit) / 5xx(서버 오류) 자동 재시도
    - 지수 백오프: 1s → 2s → 4s
    - 일반 GET/POST 모두 적용
    """
    session = requests.Session()
    retry   = Retry(
        total              = retries,
        read               = retries,
        connect            = retries,
        backoff_factor     = backoff_factor,
        status_forcelist   = status_forcelist,
        allowed_methods    = {"GET", "POST"},
        raise_on_status    = False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://",  adapter)
    session.headers.update({"User-Agent": "netpharm-research/2.0 (academic)"})
    return session


def get_cache_stats(db_path: Path = DB_PATH) -> dict:
    """캐시 DB 통계 반환 (Streamlit 내보내기용)"""
    if not Path(db_path).exists():
        return {"error": "캐시 DB 없음"}
    try:
        conn = sqlite3.connect(db_path)
        total = conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0]
        by_db = {}
        for row in conn.execute(
            "SELECT substr(key,1,instr(key,':')-1) AS db, COUNT(*) AS n "
            "FROM cache GROUP BY db ORDER BY n DESC"
        ):
            by_db[row[0]] = row[1]
        oldest = conn.execute("SELECT MIN(created_at) FROM cache").fetchone()[0]
        size_kb = Path(db_path).stat().st_size / 1024
        conn.close()
        return {
            "total_entries": total,
            "by_source":     by_db,
            "oldest_entry":  oldest,
            "db_size_kb":    round(size_kb, 1),
        }
    except Exception as e:
        return {"error": str(e)}


# ─────────────────────────────────────────────────────────────────────────────
# 1. SQLite 캐시
# ─────────────────────────────────────────────────────────────────────────────
class CacheDB:
    """API 응답 SQLite 캐시 (30일 유효)"""

    def __init__(self, db_path: Path = DB_PATH):
        self.conn = sqlite3.connect(db_path)
        self._init_tables()

    def _init_tables(self):
        self.conn.executescript("""
            CREATE TABLE IF NOT EXISTS cache (
                key        TEXT PRIMARY KEY,
                value      TEXT,
                created_at TEXT
            );
            CREATE TABLE IF NOT EXISTS compounds (
                compound_id   TEXT,
                compound_name TEXT,
                herb_name     TEXT,
                smiles        TEXT,
                mol_weight    REAL,
                source        TEXT
            );
            CREATE TABLE IF NOT EXISTS targets (
                compound_id  TEXT,
                gene_symbol  TEXT,
                uniprot_id   TEXT,
                target_name  TEXT,
                action_type  TEXT,
                source       TEXT
            );
            CREATE TABLE IF NOT EXISTS drug_targets (
                drug_name    TEXT,
                gene_symbol  TEXT,
                uniprot_id   TEXT,
                action_type  TEXT,
                source       TEXT
            );
            CREATE TABLE IF NOT EXISTS ppi (
                protein_a TEXT,
                protein_b TEXT,
                score     REAL,
                source    TEXT DEFAULT 'STRING'
            );
        """)
        self.conn.commit()

    def get(self, key: str):
        row = self.conn.execute(
            "SELECT value, created_at FROM cache WHERE key=?", (key,)
        ).fetchone()
        if row:
            age = datetime.now() - datetime.fromisoformat(row[1])
            if age < timedelta(days=30):
                return json.loads(row[0])
        return None

    def set(self, key: str, value):
        self.conn.execute(
            "INSERT OR REPLACE INTO cache VALUES (?,?,?)",
            (key, json.dumps(value, ensure_ascii=False), datetime.now().isoformat()),
        )
        self.conn.commit()

    def close(self):
        self.conn.close()


# ─────────────────────────────────────────────────────────────────────────────
# 2. ChEMBL 수집기 (화합물-타겟, 양약+한약 성분 공통)
# ─────────────────────────────────────────────────────────────────────────────
class ChEMBLCollector:
    """
    ChEMBL REST API — 화합물→타겟 활성 데이터
    - 한약 활성 성분명(영문)으로 ChEMBL ID 조회 후 인간 타겟 수집
    - 양약명으로도 동일하게 사용 가능
    """
    BASE = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self, cache: CacheDB, delay: float = 0.5):
        self.cache   = cache
        self.delay   = delay
        self.session = make_retry_session()

    def _get(self, endpoint: str, params: dict = None):
        key = f"CHEMBL:{endpoint}:{json.dumps(params or {}, sort_keys=True)}"
        cached = self.cache.get(key)
        if cached is not None:
            log.debug(f"[ChEMBL] 캐시: {endpoint}")
            return cached
        time.sleep(self.delay)
        try:
            r = self.session.get(f"{self.BASE}/{endpoint}", params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
            self.cache.set(key, data)
            return data
        except Exception as e:
            log.error(f"[ChEMBL] {endpoint} 오류: {e}")
            return {}

    def search_molecule(self, name: str) -> list:
        """
        화합물명으로 ChEMBL molecule 검색 (정확도 우선 순서)
        1. pref_name__iexact  → 정확 이름 매칭
        2. molecule/search?q  → ChEMBL fulltext search (동의어 포함)
        3. pref_name__icontains → 포함 검색 (정확 매칭 우선 정렬)
        반환: [{chembl_id, pref_name, mol_formula, mol_weight, smiles}, ...]
        """
        def _parse(mols):
            result = []
            for m in mols:
                if not m.get("molecule_chembl_id"):
                    continue
                props = m.get("molecule_properties") or {}
                result.append({
                    "chembl_id":   m.get("molecule_chembl_id", ""),
                    "pref_name":   m.get("pref_name") or name,
                    "mol_formula": props.get("full_molformula", ""),
                    "mol_weight":  props.get("full_mwt"),
                    "smiles":      (m.get("molecule_structures") or {}).get("canonical_smiles", ""),
                    "max_phase":   m.get("max_phase"),
                })
            return result

        name_upper = name.upper()

        # 1단계: pref_name 정확 매칭
        data = self._get("molecule", {
            "pref_name__iexact": name,
            "format": "json",
            "limit": 5,
        })
        mols = data.get("molecules", [])
        if mols:
            return _parse(mols)

        # 2단계: synonym 정확 매칭 (molecule_synonyms__molecule_synonym__iexact)
        # fulltext search 대신 synonym DB에서 정확히 이 이름으로 등록된 경우만 반환
        # → "ferulic acid"가 FERULATE의 synonym으로 등록돼 있으면 정확히 찾음
        # → "acetoxyvalerenic acid"가 synonym에 없으면 아무것도 반환하지 않음 (오매칭 방지)
        data2 = self._get("molecule", {
            "molecule_synonyms__molecule_synonym__iexact": name,
            "format": "json",
            "limit":  5,
        })
        mols2 = data2.get("molecules", [])
        if mols2:
            mols2.sort(key=lambda m: 0 if m.get("pref_name") else 1)
            parsed = _parse(mols2)
            if parsed:
                return parsed

        # 3단계: pref_name 포함 검색 — 정확도순으로 정렬
        # 파생물(DIMETHYL ETHER 등) 패널티: 쿼리보다 단어 수 많으면 순위 낮춤
        data3 = self._get("molecule", {
            "pref_name__icontains": name,
            "format": "json",
            "limit":  15,
        })
        mols3 = data3.get("molecules", [])
        if mols3:
            q_words = len(name_upper.split())
            def _rank(m):
                pn = (m.get("pref_name") or "").upper()
                extra = max(0, len(pn.split()) - q_words)
                if pn == name_upper:
                    return (0, extra, len(pn))
                if pn.startswith(name_upper) and extra == 0:
                    return (1, extra, len(pn))
                if pn.startswith(name_upper):
                    return (3, extra, len(pn))  # 파생물 (DIMETHYL ETHER 등)
                return (2, extra, len(pn))
            mols3.sort(key=_rank)
            return _parse(mols3[:5])

        # 3.5단계: ChEMBL 오타 보정 — 단어 길이 ±1 변형으로 재검색
        # (예: liquiritigenin → liquirtigenin, CHEMBL 오타 대응)
        for drop_idx in range(len(name) - 3, 3, -3):
            variant = name[:drop_idx] + name[drop_idx+1:]
            if variant == name:
                continue
            data35 = self._get("molecule", {
                "pref_name__iexact": variant,
                "format": "json",
                "limit":  3,
            })
            mols35 = data35.get("molecules", [])
            if mols35:
                return _parse(mols35)

        return []

    def get_herb_compounds_by_organism(
        self,
        herb_name: str,
        max_results: int = 30,
        mw_cutoff: float = 800.0,
    ) -> list[str]:
        """
        TCMSP 미보유 약재 fallback: ChEMBL에서 학명/속명 기반 천연물 성분 검색.

        검색 전략:
          1. molecule_synonyms__molecule_synonym 에 속명(genus) 포함 + natural_product=1
          2. compound_record 에 source organism 포함 (genus/species 순으로 시도)
        분자량(MW) ≤ mw_cutoff 인 것만 반환.
        """
        parts = herb_name.split()
        genus   = parts[0]
        species = parts[1] if len(parts) > 1 else ""

        found: dict[str, str] = {}  # chembl_id → pref_name

        # ── 전략 1: synonym 검색 (genus) ────────────────────────────────────
        for term in ([genus, species] if species else [genus]):
            data = self._get("molecule", {
                "molecule_synonyms__molecule_synonym__icontains": term,
                "natural_product":  1,
                "format":           "json",
                "limit":            max_results * 3,
            })
            for mol in data.get("molecules", []):
                cid  = mol.get("molecule_chembl_id", "")
                name = mol.get("pref_name", "")
                if not cid or not name:
                    continue
                props = mol.get("molecule_properties") or {}
                try:
                    mw = float(props.get("full_mwt") or 9999)
                except (TypeError, ValueError):
                    mw = 9999
                if mw <= mw_cutoff:
                    found[cid] = name

        # ── 전략 2: compound_record source organism 검색 ────────────────────
        for term in ([genus, species] if species else [genus]):
            data = self._get("compound_record", {
                "compound_source_organism__icontains": term,
                "format": "json",
                "limit":  max_results * 2,
            })
            for rec in data.get("compound_records", []):
                cid = rec.get("molecule_chembl_id", "")
                if cid and cid not in found:
                    # pref_name을 별도 조회
                    mol_data = self._get(f"molecule/{cid}", {"format": "json"})
                    name = mol_data.get("pref_name", "")
                    if not name:
                        continue
                    props = mol_data.get("molecule_properties") or {}
                    try:
                        mw = float(props.get("full_mwt") or 9999)
                    except (TypeError, ValueError):
                        mw = 9999
                    if mw <= mw_cutoff:
                        found[cid] = name

        compounds = list(found.values())[:max_results]
        log.info(f"[ChEMBL_herb] '{herb_name}' → {len(compounds)}개 성분 (fallback)")
        return compounds

    def get_human_targets(
        self,
        chembl_id:    str,
        assay_types:  list = ("B", "F"),
        max_results:  int  = 200,
    ) -> list:
        """
        ChEMBL molecule ID → 인간 타겟 리스트
        assay_types: B=binding, F=functional
        반환: [{gene_symbol, uniprot_id, target_name, action_type, pchembl_value}, ...]
        """
        all_targets = {}
        for atype in assay_types:
            data = self._get("activity", {
                "molecule_chembl_id":  chembl_id,
                "target_organism":     "Homo sapiens",
                "assay_type":          atype,
                "pchembl_value__isnull": False,
                "format":              "json",
                "limit":               max_results,
            })
            for act in data.get("activities", []):
                tid = act.get("target_chembl_id", "")
                if not tid or tid in all_targets:
                    continue
                # 타겟 상세 정보 조회 (유전자 심볼)
                tdata = self._get(f"target/{tid}", {"format": "json"})
                if not tdata:
                    continue
                comps = tdata.get("target_components", [])
                gene_symbol = ""
                uniprot_id  = ""
                for comp in comps:
                    uniprot_id = comp.get("accession", "")
                    # 컴포넌트 시노님에서 gene symbol 추출
                    for syn in comp.get("target_component_synonyms", []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            gene_symbol = syn.get("component_synonym", "")
                            break
                    if gene_symbol:
                        break

                if not gene_symbol:
                    continue

                all_targets[tid] = {
                    "gene_symbol":  gene_symbol.upper(),
                    "uniprot_id":   uniprot_id,
                    "target_name":  tdata.get("pref_name", ""),
                    "target_type":  tdata.get("target_type", ""),
                    "action_type":  act.get("action_type", ""),
                    "pchembl_value": act.get("pchembl_value"),
                }

        return list(all_targets.values())

    def collect_compound_targets(
        self,
        compound_name: str,
        source_label:  str = "ChEMBL",
        use_supplement: bool = True,
        similarity_cutoff: int = 65,
    ) -> dict:
        """
        화합물명 → {molecule, targets} 전체 수집
        한약 활성 성분 / 양약 어느 쪽이든 사용 가능

        타겟 수집 전략 (use_supplement=True 시):
          1순위: ChEMBL 실험 데이터
          2순위: PubChem BioAssay (ChEMBL 0개일 때)
          3순위: ChEMBL 유사도 검색 (PubChem도 0개이고 SMILES 있을 때)
        """
        log.info(f"[ChEMBL] '{compound_name}' 수집 시작")
        molecules = self.search_molecule(compound_name)
        if not molecules:
            log.warning(f"[ChEMBL] '{compound_name}' molecule 검색 결과 없음")
            return {}

        mol = molecules[0]  # 최상위 매칭
        chembl_id = mol["chembl_id"]
        log.info(f"[ChEMBL] 매칭: {mol['pref_name']} ({chembl_id})")

        targets = self.get_human_targets(chembl_id)
        log.info(f"[ChEMBL] '{compound_name}': 타겟 {len(targets)}개")

        # ── 보완 수집: TCMSP 사전데이터 → curated → PubChem → 유사도 → TCMSP ──
        if use_supplement and len(targets) == 0:
            # 1. TCMSP 사전 수집 데이터 (재스크래핑 없이 즉시 사용)
            tcmsp_map = getattr(self, "_tcmsp_targets_map", {})
            prefilled = tcmsp_map.get(compound_name, [])
            if prefilled:
                targets = prefilled
                log.info(f"[TCMSP캐시] '{compound_name}': 사전 타깃 {len(targets)}개 사용")
            else:
                # 2. 기타 보완 수집 (TCMSP는 사전 데이터 없을 때만)
                try:
                    from pubchem_target import supplement_targets
                    herb_name = getattr(self, "_current_herb", "")
                    supp = supplement_targets(
                        compound_name     = compound_name,
                        smiles            = mol.get("smiles", ""),
                        chembl_id         = chembl_id,
                        herb_name         = herb_name,
                        cache             = self.cache,
                        similarity_cutoff = similarity_cutoff,
                        use_tcmsp         = not bool(tcmsp_map),  # 사전 맵 있으면 TCMSP 재시도 안 함
                    )
                    if supp:
                        targets = supp
                        log.info(f"[보완] '{compound_name}': 보완 타깃 {len(targets)}개 "
                                 f"(source: {supp[0].get('source','?')})")
                except Exception as e:
                    log.warning(f"[보완] '{compound_name}' 보완 수집 실패: {e}")

        return {
            "query":    compound_name,
            "molecule": mol,
            "targets":  targets,
            "source":   source_label,
        }


# ─────────────────────────────────────────────────────────────────────────────
# 3. STRING 수집기 (PPI)
# ─────────────────────────────────────────────────────────────────────────────
class STRINGCollector:
    """STRING DB — 단백질-단백질 상호작용"""

    BASE = "https://string-db.org/api/json"

    def __init__(self, cache: CacheDB, delay: float = 1.5, species: int = 9606):
        self.cache   = cache
        self.delay   = delay
        self.species = species
        self.session = make_retry_session()

    def get_ppi(self, gene_symbols: list, min_score: int = 400) -> list:
        """
        유전자 심볼 리스트 → PPI 엣지 리스트
        min_score: 400(medium), 700(high), 900(highest)
        """
        key = f"STRING:ppi:{json.dumps(sorted(gene_symbols))}:{min_score}"
        cached = self.cache.get(key)
        if cached is not None:
            log.info("[STRING] 캐시 사용")
            return cached

        time.sleep(self.delay)
        try:
            r = self.session.post(
                f"{self.BASE}/network",
                data={
                    "identifiers":    "%0d".join(gene_symbols),
                    "species":        self.species,
                    "required_score": min_score,
                    "caller_identity":"netpharm_research",
                },
                timeout=60,
            )
            r.raise_for_status()
            data = r.json()
            self.cache.set(key, data)
            log.info(f"[STRING] PPI {len(data)}개 수집")
            return data
        except Exception as e:
            log.error(f"[STRING] 오류: {e}")
            return []

    def map_identifiers(self, gene_symbols: list) -> list:
        """유전자 심볼 → STRING ID 매핑"""
        key = f"STRING:map:{json.dumps(sorted(gene_symbols))}"
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        time.sleep(self.delay)
        try:
            r = self.session.post(
                f"{self.BASE}/get_string_ids",
                data={
                    "identifiers":    "\r".join(gene_symbols),
                    "species":        self.species,
                    "limit":          1,
                    "caller_identity":"netpharm_research",
                },
                timeout=30,
            )
            r.raise_for_status()
            data = r.json()
            self.cache.set(key, data)
            return data
        except Exception as e:
            log.error(f"[STRING] ID 매핑 오류: {e}")
            return []


# ─────────────────────────────────────────────────────────────────────────────
# 4. UniProt 보조 수집기
# ─────────────────────────────────────────────────────────────────────────────
class UniProtCollector:
    """UniProt REST API — 유전자 정보 보완"""

    BASE = "https://rest.uniprot.org/uniprotkb"

    def __init__(self, cache: CacheDB, delay: float = 0.5):
        self.cache   = cache
        self.delay   = delay
        self.session = make_retry_session()

    def get_gene_info(self, gene_symbol: str) -> dict:
        """유전자 심볼로 UniProt 기본 정보 조회"""
        key = f"UNIPROT:gene:{gene_symbol}"
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        time.sleep(self.delay)
        try:
            r = self.session.get(
                f"{self.BASE}/search",
                params={
                    "query":  f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
                    "format": "json",
                    "size":   1,
                    "fields": "accession,gene_names,protein_name,cc_function",
                },
                timeout=20,
            )
            r.raise_for_status()
            results = r.json().get("results", [])
            if not results:
                return {}
            entry = results[0]
            data = {
                "accession":   entry.get("primaryAccession", ""),
                "gene_symbol": gene_symbol.upper(),
                "protein_name": (
                    (entry.get("proteinDescription") or {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", "")
                ),
            }
            self.cache.set(key, data)
            return data
        except Exception as e:
            log.error(f"[UniProt] {gene_symbol} 오류: {e}")
            return {}


# ─────────────────────────────────────────────────────────────────────────────
# 5. TCMSP 로컬 파일 로더
# ─────────────────────────────────────────────────────────────────────────────
class TCMSPLoader:
    """
    TCMSP 로컬 TSV 파일 로더 (수동 다운로드 후 사용)

    다운로드 방법:
      1. https://tcmsp-e.com 에서 계정 생성
      2. 허브별 Ingredients + Targets TSV 다운로드
      3. data/raw/tcmsp/ 폴더에 저장
      4. load_herb(herb_name, ingredients_file, targets_file) 호출

    반환 형식은 ChEMBLCollector와 동일하므로 NetPharmCollector에서 혼용 가능
    """

    def __init__(self, data_dir: Path = RAW_DIR / "tcmsp"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def load_herb(
        self,
        herb_name:        str,
        ingredients_file: str,
        targets_file:     str,
        ob_threshold:     float = 30.0,
        dl_threshold:     float = 0.18,
    ) -> list:
        """
        TCMSP TSV → ChEMBLCollector 호환 형식으로 변환
        반환: [{"query": compound_name, "molecule": {...}, "targets": [...]}, ...]
        """
        import pandas as pd
        ing_df = pd.read_csv(ingredients_file, sep="\t")
        tgt_df = pd.read_csv(targets_file,     sep="\t")

        # ADME 필터
        ing_df = ing_df[
            (ing_df["OB (%)"] >= ob_threshold) &
            (ing_df["DL"]     >= dl_threshold)
        ].copy()
        log.info(f"[TCMSP] {herb_name}: ADME 필터 후 {len(ing_df)}개 성분")

        results = []
        for _, row in ing_df.iterrows():
            mol_name = str(row.get("Molecule Name", row.get("mol_name", ""))).strip()
            mol_id   = str(row.get("Mol ID",        row.get("mol_id",   ""))).strip()
            # 해당 성분의 타겟 추출
            comp_tgt = tgt_df[tgt_df["Molecule Name"] == mol_name] if "Molecule Name" in tgt_df else pd.DataFrame()
            targets = []
            for _, trow in comp_tgt.iterrows():
                gene = str(trow.get("Target name", "")).strip().upper()
                if gene:
                    targets.append({
                        "gene_symbol": gene,
                        "uniprot_id":  str(trow.get("UniProt ID", "")),
                        "target_name": str(trow.get("Target name", "")),
                        "action_type": "",
                    })
            results.append({
                "query":    mol_name,
                "herb_name":herb_name,
                "molecule": {
                    "chembl_id":  mol_id,
                    "pref_name":  mol_name,
                    "mol_weight": row.get("MW"),
                    "ob":         row.get("OB (%)"),
                    "dl":         row.get("DL"),
                },
                "targets": targets,
                "source":  "TCMSP",
            })
        return results


# ─────────────────────────────────────────────────────────────────────────────
# 6. 통합 수집기
# ─────────────────────────────────────────────────────────────────────────────
class NetPharmCollector:
    """
    전체 파이프라인 통합 수집기

    한약 활성 성분 기반 수집 (ChEMBL):
        results = collector.collect_herb_compounds_chembl(
            herb_name="Panax ginseng",
            compounds=["ginsenoside Rb1", "ginsenoside Rg1", "panaxadiol"]
        )

    양약 기반 수집 (ChEMBL):
        results = collector.collect_drug_targets_chembl(["aspirin", "metformin"])

    TCMSP 로컬 파일 기반:
        results = collector.collect_from_tcmsp(
            herb_name="Panax ginseng",
            ingredients_file="data/raw/tcmsp/ginseng_ingredients.tsv",
            targets_file="data/raw/tcmsp/ginseng_targets.tsv",
        )

    PPI:
        ppi = collector.collect_ppi(gene_symbols)
    """

    def __init__(self):
        self.db      = CacheDB()
        self.chembl  = ChEMBLCollector(self.db)
        self.string  = STRINGCollector(self.db)
        self.uniprot = UniProtCollector(self.db)
        self.tcmsp   = TCMSPLoader()

    # ── 한약 수집 (ChEMBL 경유) ──────────────────────────────────────────────
    def collect_herb_compounds_chembl(
        self,
        herb_name: str,
        compounds: list,
        tcmsp_targets_map: dict = None,
    ) -> dict:
        """
        한약 활성 성분 목록 → ChEMBL 경유 타겟 수집
        compounds: 한약의 주요 활성 성분명(영문) 리스트
        tcmsp_targets_map: {mol_name: [targets]} — TCMSP 사전 수집 타깃
                           (있으면 ChEMBL 실패 시 재스크래핑 없이 바로 사용)
        """
        log.info(f"[수집] {herb_name}: 성분 {len(compounds)}개 ChEMBL 조회")
        results = []
        self.chembl._current_herb = herb_name
        self.chembl._tcmsp_targets_map = tcmsp_targets_map or {}
        for comp in compounds:
            r = self.chembl.collect_compound_targets(comp, source_label="ChEMBL")
            if r:
                r["herb_name"] = herb_name
                results.append(r)
        self.chembl._current_herb = ""
        self.chembl._tcmsp_targets_map = {}

        return {herb_name: {"compounds": results, "source": "ChEMBL"}}

    # ── 양약 수집 (ChEMBL) ───────────────────────────────────────────────────
    def collect_drug_targets_chembl(self, drug_names: list) -> dict:
        """
        양약명 리스트 → ChEMBL 경유 타겟 수집
        반환: {drug_name: {"molecule", "targets"}, ...}
        """
        results = {}
        for drug in drug_names:
            r = self.chembl.collect_compound_targets(drug, source_label="ChEMBL_drug")
            if r:
                results[drug] = r
        return results

    # ── TCMSP 로컬 파일 ──────────────────────────────────────────────────────
    def collect_from_tcmsp(
        self,
        herb_name:        str,
        ingredients_file: str,
        targets_file:     str,
        ob_threshold:     float = 30.0,
        dl_threshold:     float = 0.18,
    ) -> dict:
        """TCMSP 로컬 파일 → 동일한 반환 형식"""
        compounds = self.tcmsp.load_herb(
            herb_name, ingredients_file, targets_file,
            ob_threshold, dl_threshold
        )
        return {herb_name: {"compounds": compounds, "source": "TCMSP"}}

    # ── PPI ──────────────────────────────────────────────────────────────────
    def collect_ppi(self, gene_symbols: list, min_score: int = 400) -> list:
        """유전자 목록 → STRING PPI 엣지 리스트"""
        if not gene_symbols:
            log.warning("유전자 목록이 비어 있음")
            return []
        return self.string.get_ppi(gene_symbols, min_score=min_score)

    # ── 유전자 추출 ──────────────────────────────────────────────────────────
    def extract_gene_symbols(self, herb_results: dict) -> list:
        """
        collect_herb_compounds_chembl 또는 collect_from_tcmsp 결과에서
        유전자 심볼 추출
        """
        genes = set()
        for herb_name, herb_data in herb_results.items():
            for comp_data in herb_data.get("compounds", []):
                for t in comp_data.get("targets", []):
                    gs = t.get("gene_symbol", "")
                    if gs and isinstance(gs, str):
                        genes.add(gs.strip().upper())
        return sorted(genes)

    def extract_drug_genes(self, drug_results: dict) -> list:
        """collect_drug_targets_chembl 결과에서 유전자 심볼 추출"""
        genes = set()
        for drug_name, data in drug_results.items():
            for t in data.get("targets", []):
                gs = t.get("gene_symbol", "")
                if gs:
                    genes.add(gs.strip().upper())
        return sorted(genes)

    def summary(self, herb_results: dict) -> dict:
        """수집 결과 요약"""
        stats = {}
        for herb_name, herb_data in herb_results.items():
            comps = herb_data.get("compounds", [])
            total_targets = sum(len(c.get("targets", [])) for c in comps)
            stats[herb_name] = {
                "compounds": len(comps),
                "targets":   total_targets,
                "source":    herb_data.get("source", ""),
            }
        return stats

    def close(self):
        self.db.close()


# ─────────────────────────────────────────────────────────────────────────────
# 실행 예시 (ChEMBL 기반)
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    collector = NetPharmCollector()

    # 황기(Astragalus) 주요 활성 성분 — 논문에서 자주 사용되는 성분
    herb_results = collector.collect_herb_compounds_chembl(
        herb_name = "Astragalus membranaceus",
        compounds = ["astragaloside IV", "calycosin", "formononetin"],
    )

    stats = collector.summary(herb_results)
    for herb, s in stats.items():
        log.info(f"{herb}: 성분 {s['compounds']}개, 타겟 {s['targets']}개")

    genes = collector.extract_gene_symbols(herb_results)
    log.info(f"추출 유전자: {len(genes)}개 → {genes[:15]}")

    if genes:
        ppi = collector.collect_ppi(genes[:50], min_score=400)
        log.info(f"STRING PPI 엣지: {len(ppi)}개")

    collector.close()
