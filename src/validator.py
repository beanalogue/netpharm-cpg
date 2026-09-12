"""
분석 전 DB 검증 모듈
- ChEMBL: 성분명 → 매칭 분자 존재 여부 + 정확 이름
- OpenTargets: 질환명 → 매칭 질환 존재 여부 + 타깃 수
캐시 활용, 빠른 응답 우선
"""

import time
import logging
from pathlib import Path

import requests

log = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
OT_GRAPHQL   = "https://api.platform.opentargets.org/api/v4/graphql"


def _session():
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    s = requests.Session()
    r = Retry(total=2, backoff_factor=0.5, status_forcelist=(429, 500, 502, 503))
    s.mount("https://", HTTPAdapter(max_retries=r))
    s.headers["User-Agent"] = "netpharm-validator/1.0 (academic)"
    return s

_SES = _session()


# ─────────────────────────────────────────────────────────────────────────────
# 1. 성분 / 약물 → ChEMBL 검증
# ─────────────────────────────────────────────────────────────────────────────

def validate_compound(name: str, cache=None) -> dict:
    """
    화합물명 ChEMBL 검증.
    반환:
      {
        "query":       원본 입력명,
        "found":       True/False,
        "chembl_id":   "CHEMBL12345" or "",
        "matched_name": ChEMBL pref_name,
        "smiles":      SMILES string,
        "status":      "exact" | "synonym" | "partial" | "not_found",
        "note":        사용자에게 보여줄 메시지,
      }
    """
    cache_key = f"VALIDATE:compound:{name.lower().strip()}"
    if cache:
        hit = cache.get(cache_key)
        if hit is not None:
            return hit

    name_upper = name.upper().strip()
    result = {
        "query": name, "found": False,
        "chembl_id": "", "matched_name": "", "smiles": "",
        "status": "not_found", "note": "",
    }

    def _parse_mol(m):
        props = m.get("molecule_properties") or {}
        structs = m.get("molecule_structures") or {}
        return {
            "chembl_id":    m.get("molecule_chembl_id", ""),
            "matched_name": m.get("pref_name") or "",
            "smiles":       structs.get("canonical_smiles", ""),
            "mw":           props.get("full_mwt"),
        }

    try:
        # 1) 정확 매칭
        r = _SES.get(f"{CHEMBL_BASE}/molecule",
                     params={"pref_name__iexact": name, "format": "json", "limit": 1},
                     timeout=15)
        mols = r.json().get("molecules", []) if r.ok else []
        if mols:
            m = _parse_mol(mols[0])
            result.update(found=True, status="exact",
                          note=f"정확 매칭: {m['matched_name']} ({m['chembl_id']})", **m)
            if cache: cache.set(cache_key, result)
            return result

        # 2) Synonym 정확 매칭 (molecule_synonyms__molecule_synonym__iexact)
        # molecule/search?q= 대신 synonym DB에서 정확히 등록된 이름만 검색
        r2 = _SES.get(f"{CHEMBL_BASE}/molecule",
                      params={"molecule_synonyms__molecule_synonym__iexact": name,
                              "format": "json", "limit": 3},
                      timeout=15)
        mols2 = r2.json().get("molecules", []) if r2.ok else []
        if mols2:
            mols2.sort(key=lambda x: 0 if x.get("pref_name") else 1)
            m = _parse_mol(mols2[0])
            result.update(found=True, status="synonym",
                          note=f"동의어 매칭: {m['matched_name']} ({m['chembl_id']})",
                          **m)
            if cache: cache.set(cache_key, result)
            return result

        # 3) 부분 포함 검색 (정확도순 정렬)
        # 파생물 패널티: 쿼리보다 단어 수 많으면 순위 낮춤
        r3 = _SES.get(f"{CHEMBL_BASE}/molecule",
                      params={"pref_name__icontains": name, "format": "json", "limit": 10},
                      timeout=15)
        mols3 = r3.json().get("molecules", []) if r3.ok else []
        if mols3:
            q_words = len(name_upper.split())
            def rank(x):
                pn = (x.get("pref_name") or "").upper()
                extra = max(0, len(pn.split()) - q_words)
                if pn == name_upper: return (0, extra, len(pn))
                if pn.startswith(name_upper) and extra == 0: return (1, extra, len(pn))
                if pn.startswith(name_upper): return (3, extra, len(pn))
                return (2, extra, len(pn))
            mols3.sort(key=rank)
            # 품질 필터: 쿼리가 pref_name의 완전한 단어로 포함된 경우만 사용
            # (예: ISOLIQUIRITIGENIN 같이 쿼리가 단어 중간에 묻힌 경우 제외)
            best = mols3[0]
            pn_best = (best.get("pref_name") or "").upper()
            pn_words = set(pn_best.split())
            if pn_best == name_upper or name_upper in pn_words:
                m = _parse_mol(best)
                result.update(found=True, status="partial",
                              note=f"부분 매칭: {m['matched_name']} ({m['chembl_id']}) — 이름 확인 권장",
                              **m)
                if cache: cache.set(cache_key, result)
                return result
            # 품질 미달 → stage 4(PubChem)로 계속

        # 3.5) ChEMBL 오타 보정 — 글자 하나 빠진 변형으로 재검색
        for drop_idx in range(len(name) - 3, 3, -3):
            variant = name[:drop_idx] + name[drop_idx+1:]
            r35 = _SES.get(f"{CHEMBL_BASE}/molecule",
                           params={"pref_name__iexact": variant, "format": "json", "limit": 1},
                           timeout=10)
            mols35 = r35.json().get("molecules", []) if r35.ok else []
            if mols35:
                m = _parse_mol(mols35[0])
                result.update(found=True, status="synonym",
                              note=f"오타 보정 매칭: {m['matched_name']} ({m['chembl_id']})",
                              **m)
                if cache: cache.set(cache_key, result)
                return result

    except Exception as e:
        result["note"] = f"API 오류: {e}"

    # 4) PubChem → InChIKey → ChEMBL 매칭
    suggestions = _suggest_via_pubchem(name)
    if suggestions:
        # ChEMBL에서 찾은 경우 found=True로 처리
        s = suggestions[0]
        result.update(
            found=True,
            chembl_id=s["chembl_id"],
            matched_name=s["matched_name"],
            smiles=s["smiles"],
            status="synonym",
            note=f"PubChem→ChEMBL 매칭: {s['matched_name']} ({s['chembl_id']})",
        )
    else:
        result["note"] = "ChEMBL 미등록 — 분석 시 PubChem SMILES로 자동 보완"
        result["suggestions"] = []
    if cache: cache.set(cache_key, result)
    return result


def _suggest_via_pubchem(name: str) -> list:
    """
    PubChem에서 화합물명 조회 → InChIKey로 ChEMBL 매칭 시도.
    실패하면 PubChem 동의어 중 ChEMBL에 있는 이름 반환.
    반환: [{"chembl_id", "matched_name", "smiles", "source"}, ...]
    """
    results = []
    try:
        # PubChem CID + InChIKey + synonym 조회
        r = _SES.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/InChIKey,CanonicalSMILES,IUPACName/JSON",
            timeout=30,
        )
        if not r.ok:
            return results
        props = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
        inchikey = props.get("InChIKey", "")
        smiles   = props.get("CanonicalSMILES", "")
        cid      = props.get("CID")

        # PubChem synonym에서 사용 가능한 이름 가져오기
        pubchem_name = name
        if cid:
            try:
                rs = _SES.get(
                    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON",
                    timeout=15,
                )
                syns = rs.json().get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
                # 짧고 알파벳만 있는 첫 번째 synonym 사용
                for s in syns:
                    if s and len(s) <= 50 and s.replace(" ", "").replace("-", "").isalnum():
                        pubchem_name = s
                        break
            except Exception:
                pass

        # ChEMBL InChIKey 검색
        if inchikey:
            std_key = inchikey.split("-")[0]
            r2 = _SES.get(f"{CHEMBL_BASE}/molecule",
                          params={"molecule_structures__standard_inchi_key__istartswith": std_key,
                                  "format": "json", "limit": 3},
                          timeout=30)
            for mol in (r2.json().get("molecules", []) if r2.ok else []):
                structs = mol.get("molecule_structures") or {}
                matched = mol.get("pref_name") or pubchem_name
                results.append({
                    "chembl_id":    mol.get("molecule_chembl_id", ""),
                    "matched_name": matched,
                    "smiles":       structs.get("canonical_smiles", smiles),
                    "source":       "PubChem→InChIKey→ChEMBL",
                })
            if results:
                return results
    except Exception as e:
        log.debug(f"[suggest] {name}: {e}")
    return results


def validate_drug(name: str, cache=None) -> dict:
    """약물명 ChEMBL 검증 (validate_compound와 동일 로직)"""
    res = validate_compound(name, cache=cache)
    res["type"] = "drug"
    return res


# ─────────────────────────────────────────────────────────────────────────────
# 2. 질환 → OpenTargets 검증
# ─────────────────────────────────────────────────────────────────────────────

def validate_disease(query: str, cache=None) -> dict:
    """
    질환명 OpenTargets 검증.
    반환:
      {
        "query":        원본 입력,
        "found":        True/False,
        "disease_id":   "MONDO:...",
        "disease_name": 공식 질환명,
        "n_targets":    관련 타깃 수 (추정),
        "note":         메시지,
      }
    """
    cache_key = f"VALIDATE:disease:{query.lower().strip()}"
    if cache:
        hit = cache.get(cache_key)
        if hit is not None:
            return hit

    result = {
        "query": query, "found": False,
        "disease_id": "", "disease_name": "", "n_targets": 0, "note": "",
    }

    gql = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 3}) {
        hits {
          id
          name
          entity
          object { ... on Disease { id name associatedTargets { count } } }
        }
      }
    }
    """
    try:
        r = _SES.post(OT_GRAPHQL,
                      json={"query": gql, "variables": {"q": query}},
                      timeout=20)
        if r.ok:
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            for h in hits:
                obj = h.get("object", {})
                if obj.get("id"):
                    n = obj.get("associatedTargets", {}).get("count", 0)
                    result.update(
                        found=True,
                        disease_id=obj["id"],
                        disease_name=obj.get("name", h.get("name", "")),
                        n_targets=n,
                        note=f"OpenTargets 매칭: {obj.get('name','')} ({obj['id']}) — 관련 타깃 {n}개",
                    )
                    break
    except Exception as e:
        result["note"] = f"OpenTargets API 오류: {e}"

    if not result["found"]:
        result["note"] = "OpenTargets에서 찾을 수 없음 — 영문 질환명으로 입력 권장"

    if cache: cache.set(cache_key, result)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 3. 전체 config 검증 (병렬 아닌 순차, 캐시 우선)
# ─────────────────────────────────────────────────────────────────────────────

def validate_config(cfg: dict, cache=None) -> dict:
    """
    config 전체 검증.
    반환:
      {
        "herbs":   { herb_name: { compound_name: validate_compound 결과 } },
        "drugs":   { drug_name: validate_drug 결과 },
        "disease": validate_disease 결과 or None,
        "summary": { ok, warn, fail, total },
      }
    """
    report = {"herbs": {}, "drugs": {}, "disease": None}
    ok = warn = fail = 0

    for herb, compounds in (cfg.get("herbs") or {}).items():
        report["herbs"][herb] = {}
        for comp in compounds:
            time.sleep(0.1)  # rate limit
            v = validate_compound(comp, cache=cache)
            report["herbs"][herb][comp] = v
            if v["found"] and v["status"] == "exact":
                ok += 1
            elif v["found"]:
                warn += 1
            else:
                fail += 1

    for drug in (cfg.get("drugs") or []):
        time.sleep(0.1)
        v = validate_drug(drug, cache=cache)
        report["drugs"][drug] = v
        if v["found"] and v["status"] == "exact":
            ok += 1
        elif v["found"]:
            warn += 1
        else:
            fail += 1

    dq = cfg.get("disease_query", "")
    if dq:
        time.sleep(0.1)
        report["disease"] = validate_disease(dq, cache=cache)

    total = ok + warn + fail
    report["summary"] = {"ok": ok, "warn": warn, "fail": fail, "total": total}
    return report
