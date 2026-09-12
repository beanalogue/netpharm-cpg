"""
tcmsp_scraper.py — TCMSP-E 자동 스크래퍼
requests만으로 HTML에 임베드된 Kendo Grid JSON 데이터를 추출

흐름:
  1. herb name 검색 → herb_en_name 추출
  2. 약재 상세 페이지 → 화합물 목록 (OB, DL, ADME)
  3. 각 화합물 → molecule.php → 타깃 목록 (target_name, SVM_score, RF_score)

결과는 JSON으로 캐시 저장 (data/tcmsp_cache/)
"""
import json
import logging
import re
import time
from pathlib import Path
from urllib.parse import quote, urljoin

import requests

log = logging.getLogger(__name__)

_BASE = "https://tcmsp-e.com"
_CACHE_DIR = Path(__file__).parent.parent / "data" / "tcmsp_cache"
_SESSION_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
        "Chrome/120.0.0.0 Safari/537.36"
    ),
    "Accept-Language": "en-US,en;q=0.9",
}

# 약재명 → TCMSP 검색어 별칭 (학명이 검색 안 될 때 한어병음/약전명 시도)
# 형식: "표준 학명": ["TCMSP pinyin1", "약전명", ...]
_HERB_ALIASES: dict[str, list[str]] = {
    # ── 가미소요산 ──────────────────────────────────────────────
    "Bupleurum chinense":           ["Chaihu", "Radix Bupleuri"],
    "Paeonia lactiflora":           ["Baishao", "Radix Paeoniae Alba"],
    "Angelica sinensis":            ["Danggui", "Radix Angelicae Sinensis"],
    "Poria cocos":                  ["Fuling", "Poria"],
    "Atractylodes macrocephala":    ["Baizhu", "Rhizoma Atractylodis Macrocephalae"],
    "Glycyrrhiza uralensis":        ["Gancao", "Radix Glycyrrhizae", "licorice"],
    "Gardenia jasminoides":         ["Zhizi", "Fructus Gardeniae"],
    "Moutan cortex":                ["Mudanpi", "Cortex Moutan"],
    "Paeonia suffruticosa":         ["Mudanpi", "Cortex Moutan"],
    "Mentha haplocalyx":            ["Bohe", "Herba Menthae"],
    "Zingiber officinale":          ["Shengjiang", "Ginger"],

    # ── 귀비탕/가미귀비탕 ────────────────────────────────────────
    "Panax ginseng":                ["Renshen", "Ginseng", "Radix Ginseng"],
    "Astragalus membranaceus":      ["Huangqi", "Radix Astragali"],
    # 용안육/원지는 TCMSP DB 미보유 → ChEMBL/PubChem fallback
    # "Dimocarpus longan":          # TCMSP 없음
    # "Polygala tenuifolia":        # TCMSP 없음
    "Ziziphus spinosa":             ["Suanzaoren", "Semen Ziziphi Spinosae"],
    "Ziziphus jujuba var. spinosa": ["Suanzaoren", "Semen Ziziphi Spinosae"],
    "Aucklandia lappa":             ["Muxiang", "Radix Aucklandiae"],
    "Saussurea costus":             ["Muxiang", "Radix Aucklandiae"],
    "Ziziphus jujuba":              ["Dazao", "Jujube"],

    # ── 육미지황탕/팔미지황탕 ────────────────────────────────────
    "Rehmannia glutinosa":          ["Dihuang", "Radix Rehmanniae"],
    "Cornus officinalis":           ["Shanzhuyu", "Fructus Corni"],
    "Dioscorea opposita":           ["Shanyao", "Rhizoma Dioscoreae"],
    "Dioscorea opposita Thunb.":    ["Shanyao", "Rhizoma Dioscoreae"],
    "Alisma orientale":             ["Zexie", "Rhizoma Alismatis"],
    "Alisma plantago-aquatica":     ["Zexie", "Rhizoma Alismatis"],
    "Cinnamomum cassia":            ["Rougui", "Cortex Cinnamomi"],
    "Aconitum carmichaelii":        ["Fuzi", "Radix Aconiti Lateralis"],

    # ── 보중익기탕/사군자탕 ──────────────────────────────────────
    "Codonopsis pilosula":          ["Dangshen", "Radix Codonopsis"],
    "Cimicifuga foetida":           ["Shengma", "Rhizoma Cimicifugae"],
    "Citrus reticulata":            ["Chenpi", "Pericarpium Citri Reticulatae"],
    "Citrus aurantium":             ["Zhishi", "Fructus Aurantii Immaturus"],

    # ── 사물탕 ───────────────────────────────────────────────────
    "Ligusticum chuanxiong":        ["Chuanxiong", "Rhizoma Chuanxiong"],
    "Paeonia lactiflora (Baishao)": ["Baishao", "Radix Paeoniae Alba"],

    # ── 황련해독탕/온청음 ─────────────────────────────────────────
    "Scutellaria baicalensis":      ["Huangqin", "Radix Scutellariae"],
    "Coptis chinensis":             ["Huanglian", "Rhizoma Coptidis"],
    "Phellodendron amurense":       ["Huangbai", "Cortex Phellodendri"],
    "Phellodendron chinense":       ["Huangbai", "Cortex Phellodendri"],

    # ── 반하사심탕/온담탕 ─────────────────────────────────────────
    "Pinellia ternata":             ["Banxia", "Rhizoma Pinelliae"],
    "Bambusa tuldoides":            ["Zhuru", "Caulis Bambusae In Taenia"],

    # ── 소시호탕/시호제 ───────────────────────────────────────────
    "Bupleurum scorzonerifolium":   ["Chaihu", "Radix Bupleuri"],

    # ── 마황탕/갈근탕 ─────────────────────────────────────────────
    "Ephedra sinica":               ["Mahuang", "Herba Ephedrae"],
    "Pueraria lobata":              ["Gegen", "Radix Puerariae Lobatae"],
    "Pueraria montana":             ["Gegen", "Radix Puerariae Lobatae"],
    "Prunus armeniaca":             ["Xingren", "Semen Armeniacae Amarum"],

    # ── 소풍산/형방제 ─────────────────────────────────────────────
    "Schizonepeta tenuifolia":      ["Jingjie", "Herba Schizonepetae"],
    "Saposhnikovia divaricata":     ["Fangfeng", "Radix Saposhnikoviae"],
    "Angelica dahurica":            ["Baizhi", "Radix Angelicae Dahuricae"],

    # ── 독활기생탕/관절처방 ───────────────────────────────────────
    "Angelica pubescens":           ["Duhuo", "Radix Angelicae Pubescentis"],
    "Loranthus parasiticus":        ["Sangjisheng", "Herba Taxilli"],
    "Taxillus chinensis":           ["Sangjisheng", "Herba Taxilli"],
    "Eucommia ulmoides":            ["Duzhong", "Cortex Eucommiae"],
    "Achyranthes bidentata":        ["Niuxi", "Radix Achyranthis Bidentatae"],
    "Gentiana macrophylla":         ["Qinjiao", "Radix Gentianae Macrophyllae"],
    "Asarum sieboldii":             ["Xixin", "Radix Et Rhizoma Asari"],

    # ── 혈부축어탕/어혈처방 ───────────────────────────────────────
    "Prunus persica":               ["Taoren", "Semen Persicae"],
    "Carthamus tinctorius":         ["Honghua", "Flos Carthami"],
    "Curcuma longa":                ["Jianghuang", "Rhizoma Curcumae Longae"],
    "Curcuma wenyujin":             ["Yujin", "Radix Curcumae"],

    # ── 소화/위장처방 ─────────────────────────────────────────────
    "Crataegus pinnatifida":        ["Shanzha", "Fructus Crataegi"],
    "Massa Medicata Fermentata":    ["Shenqu", "Massa Medicata Fermentata"],
    "Hordeum vulgare":              ["Maiya", "Fructus Hordei Germinatus"],
    "Amomum villosum":              ["Sharen", "Fructus Amomi"],
    "Magnolia officinalis":         ["Houpo", "Cortex Magnoliae Officinalis"],

    # ── 진해/폐처방 ───────────────────────────────────────────────
    "Platycodon grandiflorum":      ["Jiegeng", "Radix Platycodonis"],
    "Fritillaria thunbergii":       ["Zhebeimu", "Bulbus Fritillariae Thunbergii"],
    "Fritillaria cirrhosa":         ["Chuanbeimu", "Bulbus Fritillariae Cirrhosae"],
    "Stemona japonica":             ["Baibu", "Radix Stemonae"],
    "Aster tataricus":              ["Ziwan", "Radix Et Rhizoma Asteris"],
    "Tussilago farfara":            ["Kuandonghua", "Flos Farfarae"],

    # ── 이뇨/신처방 ───────────────────────────────────────────────
    "Plantago asiatica":            ["Cheqianzi", "Semen Plantaginis"],
    "Dianthus superbus":            ["Qumai", "Herba Dianthi"],
    "Lysimachia christinae":        ["Jinqiancao", "Herba Lysimachiae"],
    "Pyrrosia lingua":              ["Shiwei", "Folium Pyrrosiae"],

    # ── 보음/보양 ─────────────────────────────────────────────────
    "Ophiopogon japonicus":         ["Maidong", "Radix Ophiopogonis"],
    "Asparagus cochinchinensis":    ["Tiandong", "Radix Asparagi"],
    "Ligustrum lucidum":            ["Nvzhenzi", "Fructus Ligustri Lucidi"],
    "Eclipta prostrata":            ["Mohanlian", "Herba Ecliptae"],
    "Cuscuta chinensis":            ["Tusizi", "Semen Cuscutae"],
    "Epimedium brevicornum":        ["Yinyanghuo", "Herba Epimedii"],
    "Morinda officinalis":          ["Bajitian", "Radix Morindae Officinalis"],
    "Cistanche deserticola":        ["Roucongrong", "Herba Cistanches"],
    "Dipsacus asper":               ["Xuduan", "Radix Dipsaci"],
    "Drynaria fortunei":            ["Gusuibu", "Rhizoma Drynariae"],

    # ── 간/담처방 ─────────────────────────────────────────────────
    "Gentiana scabra":              ["Longdancao", "Radix Et Rhizoma Gentianae"],
    "Alisma orientale (Sam.) Juz.": ["Zexie", "Rhizoma Alismatis"],
    "Silybum marianum":             ["Shuifei Ji", "Fructus Silybi"],

    # ── 안신/심처방 ───────────────────────────────────────────────
    "Acorus tatarinowii":           ["Shichangpu", "Rhizoma Acori Tatarinowii"],
    "Semen Ziziphi Spinosae":       ["Suanzaoren"],
    "Biota orientalis":             ["Baiziren", "Semen Platycladi"],
    "Platycladus orientalis":       ["Baiziren", "Semen Platycladi"],
    "Albizzia julibrissin":         ["Hehuanpi", "Cortex Albiziae"],

    # ── 지혈처방 ───────────────────────────────────────────────────
    "Agrimonia pilosa":             ["Xianhecao", "Herba Agrimoniae"],
    "Sanguisorba officinalis":      ["Diyu", "Radix Sanguisorbae"],
    "Cirsium japonicum":            ["Daji", "Herba Cirsii Japonici"],

    # ── 기타 자주 쓰이는 약재 ────────────────────────────────────
    "Coptis japonica":              ["Huanglian", "Rhizoma Coptidis"],
    "Salvia miltiorrhiza":          ["Danshen", "Radix Et Rhizoma Salviae Miltiorrhizae"],
    "Panax notoginseng":            ["Sanqi", "Radix Et Rhizoma Notoginseng"],
    "Astragalus propinquus":        ["Huangqi", "Radix Astragali"],
    "Cinnamomum verum":             ["Rougui", "Cortex Cinnamomi"],
    "Alpinia oxyphylla":            ["Yizhi", "Fructus Alpiniae Oxyphyllae"],
    "Lycium chinense":              ["Gouqizi", "Fructus Lycii"],
    "Lycium barbarum":              ["Gouqizi", "Fructus Lycii"],
    "Chrysanthemum morifolium":     ["Juhua", "Flos Chrysanthemi"],
    "Schisandra chinensis":         ["Wuweizi", "Fructus Schisandrae Chinensis"],
    "Morus alba":                   ["Sangye", "Folium Mori"],
    "Coptis teeta":                 ["Huanglian", "Rhizoma Coptidis"],
}

# 모듈 레벨 세션 캐시 (토큰 재사용)
_shared_sess: requests.Session | None = None
_shared_token: str = ""


# ─────────────────────────────────────────────────────────────────
# 내부 유틸
# ─────────────────────────────────────────────────────────────────

def _session() -> requests.Session:
    s = requests.Session()
    s.headers.update(_SESSION_HEADERS)
    return s


def _get_shared_session() -> tuple[requests.Session, str]:
    """모듈 레벨 세션+토큰 공유 (요청 최소화)"""
    global _shared_sess, _shared_token
    if _shared_sess is None:
        _shared_sess = _session()
    if not _shared_token:
        _shared_token = _get_token(_shared_sess)
    return _shared_sess, _shared_token


def _get_token(sess: requests.Session) -> str:
    """메인 페이지에서 CSRF token 획득"""
    try:
        r = sess.get(f"{_BASE}/tcmsp.php", timeout=15)
        m = re.search(r'name="token"[^>]+value="([a-f0-9]+)"', r.text)
        if m:
            return m.group(1)
        # 응답에서 token 파라미터 찾기 (fallback)
        m2 = re.search(r'token=([a-f0-9]{32})', r.text)
        return m2.group(1) if m2 else ""
    except Exception as e:
        log.warning(f"[TCMSP] 토큰 획득 실패: {e}")
        return ""


_herb_list_cache: dict | None = None


def _load_herb_list() -> dict:
    """data/tcmsp_herb_list.json 로드 (모듈 레벨 캐시)"""
    global _herb_list_cache
    if _herb_list_cache is not None:
        return _herb_list_cache
    p = Path(__file__).parent.parent / "data" / "tcmsp_herb_list.json"
    if p.exists():
        try:
            _herb_list_cache = json.loads(p.read_text(encoding="utf-8"))
            return _herb_list_cache
        except Exception:
            pass
    _herb_list_cache = {}
    return _herb_list_cache


def _fuzzy_match_herb_list(herb_query: str) -> str | None:
    """
    herb_list.json의 herb_en_name(약전명)에서 genus/species 퍼지 매칭
    예: 'Angelica sinensis' → 'Angelicae Sinensis Radix' → pinyin 'Danggui'
    반환: 매칭된 pinyin 또는 None
    """
    herb_db = _load_herb_list()
    if not herb_db:
        return None

    parts = herb_query.lower().split()
    if not parts:
        return None

    genus = parts[0].rstrip("e")   # 속격 변화 보정: Angelicae→Angelica, Morindae→Morinda
    species = parts[1].rstrip("is") if len(parts) > 1 else ""  # sinensis→sinen

    best = None
    best_score = 0
    for en_name, herb_data in herb_db.items():
        en_lower = en_name.lower()
        score = 0
        if genus and genus in en_lower:
            score += 2
        if species and species in en_lower:
            score += 1
        if score > best_score:
            best_score = score
            best = herb_data.get("herb_pinyin", "")

    if best_score >= 2:  # genus 매칭은 필수
        return best
    return None


def _resolve_herb_query(herb_query: str, sess: requests.Session, token: str) -> list[dict]:
    """
    약재명 검색 우선순위:
    1. 입력값 직접 검색
    2. _HERB_ALIASES 명시적 별칭
    3. herb_list.json 퍼지 매칭 (genus/species → 약전명 → pinyin)
    """
    # 1. 직접 검색
    results = search_herb(herb_query, sess, token)
    if results:
        return results

    # 2. 명시적 별칭
    for alias in _HERB_ALIASES.get(herb_query, []):
        results = search_herb(alias, sess, token)
        if results:
            log.info(f"[TCMSP] '{herb_query}' → 별칭 '{alias}'으로 발견")
            return results

    # 3. herb_list.json 퍼지 매칭
    pinyin = _fuzzy_match_herb_list(herb_query)
    if pinyin:
        results = search_herb(pinyin, sess, token)
        if results:
            log.info(f"[TCMSP] '{herb_query}' → 퍼지 매칭 '{pinyin}'으로 발견")
            return results

    log.warning(f"[TCMSP] '{herb_query}': 검색 실패 (TCMSP 미보유 약재)")
    return []


def _extract_kendo_data(html: str) -> list[dict]:
    """HTML에서 kendoGrid data 배열 추출"""
    m = re.search(r'data:\s*(\[.*?\]),\s*\n?\s*pageSize', html, re.DOTALL)
    if not m:
        return []
    try:
        return json.loads(m.group(1))
    except json.JSONDecodeError as e:
        log.warning(f"[TCMSP] JSON 파싱 실패: {e}")
        return []


def _extract_var(html: str, var_name: str):
    """var <name> = <value>; 형식의 JS 변수 값 추출"""
    m = re.search(rf'var {re.escape(var_name)}\s*=\s*(.*?);', html, re.DOTALL)
    if not m:
        return None
    val = m.group(1).strip()
    if val == "null":
        return None
    try:
        return json.loads(val)
    except json.JSONDecodeError:
        return None


def _safe_get(sess, url, referer="", retries=2, delay=2) -> requests.Response | None:
    """재시도 포함 GET"""
    for attempt in range(retries + 1):
        try:
            headers = {"Referer": referer or _BASE} if referer else {}
            r = sess.get(url, headers=headers, timeout=20)
            r.raise_for_status()
            return r
        except Exception as e:
            if attempt < retries:
                log.debug(f"[TCMSP] 재시도 {attempt+1}/{retries}: {e}")
                time.sleep(delay)
            else:
                log.warning(f"[TCMSP] 요청 실패: {url} → {e}")
    return None


# ─────────────────────────────────────────────────────────────────
# 캐시
# ─────────────────────────────────────────────────────────────────

def _cache_path(herb_name: str) -> Path:
    safe = re.sub(r'[^\w\s-]', '', herb_name).strip().replace(' ', '_')
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return _CACHE_DIR / f"{safe}.json"


def _load_cache(herb_name: str) -> dict | None:
    p = _cache_path(herb_name)
    if p.exists():
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
            # OB 필터가 적용된 구버전 캐시는 무효화 → 전체 성분으로 재수집
            if data.get("ob_cutoff", 0) > 0:
                log.info(f"[TCMSP] 구버전 캐시 무효화(OB필터 제거): {herb_name}")
                return None
            return data
        except Exception:
            pass
    return None


def _save_cache(herb_name: str, data: dict) -> None:
    p = _cache_path(herb_name)
    try:
        p.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
        log.info(f"[TCMSP] 캐시 저장: {p}")
    except Exception as e:
        log.warning(f"[TCMSP] 캐시 저장 실패: {e}")


# ─────────────────────────────────────────────────────────────────
# 핵심 스크래핑 함수
# ─────────────────────────────────────────────────────────────────

def search_herb(herb_query: str, sess: requests.Session, token: str) -> list[dict]:
    """
    약재명(영문/중문/한어병음)으로 검색
    반환: [{"herb_cn_name", "herb_en_name", "herb_pinyin"}, ...]
    """
    url = f"{_BASE}/tcmspsearch.php?qs=herb_all_name&q={quote(herb_query)}&token={token}"
    r = _safe_get(sess, url, referer=f"{_BASE}/tcmsp.php")
    if r is None:
        return []
    return _extract_kendo_data(r.text)


def get_herb_compounds(herb_en_name: str, sess: requests.Session, token: str,
                       ob_cutoff: float = 30.0, dl_cutoff: float = 0.18) -> list[dict]:
    """
    약재 화합물 목록 가져오기 (OB/DL 필터링 적용)
    반환: [{molecule_ID, MOL_ID, molecule_name, ob, dl, mw, alogp, ...}, ...]
    """
    url = f"{_BASE}/tcmspsearch.php?qr={quote(herb_en_name)}&qsr=herb_en_name&token={token}"
    r = _safe_get(sess, url, referer=f"{_BASE}/tcmsp.php")
    if r is None:
        return []

    all_compounds = _extract_kendo_data(r.text)
    log.info(f"[TCMSP] {herb_en_name}: 전체 {len(all_compounds)}개 화합물")

    filtered = []
    for c in all_compounds:
        try:
            ob = float(c.get("ob", 0))
            dl = float(c.get("dl", 0))
            if ob >= ob_cutoff and dl >= dl_cutoff:
                filtered.append(c)
        except (ValueError, TypeError):
            pass

    log.info(f"[TCMSP] {herb_en_name}: OB≥{ob_cutoff}, DL≥{dl_cutoff} → {len(filtered)}개")
    return filtered


def get_compound_targets(molecule_id: str, sess: requests.Session,
                         score_cutoff: float = 0.0) -> list[dict]:
    """
    화합물 타깃 가져오기
    반환: [{target_ID, target_name, drugbank_ID, SVM_score, RF_score}, ...]
    score_cutoff: SVM_score 필터 (기본값 0 = 필터 없음)
    """
    url = f"{_BASE}/molecule.php?qn={molecule_id}"
    r = _safe_get(sess, url)
    if r is None:
        return []

    targets = _extract_var(r.text, "tar")
    if not targets:
        return []

    if score_cutoff > 0:
        targets = [t for t in targets
                   if float(t.get("SVM_score", 0)) >= score_cutoff]

    return targets


# ─────────────────────────────────────────────────────────────────
# 타깃명 → 유전자 심볼 변환
# ─────────────────────────────────────────────────────────────────

# 간단한 로컬 매핑 (자주 등장하는 TCMSP 타깃들)
_TARGET_NAME_TO_GENE: dict[str, str] = {
    "Gamma-aminobutyric-acid receptor alpha-2 subunit": "GABRA2",
    "Gamma-aminobutyric-acid receptor alpha-1 subunit": "GABRA1",
    "Gamma-aminobutyric-acid receptor beta-2 subunit": "GABRB2",
    "Prostaglandin G/H synthase 2": "PTGS2",
    "Prostaglandin G/H synthase 1": "PTGS1",
    "Nitric oxide synthase, inducible": "NOS2",
    "Nitric oxide synthase, endothelial": "NOS3",
    "RAC-alpha serine/threonine-protein kinase": "AKT1",
    "RAC-beta serine/threonine-protein kinase": "AKT2",
    "Tumor necrosis factor": "TNF",
    "Interleukin-6": "IL6",
    "Interleukin-1 beta": "IL1B",
    "Vascular endothelial growth factor A": "VEGFA",
    "Epidermal growth factor receptor": "EGFR",
    "Serine/threonine-protein kinase mTOR": "MTOR",
    "Phosphatidylinositol 3-kinase regulatory subunit alpha": "PIK3R1",
    "Caspase-3": "CASP3",
    "Caspase-9": "CASP9",
    "Tumor protein p53": "TP53",
    "Cellular tumor antigen p53": "TP53",
    "Androgen receptor": "AR",
    "Estrogen receptor alpha": "ESR1",
    "Estrogen receptor": "ESR1",
    "Peroxisome proliferator-activated receptor gamma": "PPARG",
    "Peroxisome proliferator-activated receptor alpha": "PPARA",
    "Nuclear receptor subfamily 3 group C member 1": "NR3C1",
    "Glucocorticoid receptor": "NR3C1",
    "Matrix metalloproteinase-9": "MMP9",
    "Matrix metalloproteinase-2": "MMP2",
    "Mitogen-activated protein kinase 1": "MAPK1",
    "Mitogen-activated protein kinase 3": "MAPK3",
    "Mitogen-activated protein kinase 8": "MAPK8",
    "Mitogen-activated protein kinase 14": "MAPK14",
    "Signal transducer and activator of transcription 3": "STAT3",
    "Signal transducer and activator of transcription 1": "STAT1",
    "Proto-oncogene tyrosine-protein kinase Src": "SRC",
    "Glycogen synthase kinase-3 beta": "GSK3B",
    "Bcl-2 homologous antagonist/killer": "BAK1",
    "Apoptosis regulator Bcl-2": "BCL2",
    "Apoptosis regulator BAX": "BAX",
    "Heat shock protein HSP 90-alpha": "HSP90AA1",
    "Cyclooxygenase-2": "PTGS2",
    "Acetylcholinesterase": "ACHE",
    "Butyrylcholinesterase": "BCHE",
    "Beta-secretase 1": "BACE1",
    "Arachidonate 5-lipoxygenase": "ALOX5",
    "Transforming growth factor beta-1": "TGFB1",
    "Mothers against decapentaplegic homolog 3": "SMAD3",
    "Retinoic acid receptor ROR-gamma": "RORC",
    "Vitamin D3 receptor": "VDR",
    "Androgen-binding protein": "AR",
}


def target_name_to_gene(target_name: str) -> str:
    """타깃명 → 유전자 심볼 (로컬 매핑 우선, 실패 시 UniProt 검색)"""
    # 정확 매칭
    gene = _TARGET_NAME_TO_GENE.get(target_name)
    if gene:
        return gene

    # 부분 매칭
    tn_lower = target_name.lower()
    for k, v in _TARGET_NAME_TO_GENE.items():
        if k.lower() in tn_lower or tn_lower in k.lower():
            return v

    # UniProt 검색 (rate limit 고려하여 마지막 수단)
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": f'protein_name:"{target_name}" AND organism_id:9606',
                    "fields": "gene_names", "format": "json", "size": 1},
            timeout=10,
        )
        results = r.json().get("results", [])
        if results:
            genes = results[0].get("genes", [])
            if genes:
                gn = genes[0].get("geneName", {}).get("value", "")
                if gn:
                    log.debug(f"[UniProt] {target_name} → {gn}")
                    _TARGET_NAME_TO_GENE[target_name] = gn  # 런타임 캐시
                    return gn
    except Exception:
        pass

    return ""


# ─────────────────────────────────────────────────────────────────
# 공개 API
# ─────────────────────────────────────────────────────────────────

def scrape_herb(
    herb_query: str,
    ob_cutoff: float = 0.0,
    dl_cutoff: float = 0.0,
    svm_cutoff: float = 0.0,
    use_cache: bool = True,
    inter_request_delay: float = 0.5,
) -> dict:
    """
    약재 이름으로 TCMSP-E 스크래핑

    Parameters
    ----------
    herb_query : str
        검색어 (영문/한어병음/한자 모두 가능)
        예: "Poria", "Fuling", "茯苓", "Atractylodes macrocephala"
    ob_cutoff : float
        Oral Bioavailability 최소값 (기본 0 = 필터 없음, 전체 수집)
    dl_cutoff : float
        Drug-likeness 최소값 (기본 0 = 필터 없음, 전체 수집)
    svm_cutoff : float
        타깃 SVM score 필터 (기본 0 = 필터 없음)
    use_cache : bool
        캐시 사용 여부 (True 권장)
    inter_request_delay : float
        요청 간 지연 (초, 기본 0.5)

    Returns
    -------
    dict with keys:
        herb_query    : 검색어
        herb_cn_name  : 한자명
        herb_en_name  : 라틴명
        herb_pinyin   : 한어병음
        compounds     : [{MOL_ID, molecule_name, ob, dl, mw, targets: [{...}]}, ...]
        total_compounds : 전체 화합물 수
        active_compounds: ob_cutoff/dl_cutoff 조건 충족 화합물 수 (참고용)
        source        : "tcmsp"
    """
    # 캐시 확인
    if use_cache:
        cached = _load_cache(herb_query)
        if cached:
            log.info(f"[TCMSP] 캐시 로드: {herb_query}")
            return cached

    # 공유 세션 + 토큰 사용 (요청 최소화)
    sess, token = _get_shared_session()
    if not token:
        log.warning("[TCMSP] 토큰 없음 — 재시도")
        global _shared_token
        _shared_token = _get_token(sess)
        token = _shared_token
    if not token:
        return {"herb_query": herb_query, "compounds": [], "error": "token_failed"}

    # 2. 약재 검색 (별칭 포함)
    log.info(f"[TCMSP] 약재 검색: '{herb_query}'")
    herbs = _resolve_herb_query(herb_query, sess, token)
    if not herbs:
        log.warning(f"[TCMSP] '{herb_query}' 검색 결과 없음")
        return {"herb_query": herb_query, "compounds": [], "error": "herb_not_found"}

    herb = herbs[0]
    herb_en = herb.get("herb_en_name", "")
    log.info(f"[TCMSP] 발견: {herb.get('herb_cn_name')} / {herb.get('herb_pinyin')} / {herb_en}")

    time.sleep(inter_request_delay)

    # 3. 화합물 목록 (필터 적용)
    # 전체 화합물 수도 필요하므로 별도 요청
    url = f"{_BASE}/tcmspsearch.php?qr={quote(herb_en)}&qsr=herb_en_name&token={token}"
    r = _safe_get(sess, url, referer=f"{_BASE}/tcmsp.php")
    all_comps_raw = _extract_kendo_data(r.text) if r else []
    total_count = len(all_comps_raw)

    # ob_cutoff/dl_cutoff가 0이면 전체 수집 (필터 없음)
    if ob_cutoff > 0 or dl_cutoff > 0:
        comps_to_fetch = [c for c in all_comps_raw if _passes_filter(c, ob_cutoff, dl_cutoff)]
    else:
        comps_to_fetch = all_comps_raw
    # 참고용 active_compounds: OB≥30, DL≥0.18 기준 (논문 기준)
    n_active_ref = sum(1 for c in all_comps_raw if _passes_filter(c, 30.0, 0.18))
    log.info(f"[TCMSP] 전체 {total_count}개 수집 (논문기준 활성 {n_active_ref}개)")

    # 4. 각 화합물 타깃 조회
    compounds_out = []
    for i, comp in enumerate(comps_to_fetch):
        if not isinstance(comp, dict):
            continue
        mol_id   = comp.get("molecule_ID") or ""
        mol_name = comp.get("molecule_name") or ""
        log.info(f"[TCMSP] 타깃 조회 [{i+1}/{len(comps_to_fetch)}]: {mol_name[:50]}")

        time.sleep(inter_request_delay)
        if not mol_id:
            raw_targets = []
        else:
            try:
                raw_targets = get_compound_targets(mol_id, sess, svm_cutoff) or []
            except Exception as te:
                log.warning(f"[TCMSP] 타깃 조회 실패 ({mol_name[:30]}): {te}")
                raw_targets = []

        # 유전자 심볼 변환
        targets_out = []
        for t in (raw_targets or []):
            if not isinstance(t, dict):
                continue
            try:
                gene = target_name_to_gene(t.get("target_name", ""))
            except Exception:
                gene = ""
            targets_out.append({
                "target_name": t.get("target_name", ""),
                "gene_symbol": gene,
                "TAR_ID":      t.get("TAR_ID", ""),
                "target_ID":   t.get("target_ID", ""),
                "drugbank_ID": t.get("drugbank_ID", ""),
                "SVM_score":   t.get("SVM_score", ""),
                "RF_score":    t.get("RF_score", ""),
                "source":      "tcmsp",
            })

        compounds_out.append({
            "MOL_ID":        comp.get("MOL_ID") or "",
            "molecule_ID":   mol_id,
            "molecule_name": mol_name,
            "ob":       comp.get("ob") or "",
            "dl":       comp.get("dl") or "",
            "mw":       comp.get("mw") or "",
            "alogp":    comp.get("alogp") or "",
            "hdon":     comp.get("hdon") or "",
            "hacc":     comp.get("hacc") or "",
            "bbb":      comp.get("bbb") or "",
            "caco2":    comp.get("caco2") or "",
            "halflife": comp.get("halflife") or "",
            "targets":  targets_out,
        })

    result = {
        "herb_query": herb_query,
        "herb_cn_name": herb.get("herb_cn_name", ""),
        "herb_en_name": herb_en,
        "herb_pinyin": herb.get("herb_pinyin", ""),
        "compounds": compounds_out,
        "total_compounds": total_count,
        "active_compounds": n_active_ref,   # 참고용: OB≥30, DL≥0.18 기준
        "ob_cutoff": 0.0,                   # 필터 없음 표시
        "dl_cutoff": 0.0,
        "source": "tcmsp",
    }

    if use_cache:
        _save_cache(herb_query, result)

    return result


def _passes_filter(comp: dict, ob_cutoff: float, dl_cutoff: float) -> bool:
    try:
        return float(comp.get("ob", 0)) >= ob_cutoff and float(comp.get("dl", 0)) >= dl_cutoff
    except (ValueError, TypeError):
        return False


# ─────────────────────────────────────────────────────────────────
# collect.py 연동용 함수
# ─────────────────────────────────────────────────────────────────

def get_tcmsp_targets_for_compound(
    compound_name: str,
    herb_query: str,
    ob_cutoff: float = 0.0,
    dl_cutoff: float = 0.0,
) -> list[dict]:
    """
    특정 화합물의 TCMSP 타깃 반환 (pubchem_target.supplement_targets에서 호출)

    Parameters
    ----------
    compound_name : str
        화합물 이름 (부분 매칭)
    herb_query : str
        소속 약재명 (캐시 키)

    Returns
    -------
    list of {gene_symbol, target_name, source, ...}
    """
    data = scrape_herb(herb_query, ob_cutoff=ob_cutoff, dl_cutoff=dl_cutoff)
    compounds = data.get("compounds", [])

    # 화합물 이름 매칭
    comp_lower = compound_name.lower()
    for c in compounds:
        cname = c.get("molecule_name", "").lower()
        if comp_lower in cname or cname in comp_lower:
            targets = c.get("targets", [])
            log.info(f"[TCMSP] {compound_name}: {len(targets)}개 타깃")
            # gene_symbol이 있는 것만 반환
            return [t for t in targets if t.get("gene_symbol")]

    log.info(f"[TCMSP] '{compound_name}' → 약재 '{herb_query}'에서 찾지 못함")
    return []


# ─────────────────────────────────────────────────────────────────
# CLI 테스트용
# ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    query = sys.argv[1] if len(sys.argv) > 1 else "Fuling"
    print(f"\n=== TCMSP 스크래핑: {query} ===\n")

    result = scrape_herb(query, use_cache=False)

    print(f"약재: {result.get('herb_cn_name')} / {result.get('herb_pinyin')} / {result.get('herb_en_name')}")
    print(f"전체 화합물: {result.get('total_compounds')}")
    print(f"활성 화합물: {result.get('active_compounds')}")
    print()

    for c in result.get("compounds", []):
        tcount = len(c.get("targets", []))
        gene_count = sum(1 for t in c.get("targets", []) if t.get("gene_symbol"))
        print(f"  [{c['MOL_ID']}] {c['molecule_name'][:55]:55s} OB={c['ob']:6} DL={c['dl']:5} → {tcount}타깃 ({gene_count}개 gene)")
        for t in c.get("targets", [])[:3]:
            print(f"    - {t.get('gene_symbol') or '?':10s} | {t['target_name'][:60]}")

    if result.get("error"):
        print(f"\nError: {result['error']}")
