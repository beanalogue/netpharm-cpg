"""
약재 발굴 모드 (Reverse Network Pharmacology)
질환명 입력 → TCMSP 캐시 전체 스코어링 → 유효 가능 약재 랭킹

알고리즘:
  1. TCMSP 캐시에서 모든 약재의 활성 성분 타깃(OB/DL 필터) 로드
  2. OpenTargets에서 질환 관련 타깃 유전자 수집
  3. 단일 스코어링 (hypergeometric p-value + 정규화 점수)
  4. 조합 스코어링 (Top N 내 2~4개 조합)
  5. Greedy Set Cover (커버리지 최대화 약재 순차 선택)
"""

import json
import math
import logging
from itertools import combinations
from pathlib import Path
from typing import Callable, Optional

log = logging.getLogger(__name__)

BASE_DIR  = Path(__file__).parent.parent
CACHE_DIR = BASE_DIR / "data" / "tcmsp_cache"

# hypergeometric test 배경: 인간 단백질코딩 유전자 수
HUMAN_GENE_BACKGROUND = 20_000


def _safe_float(v, default: float = 0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _hypergeom_pvalue(k: int, H: int, D: int, M: int = HUMAN_GENE_BACKGROUND) -> float:
    """
    hypergeometric test p-value
    k=오버랩, H=약재타깃수, D=질환타깃수, M=배경유전자수
    """
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(k - 1, M, D, H))
    except Exception:
        return 1.0


# ─────────────────────────────────────────────────────────────────────────────
# 1. TCMSP 캐시 로드
# ─────────────────────────────────────────────────────────────────────────────
def load_all_herb_targets(
    ob_cutoff: float = 30.0,
    dl_cutoff: float = 0.18,
) -> dict:
    """
    TCMSP 캐시 로드 → 약재별 활성 성분 타깃 집합
    tcmsp_herb_list.json 등록 500종만 사용 (재현성 보장)
    반환: {pinyin: {herb_info, targets(set), active_compounds, ...}}
    """
    # 공식 500종 병음 목록
    herb_list_path = BASE_DIR / "data" / "tcmsp_herb_list.json"
    if herb_list_path.exists():
        _hl = json.loads(herb_list_path.read_text(encoding="utf-8"))
        official_pinyins = {v.get("herb_pinyin", "") for v in _hl.values()}
    else:
        official_pinyins = None  # 파일 없으면 필터 없이 전체 로드

    herb_map = {}
    cache_files = list(CACHE_DIR.glob("*.json"))
    log.info(f"[발굴] TCMSP 캐시 로드: {len(cache_files)}개 파일 (공식 500종 필터 적용)")

    for path in cache_files:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if data.get("error"):
            continue

        pinyin = data.get("herb_pinyin", path.stem)

        # 공식 500종 외 제외 (JSON 내부 pinyin 또는 파일명으로 매칭)
        if official_pinyins is not None:
            if pinyin not in official_pinyins and path.stem not in official_pinyins:
                continue
            # 파일명이 매칭되면 pinyin을 파일명으로 통일 (herb_list 기준 키 사용)
            if pinyin not in official_pinyins and path.stem in official_pinyins:
                pinyin = path.stem
        compounds = data.get("compounds", [])

        active = [
            c for c in compounds
            if _safe_float(c.get("ob", 0)) >= ob_cutoff
            and _safe_float(c.get("dl", 0)) >= dl_cutoff
        ]

        targets = set()
        for comp in active:
            for t in comp.get("targets", []):
                gs = t.get("gene_symbol", "").strip().upper()
                if gs:
                    targets.add(gs)

        if not targets:
            continue

        herb_map[pinyin] = {
            "pinyin":           pinyin,
            "cn_name":          data.get("herb_cn_name", ""),
            "en_name":          data.get("herb_en_name", ""),
            "targets":          targets,
            "target_count":     len(targets),
            "active_count":     len(active),
            "total_count":      len(compounds),
            "active_compounds": [c.get("molecule_name", "") for c in active],
        }

    log.info(f"[발굴] 타깃 보유 약재: {len(herb_map)}개")
    return herb_map


# ─────────────────────────────────────────────────────────────────────────────
# 2. 단일 약재 스코어링
# ─────────────────────────────────────────────────────────────────────────────
def score_herbs(
    herb_map: dict,
    disease_genes: set,
    min_overlap: int = 2,
    top_n: int = 20,
) -> list:
    """
    약재-질환 타깃 오버랩 스코어링 → 상위 top_n 반환
    정렬 기준: hypergeometric p-value (오름차순) → overlap_count (내림차순)
    """
    disease_set = {g.upper() for g in disease_genes}
    D = len(disease_set)
    if D == 0:
        return []

    results = []
    for pinyin, info in herb_map.items():
        herb_set = info["targets"]
        overlap  = herb_set & disease_set
        k = len(overlap)
        if k < min_overlap:
            continue

        H        = len(herb_set)
        jaccard  = k / (H + D - k)
        coverage = k / D
        score    = k / math.sqrt(H * D) if H > 0 else 0.0
        pval     = _hypergeom_pvalue(k, H, D)

        results.append({
            **info,
            "overlap_genes": sorted(overlap),
            "overlap_count": k,
            "jaccard":       round(jaccard, 4),
            "coverage":      round(coverage, 4),
            "score":         round(score, 4),
            "pvalue":        pval,
        })

    results.sort(key=lambda x: (x["pvalue"], -x["overlap_count"]))
    return results[:top_n]


# ─────────────────────────────────────────────────────────────────────────────
# 3. 조합 스코어링 (Top N 약재 내 k-조합)
# ─────────────────────────────────────────────────────────────────────────────
def score_herb_combinations(
    herb_map: dict,
    disease_genes: set,
    top_n_single: int = 30,
    combo_size: int = 2,
    top_n_combo: int = 10,
    min_overlap: int = 3,
) -> list:
    """
    상위 top_n_single개 약재 내에서 combo_size 조합 전수 스코어링
    synergy = 조합 오버랩 - 단일 최대 오버랩  (추가 커버리지)
    """
    disease_set = {g.upper() for g in disease_genes}
    D = len(disease_set)
    if D == 0:
        return []

    # 단일 상위 약재 추출 (기준 완화해서 후보풀 확보)
    single = score_herbs(herb_map, disease_genes, min_overlap=1, top_n=top_n_single)
    if len(single) < combo_size:
        return []

    results = []
    for combo in combinations(single, combo_size):
        union_targets = set().union(*(h["targets"] for h in combo))
        overlap = union_targets & disease_set
        k = len(overlap)
        if k < min_overlap:
            continue

        H        = len(union_targets)
        coverage = k / D
        score    = k / math.sqrt(H * D) if H > 0 else 0.0
        pval     = _hypergeom_pvalue(k, H, D)
        synergy  = k - max(h["overlap_count"] for h in combo)

        results.append({
            "herbs":             [h["pinyin"] for h in combo],
            "cn_names":          [h.get("cn_name", "") for h in combo],
            "union_target_count": H,
            "overlap_genes":     sorted(overlap),
            "overlap_count":     k,
            "coverage":          round(coverage, 4),
            "score":             round(score, 4),
            "pvalue":            pval,
            "synergy":           synergy,
        })

    results.sort(key=lambda x: (x["pvalue"], -x["overlap_count"]))
    return results[:top_n_combo]


# ─────────────────────────────────────────────────────────────────────────────
# 4. Greedy Set Cover
# ─────────────────────────────────────────────────────────────────────────────
def greedy_set_cover(
    herb_map: dict,
    disease_genes: set,
    n_herbs: int = 5,
) -> dict:
    """
    Greedy Set Cover: 질환 타깃 커버리지를 최대화하는 약재를 순차 선택
    각 단계에서 미커버 타깃을 가장 많이 추가하는 약재 선택
    전통 복합처방(방제) 구성 원리와 유사
    """
    disease_set = {g.upper() for g in disease_genes}
    remaining   = disease_set.copy()
    covered     = set()
    selected    = []
    used        = set()

    for step in range(n_herbs):
        best_herb    = None
        best_new     = set()

        for pinyin, info in herb_map.items():
            if pinyin in used:
                continue
            new = info["targets"] & remaining
            if len(new) > len(best_new):
                best_new  = new
                best_herb = info

        if best_herb is None or not best_new:
            break

        covered  |= best_new
        remaining -= best_new
        used.add(best_herb["pinyin"])

        selected.append({
            **best_herb,
            "new_targets":         sorted(best_new),
            "new_count":           len(best_new),
            "cumulative_covered":  len(covered),
            "cumulative_pct":      round(len(covered) / len(disease_set) * 100, 1),
        })

    return {
        "selected_herbs":        selected,
        "total_covered":         len(covered),
        "total_disease_targets": len(disease_set),
        "coverage_pct":          round(len(covered) / len(disease_set) * 100, 1),
        "covered_genes":         sorted(covered),
        "uncovered_genes":       sorted(remaining),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 5. 통합 실행
# ─────────────────────────────────────────────────────────────────────────────
def run_discovery(
    disease_name: str,
    top_n: int = 20,
    ob_cutoff: float = 30.0,
    dl_cutoff: float = 0.18,
    min_overlap: int = 2,
    combo_sizes: list = None,       # [2, 3] 이면 2조합+3조합 실행, None이면 스킵
    top_n_single_for_combo: int = 30,
    greedy_n: int = 5,              # Greedy Set Cover 약재 수
    use_ppi: bool = False,          # STRING DB PPI 시너지 보정 여부
    ppi_min_score: int = 400,       # STRING 최소 상호작용 점수 (0~1000)
    ppi_top_n: int = 10,            # PPI 계산할 상위 조합 수
    ot_min_score: float = 0.1,      # OpenTargets 질환 타깃 최소 score
    on_progress: Optional[Callable[[str], None]] = None,
) -> dict:
    """
    약재 발굴 메인 함수
    반환: {disease_name, disease_genes, herb_results, combinations, greedy, ...}
    """
    def _prog(msg: str):
        log.info(msg)
        if on_progress:
            on_progress(msg)

    # 1. 캐시 로드
    _prog("TCMSP 캐시 로드 중...")
    herb_map = load_all_herb_targets(ob_cutoff=ob_cutoff, dl_cutoff=dl_cutoff)
    _prog(f"캐시 완료: {len(herb_map)}개 약재 로드됨")

    # 2. 질환 타깃 수집
    _prog(f"OpenTargets에서 '{disease_name}' 타깃 수집 중...")
    from collect import CacheDB
    cache = CacheDB(BASE_DIR / "db" / "cache.db")
    from disease import DiseaseTargetCollector
    collector = DiseaseTargetCollector(cache)
    disease_id, disease_targets = collector.get_targets_by_name(
        disease_name, min_score=ot_min_score
    )

    if not disease_targets:
        return {
            "error":        f"'{disease_name}'에 대한 질환 타깃을 찾을 수 없습니다.",
            "disease_name": disease_name,
        }

    disease_genes = {t["gene_symbol"] for t in disease_targets}
    _prog(f"질환 타깃 수집 완료: {len(disease_genes)}개 유전자")

    # 3. 단일 약재 스코어링
    _prog(f"{len(herb_map)}개 약재 단일 스코어링 중...")
    scored = score_herbs(herb_map, disease_genes, min_overlap=min_overlap, top_n=top_n)
    _prog(f"단일 스코어링 완료 — 상위 {len(scored)}개 선별")

    # 4. 조합 스코어링
    combo_results = {}
    if combo_sizes:
        for cs in combo_sizes:
            n_candidates = math.comb(min(top_n_single_for_combo, len(scored) + 20), cs)
            _prog(f"{cs}개 조합 스코어링 중... (후보 최대 {n_candidates:,}개)")
            combos = score_herb_combinations(
                herb_map, disease_genes,
                top_n_single=top_n_single_for_combo,
                combo_size=cs,
                top_n_combo=10,
                min_overlap=min_overlap,
            )
            # PPI 시너지 보정 (2-herb 조합 + use_ppi 옵션)
            if use_ppi and cs == 2 and combos:
                _prog(f"STRING PPI 시너지 계산 중 (상위 {ppi_top_n}개 조합)...")
                from ppi_synergy import enrich_combos_with_ppi
                herb_disease_targets = {
                    h: herb_map[h]["targets"] & disease_genes
                    for h in herb_map if herb_map[h]["targets"] & disease_genes
                }
                combos = enrich_combos_with_ppi(
                    combos,
                    herb_disease_targets,
                    disease_genes=disease_genes,
                    min_score=ppi_min_score,
                    top_n=ppi_top_n,
                    on_progress=_prog,
                )
                _prog(f"PPI 시너지 보정 완료")
            combo_results[cs] = combos
            _prog(f"{cs}개 조합 완료 — 상위 {len(combo_results[cs])}개")

    # 5. Greedy Set Cover
    _prog(f"Greedy Set Cover 실행 중 (최대 {greedy_n}개 약재)...")
    greedy = greedy_set_cover(herb_map, disease_genes, n_herbs=greedy_n)
    _prog(
        f"Greedy 완료 — {greedy_n}개 약재로 질환 타깃 "
        f"{greedy['coverage_pct']}% ({greedy['total_covered']}/{greedy['total_disease_targets']}개) 커버"
    )

    return {
        "disease_name":  disease_name,
        "disease_id":    disease_id,
        "disease_genes": sorted(disease_genes),
        "disease_count": len(disease_genes),
        "herb_results":  scored,
        "combinations":  combo_results,
        "greedy":        greedy,
        "cache_count":   len(herb_map),
        "ob_cutoff":     ob_cutoff,
        "dl_cutoff":     dl_cutoff,
        "min_overlap":   min_overlap,
    }
