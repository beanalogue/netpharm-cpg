"""
전통 처방 기반 약재 발굴 결과 검증 모듈
curated TCM formula dataset (Chinese Pharmacopoeia 2020 + 방제학 교과서) 사용

분석 단위 (2025 개정):
  ─ 1차 outcome (primary): 질환 pool AUROC
      전체 CPG 처방 약재 합집합 vs 전체 랭킹 리스트
      Wilcoxon rank-sum + permutation test (n=1000)
      다중 검정: 질환 수(3개) 기준 Bonferroni α=0.0167
  ─ 2차 outcome (exploratory): 처방별 AUROC
      유의성 별표 없이 기술 통계로만 보고
      (n이 작고 다중 검정 후 유의 처방 없음 → 탐색적 목적)
"""

import json
import logging
import random
from pathlib import Path

log = logging.getLogger(__name__)

BASE_DIR     = Path(__file__).parent.parent
FORMULA_PATH = BASE_DIR / "data" / "tcm_validation_formulas.json"

# ── 다중 검정 기준 ────────────────────────────────────────────────
N_DISEASES           = 3       # 검증 대상 질환 수
BONFERRONI_ALPHA     = 0.05 / N_DISEASES   # = 0.0167
PERMUTATION_N        = 1_000  # permutation test 반복 수


def load_formula_data() -> dict:
    if not FORMULA_PATH.exists():
        raise FileNotFoundError(f"처방 데이터 없음: {FORMULA_PATH}")
    return json.loads(FORMULA_PATH.read_text(encoding="utf-8"))


def get_formula_herb_pool(disease_name: str) -> dict:
    """
    질환명 → 전통 처방 약재 풀
    반환: {formula_herbs(set), formulas(list), disease_cn, matched_disease, excluded_herbs}
    """
    data = load_formula_data()

    matched_key = None
    q = disease_name.lower().strip()
    for key in data:
        if key == "_meta":
            continue
        if q in key.lower() or key.lower() in q:
            matched_key = key
            break
        en = data[key].get("disease_en", "").lower()
        if q in en or en in q:
            matched_key = key
            break

    if not matched_key:
        log.warning(f"[검증] '{disease_name}'에 해당하는 처방 데이터 없음")
        return {"formula_herbs": set(), "formulas": [], "matched_disease": None,
                "excluded_herbs": {}}

    disease_data = data[matched_key]
    formulas     = disease_data.get("formulas", {})

    all_herbs = set()
    excluded_herbs = {}   # {formula_name: [excluded]}
    formula_list = []
    for fname, fdata in formulas.items():
        if fname.startswith("_") or not isinstance(fdata, dict):
            continue
        herbs_used = set(fdata.get("herbs", []))
        herbs_all  = set(fdata.get("herbs_all", fdata.get("herbs", [])))
        excluded   = sorted(herbs_all - herbs_used)
        all_herbs |= herbs_used
        if excluded:
            excluded_herbs[fdata.get("name_cn", fname)] = excluded
        formula_list.append({
            "formula_key":  fname,
            "name_cn":      fdata.get("name_cn", ""),
            "name_en":      fdata.get("name_en", ""),
            "source":       fdata.get("source", ""),
            "herbs":        sorted(herbs_used),
            "herb_count":   len(herbs_used),
            "excluded":     excluded,
        })

    log.info(
        f"[검증] '{disease_name}' → {disease_data.get('disease_cn')} "
        f"처방 {len(formula_list)}개 / 총 약재 {len(all_herbs)}개"
    )
    return {
        "matched_disease":    matched_key,
        "disease_cn":         disease_data.get("disease_cn", ""),
        "disease_en":         disease_data.get("disease_en", ""),
        "formula_herbs":      all_herbs,
        "formulas":           formula_list,
        "total_unique_herbs": len(all_herbs),
        "excluded_herbs":     excluded_herbs,
    }


# ── AUROC 계산 헬퍼 ────────────────────────────────────────────────
def _compute_auroc(formula_ranks: list, non_formula_ranks: list) -> float:
    """Mann-Whitney U 기반 AUROC"""
    pairs = sum(1 for fr in formula_ranks for nr in non_formula_ranks if fr < nr)
    total = len(formula_ranks) * len(non_formula_ranks)
    return pairs / total if total > 0 else 0.5


def _wilcoxon_pvalue(non_formula_ranks: list, formula_ranks: list) -> float:
    try:
        from scipy.stats import ranksums
        _, p = ranksums(non_formula_ranks, formula_ranks, alternative="greater")
        return float(p)
    except Exception:
        return 1.0


def _permutation_pvalue(
    all_herbs_ranked: list,
    formula_herbs: set,
    observed_auroc: float,
    n_perm: int = PERMUTATION_N,
    seed: int = 42,
) -> float:
    """
    Permutation test:
    랜덤하게 동일 크기 약재 집합을 n_perm회 추출 →
    관찰 AUROC ≥ permuted AUROC 비율로 p-value 추정
    """
    rng = random.Random(seed)
    n_formula = len([h for h in formula_herbs if h in all_herbs_ranked])
    if n_formula == 0:
        return 1.0

    total = len(all_herbs_ranked)
    count_ge = 0
    for _ in range(n_perm):
        rand_idx = set(rng.sample(range(total), n_formula))
        rand_ranks      = [i + 1 for i in rand_idx]
        non_rand_ranks  = [i + 1 for i in range(total) if i not in rand_idx]
        rand_auroc = _compute_auroc(rand_ranks, non_rand_ranks)
        if rand_auroc >= observed_auroc:
            count_ge += 1
    return count_ge / n_perm


# ── 1차 outcome: 질환 pool AUROC ──────────────────────────────────
def validate_pool(
    discovery_result: dict,
    n_permutations: int = PERMUTATION_N,
) -> dict:
    """
    PRIMARY analysis: 질환 전체 CPG 처방 약재 풀 기반 AUROC
    다중 검정 기준: Bonferroni α = 0.05 / N_DISEASES = 0.0167

    반환 dict 주요 키:
      pool_auroc, wilcoxon_p, permutation_p, significant (bool),
      bonferroni_threshold, n_formula_herbs_ranked, n_total_herbs
    """
    disease_name = discovery_result.get("disease_name", "")
    formula_info = get_formula_herb_pool(disease_name)
    formula_herbs = formula_info["formula_herbs"]

    if not formula_herbs:
        return {"error": f"'{disease_name}' 처방 데이터 없음"}

    all_herbs_ranked = [h["pinyin"] for h in discovery_result.get("herb_results", [])]
    total = len(all_herbs_ranked)

    formula_ranks = [
        all_herbs_ranked.index(h) + 1
        for h in formula_herbs if h in all_herbs_ranked
    ]
    non_formula_ranks = [
        i + 1 for i, h in enumerate(all_herbs_ranked) if h not in formula_herbs
    ]

    if not formula_ranks:
        return {"error": "처방 약재 중 순위 계산 가능한 약재 없음"}

    auroc   = round(_compute_auroc(formula_ranks, non_formula_ranks), 4)
    wilc_p  = round(_wilcoxon_pvalue(non_formula_ranks, formula_ranks), 6)
    perm_p  = round(_permutation_pvalue(all_herbs_ranked, formula_herbs, auroc,
                                        n_perm=n_permutations), 4)
    mean_rank = round(sum(formula_ranks) / len(formula_ranks), 1)

    # 사용 p-value: Wilcoxon (이론적) + permutation (경험적) 둘 다 보고
    # 유의성 판단 기준: permutation p < Bonferroni threshold
    significant = perm_p < BONFERRONI_ALPHA

    log.info(
        f"[1차검증] {disease_name} pool AUROC={auroc:.4f} "
        f"Wilcoxon p={wilc_p:.4f} Perm p={perm_p:.4f} "
        f"({'유의' if significant else '비유의'}, α={BONFERRONI_ALPHA:.4f})"
    )

    return {
        "disease_name":           disease_name,
        "disease_cn":             formula_info["disease_cn"],
        "disease_en":             formula_info["disease_en"],
        "analysis_type":          "primary",
        "n_total_herbs":          total,
        "n_formula_herbs_pool":   len(formula_herbs),
        "n_formula_herbs_ranked": len(formula_ranks),
        "n_excluded_from_tcmsp":  len(formula_herbs) - len(formula_ranks),
        "formula_ranks":          sorted(formula_ranks),
        "mean_rank":              mean_rank,
        "expected_random_rank":   round((total + 1) / 2, 1),
        "pool_auroc":             auroc,
        "wilcoxon_p":             wilc_p,
        "permutation_p":          perm_p,
        "n_permutations":         n_permutations,
        "bonferroni_threshold":   round(BONFERRONI_ALPHA, 4),
        "significant":            significant,
        "excluded_herbs":         formula_info["excluded_herbs"],
        "n_formulas":             len(formula_info["formulas"]),
        "formula_herb_pool":      sorted(formula_herbs),
    }


# ── 2차 outcome: 처방별 AUROC (exploratory) ───────────────────────
def validate_per_formula(
    discovery_result: dict,
) -> list:
    """
    EXPLORATORY analysis: 처방별 개별 AUROC
    유의성 별표 없이 순수 기술 통계로 보고
    (다중 검정 보정 후 유의 처방 없음 — 가설 생성 목적으로만 사용)
    """
    disease_name = discovery_result.get("disease_name", "")
    formula_info = get_formula_herb_pool(disease_name)
    all_herbs_ranked = [h["pinyin"] for h in discovery_result.get("herb_results", [])]
    total = len(all_herbs_ranked)

    results = []
    for f in formula_info["formulas"]:
        herbs = set(f["herbs"])
        f_ranks = [all_herbs_ranked.index(h) + 1 for h in herbs if h in all_herbs_ranked]
        non_ranks = [i + 1 for i, h in enumerate(all_herbs_ranked) if h not in herbs]

        if len(f_ranks) < 2:
            auroc, wilc_p = None, None
        else:
            auroc  = round(_compute_auroc(f_ranks, non_ranks), 4)
            wilc_p = round(_wilcoxon_pvalue(non_ranks, f_ranks), 4)

        results.append({
            "formula_key":        f["formula_key"],
            "name_cn":            f["name_cn"],
            "name_en":            f["name_en"],
            "n_herbs_total":      f["herb_count"],
            "n_herbs_ranked":     len(f_ranks),
            "n_herbs_excluded":   f["herb_count"] - len(f_ranks),
            "excluded_herbs":     f["excluded"],
            "auroc":              auroc,
            "wilcoxon_p":         wilc_p,
            "analysis_type":      "exploratory",
            "note":               "No significance stars — descriptive only (post-hoc MTC fails)",
        })

    results.sort(key=lambda x: (x["auroc"] or 0), reverse=True)
    return results


# ── 통합 래퍼 ─────────────────────────────────────────────────────
def validate_discovery(
    discovery_result: dict,
    top_n: int = 10,              # legacy 호환 (pool 분석에선 미사용)
    n_permutations: int = PERMUTATION_N,
) -> dict:
    """
    통합 검증 함수 (이전 버전 호환 유지)
    반환:
      primary   → pool AUROC 결과
      exploratory → 처방별 AUROC 목록 (기술 통계)
    """
    primary     = validate_pool(discovery_result, n_permutations=n_permutations)
    exploratory = validate_per_formula(discovery_result)

    # legacy 키 유지 (app.py 호환)
    primary["per_formula"]    = exploratory
    primary["auroc"]          = primary.get("pool_auroc")
    primary["auroc_pvalue"]   = primary.get("permutation_p")
    primary["total_herbs_analyzed"] = primary.get("n_total_herbs")

    return primary


# ── 리포트 포매터 ─────────────────────────────────────────────────
def format_validation_report(val: dict) -> str:
    sig_str = (
        f"**유의** (α={val['bonferroni_threshold']}, Bonferroni 보정)"
        if val.get("significant")
        else f"비유의 (α={val['bonferroni_threshold']}, Bonferroni 보정)"
    )
    lines = [
        f"### [1차 검증] 질환 Pool AUROC — {val.get('disease_cn','')} ({val.get('disease_name','')})",
        f"",
        f"| 지표 | 값 |",
        f"|------|----|",
        f"| Pool AUROC | **{val['pool_auroc']:.4f}** |",
        f"| Wilcoxon p | {val['wilcoxon_p']:.4f} |",
        f"| Permutation p (n={val['n_permutations']}) | {val['permutation_p']:.4f} |",
        f"| 유의성 (Bonferroni α={val['bonferroni_threshold']}) | {sig_str} |",
        f"| 처방 약재 수 (pool) | {val['n_formula_herbs_pool']}개 |",
        f"| 순위 계산 가능 | {val['n_formula_herbs_ranked']}개 |",
        f"| TCMSP 미수록 제외 | {val['n_excluded_from_tcmsp']}개 |",
        f"| 평균 순위 / 기대 순위 | {val['mean_rank']} / {val['expected_random_rank']} |",
        f"",
    ]

    # 제외 약재 목록
    excl = val.get("excluded_herbs", {})
    if excl:
        lines.append("**TCMSP 미수록 제외 약재** (처방별):")
        for fname, herbs in excl.items():
            lines.append(f"- {fname}: {', '.join(herbs)}")
        lines.append("")

    # 처방별 exploratory
    lines += [
        f"### [2차 검증] 처방별 AUROC (탐색적 — 유의성 별표 없음)",
        f"",
        f"| 처방 | n약재 (사용/전체) | AUROC | Wilcoxon p |",
        f"|------|-----------------|-------|-----------|",
    ]
    for f in val.get("per_formula", []):
        auroc_str = f"{f['auroc']:.4f}" if f["auroc"] is not None else "N/A"
        p_str     = f"{f['wilcoxon_p']:.4f}" if f["wilcoxon_p"] is not None else "N/A"
        lines.append(
            f"| {f['name_cn']} | {f['n_herbs_ranked']}/{f['n_herbs_total']} "
            f"| {auroc_str} | {p_str} |"
        )
    lines += [
        f"",
        f"> ⚠️ 처방별 분석은 탐색적 목적으로만 사용. 다중 검정 보정(Bonferroni/BH-FDR) 후 "
        f"유의한 처방 없음. 향후 가설 생성 용도로만 해석 권장.",
        f"",
        f"*출처: Chinese Pharmacopoeia 2020; 方剂学 2nd ed., 2011*",
    ]
    return "\n".join(lines)
