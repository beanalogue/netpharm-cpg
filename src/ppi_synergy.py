"""
ppi_synergy.py
STRING DB 기반 약재 간 PPI 시너지 스코어 계산

약재 A의 타깃 집합과 약재 B의 타깃 집합 사이의
단백질-단백질 상호작용(PPI) 강도를 정량화하여
단순 합집합 오버랩을 넘어선 생물학적 시너지를 측정한다.
"""

import json
import time
import logging
import requests
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

STRING_API   = "https://string-db.org/api/json/network"
SPECIES      = 9606          # Homo sapiens
CALLER_ID    = "netpharm_research"
_CACHE_PATH  = Path("db/ppi_cache.json")
_MAX_GENES   = 100           # STRING API 실용적 한계
_RATE_SLEEP  = 0.35          # 초당 ~3 요청


# ─────────────────────────────────────────────────────────────────────────────
# 캐시 관리
# ─────────────────────────────────────────────────────────────────────────────
def _load_cache() -> dict:
    if _CACHE_PATH.exists():
        try:
            return json.loads(_CACHE_PATH.read_text())
        except Exception:
            pass
    return {}


def _save_cache(key: str, data: list) -> None:
    _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    cache = _load_cache()
    cache[key] = data
    _CACHE_PATH.write_text(json.dumps(cache, ensure_ascii=False))


# ─────────────────────────────────────────────────────────────────────────────
# STRING API 호출
# ─────────────────────────────────────────────────────────────────────────────
def get_string_interactions(
    genes: list[str],
    min_score: int = 400,
    use_cache: bool = True,
) -> list[dict]:
    """
    유전자 목록에 대한 STRING 상호작용 조회.
    반환: [{"gene_a", "gene_b", "score"}, ...]  score: 0~1
    """
    if not genes:
        return []

    genes = [g.upper() for g in genes if g]
    if len(genes) > _MAX_GENES:
        genes = genes[:_MAX_GENES]

    cache_key = "|".join(sorted(genes)) + f"@{min_score}"

    if use_cache:
        cache = _load_cache()
        if cache_key in cache:
            return cache[cache_key]

    try:
        resp = requests.post(
            STRING_API,
            data={
                "identifiers":    "\n".join(genes),
                "species":        SPECIES,
                "caller_identity": CALLER_ID,
                "required_score": min_score,
            },
            timeout=15,
        )
        raw = resp.json()
        result = [
            {
                "gene_a": item.get("preferredName_A", "").upper(),
                "gene_b": item.get("preferredName_B", "").upper(),
                "score":  round(float(item.get("score", 0)), 4),
            }
            for item in raw
            if isinstance(item, dict)
        ]
        if use_cache:
            _save_cache(cache_key, result)
        time.sleep(_RATE_SLEEP)
        return result

    except Exception as e:
        log.warning(f"[STRING] API 오류: {e}")
        return []


# ─────────────────────────────────────────────────────────────────────────────
# PPI 시너지 계산
# ─────────────────────────────────────────────────────────────────────────────
def compute_ppi_synergy(
    genes_a: set[str],
    genes_b: set[str],
    disease_genes: Optional[set[str]] = None,
    min_score: int = 400,
) -> dict:
    """
    약재 A 타깃 vs 약재 B 타깃 간 cross-herb PPI 시너지 계산.

    반환:
      ppi_cross_n     : 두 약재 타깃 사이 상호작용 수
      ppi_cross_score : cross-herb 상호작용 점수 합 (0~N)
      ppi_mean_score  : cross-herb 평균 상호작용 강도 (0~1)
      ppi_disease_n   : 질환 타깃 관련 cross-herb 상호작용 수
      top_pairs       : 상위 5 쌍 [(gene_a, gene_b, score), ...]
    """
    ga = {g.upper() for g in genes_a if g}
    gb = {g.upper() for g in genes_b if g}
    dg = {g.upper() for g in disease_genes} if disease_genes else set()

    all_genes = list(ga | gb)
    if len(all_genes) < 2:
        return _empty_ppi()

    interactions = get_string_interactions(all_genes, min_score=min_score)

    cross, disease_cross = [], []
    for itr in interactions:
        a, b, sc = itr["gene_a"], itr["gene_b"], itr["score"]
        is_cross = (a in ga and b in gb) or (a in gb and b in ga)
        if is_cross:
            cross.append((a, b, sc))
            if dg and (a in dg or b in dg):
                disease_cross.append((a, b, sc))

    cross.sort(key=lambda x: -x[2])
    total_score = sum(x[2] for x in cross)
    mean_score  = total_score / len(cross) if cross else 0.0

    return {
        "ppi_cross_n":     len(cross),
        "ppi_cross_score": round(total_score, 3),
        "ppi_mean_score":  round(mean_score, 3),
        "ppi_disease_n":   len(disease_cross),
        "top_pairs":       [(a, b, round(sc, 3)) for a, b, sc in cross[:5]],
    }


def _empty_ppi() -> dict:
    return {
        "ppi_cross_n":     0,
        "ppi_cross_score": 0.0,
        "ppi_mean_score":  0.0,
        "ppi_disease_n":   0,
        "top_pairs":       [],
    }


# ─────────────────────────────────────────────────────────────────────────────
# 조합 목록에 PPI 시너지 일괄 추가
# ─────────────────────────────────────────────────────────────────────────────
def enrich_combos_with_ppi(
    combos: list[dict],
    herb_targets: dict[str, set[str]],
    disease_genes: Optional[set[str]] = None,
    min_score: int = 400,
    top_n: int = 10,
    on_progress=None,
) -> list[dict]:
    """
    조합 결과 리스트에 PPI 시너지 필드를 추가하여 반환.
    top_n개 조합에만 적용 (API 호출 절약).

    입력 combo dict에 추가되는 필드:
      ppi_cross_n, ppi_cross_score, ppi_mean_score,
      ppi_disease_n, top_pairs, combined_score
    """
    enriched = []
    for i, combo in enumerate(combos[:top_n]):
        herbs = combo.get("herbs", [])
        if len(herbs) != 2:
            # 현재 PPI는 2-herb 조합에 최적화
            combo["ppi_cross_n"]     = 0
            combo["ppi_cross_score"] = 0.0
            combo["ppi_mean_score"]  = 0.0
            combo["ppi_disease_n"]   = 0
            combo["top_pairs"]       = []
            combo["combined_score"]  = combo.get("score", 0)
            enriched.append(combo)
            continue

        h1, h2 = herbs[0], herbs[1]
        t1 = herb_targets.get(h1, set())
        t2 = herb_targets.get(h2, set())

        if on_progress:
            on_progress(f"PPI 조회 중: {h1} × {h2} ({i+1}/{min(top_n, len(combos))})")

        ppi = compute_ppi_synergy(t1, t2, disease_genes=disease_genes, min_score=min_score)

        # combined_score = 타깃 오버랩 점수 + PPI 시너지 보정
        base_score  = combo.get("score", 0)
        ppi_bonus   = ppi["ppi_cross_score"] * 0.1   # 가중치 조정 가능
        combined    = round(base_score + ppi_bonus, 4)

        combo.update(ppi)
        combo["combined_score"] = combined
        enriched.append(combo)

    # PPI 정보 없는 나머지 추가
    for combo in combos[top_n:]:
        combo.update(_empty_ppi())
        combo["combined_score"] = combo.get("score", 0)
        enriched.append(combo)

    # combined_score 기준 재정렬
    enriched.sort(key=lambda x: (-x["combined_score"], -x.get("ppi_cross_n", 0)))
    return enriched
