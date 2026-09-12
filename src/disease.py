"""
질환-타깃 교차 분석
OpenTargets API (무료, 인증 불필요)로 질환 타깃 수집
+ Venn Diagram / UpSet Plot 생성
"""

import time
import json
import logging
from pathlib import Path

import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BASE_DIR     = Path(__file__).parent.parent
RESULTS_DIR  = BASE_DIR / "results"

log = logging.getLogger(__name__)

OPENTARGETS_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

try:
    from matplotlib_venn import venn2, venn3
    VENN_OK = True
except ImportError:
    VENN_OK = False
    log.warning("[Disease] matplotlib-venn 미설치")

try:
    from upsetplot import UpSet, from_memberships
    UPSET_OK = True
except ImportError:
    UPSET_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# 1. OpenTargets 수집기
# ─────────────────────────────────────────────────────────────────────────────
class DiseaseTargetCollector:

    def __init__(self, cache, delay: float = 1.0):
        self.cache   = cache
        self.delay   = delay
        from collect import make_retry_session
        self.session = make_retry_session()

    def search_disease(self, query: str, limit: int = 10) -> list:
        """질환명 검색 → [{id, name}, ...]"""
        key = f"OT:search:{query}"
        cached = self.cache.get(key)
        if cached is not None:
            return cached

        gql = """
        query($q: String!) {
          search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: %d}) {
            hits { id name }
          }
        }
        """ % limit

        time.sleep(self.delay)
        try:
            r = self.session.post(
                OPENTARGETS_GQL,
                json={"query": gql, "variables": {"q": query}},
                timeout=20,
            )
            r.raise_for_status()
            hits = r.json()["data"]["search"]["hits"]
            self.cache.set(key, hits)
            return hits
        except Exception as e:
            log.error(f"[OpenTargets] 검색 오류: {e}")
            return []

    def get_disease_targets(
        self,
        disease_id:  str,
        min_score:   float = 0.1,
        max_targets: int   = 500,
    ) -> list:
        """
        질환 EFO/MONDO ID → 연관 타겟 유전자 목록
        반환: [{"gene_symbol", "score", "disease_id"}, ...]
        """
        key = f"OT:targets:{disease_id}:{min_score}"
        cached = self.cache.get(key)
        if cached is not None:
            log.info(f"[OpenTargets] 캐시: {disease_id}")
            return cached

        gql = """
        query($id: String!, $size: Int!) {
          disease(efoId: $id) {
            name
            associatedTargets(page: {index: 0, size: $size}) {
              rows {
                target { approvedSymbol }
                score
              }
            }
          }
        }
        """
        time.sleep(self.delay)
        try:
            r = self.session.post(
                OPENTARGETS_GQL,
                json={"query": gql, "variables": {"id": disease_id, "size": max_targets}},
                timeout=30,
            )
            r.raise_for_status()
            data = r.json()["data"]["disease"]
            if not data:
                log.warning(f"[OpenTargets] {disease_id}: 결과 없음")
                return []

            rows = data["associatedTargets"]["rows"]
            targets = [
                {
                    "gene_symbol": row["target"]["approvedSymbol"].upper(),
                    "score":       round(row["score"], 4),
                    "disease_id":  disease_id,
                }
                for row in rows
                if row["score"] >= min_score
            ]
            self.cache.set(key, targets)
            log.info(f"[OpenTargets] {data['name']}: 타겟 {len(targets)}개 (score≥{min_score})")
            return targets

        except Exception as e:
            log.error(f"[OpenTargets] {disease_id} 타겟 조회 오류: {e}")
            return []

    def get_targets_by_name(self, disease_name: str, min_score: float = 0.1) -> tuple:
        """질환명 → (disease_id, gene_list) 자동 처리"""
        hits = self.search_disease(disease_name, limit=5)
        if not hits:
            return "", []
        best = hits[0]
        targets = self.get_disease_targets(best["id"], min_score=min_score)
        log.info(f"[OpenTargets] '{disease_name}' → '{best['name']}' ({best['id']})")
        return best["id"], targets


# ─────────────────────────────────────────────────────────────────────────────
# 2. 교차 분석
# ─────────────────────────────────────────────────────────────────────────────
class IntersectionAnalyzer:

    @staticmethod
    def intersect(
        herb_genes:    set,
        disease_genes: set,
        herb_label:    str = "Herb",
        disease_label: str = "Disease",
    ) -> dict:
        """두 유전자 집합의 교집합 분석"""
        herb_set    = {g.upper() for g in herb_genes}
        disease_set = {g.upper() for g in disease_genes}
        shared      = herb_set & disease_set

        result = {
            "herb_only":    sorted(herb_set    - shared),
            "disease_only": sorted(disease_set - shared),
            "shared":       sorted(shared),
            "herb_total":   len(herb_set),
            "disease_total":len(disease_set),
            "shared_count": len(shared),
            "herb_label":   herb_label,
            "disease_label":disease_label,
        }
        log.info(
            f"[교차분석] {herb_label}({len(herb_set)}) ∩ {disease_label}({len(disease_set)}) "
            f"= {len(shared)}개"
        )
        return result

    @staticmethod
    def multi_intersect(gene_sets: dict) -> dict:
        """
        3개 이상 집합 교차 분석
        gene_sets: {label: set_of_genes}
        반환: 모든 집합의 교집합 포함 dict
        """
        labels = list(gene_sets.keys())
        sets   = [set(g.upper() for g in gene_sets[l]) for l in labels]

        result = {}
        # 전체 공통
        all_shared = sets[0].copy()
        for s in sets[1:]:
            all_shared &= s
        result["all_shared"] = sorted(all_shared)

        # 쌍별 교집합
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                key = f"{labels[i]}_x_{labels[j]}"
                result[key] = sorted(sets[i] & sets[j])

        return result


# ─────────────────────────────────────────────────────────────────────────────
# 3. 시각화
# ─────────────────────────────────────────────────────────────────────────────
class IntersectionPlotter:

    def __init__(self, output_dir: Path = RESULTS_DIR, dpi: int = 200):
        self.out = Path(output_dir)
        self.out.mkdir(exist_ok=True)
        self.dpi = dpi

    def venn_diagram(
        self,
        intersection: dict,
        prefix: str = "analysis",
        filename: str = None,
    ) -> Path:
        """Venn Diagram (2집합)"""
        if not VENN_OK:
            log.warning("[Venn] matplotlib-venn 미설치")
            return None

        herb_only    = len(intersection["herb_only"])
        disease_only = len(intersection["disease_only"])
        shared       = intersection["shared_count"]

        fig, ax = plt.subplots(figsize=(7, 5))
        v = venn2(
            subsets=(herb_only, disease_only, shared),
            set_labels=(intersection["herb_label"], intersection["disease_label"]),
            ax=ax,
        )
        # 색상
        if v.get_patch_by_id("10"):
            v.get_patch_by_id("10").set_color("#4C72B0")
            v.get_patch_by_id("10").set_alpha(0.6)
        if v.get_patch_by_id("01"):
            v.get_patch_by_id("01").set_color("#DD8452")
            v.get_patch_by_id("01").set_alpha(0.6)
        if v.get_patch_by_id("11"):
            v.get_patch_by_id("11").set_color("#55A868")
            v.get_patch_by_id("11").set_alpha(0.8)

        ax.set_title(
            f"Target Gene Intersection\n"
            f"Shared: {shared} genes",
            fontsize=12, fontweight="bold",
        )
        plt.tight_layout()

        fname = filename or f"{prefix}_venn.png"
        path  = self.out / fname
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] Venn Diagram → {path}")
        return path

    def venn3_diagram(
        self,
        sets:   dict,
        prefix: str = "analysis",
    ) -> Path:
        """Venn Diagram (3집합)"""
        if not VENN_OK or len(sets) < 3:
            return None

        labels = list(sets.keys())[:3]
        gene_sets = [set(g.upper() for g in sets[l]) for l in labels]

        fig, ax = plt.subplots(figsize=(8, 6))
        venn3(gene_sets, set_labels=labels, ax=ax)
        ax.set_title("Target Gene Intersection (3 sets)", fontsize=12, fontweight="bold")
        plt.tight_layout()

        path = self.out / f"{prefix}_venn3.png"
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] Venn3 → {path}")
        return path

    def upset_plot(
        self,
        gene_sets: dict,
        prefix:    str = "analysis",
    ) -> Path:
        """UpSet Plot (4개 이상 집합)"""
        if not UPSET_OK:
            log.warning("[UpSet] upsetplot 미설치")
            return None

        all_genes = sorted(set().union(*[set(g.upper() for g in v) for v in gene_sets.values()]))
        memberships = []
        for gene in all_genes:
            mem = tuple(label for label, genes in gene_sets.items()
                        if gene.upper() in {g.upper() for g in genes})
            memberships.append(mem)

        data = from_memberships(memberships)
        fig  = plt.figure(figsize=(12, 6))
        UpSet(data, show_counts=True).plot(fig)
        plt.suptitle("Target Gene Overlap (UpSet Plot)", fontsize=12, fontweight="bold")

        path = self.out / f"{prefix}_upset.png"
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] UpSet Plot → {path}")
        return path

    def shared_gene_barplot(
        self,
        intersection: dict,
        disease_targets: list,
        prefix: str = "analysis",
    ) -> Path:
        """공유 유전자를 질환 연관성 점수로 정렬한 bar chart"""
        shared = intersection["shared"]
        if not shared:
            return None

        score_map = {t["gene_symbol"]: t.get("score", 0) for t in disease_targets}
        df = pd.DataFrame([
            {"gene": g, "score": score_map.get(g, 0)}
            for g in shared
        ]).sort_values("score", ascending=True)

        fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.35)))
        ax.barh(df["gene"], df["score"], color="#55A868", alpha=0.8, edgecolor="white")
        ax.set_xlabel("OpenTargets Association Score", fontsize=10)
        ax.set_title(
            f"Shared Target Genes — Disease Association Score\n"
            f"({intersection['herb_label']} ∩ {intersection['disease_label']})",
            fontsize=11, fontweight="bold",
        )
        ax.tick_params(axis="y", labelsize=8)
        plt.tight_layout()

        path = self.out / f"{prefix}_shared_genes.png"
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] 공유 유전자 bar chart → {path}")
        return path
