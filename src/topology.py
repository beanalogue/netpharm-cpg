"""
허브 유전자 다중 토폴로지 분석
Degree / Betweenness / Closeness / ECC / MCC 복합 랭킹

MCC (Maximal Clique Centrality) 근사:
  networkx에 내장 함수 없으므로 clique 기반으로 직접 구현
  MCC(v) = sum of (k*(k-1)) for each maximal clique containing v, where k = clique size
"""

import logging
from itertools import combinations
from pathlib import Path

import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger(__name__)

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# MCC 계산
# ─────────────────────────────────────────────────────────────────────────────
def calc_mcc(G: nx.Graph) -> dict:
    """
    Maximal Clique Centrality 계산
    MCC(v) = Σ C(k,2) for each maximal clique of size k containing v
    = Σ k*(k-1)/2
    대형 그래프에서 느릴 수 있으므로 타겟 서브그래프에만 적용 권장
    """
    mcc = {n: 0 for n in G.nodes()}
    for clique in nx.find_cliques(G):
        k    = len(clique)
        score = k * (k - 1) // 2
        for node in clique:
            mcc[node] += score
    return mcc


# ─────────────────────────────────────────────────────────────────────────────
# 복합 토폴로지 분석기
# ─────────────────────────────────────────────────────────────────────────────
class TopologyAnalyzerExtended:
    """
    5가지 중심성 지표 + 복합 순위 점수

    지표:
      1. Degree Centrality
      2. Betweenness Centrality
      3. Closeness Centrality
      4. ECC (Eigenvector Centrality) — 연결 이웃의 중요성 반영
      5. MCC (Maximal Clique Centrality)

    복합 점수: Z-score 표준화 후 합산 → 정규화 (0~1)
    """

    def __init__(self, G: nx.Graph | nx.DiGraph):
        self.G_orig       = G
        self.G            = G.to_undirected() if G.is_directed() else G

    def compute(self, node_type_filter: str = "target") -> pd.DataFrame:
        """
        전체 지표 계산
        node_type_filter: 해당 node_type 노드만 서브그래프로 분석 (None=전체)
        """
        if node_type_filter:
            nodes = [n for n, d in self.G.nodes(data=True)
                     if d.get("node_type") == node_type_filter]
            G = self.G.subgraph(nodes).copy()
        else:
            G    = self.G
            nodes = list(G.nodes())

        if len(nodes) < 2:
            log.warning("[토폴로지] 노드 수 부족 — 계산 건너뜀")
            return pd.DataFrame()

        log.info(f"[토폴로지] {len(nodes)}개 노드 분석 중...")

        degree      = dict(G.degree())
        betweenness = nx.betweenness_centrality(G, normalized=True)
        closeness   = nx.closeness_centrality(G)
        clustering  = nx.clustering(G)

        # Eigenvector (연결 끊긴 그래프에서 오류 방지)
        try:
            eigenvector = nx.eigenvector_centrality(G, max_iter=500)
        except nx.PowerIterationFailedConvergence:
            eigenvector = {n: 0.0 for n in G.nodes()}

        # MCC (작은 서브그래프에서만)
        if len(nodes) <= 300:
            mcc = calc_mcc(G)
        else:
            log.warning("[토폴로지] MCC: 노드 >300, 근사값(degree²) 사용")
            mcc = {n: degree.get(n, 0) ** 2 for n in G.nodes()}

        rows = []
        for n in nodes:
            attr = self.G_orig.nodes.get(n, {})
            rows.append({
                "node":         n,
                "node_type":    attr.get("node_type", "unknown"),
                "label":        attr.get("label", n),
                "degree":       degree.get(n, 0),
                "betweenness":  round(betweenness.get(n, 0), 6),
                "closeness":    round(closeness.get(n,   0), 6),
                "eigenvector":  round(eigenvector.get(n, 0), 6),
                "clustering":   round(clustering.get(n,  0), 6),
                "mcc":          mcc.get(n, 0),
            })

        df = pd.DataFrame(rows)
        df = self._add_composite_score(df)
        df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
        df.insert(0, "rank", range(1, len(df) + 1))

        log.info(f"[토폴로지] 완료. 상위 허브: {df['label'].head(5).tolist()}")
        return df

    @staticmethod
    def _add_composite_score(df: pd.DataFrame) -> pd.DataFrame:
        """Z-score 표준화 후 5개 지표 합산 → 0~1 정규화"""
        metrics = ["degree", "betweenness", "closeness", "eigenvector", "mcc"]
        df = df.copy()
        z_sum = pd.Series(0.0, index=df.index)
        for m in metrics:
            col = df[m].astype(float)
            std = col.std()
            if std > 0:
                z_sum += (col - col.mean()) / std
            else:
                z_sum += 0
        mn, mx = z_sum.min(), z_sum.max()
        df["composite_score"] = ((z_sum - mn) / (mx - mn)).round(4) if mx > mn else 0.0
        return df

    def top_hubs(self, df: pd.DataFrame, n: int = 10) -> pd.DataFrame:
        return df.head(n)

    def plot_hub_ranking(
        self,
        df:       pd.DataFrame,
        top_n:    int  = 10,
        prefix:   str  = "analysis",
        filename: str  = None,
    ) -> Path:
        """허브 랭킹 레이더 차트 + 순위 막대 복합 Figure"""
        top = df.head(top_n).copy()
        metrics = ["degree", "betweenness", "closeness", "eigenvector", "mcc"]

        fig, axes = plt.subplots(1, 2, figsize=(14, max(5, top_n * 0.5)))

        # ── 왼쪽: 복합 점수 막대 ──────────────────────────────────────────
        ax = axes[0]
        colors = plt.cm.RdYlGn(np.linspace(0.3, 0.9, len(top)))[::-1]
        bars = ax.barh(top["label"][::-1], top["composite_score"][::-1],
                       color=colors, edgecolor="white", linewidth=0.5)
        for bar, score in zip(bars, top["composite_score"][::-1]):
            ax.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height() / 2,
                    f"{score:.3f}", va="center", ha="left", fontsize=8)
        ax.set_xlabel("Composite Score (normalized)", fontsize=10)
        ax.set_title("Hub Gene Ranking", fontsize=11, fontweight="bold")
        ax.set_xlim(0, 1.15)
        ax.tick_params(axis="y", labelsize=9)

        # ── 오른쪽: 지표별 히트맵 ────────────────────────────────────────
        ax2 = axes[1]
        heatmap_data = top[metrics].copy()
        # 각 지표 0~1 정규화
        for m in metrics:
            col = heatmap_data[m].astype(float)
            rng = col.max() - col.min()
            heatmap_data[m] = (col - col.min()) / rng if rng > 0 else 0

        im = ax2.imshow(heatmap_data.values, cmap="YlOrRd", aspect="auto",
                        vmin=0, vmax=1)
        ax2.set_xticks(range(len(metrics)))
        ax2.set_xticklabels(["Degree", "Betweenness", "Closeness", "Eigenvector", "MCC"],
                             rotation=30, ha="right", fontsize=9)
        ax2.set_yticks(range(len(top)))
        ax2.set_yticklabels(top["label"], fontsize=9)
        ax2.set_title("Centrality Metrics Heatmap", fontsize=11, fontweight="bold")
        plt.colorbar(im, ax=ax2, shrink=0.6, label="Normalized score")

        for i in range(len(top)):
            for j in range(len(metrics)):
                ax2.text(j, i, f"{heatmap_data.iloc[i, j]:.2f}",
                         ha="center", va="center", fontsize=7, color="black")

        plt.tight_layout()
        fname = filename or f"{prefix}_hub_ranking.png"
        path  = RESULTS_DIR / fname
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] 허브 랭킹 → {path}")
        return path
