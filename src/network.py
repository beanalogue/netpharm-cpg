"""
Network Pharmacology — Network Builder
herb-compound-target 네트워크 + PPI 통합 + 위상 분석 + Cytoscape 내보내기

주요 기능:
  1. HCT 네트워크   : Herb → Compound → Target 다층 네트워크
  2. PPI 네트워크   : STRING 기반 단백질 상호작용
  3. 통합 네트워크  : HCT + PPI 병합 (한약-양약 공유 타겟 분석)
  4. 위상 지표      : degree, betweenness, closeness, hub 판별
  5. 내보내기       : GraphML (Cytoscape), edge/node CSV (논문용)
"""

import json
import logging
from pathlib import Path
from collections import defaultdict

import networkx as nx
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR    = Path(__file__).parent.parent
NETWORK_DIR = BASE_DIR / "networks"
NETWORK_DIR.mkdir(exist_ok=True)

log = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)


# ─────────────────────────────────────────────────────────────────────────────
# 1. HCT 네트워크 구축
# ─────────────────────────────────────────────────────────────────────────────
class HCTNetwork:
    """
    Herb-Compound-Target (HCT) 다층 네트워크

    노드 유형:
      herb     : 한약 (예: Panax ginseng)
      compound : 활성 성분 (예: Ginsenoside Rb1)
      target   : 타겟 단백질/유전자 (예: TP53)

    엣지:
      herb → compound  : 'herb_compound'
      compound → target: 'compound_target'
    """

    def __init__(self):
        self.G = nx.DiGraph()  # 방향 그래프 (herb→compound→target)

    def build_from_herb_results(self, herb_results: dict, source_label: str = "herb"):
        """
        NetPharmCollector.collect_herb_compounds_chembl() 또는
        collect_from_tcmsp() 결과로 HCT 네트워크 구축

        herb_results:
          {herb_name: {"compounds": [{"query", "molecule", "targets"}, ...], "source": ...}}
        """
        for herb_name, herb_data in herb_results.items():
            src = herb_data.get("source", source_label)

            # 허브 노드
            self.G.add_node(herb_name, node_type="herb", label=herb_name, source=src)

            for comp_data in herb_data.get("compounds", []):
                mol   = comp_data.get("molecule", {})
                cid   = mol.get("chembl_id") or comp_data.get("query", "")
                cname = mol.get("pref_name")  or comp_data.get("query", cid)
                if not cid:
                    continue

                # 성분 노드
                self.G.add_node(
                    cid,
                    node_type="compound",
                    label=cname,
                    mol_weight=mol.get("mol_weight"),
                    ob=mol.get("ob"),
                    dl=mol.get("dl"),
                    smiles=mol.get("smiles", ""),
                    source=src,
                )
                self.G.add_edge(herb_name, cid, edge_type="herb_compound")

                for t in comp_data.get("targets", []):
                    gs = t.get("gene_symbol", "").strip().upper()
                    if not gs:
                        continue
                    if gs not in self.G:
                        self.G.add_node(
                            gs,
                            node_type="target",
                            label=gs,
                            uniprot=t.get("uniprot_id", ""),
                            target_name=t.get("target_name", ""),
                            source=src,
                        )
                    self.G.add_edge(cid, gs, edge_type="compound_target")

        log.info(
            f"[HCT] 네트워크 구축 완료: "
            f"노드 {self.G.number_of_nodes()}개 / 엣지 {self.G.number_of_edges()}개"
        )
        return self.G

    def build_from_tcmsp(self, compounds_df: pd.DataFrame, targets_df: pd.DataFrame,
                          herb_col="Herb_name", comp_col="mol_name",
                          gene_col="Target_name"):
        """TCMSP 로컬 데이터프레임으로 HCT 구축"""
        merged = pd.merge(compounds_df, targets_df, on=comp_col, how="inner")
        for _, row in merged.iterrows():
            herb  = str(row[herb_col]).strip()
            comp  = str(row[comp_col]).strip()
            gene  = str(row[gene_col]).strip().upper()

            self.G.add_node(herb, node_type="herb",     label=herb)
            self.G.add_node(comp, node_type="compound", label=comp,
                            ob=row.get("OB (%)"), dl=row.get("DL"))
            self.G.add_node(gene, node_type="target",   label=gene)
            self.G.add_edge(herb, comp, edge_type="herb_compound")
            self.G.add_edge(comp, gene, edge_type="compound_target")

        log.info(f"[HCT/TCMSP] 노드 {self.G.number_of_nodes()} / 엣지 {self.G.number_of_edges()}")
        return self.G

    def get_targets(self) -> list:
        """타겟(유전자) 노드 목록 반환"""
        return [n for n, d in self.G.nodes(data=True) if d.get("node_type") == "target"]

    def get_herb_target_matrix(self) -> pd.DataFrame:
        """herb × target 공유 행렬 (논문 Table 작성용)"""
        herbs   = [n for n, d in self.G.nodes(data=True) if d.get("node_type") == "herb"]
        targets = self.get_targets()
        if not herbs or not targets:
            return pd.DataFrame()
        mat = pd.DataFrame(0, index=herbs, columns=targets)
        for herb in herbs:
            for comp in self.G.successors(herb):
                for target in self.G.successors(comp):
                    if target in mat.columns:
                        mat.loc[herb, target] = 1
        return mat


# ─────────────────────────────────────────────────────────────────────────────
# 2. PPI 네트워크 구축
# ─────────────────────────────────────────────────────────────────────────────
class PPINetwork:
    """STRING PPI 네트워크"""

    def __init__(self):
        self.G = nx.Graph()

    def build_from_string(self, ppi_edges: list, min_score: int = 400):
        """
        STRINGCollector.get_ppi() 결과로 PPI 네트워크 구축

        ppi_edges: [{stringId_A, stringId_B, preferredName_A, preferredName_B, score}, ...]
        """
        for edge in ppi_edges:
            a     = edge.get("preferredName_A", edge.get("stringId_A", ""))
            b     = edge.get("preferredName_B", edge.get("stringId_B", ""))
            score = float(edge.get("score", 0))

            if score < min_score / 1000:  # STRING score는 0~1 스케일
                continue

            self.G.add_node(a.upper(), node_type="protein", label=a.upper())
            self.G.add_node(b.upper(), node_type="protein", label=b.upper())
            self.G.add_edge(a.upper(), b.upper(), weight=score, source="STRING")

        log.info(
            f"[PPI] 네트워크: 노드 {self.G.number_of_nodes()} / "
            f"엣지 {self.G.number_of_edges()}"
        )
        return self.G


# ─────────────────────────────────────────────────────────────────────────────
# 3. 통합 네트워크
# ─────────────────────────────────────────────────────────────────────────────
class IntegratedNetwork:
    """
    HCT + PPI 통합 네트워크
    - 공유 타겟(허브 유전자) 식별
    - 한약-양약 상호작용 교집합 분석
    """

    def __init__(self, hct: HCTNetwork, ppi: PPINetwork):
        self.hct = hct
        self.ppi = ppi
        # HCT(DiGraph) + PPI(Graph) → 무방향 통합 그래프
        # PPI 먼저, HCT 나중 → HCT node_type("target") 우선 보존
        self.G = nx.Graph()
        self.G.update(ppi.G)
        self.G.update(hct.G)
        log.info(
            f"[통합] 노드 {self.G.number_of_nodes()} / 엣지 {self.G.number_of_edges()}"
        )

    def get_shared_targets(self, herb_a: str, herb_b: str) -> list:
        """두 허브가 공유하는 타겟 유전자 목록"""
        def targets_of(herb):
            result = set()
            for comp in self.hct.G.successors(herb):
                for t in self.hct.G.successors(comp):
                    result.add(t)
            return result

        a_targets = targets_of(herb_a)
        b_targets = targets_of(herb_b)
        shared    = a_targets & b_targets
        log.info(f"[공유타겟] {herb_a} ∩ {herb_b} = {len(shared)}개")
        return sorted(shared)

    def get_shared_targets_herb_drug(self, herb_name: str, drug_targets: list) -> list:
        """한약 타겟 ∩ 양약 타겟"""
        herb_t = set()
        for comp in self.hct.G.successors(herb_name):
            for t in self.hct.G.successors(comp):
                herb_t.add(t)
        drug_t  = set(g.strip().upper() for g in drug_targets)
        shared  = herb_t & drug_t
        log.info(f"[한약-양약] {herb_name} ∩ 양약 = {len(shared)}개")
        return sorted(shared)


# ─────────────────────────────────────────────────────────────────────────────
# 4. 위상 분석
# ─────────────────────────────────────────────────────────────────────────────
class TopologyAnalyzer:
    """
    네트워크 위상 지표 계산
    논문에 필수로 보고되는 지표:
      - degree (차수)
      - betweenness centrality (매개 중심성)
      - closeness centrality (근접 중심성)
      - clustering coefficient
    """

    def __init__(self, G: nx.Graph | nx.DiGraph):
        self.G = G
        # 무방향 처리 (중심성 계산용)
        self.G_undirected = G.to_undirected() if G.is_directed() else G

    def compute_all(self, target_nodes: list = None) -> pd.DataFrame:
        """
        모든 위상 지표 계산 → DataFrame 반환
        target_nodes: 계산할 노드 한정 (None이면 전체)
        """
        log.info("[위상분석] 지표 계산 시작...")

        nodes = target_nodes or list(self.G.nodes())

        degree       = dict(self.G_undirected.degree())
        betweenness  = nx.betweenness_centrality(self.G_undirected, normalized=True)
        closeness    = nx.closeness_centrality(self.G_undirected)
        clustering   = nx.clustering(self.G_undirected)

        rows = []
        for n in nodes:
            attr = self.G.nodes[n]
            rows.append({
                "node":          n,
                "node_type":     attr.get("node_type", "unknown"),
                "label":         attr.get("label", n),
                "degree":        degree.get(n, 0),
                "betweenness":   round(betweenness.get(n, 0), 6),
                "closeness":     round(closeness.get(n,   0), 6),
                "clustering":    round(clustering.get(n,  0), 6),
            })

        df = pd.DataFrame(rows).sort_values("degree", ascending=False)
        log.info(f"[위상분석] {len(df)}개 노드 지표 계산 완료")
        return df

    def get_hub_nodes(self, df: pd.DataFrame, top_n: int = 20,
                      node_type: str = "target") -> pd.DataFrame:
        """
        허브 노드 식별 (degree 상위 노드)
        node_type: 'target', 'compound', 'herb', None(전체)
        """
        filtered = df[df["node_type"] == node_type] if node_type else df
        hubs = filtered.nlargest(top_n, "degree")
        log.info(f"[허브] 상위 {len(hubs)}개 {node_type or '전체'} 노드:")
        for _, row in hubs.iterrows():
            log.info(f"  {row['label']:25s}  degree={row['degree']:4d}  "
                     f"betweenness={row['betweenness']:.4f}")
        return hubs

    def compute_target_only(self, target_nodes: list) -> pd.DataFrame:
        """타겟 유전자만 대상으로 PPI 서브그래프 위상 분석"""
        sub = self.G_undirected.subgraph(target_nodes)
        analyzer = TopologyAnalyzer(sub)
        return analyzer.compute_all()


# ─────────────────────────────────────────────────────────────────────────────
# 5. 내보내기
# ─────────────────────────────────────────────────────────────────────────────
class NetworkExporter:
    """Cytoscape / 논문용 파일 내보내기"""

    def __init__(self, output_dir: Path = NETWORK_DIR):
        self.out = Path(output_dir)
        self.out.mkdir(exist_ok=True)

    def to_graphml(self, G: nx.Graph | nx.DiGraph, filename: str):
        """GraphML 내보내기 (Cytoscape에서 직접 열기 가능)"""
        # GraphML은 bool/None 값을 지원하지 않으므로 문자열 변환
        G2 = G.copy()
        for n, data in G2.nodes(data=True):
            for k, v in list(data.items()):
                if v is None:
                    data[k] = ""
                elif not isinstance(v, (str, int, float)):
                    data[k] = str(v)

        path = self.out / filename
        nx.write_graphml(G2, path)
        log.info(f"[내보내기] GraphML → {path}")
        return path

    def to_edge_csv(self, G: nx.Graph | nx.DiGraph, filename: str):
        """엣지 리스트 CSV (Cytoscape import / 논문 보충 자료)"""
        rows = []
        for u, v, data in G.edges(data=True):
            rows.append({
                "source":    u,
                "target":    v,
                "edge_type": data.get("edge_type", data.get("source", "")),
                "weight":    data.get("weight", 1.0),
            })
        df = pd.DataFrame(rows)
        path = self.out / filename
        df.to_csv(path, index=False)
        log.info(f"[내보내기] 엣지 CSV ({len(df)}행) → {path}")
        return path

    def to_node_csv(self, G: nx.Graph | nx.DiGraph, filename: str):
        """노드 속성 CSV (Cytoscape node table)"""
        rows = []
        for n, data in G.nodes(data=True):
            rows.append({"id": n, **data})
        df = pd.DataFrame(rows)
        path = self.out / filename
        df.to_csv(path, index=False)
        log.info(f"[내보내기] 노드 CSV ({len(df)}행) → {path}")
        return path

    def to_cytoscape_json(self, G: nx.Graph | nx.DiGraph, filename: str):
        """Cytoscape.js JSON 포맷 (웹 시각화용)"""
        data = nx.cytoscape_data(G)
        path = self.out / filename
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        log.info(f"[내보내기] Cytoscape JSON → {path}")
        return path

    def topology_to_csv(self, df: pd.DataFrame, filename: str):
        """위상 분석 결과 CSV"""
        path = self.out / filename
        df.to_csv(path, index=False)
        log.info(f"[내보내기] 위상 분석 CSV ({len(df)}행) → {path}")
        return path


# ─────────────────────────────────────────────────────────────────────────────
# 6. 파이프라인 통합 함수
# ─────────────────────────────────────────────────────────────────────────────
def draw_network_png(
    G:           nx.Graph,
    hub_genes:   list,
    prefix:      str  = "analysis",
    top_n_label: int  = 15,
    dpi:         int  = 300,
) -> Path:
    """
    통합 네트워크 시각화 PNG 생성 (논문용 300 DPI)

    노드 색:
      herb     : #2E7D32 (진녹색)
      compound : #1565C0 (진파랑)
      target   : #C62828 (진빨강, hub은 주황)
      drug     : #6A1B9A (보라)
    """
    RESULTS_DIR = BASE_DIR / "results"
    RESULTS_DIR.mkdir(exist_ok=True)

    # ── 레이아웃 ─────────────────────────────────────────────────
    # 노드 수가 많으면 spring, 적으면 kamada_kawai
    n = G.number_of_nodes()
    if n > 150:
        pos = nx.spring_layout(G, k=1.8 / (n ** 0.5), seed=42, iterations=50)
    else:
        try:
            pos = nx.kamada_kawai_layout(G)
        except Exception:
            pos = nx.spring_layout(G, seed=42)

    hub_set = set(hub_genes)

    # ── 노드 속성 분리 ────────────────────────────────────────────
    node_color, node_size, node_alpha = [], [], []
    for node in G.nodes():
        ntype = G.nodes[node].get("node_type", "")
        if ntype == "herb":
            node_color.append("#2E7D32")
            node_size.append(400)
            node_alpha.append(0.9)
        elif ntype == "compound":
            node_color.append("#1565C0")
            node_size.append(120)
            node_alpha.append(0.6)
        elif ntype == "drug":
            node_color.append("#6A1B9A")
            node_size.append(350)
            node_alpha.append(0.9)
        elif node in hub_set:          # hub target
            node_color.append("#E65100")
            node_size.append(280)
            node_alpha.append(0.95)
        else:                          # non-hub target
            node_color.append("#C62828")
            node_size.append(60)
            node_alpha.append(0.45)

    # ── 엣지 색 ────────────────────────────────────────────────────
    edge_color, edge_alpha, edge_width = [], [], []
    for u, v, d in G.edges(data=True):
        etype = d.get("edge_type", "")
        if etype == "herb_compound":
            edge_color.append("#43A047"); edge_alpha.append(0.5); edge_width.append(0.8)
        elif etype == "compound_target":
            edge_color.append("#42A5F5"); edge_alpha.append(0.4); edge_width.append(0.5)
        else:                                # PPI
            edge_color.append("#BDBDBD"); edge_alpha.append(0.25); edge_width.append(0.3)

    # ── 그리기 ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 11), facecolor="white")

    nx.draw_networkx_edges(
        G, pos, ax=ax,
        edge_color=edge_color,
        alpha=0.35,
        width=[w * 0.7 for w in edge_width],
    )
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color=node_color,
        node_size=node_size,
        alpha=0.85,
        linewidths=0.3,
        edgecolors="white",
    )

    # Hub 유전자 라벨 (상위 top_n_label 개)
    hub_nodes_sorted = [n for n in G.nodes() if n in hub_set][:top_n_label]
    label_dict = {n: n for n in hub_nodes_sorted}
    # Herb 라벨도 추가
    for node in G.nodes():
        if G.nodes[node].get("node_type") in ("herb", "drug"):
            label_dict[node] = node.split()[0]  # 첫 단어만

    nx.draw_networkx_labels(
        G, pos, labels=label_dict, ax=ax,
        font_size=6.5, font_color="#212121",
        font_weight="bold",
        bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.55, lw=0),
    )

    # ── 범례 ────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(color="#2E7D32", label="Herb"),
        mpatches.Patch(color="#6A1B9A", label="Drug"),
        mpatches.Patch(color="#1565C0", label="Compound"),
        mpatches.Patch(color="#E65100", label="Hub Target"),
        mpatches.Patch(color="#C62828", label="Target"),
    ]
    ax.legend(handles=legend_handles, loc="upper left",
              fontsize=9, framealpha=0.8, edgecolor="#BDBDBD")

    ax.set_title(
        f"Integrated Herb-Compound-Target Network\n"
        f"({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)",
        fontsize=13, fontweight="bold", pad=10,
    )
    ax.axis("off")
    plt.tight_layout(pad=0.5)

    out = RESULTS_DIR / f"{prefix}_network.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"[네트워크 PNG] → {out}")

    # ── 허브 유전자 막대 그래프 ────────────────────────────────────
    _draw_hub_ranking(G, hub_genes, prefix, RESULTS_DIR, dpi)

    return out


def _draw_hub_ranking(
    G: nx.Graph, hub_genes: list, prefix: str,
    out_dir: Path, dpi: int = 300, top_n: int = 20,
) -> Path:
    """Hub 유전자 Degree 막대 그래프"""
    degrees = dict(G.degree())
    hub_deg = [(g, degrees.get(g, 0)) for g in hub_genes if g in degrees]
    hub_deg.sort(key=lambda x: -x[1])
    hub_deg = hub_deg[:top_n]

    if not hub_deg:
        return None

    genes  = [x[0] for x in hub_deg]
    values = [x[1] for x in hub_deg]

    # 색: 상위 5개 강조
    colors = ["#E65100" if i < 5 else "#EF9A9A" for i in range(len(genes))]

    fig, ax = plt.subplots(figsize=(8, 6), facecolor="white")
    bars = ax.barh(range(len(genes)), values, color=colors, edgecolor="white", height=0.7)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Degree", fontsize=11)
    ax.set_title(f"Hub Target Genes (Top {top_n} by Degree)", fontsize=12, fontweight="bold")

    for bar, val in zip(bars, values):
        ax.text(val + 0.5, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=9, color="#424242")

    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()

    out = out_dir / f"{prefix}_hub_ranking.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"[허브 랭킹 PNG] → {out}")
    return out


def build_full_network(
    herb_results: dict,
    ppi_edges:    list,
    prefix:       str = "analysis",
    min_ppi_score: int = 400,
    make_plots:    bool = True,
) -> dict:
    """
    수집 결과 → 전체 네트워크 구축 + 분석 + 내보내기 원스텝 함수

    반환:
        {
            "hct":        HCTNetwork,
            "ppi":        PPINetwork,
            "integrated": IntegratedNetwork,
            "topology":   pd.DataFrame,
            "hubs":       pd.DataFrame,
            "files":      {graphml, edge_csv, node_csv, topology_csv}
        }
    """
    exporter = NetworkExporter()

    # 1. HCT 네트워크
    hct = HCTNetwork()
    hct.build_from_herb_results(herb_results)

    # 2. PPI 네트워크
    ppi = PPINetwork()
    ppi.build_from_string(ppi_edges, min_score=min_ppi_score)

    # 3. 통합
    integrated = IntegratedNetwork(hct, ppi)

    # 4. 위상 분석 (통합 네트워크 전체)
    analyzer  = TopologyAnalyzer(integrated.G)
    topo_df   = analyzer.compute_all()
    hub_df    = analyzer.get_hub_nodes(topo_df, top_n=20, node_type="target")

    # 5. 내보내기
    files = {
        "hct_graphml":    exporter.to_graphml(hct.G,         f"{prefix}_hct.graphml"),
        "ppi_graphml":    exporter.to_graphml(ppi.G,         f"{prefix}_ppi.graphml"),
        "full_graphml":   exporter.to_graphml(integrated.G,  f"{prefix}_full.graphml"),
        "edge_csv":       exporter.to_edge_csv(integrated.G, f"{prefix}_edges.csv"),
        "node_csv":       exporter.to_node_csv(integrated.G, f"{prefix}_nodes.csv"),
        "topology_csv":   exporter.topology_to_csv(topo_df,  f"{prefix}_topology.csv"),
        "cytoscape_json": exporter.to_cytoscape_json(integrated.G, f"{prefix}_cytoscape.json"),
    }

    # 6. 네트워크 시각화 PNG
    if make_plots:
        hub_gene_list = hub_df["label"].tolist() if not hub_df.empty else []
        network_png   = draw_network_png(integrated.G, hub_gene_list, prefix=prefix)
        files["network_png"]     = network_png
        files["hub_ranking_png"] = BASE_DIR / "results" / f"{prefix}_hub_ranking.png"

    return {
        "hct":        hct,
        "ppi":        ppi,
        "integrated": integrated,
        "topology":   topo_df,
        "hubs":       hub_df,
        "files":      files,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 실행 예시 (더미 데이터)
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # 더미 herb_results (실제 사용 시 collect.py NetPharmCollector 결과 전달)
    dummy_herb_results = {
        "Panax ginseng": {
            "herb_001": {
                "herb": {"herb_name_cn": "人参", "herb_id": "herb_001"},
                "compounds": [
                    {"compound_id": "CID_001", "compound_name": "Ginsenoside Rb1",
                     "OB": 36.19, "DL": 0.75},
                    {"compound_id": "CID_002", "compound_name": "Ginsenoside Rg1",
                     "OB": 36.19, "DL": 0.30},
                ],
                "targets": [
                    {"compound_id": "CID_001", "gene_symbol": "TP53"},
                    {"compound_id": "CID_001", "gene_symbol": "AKT1"},
                    {"compound_id": "CID_002", "gene_symbol": "EGFR"},
                    {"compound_id": "CID_002", "gene_symbol": "TP53"},
                ],
            }
        },
        "Astragalus membranaceus": {
            "herb_002": {
                "herb": {"herb_name_cn": "黄芪", "herb_id": "herb_002"},
                "compounds": [
                    {"compound_id": "CID_003", "compound_name": "Astragaloside IV",
                     "OB": 38.51, "DL": 0.20},
                ],
                "targets": [
                    {"compound_id": "CID_003", "gene_symbol": "AKT1"},
                    {"compound_id": "CID_003", "gene_symbol": "TNF"},
                ],
            }
        },
    }

    dummy_ppi = [
        {"preferredName_A": "TP53", "preferredName_B": "AKT1", "score": 0.9},
        {"preferredName_A": "AKT1", "preferredName_B": "EGFR", "score": 0.85},
        {"preferredName_A": "TP53", "preferredName_B": "TNF",  "score": 0.7},
    ]

    result = build_full_network(dummy_herb_results, dummy_ppi, prefix="test")

    print("\n=== 허브 유전자 TOP 20 ===")
    print(result["hubs"][["label", "degree", "betweenness"]].to_string(index=False))

    print("\n=== 생성된 파일 ===")
    for k, v in result["files"].items():
        print(f"  {k:20s}: {v}")

    shared = result["integrated"].get_shared_targets(
        "Panax ginseng", "Astragalus membranaceus"
    )
    print(f"\n공유 타겟: {shared}")
