"""
투고용 데이터 패키징 + 재현성 메타데이터 로깅
Feature 5: Excel 다중 시트 + SVG/PDF 고해상도
Feature 6: DB 버전·쿼리 일시·UniProt 매핑 로그
"""

import json
import logging
import platform
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR    = Path(__file__).parent.parent
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# 메타데이터 로거 (Feature 6)
# ─────────────────────────────────────────────────────────────────────────────
class MetadataLogger:
    """
    분석 재현성 보장 메타데이터 자동 수집
    - 외부 DB 쿼리 기록 (DB명, 접근 일시, 파라미터)
    - Python/패키지 버전
    - UniProt ID 매핑 로그
    """

    def __init__(self):
        self.records: list[dict] = []
        self.start_time = datetime.now().isoformat()

    def log_db_access(
        self,
        db_name:    str,
        endpoint:   str,
        query:      dict,
        n_results:  int,
        version:    str = "N/A",
    ):
        self.records.append({
            "timestamp":  datetime.now().isoformat(),
            "db_name":    db_name,
            "endpoint":   endpoint,
            "query":      json.dumps(query, ensure_ascii=False),
            "n_results":  n_results,
            "db_version": version,
        })

    def log_uniprot_mapping(self, gene_symbol: str, uniprot_id: str, source: str):
        self.records.append({
            "timestamp":  datetime.now().isoformat(),
            "db_name":    "UniProt",
            "endpoint":   "gene_mapping",
            "query":      gene_symbol,
            "n_results":  1 if uniprot_id else 0,
            "db_version": uniprot_id,  # UniProt accession을 버전 필드에 기록
        })

    def get_env_info(self) -> dict:
        """실행 환경 정보"""
        import importlib.metadata
        packages = ["networkx", "pandas", "numpy", "scipy",
                    "rdkit", "requests", "google-genai"]
        versions = {}
        for pkg in packages:
            try:
                versions[pkg] = importlib.metadata.version(pkg)
            except Exception:
                versions[pkg] = "N/A"
        return {
            "python_version": sys.version,
            "platform":       platform.platform(),
            "analysis_start": self.start_time,
            "package_versions": versions,
        }

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.records) if self.records else pd.DataFrame(
            columns=["timestamp", "db_name", "endpoint", "query", "n_results", "db_version"]
        )

    def to_text_report(self) -> str:
        env = self.get_env_info()
        lines = [
            "# Reproducibility Metadata Report",
            f"Generated: {datetime.now().isoformat()}",
            "",
            "## Environment",
            f"Python: {env['python_version']}",
            f"Platform: {env['platform']}",
            f"Analysis started: {env['analysis_start']}",
            "",
            "## Package Versions",
        ]
        for pkg, ver in env["package_versions"].items():
            lines.append(f"  {pkg}: {ver}")
        lines += [
            "",
            "## Database Access Log",
            f"{'Timestamp':<26} {'Database':<15} {'Endpoint':<25} {'Results':>8}",
            "-" * 80,
        ]
        for r in self.records:
            lines.append(
                f"{r['timestamp']:<26} {r['db_name']:<15} "
                f"{r['endpoint']:<25} {r['n_results']:>8}"
            )
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Excel 패키져 (Feature 5)
# ─────────────────────────────────────────────────────────────────────────────
class SupplementaryExporter:
    """
    투고용 Supplementary Data Excel (.xlsx) 자동 생성
    각 분석 결과를 별도 시트로 구성
    """

    def __init__(self, output_dir: Path = RESULTS_DIR):
        self.out = Path(output_dir)
        self.out.mkdir(exist_ok=True)

    def build_excel(
        self,
        prefix:        str,
        herb_results:  dict    = None,
        admet_df:      pd.DataFrame = None,
        topology_df:   pd.DataFrame = None,
        enrich_results: dict   = None,
        intersection:  dict    = None,
        disease_targets: list  = None,
        metadata:      MetadataLogger = None,
    ) -> Path:
        """
        모든 분석 결과를 하나의 xlsx로 패키징
        시트 구성:
          Sheet1: Compound_ADMET
          Sheet2: Network_Topology
          Sheet3: GO_BP / Sheet4: GO_MF / Sheet5: GO_CC
          Sheet6: KEGG_Pathways
          Sheet7: Reactome
          Sheet8: Disease_Targets
          Sheet9: Shared_Genes
          Sheet10: Metadata
        """
        path = self.out / f"{prefix}_supplementary.xlsx"

        with pd.ExcelWriter(path, engine="openpyxl") as writer:

            # ── Sheet 1: 화합물 ADMET ────────────────────────────────────
            if admet_df is not None and not admet_df.empty:
                cols = [c for c in [
                    "compound_name", "compound_id", "ob", "dl",
                    "mw", "logp", "hbd", "hba", "tpsa", "rotbonds",
                    "pass_filter", "fail_reasons"
                ] if c in admet_df.columns]
                admet_df[cols].to_excel(writer, sheet_name="Compound_ADMET", index=False)
                self._format_sheet(writer, "Compound_ADMET")

            # ── Sheet 2: 네트워크 토폴로지 ───────────────────────────────
            if topology_df is not None and not topology_df.empty:
                cols = [c for c in [
                    "rank", "label", "node_type", "degree",
                    "betweenness", "closeness", "eigenvector", "mcc", "composite_score"
                ] if c in topology_df.columns]
                topology_df[cols].to_excel(writer, sheet_name="Network_Topology", index=False)
                self._format_sheet(writer, "Network_Topology")

            # ── Sheet 3-7: 농축분석 결과 ─────────────────────────────────
            if enrich_results:
                sheet_map = {
                    "GO_BP":    "GO_Biological_Process",
                    "GO_MF":    "GO_Molecular_Function",
                    "GO_CC":    "GO_Cellular_Component",
                    "KEGG":     "KEGG_Pathways",
                    "Reactome": "Reactome",
                }
                for alias, sheet_name in sheet_map.items():
                    df = enrich_results.get(alias)
                    if df is not None and not df.empty:
                        cols = [c for c in [
                            "term", "gene_count", "p_value", "adj_p_value",
                            "combined_score", "genes"
                        ] if c in df.columns]
                        df[cols].to_excel(writer, sheet_name=sheet_name, index=False)
                        self._format_sheet(writer, sheet_name)

            # ── Sheet 8: 질환 타겟 ───────────────────────────────────────
            if disease_targets:
                pd.DataFrame(disease_targets).to_excel(
                    writer, sheet_name="Disease_Targets", index=False
                )

            # ── Sheet 9: 공유 유전자 ─────────────────────────────────────
            if intersection and intersection.get("shared"):
                shared_df = pd.DataFrame({
                    "shared_gene":     intersection["shared"],
                    "herb_label":      intersection.get("herb_label", ""),
                    "disease_label":   intersection.get("disease_label", ""),
                })
                shared_df.to_excel(writer, sheet_name="Shared_Genes", index=False)

            # ── Sheet 10: 메타데이터 ─────────────────────────────────────
            if metadata:
                meta_df = metadata.to_dataframe()
                if not meta_df.empty:
                    meta_df.to_excel(writer, sheet_name="Metadata_Log", index=False)
                env = metadata.get_env_info()
                env_df = pd.DataFrame([
                    {"key": "python_version", "value": env["python_version"]},
                    {"key": "platform",       "value": env["platform"]},
                    {"key": "analysis_start", "value": env["analysis_start"]},
                ] + [
                    {"key": f"pkg_{k}", "value": v}
                    for k, v in env["package_versions"].items()
                ])
                env_df.to_excel(writer, sheet_name="Environment", index=False)

        log.info(f"[Excel] Supplementary Data → {path}")
        return path

    @staticmethod
    def _format_sheet(writer, sheet_name: str):
        """열 너비 자동 조정"""
        try:
            ws = writer.sheets[sheet_name]
            for col in ws.columns:
                max_len = max(
                    (len(str(cell.value)) for cell in col if cell.value), default=10
                )
                ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 50)
        except Exception:
            pass


# ─────────────────────────────────────────────────────────────────────────────
# 고해상도 이미지 내보내기 (SVG / PDF)
# ─────────────────────────────────────────────────────────────────────────────
class FigureExporter:
    """PNG → SVG / PDF 변환 (논문 투고용 벡터 이미지)"""

    def __init__(self, output_dir: Path = RESULTS_DIR):
        self.out = Path(output_dir)
        self.out.mkdir(exist_ok=True)

    def save_svg(self, fig: plt.Figure, filename: str) -> Path:
        path = self.out / filename
        fig.savefig(path, format="svg", bbox_inches="tight")
        log.info(f"[SVG] → {path}")
        return path

    def save_pdf(self, fig: plt.Figure, filename: str) -> Path:
        path = self.out / filename
        fig.savefig(path, format="pdf", bbox_inches="tight", dpi=300)
        log.info(f"[PDF] → {path}")
        return path

    def png_to_svg(self, png_path: Path) -> Path:
        """기존 PNG를 SVG로 재저장 (matplotlib 재렌더링)"""
        # PNG를 직접 SVG로 변환하는 대신, 동일 경로에 SVG 파일 생성 안내
        svg_path = png_path.with_suffix(".svg")
        log.info(f"[SVG] {png_path.name} → {svg_path.name} (matplotlib 재렌더링 필요)")
        return svg_path

    def export_enrichment_svg(
        self,
        enrich_results: dict,
        prefix:         str,
        top_n:          int = 15,
    ) -> list:
        """농축분석 결과를 SVG로 재생성"""
        from enrichment import EnrichmentPlotter, EnrichmentProcessor
        plotter   = EnrichmentPlotter(output_dir=self.out, dpi=300)
        processor = EnrichmentProcessor()
        paths     = []

        for alias, df in enrich_results.items():
            if df is None or df.empty or len(df) < 2:
                continue
            # Dot plot SVG
            sub = processor.add_neg_log10p(df).head(top_n).copy()
            sub["term_clean"] = sub["term"].apply(processor.clean_term)
            sub = sub.sort_values("-log10(adj_p)", ascending=True)

            fig, ax = plt.subplots(figsize=(9, max(4, len(sub) * 0.38)))
            sc = ax.scatter(
                sub["-log10(adj_p)"], sub["term_clean"],
                s=sub["gene_count"] * 12,
                c=sub["adj_p_value"], cmap="RdYlBu_r",
                vmin=0, vmax=0.05, alpha=0.85,
                linewidths=0.4, edgecolors="gray",
            )
            plt.colorbar(sc, ax=ax, shrink=0.5).set_label("Adj. p-value", fontsize=9)
            ax.set_xlabel("-log₁₀(Adjusted p-value)", fontsize=10)
            ax.set_title(f"{alias} Enrichment Analysis", fontsize=11, fontweight="bold")
            ax.tick_params(axis="y", labelsize=8)
            plt.tight_layout()

            svg_path = self.out / f"{prefix}_{alias}_dotplot.svg"
            fig.savefig(svg_path, format="svg", bbox_inches="tight")
            plt.close(fig)
            paths.append(svg_path)
            log.info(f"[SVG] {alias} dotplot → {svg_path}")

        return paths


# ─────────────────────────────────────────────────────────────────────────────
# 통합 내보내기
# ─────────────────────────────────────────────────────────────────────────────
def export_all(
    prefix:          str,
    pipeline_output: dict,
    admet_results:   dict   = None,
    topology_df:     pd.DataFrame = None,
    intersection:    dict   = None,
    disease_targets: list   = None,
    metadata:        MetadataLogger = None,
    export_svg:      bool   = True,
) -> dict:
    """
    모든 결과를 한 번에 내보내기
    반환: {excel, svgs, metadata_report}
    """
    supp  = SupplementaryExporter()
    fig_e = FigureExporter()

    enrich = pipeline_output.get("enrichment", {})
    enrich_results = (enrich or {}).get("results", {})

    admet_df = None
    if admet_results:
        admet_df = admet_results.get("summary")

    # Excel
    excel_path = supp.build_excel(
        prefix         = prefix,
        admet_df       = admet_df,
        topology_df    = topology_df,
        enrich_results = enrich_results,
        intersection   = intersection,
        disease_targets= disease_targets,
        metadata       = metadata,
    )

    # SVG
    svg_paths = []
    if export_svg and enrich_results:
        svg_paths = fig_e.export_enrichment_svg(enrich_results, prefix)

    # 메타데이터 텍스트 보고서
    meta_report_path = None
    if metadata:
        report_text = metadata.to_text_report()
        meta_report_path = RESULTS_DIR / f"{prefix}_metadata.txt"
        meta_report_path.write_text(report_text, encoding="utf-8")
        log.info(f"[Metadata] → {meta_report_path}")

    return {
        "excel":           excel_path,
        "svgs":            svg_paths,
        "metadata_report": meta_report_path,
    }
