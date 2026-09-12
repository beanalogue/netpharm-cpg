"""
Network Pharmacology — Enrichment Analysis
Enrichr API 기반 GO / KEGG / Reactome 농축분석 + 시각화

지원 데이터베이스:
  GO_Biological_Process_2023
  GO_Molecular_Function_2023
  GO_Cellular_Component_2023
  KEGG_2021_Human
  Reactome_2022
  WikiPathways_2023_Human

출력:
  - enrichment/  폴더에 DB별 CSV
  - enrichment/  폴더에 dotplot, barplot (PNG/SVG)
  - 통합 요약 CSV (논문 보충 자료용)
"""

import time
import logging
from pathlib import Path

import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # 헤드리스 서버 — GUI 없이 렌더링
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE_DIR        = Path(__file__).parent.parent
ENRICHMENT_DIR  = BASE_DIR / "enrichment"
ENRICHMENT_DIR.mkdir(exist_ok=True)

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

log = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)

# ── 분석할 데이터베이스 목록 ─────────────────────────────────────────────────
DATABASES = {
    "GO_BP":       "GO_Biological_Process_2023",
    "GO_MF":       "GO_Molecular_Function_2023",
    "GO_CC":       "GO_Cellular_Component_2023",
    "KEGG":        "KEGG_2021_Human",
    "Reactome":    "Reactome_2022",
    "WikiPathways":"WikiPathways_2023_Human",
}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Enrichr API 클라이언트
# ─────────────────────────────────────────────────────────────────────────────
class EnrichrClient:
    """Enrichr REST API 래퍼"""

    COLUMNS = [
        "rank", "term", "p_value", "z_score",
        "combined_score", "genes", "adj_p_value",
        "old_p_value", "old_adj_p_value",
    ]

    def __init__(self, delay: float = 0.5):
        self.delay   = delay
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "netpharm-research/1.0"})

    def submit_gene_list(self, gene_symbols: list, description: str = "netpharm") -> int:
        """
        유전자 리스트 제출 → userListId 반환
        이후 enrich() 호출에 사용
        """
        payload = {
            "list":        (None, "\n".join(gene_symbols)),
            "description": (None, description),
        }
        r = self.session.post(f"{ENRICHR_BASE}/addList", files=payload, timeout=30)
        r.raise_for_status()
        user_list_id = r.json()["userListId"]
        log.info(f"[Enrichr] 유전자 {len(gene_symbols)}개 제출 → userListId={user_list_id}")
        return user_list_id

    def enrich(self, user_list_id: int, database: str) -> pd.DataFrame:
        """
        단일 데이터베이스 농축분석 수행
        반환: DataFrame (rank, term, p_value, adj_p_value, genes, ...)
        """
        time.sleep(self.delay)
        r = self.session.get(
            f"{ENRICHR_BASE}/enrich",
            params={"userListId": user_list_id, "backgroundType": database},
            timeout=60,
        )
        r.raise_for_status()
        raw = r.json().get(database, [])

        if not raw:
            log.warning(f"[Enrichr] {database}: 결과 없음")
            return pd.DataFrame(columns=self.COLUMNS)

        df = pd.DataFrame(raw, columns=self.COLUMNS)
        df["database"] = database
        df["genes"]    = df["genes"].apply(lambda x: ";".join(x) if isinstance(x, list) else x)
        df["gene_count"] = df["genes"].apply(lambda x: len(x.split(";")) if x else 0)
        log.info(f"[Enrichr] {database}: {len(df)}개 term 반환")
        return df

    def enrich_all(
        self,
        gene_symbols: list,
        databases:    dict  = None,
        adj_p_cutoff: float = 0.05,
        description:  str   = "netpharm",
    ) -> dict:
        """
        전체 DB 농축분석 수행
        반환: {db_alias: DataFrame}
        """
        if not gene_symbols:
            log.error("[Enrichr] 유전자 목록이 비어 있음")
            return {}

        dbs = databases or DATABASES
        user_list_id = self.submit_gene_list(gene_symbols, description)

        results = {}
        for alias, db_name in dbs.items():
            df = self.enrich(user_list_id, db_name)
            if df.empty:
                results[alias] = df
                continue
            # FDR 필터
            filtered = df[df["adj_p_value"] <= adj_p_cutoff].copy()
            filtered = filtered.sort_values("adj_p_value")
            results[alias] = filtered
            log.info(
                f"[Enrichr] {alias}: FDR≤{adj_p_cutoff} 통과 {len(filtered)}/{len(df)}개"
            )

        return results


# ─────────────────────────────────────────────────────────────────────────────
# 2. 결과 후처리
# ─────────────────────────────────────────────────────────────────────────────
class EnrichmentProcessor:
    """농축분석 결과 정제 및 요약"""

    @staticmethod
    def clean_term(term: str) -> str:
        """GO/KEGG term 이름에서 ID 부분 분리"""
        # 예: "apoptotic process (GO:0006915)" → "apoptotic process"
        if "(" in term:
            return term[:term.rfind("(")].strip()
        return term.strip()

    @staticmethod
    def add_neg_log10p(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["-log10(adj_p)"] = -np.log10(df["adj_p_value"].clip(lower=1e-300))
        return df

    @staticmethod
    def extract_go_id(term: str) -> str:
        """GO:XXXXXXX ID 추출"""
        if "GO:" in term:
            start = term.index("GO:")
            return term[start:start+10].rstrip(")")
        return ""

    def summarize(self, results: dict, top_n: int = 10) -> pd.DataFrame:
        """
        전체 DB 결과를 하나의 요약 DataFrame으로 통합
        각 DB에서 상위 top_n개 term 포함
        """
        frames = []
        for alias, df in results.items():
            if df.empty:
                continue
            top = df.head(top_n).copy()
            top["db_alias"]  = alias
            top["term_clean"] = top["term"].apply(self.clean_term)
            top = self.add_neg_log10p(top)
            frames.append(top)

        if not frames:
            return pd.DataFrame()

        summary = pd.concat(frames, ignore_index=True)
        cols = ["db_alias", "term_clean", "gene_count", "adj_p_value",
                "-log10(adj_p)", "combined_score", "genes"]
        return summary[[c for c in cols if c in summary.columns]]


# ─────────────────────────────────────────────────────────────────────────────
# 3. 시각화
# ─────────────────────────────────────────────────────────────────────────────
class EnrichmentPlotter:
    """
    논문 수준 시각화
      - Dot plot  (x=gene ratio 또는 -log10p, y=term, size=gene count, color=adj_p)
      - Bar plot  (x=-log10(adj_p), y=term)
    """

    PALETTE = {
        "GO_BP":       "#4C72B0",
        "GO_MF":       "#55A868",
        "GO_CC":       "#C44E52",
        "KEGG":        "#DD8452",
        "Reactome":    "#8172B2",
        "WikiPathways":"#937860",
    }

    def __init__(self, output_dir: Path = ENRICHMENT_DIR, dpi: int = 200):
        self.out = Path(output_dir)
        self.dpi = dpi

    def dotplot(
        self,
        df:        pd.DataFrame,
        db_alias:  str,
        top_n:     int   = 20,
        filename:  str   = None,
    ) -> Path:
        """
        Dot plot
        x축: -log10(adjusted p-value)
        y축: term (상위 top_n)
        점 크기: gene count
        점 색상: adjusted p-value (낮을수록 진한 색)
        """
        proc = EnrichmentProcessor()
        sub  = proc.add_neg_log10p(df).head(top_n).copy()
        sub["term_clean"] = sub["term"].apply(proc.clean_term)
        sub = sub.sort_values("-log10(adj_p)", ascending=True)  # 아래→위 정렬

        fig, ax = plt.subplots(figsize=(9, max(4, len(sub) * 0.38)))

        sc = ax.scatter(
            sub["-log10(adj_p)"],
            sub["term_clean"],
            s=sub["gene_count"] * 12,
            c=sub["adj_p_value"],
            cmap="RdYlBu_r",
            vmin=0, vmax=0.05,
            alpha=0.85,
            linewidths=0.4,
            edgecolors="gray",
        )

        cbar = plt.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
        cbar.set_label("Adjusted p-value", fontsize=9)

        # 범례 (점 크기)
        for size_val in [1, 5, 10, 20]:
            ax.scatter([], [], s=size_val * 12, c="gray", alpha=0.6,
                       label=f"{size_val} genes")
        ax.legend(
            title="Gene count", loc="lower right",
            fontsize=8, title_fontsize=8,
            framealpha=0.7,
        )

        ax.axvline(-np.log10(0.05), color="red", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.set_xlabel("-log₁₀(Adjusted p-value)", fontsize=10)
        ax.set_ylabel("")
        ax.set_title(f"{db_alias} Enrichment Analysis", fontsize=11, fontweight="bold")
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=9)
        plt.tight_layout()

        fname = filename or f"{db_alias}_dotplot.png"
        path  = self.out / fname
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] Dot plot → {path}")
        return path

    def barplot(
        self,
        df:       pd.DataFrame,
        db_alias: str,
        top_n:    int  = 20,
        filename: str  = None,
    ) -> Path:
        """
        수평 Bar plot (-log10 adj_p)
        """
        proc = EnrichmentProcessor()
        sub  = proc.add_neg_log10p(df).head(top_n).copy()
        sub["term_clean"] = sub["term"].apply(proc.clean_term)
        sub = sub.sort_values("-log10(adj_p)", ascending=True)

        color = self.PALETTE.get(db_alias, "#4C72B0")
        fig, ax = plt.subplots(figsize=(9, max(4, len(sub) * 0.38)))

        bars = ax.barh(
            sub["term_clean"],
            sub["-log10(adj_p)"],
            color=color, alpha=0.82,
            edgecolor="white", linewidth=0.5,
        )
        # gene count 레이블
        for bar, count in zip(bars, sub["gene_count"]):
            ax.text(
                bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                f"n={count}", va="center", ha="left", fontsize=7, color="#555555",
            )

        ax.axvline(-np.log10(0.05), color="red", linestyle="--",
                   linewidth=0.8, alpha=0.6, label="FDR=0.05")
        ax.set_xlabel("-log₁₀(Adjusted p-value)", fontsize=10)
        ax.set_title(f"{db_alias} Top {top_n} Enriched Terms", fontsize=11, fontweight="bold")
        ax.tick_params(axis="y", labelsize=8)
        ax.legend(fontsize=8)
        plt.tight_layout()

        fname = filename or f"{db_alias}_barplot.png"
        path  = self.out / fname
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] Bar plot → {path}")
        return path

    def combined_dotplot(
        self,
        results: dict,
        top_n_per_db: int = 5,
        filename: str = "combined_dotplot.png",
    ) -> Path:
        """
        여러 DB 결과를 하나의 dot plot에 통합 (논문 Figure용)
        """
        proc   = EnrichmentProcessor()
        frames = []
        for alias, df in results.items():
            if df.empty:
                continue
            sub = proc.add_neg_log10p(df).head(top_n_per_db).copy()
            sub["term_clean"] = sub["term"].apply(proc.clean_term)
            sub["db_alias"]   = alias
            frames.append(sub)

        if not frames:
            log.warning("[플롯] 통합 dotplot: 데이터 없음")
            return None

        all_df = pd.concat(frames, ignore_index=True)
        all_df = all_df.sort_values(["db_alias", "-log10(adj_p)"], ascending=[True, False])

        # y축: DB별로 구분되는 term 레이블
        all_df["y_label"] = all_df["term_clean"].str[:55]

        fig, ax = plt.subplots(figsize=(11, max(6, len(all_df) * 0.38)))

        for i, (_, row) in enumerate(all_df.iterrows()):
            color = self.PALETTE.get(row["db_alias"], "#888888")
            ax.scatter(
                row["-log10(adj_p)"], i,
                s=row["gene_count"] * 14,
                c=[color], alpha=0.85,
                linewidths=0.4, edgecolors="gray",
                zorder=3,
            )

        ax.set_yticks(range(len(all_df)))
        ax.set_yticklabels(all_df["y_label"], fontsize=7.5)
        ax.axvline(-np.log10(0.05), color="red", linestyle="--",
                   linewidth=0.8, alpha=0.6)
        ax.set_xlabel("-log₁₀(Adjusted p-value)", fontsize=10)
        ax.set_title("GO / KEGG Enrichment Analysis", fontsize=12, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)

        # DB별 색상 범례
        patches = [
            mpatches.Patch(color=c, label=alias)
            for alias, c in self.PALETTE.items()
            if alias in results and not results[alias].empty
        ]
        ax.legend(handles=patches, loc="lower right", fontsize=8,
                  title="Database", title_fontsize=8)

        plt.tight_layout()
        path = self.out / filename
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        log.info(f"[플롯] 통합 Dot plot → {path}")
        return path


# ─────────────────────────────────────────────────────────────────────────────
# 4. 파일 내보내기
# ─────────────────────────────────────────────────────────────────────────────
class EnrichmentExporter:
    """농축분석 결과 CSV 내보내기"""

    def __init__(self, output_dir: Path = ENRICHMENT_DIR):
        self.out = Path(output_dir)
        self.out.mkdir(exist_ok=True)

    def save_per_db(self, results: dict, prefix: str = "analysis") -> dict:
        """DB별 결과를 개별 CSV로 저장"""
        paths = {}
        for alias, df in results.items():
            if df.empty:
                continue
            path = self.out / f"{prefix}_{alias}.csv"
            df.to_csv(path, index=False)
            log.info(f"[저장] {alias} → {path} ({len(df)}개 term)")
            paths[alias] = path
        return paths

    def save_summary(self, summary_df: pd.DataFrame, prefix: str = "analysis") -> Path:
        """통합 요약 CSV 저장 (논문 보충 자료용)"""
        path = self.out / f"{prefix}_enrichment_summary.csv"
        summary_df.to_csv(path, index=False)
        log.info(f"[저장] 통합 요약 → {path} ({len(summary_df)}행)")
        return path


# ─────────────────────────────────────────────────────────────────────────────
# 5. 파이프라인 통합 함수
# ─────────────────────────────────────────────────────────────────────────────
def run_enrichment(
    gene_symbols:  list,
    prefix:        str   = "analysis",
    adj_p_cutoff:  float = 0.05,
    top_n:         int   = 20,
    databases:     dict  = None,
    make_plots:    bool  = True,
) -> dict:
    """
    유전자 리스트 → 전체 농축분석 + 시각화 + 저장 원스텝 함수

    반환:
        {
            "results":   {db_alias: DataFrame},
            "summary":   DataFrame,
            "files":     {csv_paths, plot_paths},
        }
    """
    if not gene_symbols:
        log.error("유전자 목록이 비어 있습니다.")
        return {}

    log.info(f"[농축분석] 시작: {len(gene_symbols)}개 유전자, prefix='{prefix}'")

    client    = EnrichrClient()
    processor = EnrichmentProcessor()
    exporter  = EnrichmentExporter()
    plotter   = EnrichmentPlotter()

    # 1. 농축분석
    results = client.enrich_all(
        gene_symbols,
        databases=databases or DATABASES,
        adj_p_cutoff=adj_p_cutoff,
        description=prefix,
    )

    # 2. 요약
    summary = processor.summarize(results, top_n=top_n)

    # 3. 저장
    csv_paths  = exporter.save_per_db(results, prefix=prefix)
    summary_path = None
    if not summary.empty:
        summary_path = exporter.save_summary(summary, prefix=prefix)

    # 4. 시각화
    plot_paths = {}
    if make_plots:
        for alias, df in results.items():
            if df.empty or len(df) < 2:
                continue
            plot_paths[f"{alias}_dot"] = plotter.dotplot(df, alias, top_n=top_n,
                                                          filename=f"{prefix}_{alias}_dotplot.png")
            plot_paths[f"{alias}_bar"] = plotter.barplot(df, alias, top_n=top_n,
                                                          filename=f"{prefix}_{alias}_barplot.png")

        non_empty = {k: v for k, v in results.items() if not v.empty}
        if len(non_empty) >= 2:
            plot_paths["combined"] = plotter.combined_dotplot(
                non_empty,
                top_n_per_db=5,
                filename=f"{prefix}_combined_dotplot.png",
            )

    return {
        "results":      results,
        "summary":      summary,
        "files": {
            "csv":     csv_paths,
            "summary": summary_path,
            "plots":   plot_paths,
        },
    }


# ─────────────────────────────────────────────────────────────────────────────
# 실행 예시
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # 예시: 네트워크약리학 논문에서 자주 등장하는 타겟 유전자
    TEST_GENES = [
        "TP53", "AKT1", "EGFR", "TNF", "IL6", "VEGFA", "MYC",
        "CASP3", "BCL2", "MAPK1", "MAPK3", "PIK3CA", "STAT3",
        "JUN", "FOS", "NF1", "PTEN", "CDK2", "CCND1", "HSP90AA1",
    ]

    output = run_enrichment(
        gene_symbols=TEST_GENES,
        prefix="test",
        adj_p_cutoff=0.05,
        top_n=20,
        make_plots=True,
    )

    print("\n=== 농축분석 요약 ===")
    if not output["summary"].empty:
        print(output["summary"][["db_alias", "term_clean", "gene_count", "adj_p_value"]]
              .head(15).to_string(index=False))

    print("\n=== 생성된 파일 ===")
    for category, items in output["files"].items():
        if isinstance(items, dict):
            for k, v in items.items():
                if v:
                    print(f"  [{category}] {k}: {v}")
        elif items:
            print(f"  [{category}]: {items}")
