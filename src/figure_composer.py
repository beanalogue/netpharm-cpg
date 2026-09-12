"""
논문 투고용 멀티패널 Figure 합성기
최대 4개 패널을 2×2 그리드로 배치, A/B/C/D 라벨링, 300+ DPI 출력

사용:
    from figure_composer import PanelComposer
    composer = PanelComposer()
    out = composer.compose(
        panels=[
            ("results/venn.png",            "Venn Diagram"),
            ("results/hub_ranking.png",     "Hub Gene Ranking (MCC)"),
            ("enrichment/combined_dot.png", "GO/KEGG Enrichment"),
            ("results/dock_wheel.png",      "Molecular Docking"),
        ],
        title="Network Pharmacology Analysis",
        prefix="my_study",
        dpi=300,
    )
    # out: {png: Path, svg: Path}
"""

import logging
from pathlib import Path
from typing import Optional
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

log = logging.getLogger(__name__)

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

try:
    from PIL import Image
    PIL_OK = True
except ImportError:
    PIL_OK = False


# ─────────────────────────────────────────────────────────────────────────────
# 유틸리티
# ─────────────────────────────────────────────────────────────────────────────
def _load_image(path: str | Path) -> Optional[np.ndarray]:
    """이미지 파일 → numpy array (실패 시 None)"""
    p = Path(path)
    if not p.exists():
        log.warning(f"[Composer] 파일 없음: {p}")
        return None
    try:
        if PIL_OK:
            img = Image.open(p).convert("RGB")
            return np.array(img)
        else:
            # PIL 없으면 matplotlib으로 로드
            return plt.imread(str(p))
    except Exception as e:
        log.warning(f"[Composer] 로드 실패 {p.name}: {e}")
        return None


def _placeholder(label: str, title: str) -> np.ndarray:
    """이미지 없을 때 빈 패널 생성"""
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.set_facecolor("#F5F5F5")
    ax.text(0.5, 0.5, f"Panel {label}\n({title})\n[Not generated]",
            ha="center", va="center", fontsize=11,
            color="#9E9E9E", style="italic",
            transform=ax.transAxes)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#BDBDBD")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=72, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    if PIL_OK:
        return np.array(Image.open(buf).convert("RGB"))
    return plt.imread(buf)


# ─────────────────────────────────────────────────────────────────────────────
# 메인 합성기
# ─────────────────────────────────────────────────────────────────────────────
class PanelComposer:
    """
    최대 4개 PNG/SVG 이미지를 학술지 규격 멀티패널 Figure로 합성

    패널 배치 (2×2):
      A (0,0) │ B (0,1)
      ────────┼────────
      C (1,0) │ D (1,1)
    """

    LABEL_FONTSIZE  = 22
    TITLE_FONTSIZE  = 13
    LABEL_COLOR     = "black"
    PANEL_PADDING   = 0.02   # 패널 간격 비율

    def compose(
        self,
        panels:   list[tuple],    # [(path_or_None, subtitle), ...]  최대 4개
        title:    str  = "",
        prefix:   str  = "figure",
        dpi:      int  = 300,
        figsize:  tuple = (16, 14),
        fmt:      list  = ("png", "svg"),
        output_dir: Path = None,
    ) -> dict:
        """
        panels: [(path_or_None, subtitle), ...]
          - path_or_None : PNG/SVG 경로 or None (빈 패널)
          - subtitle      : 패널 제목 (하단 표시)
        반환: {"png": Path, "svg": Path, "panels_used": int}
        """
        if output_dir is None:
            output_dir = RESULTS_DIR
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        n = min(len(panels), 4)
        labels = ["A","B","C","D"]

        # ── 그리드 설정 ────────────────────────────────────────────────
        nrows = 1 if n <= 2 else 2
        ncols = min(n, 2) if n <= 2 else 2
        if n == 1:
            nrows, ncols = 1, 1
        elif n == 3:
            nrows, ncols = 2, 2

        # 제목 여백
        title_space = 0.06 if title else 0.0
        fig = plt.figure(figsize=figsize, facecolor="white",
                         constrained_layout=False)
        if title:
            fig.suptitle(title, fontsize=16, fontweight="bold",
                         y=0.98, va="top")

        top = 1.0 - title_space - 0.01
        gs = gridspec.GridSpec(
            nrows, ncols,
            figure=fig,
            hspace=0.20, wspace=0.12,
            top=top, bottom=0.04, left=0.04, right=0.97,
        )

        used = 0
        for i, (path_or_none, subtitle) in enumerate(panels[:4]):
            row = i // 2
            col = i % 2
            ax  = fig.add_subplot(gs[row, col])

            # 이미지 로드
            img = None
            if path_or_none:
                img = _load_image(path_or_none)
            if img is None:
                img = _placeholder(labels[i], subtitle)
            else:
                used += 1

            ax.imshow(img, aspect="auto", interpolation="lanczos")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_edgecolor("#BDBDBD")

            # A/B/C/D 라벨 (좌상단 외부)
            ax.text(
                -0.02, 1.02, labels[i],
                transform=ax.transAxes,
                fontsize=self.LABEL_FONTSIZE,
                fontweight="bold",
                color=self.LABEL_COLOR,
                va="bottom", ha="right",
            )

            # 패널 하단 부제목
            if subtitle:
                ax.set_title(subtitle, fontsize=self.TITLE_FONTSIZE,
                             pad=4, color="#333333")

        # n == 3: 마지막 셀 (1,1) 비우기
        if n == 3 and nrows == 2 and ncols == 2:
            ax_empty = fig.add_subplot(gs[1, 1])
            ax_empty.axis("off")

        # GridSpec이 이미 여백을 관리하므로 tight_layout은 건너뜀

        # ── 저장 ──────────────────────────────────────────────────────
        out = {}
        if "png" in fmt:
            p = output_dir / f"{prefix}_figure_multipanel.png"
            fig.savefig(p, dpi=dpi, bbox_inches="tight",
                        facecolor="white", format="png")
            out["png"] = p
            log.info(f"[Composer] PNG → {p}")
        if "svg" in fmt:
            p = output_dir / f"{prefix}_figure_multipanel.svg"
            fig.savefig(p, format="svg", bbox_inches="tight")
            out["svg"] = p
            log.info(f"[Composer] SVG → {p}")

        plt.close(fig)
        out["panels_used"] = used
        return out

    def compose_custom(
        self,
        panel_paths: dict,   # {"A": path, "B": path, ...}
        subtitles:   dict,   # {"A": "Venn", ...}
        **kwargs,
    ) -> dict:
        """
        라벨 지정 방식 인터페이스
        panel_paths: {"A": path_or_None, "B": ..., "C": ..., "D": ...}
        """
        ordered = [
            (panel_paths.get(k), subtitles.get(k, ""))
            for k in ["A","B","C","D"]
            if k in panel_paths or k in subtitles
        ]
        return self.compose(ordered, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# 도킹 결과 히트맵
# ─────────────────────────────────────────────────────────────────────────────
def make_docking_heatmap(
    dock_df,
    prefix:     str  = "analysis",
    output_dir: Path = None,
    dpi:        int  = 150,
) -> Optional[Path]:
    """
    도킹 결과 DataFrame → 히트맵 PNG
    컬럼: compound, target, affinity
    반환: 저장된 PNG Path (실패 시 None)
    """
    try:
        import pandas as pd
        output_dir = Path(output_dir or RESULTS_DIR)
        output_dir.mkdir(exist_ok=True)

        df = dock_df[dock_df["success"] == True].copy() if "success" in dock_df.columns else dock_df.copy()
        if df.empty:
            return None

        # compound 이름 정리 (CHEMBL ID → 짧게 표시)
        df["compound_label"] = df["compound"].astype(str).str[:12]

        pivot = df.pivot_table(index="compound_label", columns="target",
                               values="affinity", aggfunc="min")
        pivot = pivot.sort_index(axis=0)  # compound 정렬
        # target을 평균 친화력 기준으로 정렬 (낮을수록 강함)
        pivot = pivot[pivot.mean().sort_values().index]

        fig, ax = plt.subplots(figsize=(max(5, len(pivot.columns)*1.4),
                                        max(3, len(pivot)*0.7) + 1))

        # 색상: 음수일수록 진파랑 (강한 결합)
        vmin = pivot.values[~np.isnan(pivot.values)].min() if pivot.notna().any().any() else -10
        vmax = min(0, pivot.values[~np.isnan(pivot.values)].max()) if pivot.notna().any().any() else 0
        im = ax.imshow(pivot.values, cmap="Blues_r", aspect="auto",
                       vmin=vmin, vmax=vmax)

        # 축 라벨
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=35, ha="right", fontsize=9)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index, fontsize=9)

        # 셀 값 표시
        for r in range(len(pivot.index)):
            for c in range(len(pivot.columns)):
                val = pivot.iloc[r, c]
                if not np.isnan(val):
                    ax.text(c, r, f"{val:.2f}", ha="center", va="center",
                            fontsize=8, color="white" if val < (vmin + vmax) / 2 else "black")

        cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label("Binding Affinity (kcal/mol)", fontsize=9)
        ax.set_title("Molecular Docking — Binding Affinity Heatmap", fontsize=11, pad=8)
        fig.tight_layout()

        out_path = output_dir / f"{prefix}_docking_heatmap.png"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        log.info(f"[Composer] Docking heatmap → {out_path}")
        return out_path
    except Exception as e:
        log.warning(f"[Composer] 도킹 히트맵 생성 실패: {e}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 단축 함수
# ─────────────────────────────────────────────────────────────────────────────
def make_publication_figure(
    venn_path:     Optional[str | Path] = None,
    hub_path:      Optional[str | Path] = None,
    dotplot_path:  Optional[str | Path] = None,
    docking_path:  Optional[str | Path] = None,
    title:         str  = "Network Pharmacology Analysis",
    prefix:        str  = "analysis",
    dpi:           int  = 300,
    output_dir:    Path = None,
) -> dict:
    """
    네트워크약리학 표준 4-패널 Figure 생성

    A: Venn/UpSet Diagram (Target Gene Intersection)
    B: Hub Gene Ranking / PPI Network (MCC)
    C: GO/KEGG Dot Plot
    D: Molecular Docking Interaction
    """
    panels = [
        (venn_path,    "A  Target Gene Intersection"),
        (hub_path,     "B  Hub Gene Analysis (MCC)"),
        (dotplot_path, "C  GO/KEGG Enrichment"),
        (docking_path, "D  Molecular Docking"),
    ]
    composer = PanelComposer()
    return composer.compose(
        panels     = panels,
        title      = title,
        prefix     = prefix,
        dpi        = dpi,
        output_dir = output_dir or RESULTS_DIR,
    )
