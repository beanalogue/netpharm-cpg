"""
Network Pharmacology Pipeline — 메인 진입점
수집 → 네트워크 → 농축분석 → 요약 보고서 자동 실행

사용법:
  # 기본 실행 (config 파일)
  python src/pipeline.py --config configs/my_analysis.json

  # 커맨드라인 직접 입력
  python src/pipeline.py \\
      --herbs "Panax ginseng" "Astragalus membranaceus" \\
      --drugs "metformin" "aspirin" \\
      --name ginseng_vs_metformin

  # 도움말
  python src/pipeline.py --help
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(BASE_DIR / "src"))

from collect     import NetPharmCollector
from network     import build_full_network
from enrichment  import run_enrichment

RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

CONFIG_DIR = BASE_DIR / "configs"
CONFIG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("pipeline")


# ─────────────────────────────────────────────────────────────────────────────
# 설정 구조
# ─────────────────────────────────────────────────────────────────────────────
DEFAULT_CONFIG = {
    "name":             "analysis",
    # 한약: {herb_name: [compound1, compound2, ...]} 형식
    # 예: {"Panax ginseng": ["ginsenoside Rb1", "ginsenoside Rg1"]}
    "herbs":            {},
    # 양약: 약물명 리스트
    # 예: ["metformin", "aspirin"]
    "drugs":            [],
    # TCMSP 로컬 파일: {herb_name: {ingredients: path, targets: path}}
    "tcmsp_files":      {},
    "ppi_min_score":    400,
    "adj_p_cutoff":     0.05,
    "top_n_enrichment": 20,
    "adme_ob":          30.0,
    "adme_dl":          0.18,
    "string_species":   9606,
    "make_plots":       True,
    "enrichment_dbs":   None,  # None → 기본 6종 전체
}


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return {**DEFAULT_CONFIG, **cfg}


def save_config_template(name: str = "example"):
    """예시 config 파일 생성"""
    template = {
        "name": name,
        "_comment_herbs": "herb_name: [ChEMBL에서 검색할 활성 성분명(영문) 목록]",
        "herbs": {
            "Panax ginseng": [
                "ginsenoside Rb1", "ginsenoside Rg1", "panaxadiol"
            ],
            "Astragalus membranaceus": [
                "astragaloside IV", "calycosin", "formononetin"
            ]
        },
        "_comment_drugs": "양약명 리스트 (ChEMBL pref_name 기준)",
        "drugs": ["metformin", "aspirin"],
        "_comment_tcmsp": "TCMSP 로컬 파일 사용 시 (선택사항)",
        "tcmsp_files": {},
        "ppi_min_score":    400,
        "adj_p_cutoff":     0.05,
        "top_n_enrichment": 20,
        "adme_ob":          30.0,
        "adme_dl":          0.18,
        "make_plots":       True,
        "enrichment_dbs":   None
    }
    path = CONFIG_DIR / f"{name}.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(template, f, ensure_ascii=False, indent=2)
    print(f"Config 템플릿 생성: {path}")
    return path


# ─────────────────────────────────────────────────────────────────────────────
# 단계별 실행 함수
# ─────────────────────────────────────────────────────────────────────────────
def step_collect(cfg: dict, collector: NetPharmCollector) -> tuple[dict, dict, list]:
    """Step 1: 데이터 수집"""
    log.info("=" * 60)
    log.info("STEP 1 — 데이터 수집")
    log.info("=" * 60)

    herb_results = {}

    # ChEMBL 경유 한약 수집 (herbs: {herb_name: [compound, ...]})
    herbs_cfg = cfg.get("herbs", {})
    if isinstance(herbs_cfg, dict) and herbs_cfg:
        for herb_name, compounds in herbs_cfg.items():
            tcmsp_targets_map: dict = {}
            # TCMSP 항상 먼저 시도 (성분 직접 입력 여부 무관)
            try:
                from tcmsp_scraper import scrape_herb
                tcmsp_data = scrape_herb(herb_name, use_cache=True)
                tcmsp_compounds = tcmsp_data.get("compounds", [])
                if tcmsp_compounds:
                    compounds = [c["molecule_name"] for c in tcmsp_compounds]
                    # TCMSP 타깃 사전 구성 (재스크래핑 방지)
                    for c in tcmsp_compounds:
                        tgts = [
                            {
                                "gene_symbol":  t.get("gene_symbol", ""),
                                "uniprot_id":   "",
                                "target_name":  t.get("target_name", ""),
                                "action_type":  "",
                                "pchembl_value": None,
                                "source":       "TCMSP",
                            }
                            for t in c.get("targets", [])
                            if t.get("gene_symbol")
                        ]
                        if tgts:
                            tcmsp_targets_map[c["molecule_name"]] = tgts
                    log.info(f"[TCMSP 자동] {herb_name} → {len(compounds)}개 성분, "
                             f"타깃 사전 {len(tcmsp_targets_map)}개 성분 보유")
                elif compounds:
                    log.info(f"[TCMSP 없음] {herb_name} — 수동 입력 성분 {len(compounds)}개 사용")
            except Exception as e:
                log.warning(f"[TCMSP 자동] {herb_name} 조회 실패: {e}")

            # TCMSP·수동 입력 모두 없으면 → ChEMBL organism 검색 fallback
            if not compounds:
                log.warning(f"[TCMSP 없음] '{herb_name}' — ChEMBL organism 검색으로 fallback")
                try:
                    compounds = collector.chembl.get_herb_compounds_by_organism(herb_name)
                    if compounds:
                        log.info(f"[ChEMBL fallback] {herb_name} → {len(compounds)}개 성분")
                    else:
                        log.warning(f"[ChEMBL fallback] {herb_name}: 성분 없음 — config에 직접 입력 필요")
                except Exception as e:
                    log.warning(f"[ChEMBL fallback] {herb_name} 조회 실패: {e}")

            log.info(f"한약 수집 (ChEMBL): {herb_name} — {len(compounds)}개 성분")
            r = collector.collect_herb_compounds_chembl(
                herb_name, compounds, tcmsp_targets_map=tcmsp_targets_map
            )
            herb_results.update(r)
    elif isinstance(herbs_cfg, list) and herbs_cfg:
        # 하위 호환: 리스트 형식이면 각 항목을 herb_name으로만 처리 (경고)
        log.warning("herbs는 {herb_name: [compounds]} 형식 권장. 현재 리스트는 성분 없이 처리됩니다.")

    # TCMSP 로컬 파일 추가
    for herb_name, files in cfg.get("tcmsp_files", {}).items():
        log.info(f"한약 수집 (TCMSP 로컬): {herb_name}")
        r = collector.collect_from_tcmsp(
            herb_name,
            files["ingredients"],
            files["targets"],
            cfg.get("adme_ob", 30.0),
            cfg.get("adme_dl", 0.18),
        )
        herb_results.update(r)

    # 양약 수집 (ChEMBL)
    drug_results = {}
    if cfg.get("drugs"):
        log.info(f"양약 수집 (ChEMBL): {cfg['drugs']}")
        drug_results = collector.collect_drug_targets_chembl(cfg["drugs"])

    # 유전자 심볼 추출
    herb_genes = set(collector.extract_gene_symbols(herb_results))
    drug_genes = set(collector.extract_drug_genes(drug_results))
    all_genes  = sorted(herb_genes | drug_genes)

    log.info(
        f"수집 완료 — 한약 유전자: {len(herb_genes)}개 / "
        f"양약 유전자: {len(drug_genes)}개 / 합계: {len(all_genes)}개"
    )

    # ── 사전 QC 체크 ──────────────────────────────────────────────
    _qc_check(herb_results)

    return herb_results, drug_results, all_genes


def _qc_check(herb_results: dict):
    """수집된 약재별 성분·타겟 수를 검사해 부족한 경우 경고"""
    issues = []
    for herb, data in herb_results.items():
        comps = data.get("compounds", [])
        n_comp = len(comps)
        n_tgt  = sum(len(c.get("targets", [])) for c in comps)
        if n_comp == 0:
            issues.append(f"  ❌ {herb}: 성분 0개 — DB 미보유 또는 검색 실패")
        elif n_comp < 3:
            issues.append(f"  ⚠️  {herb}: 성분 {n_comp}개, 타겟 {n_tgt}개 — 데이터 부족")
        elif n_tgt == 0:
            issues.append(f"  ⚠️  {herb}: 성분 {n_comp}개이나 타겟 0개 — 타겟 수집 실패")

    if issues:
        log.warning("=" * 60)
        log.warning("QC 경고 — 데이터 부족 약재:")
        for msg in issues:
            log.warning(msg)
        log.warning("대응 방법:")
        log.warning("  1. config에 성분명 직접 지정")
        log.warning("  2. HERB/BATMAN-TCM 등 대체 DB 활용")
        log.warning("  3. 해당 약재를 분석에서 제외")
        log.warning("=" * 60)
        try:
            from telegram_notify import send as tg_send
            warn_lines = "\n".join(issues)
            tg_send(f"⚠️ <b>QC 경고</b>\n{warn_lines}")
        except Exception:
            pass


def step_network(cfg: dict, herb_results: dict, ppi_edges: list) -> dict:
    """Step 2: 네트워크 구축"""
    log.info("=" * 60)
    log.info("STEP 2 — 네트워크 구축")
    log.info("=" * 60)

    net = build_full_network(
        herb_results   = herb_results,
        ppi_edges      = ppi_edges,
        prefix         = cfg["name"],
        min_ppi_score  = cfg["ppi_min_score"],
        make_plots     = cfg.get("make_plots", True),
    )
    return net


def step_enrichment(cfg: dict, genes: list) -> dict:
    """Step 3: 농축분석"""
    log.info("=" * 60)
    log.info("STEP 3 — 농축분석 (GO / KEGG)")
    log.info("=" * 60)

    if not genes:
        log.warning("유전자 목록이 비어 있어 농축분석을 건너뜁니다.")
        return {}

    return run_enrichment(
        gene_symbols  = genes,
        prefix        = cfg["name"],
        adj_p_cutoff  = cfg["adj_p_cutoff"],
        top_n         = cfg["top_n_enrichment"],
        databases     = cfg.get("enrichment_dbs"),
        make_plots    = cfg["make_plots"],
    )


# ─────────────────────────────────────────────────────────────────────────────
# 보고서 생성
# ─────────────────────────────────────────────────────────────────────────────
def generate_report(
    cfg:          dict,
    herb_results: dict,
    drug_results: dict,
    all_genes:    list,
    net:          dict,
    enrich:       dict,
    elapsed:      float,
) -> Path:
    """분석 요약 보고서 (Markdown)"""
    now   = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    name  = cfg["name"]
    lines = [
        f"# Network Pharmacology Analysis Report",
        f"",
        f"**Analysis name**: `{name}`  ",
        f"**Generated**: {now}  ",
        f"**Elapsed**: {elapsed:.1f} seconds  ",
        f"",
        f"---",
        f"",
        f"## 1. Input",
        f"",
        f"| Category | Items |",
        f"|----------|-------|",
        f"| Herbs  | {', '.join(cfg['herbs']) if cfg['herbs'] else '—'} |",
        f"| Drugs  | {', '.join(cfg['drugs']) if cfg['drugs'] else '—'} |",
        f"",
        f"## 2. Data Collection",
        f"",
    ]

    # 한약 수집 요약
    if herb_results:
        lines += ["| Herb | Compounds | Targets | Source |",
                  "|------|-----------|---------|--------|"]
        for herb, herb_data in herb_results.items():
            comps  = herb_data.get("compounds", [])
            n_comp = len(comps)
            n_tgt  = sum(len(c.get("targets", [])) for c in comps)
            src    = herb_data.get("source", "")
            lines.append(f"| {herb} | {n_comp} | {n_tgt} | {src} |")
        lines.append("")

    lines += [
        f"**Total unique target genes**: {len(all_genes)}",
        f"",
        f"## 3. Network",
        f"",
    ]

    if net:
        hct_g = net["hct"].G
        ppi_g = net["ppi"].G
        int_g = net["integrated"].G
        lines += [
            f"| Network | Nodes | Edges |",
            f"|---------|-------|-------|",
            f"| HCT (Herb-Compound-Target) | {hct_g.number_of_nodes()} | {hct_g.number_of_edges()} |",
            f"| PPI (STRING)               | {ppi_g.number_of_nodes()} | {ppi_g.number_of_edges()} |",
            f"| Integrated                 | {int_g.number_of_nodes()} | {int_g.number_of_edges()} |",
            f"",
        ]

        # 허브 유전자
        hubs = net.get("hubs")
        if hubs is not None and not hubs.empty:
            lines += ["### Hub Targets (Top 10 by degree)", ""]
            lines += ["| Gene | Degree | Betweenness |",
                      "|------|--------|-------------|"]
            for _, row in hubs.head(10).iterrows():
                lines.append(
                    f"| {row['label']} | {row['degree']} | {row['betweenness']:.4f} |"
                )
            lines.append("")

        # 공유 타겟 (한약 간)
        herb_names = list(cfg["herbs"].keys()) if isinstance(cfg["herbs"], dict) else list(cfg["herbs"])
        if len(herb_names) >= 2:
            lines += ["### Shared Targets Between Herbs", ""]
            for i in range(len(herb_names)):
                for j in range(i + 1, len(herb_names)):
                    try:
                        shared = net["integrated"].get_shared_targets(herb_names[i], herb_names[j])
                    except Exception:
                        shared = []
                    lines.append(
                        f"- **{herb_names[i]}** ∩ **{herb_names[j]}**: "
                        f"{len(shared)} genes — {', '.join(shared[:10])}"
                        f"{'...' if len(shared) > 10 else ''}"
                    )
            lines.append("")

    lines += ["## 4. Enrichment Analysis", ""]

    if enrich and enrich.get("results"):
        lines += ["| Database | Significant terms (FDR≤{:.2f}) |".format(cfg["adj_p_cutoff"]),
                  "|----------|-------------------------------|"]
        for alias, df in enrich["results"].items():
            lines.append(f"| {alias} | {len(df)} |")
        lines.append("")

        # KEGG 상위 10
        kegg_df = enrich["results"].get("KEGG")
        if kegg_df is not None and not kegg_df.empty:
            lines += ["### Top 10 KEGG Pathways", ""]
            lines += ["| Pathway | Gene count | Adj. p-value |",
                      "|---------|-----------|-------------|"]
            for _, row in kegg_df.head(10).iterrows():
                term = row["term"][:60] if len(row["term"]) > 60 else row["term"]
                lines.append(f"| {term} | {row['gene_count']} | {row['adj_p_value']:.2e} |")
            lines.append("")

    # ── Publication Guide ─────────────────────────────────────────
    re_dir  = RESULTS_DIR
    en_dir  = BASE_DIR / "enrichment"
    net_dir = BASE_DIR / "networks"

    def _exists(folder, pattern):
        import glob as _gl
        return bool(_gl.glob(str(folder / pattern)))

    main_figs   = []
    supp_figs   = []
    main_tables = []
    supp_tables = []

    # Main Figures
    if _exists(re_dir, f"{name}_figure_multipanel*"):
        main_figs.append(f"**Fig 1** — Integrated multipanel figure: `results/{name}_figure_multipanel.png`")
    else:
        if _exists(re_dir, f"{name}_network*") or _exists(re_dir, f"{name}_*network*"):
            main_figs.append(f"**Fig 1** — PPI / Integrated network: `results/{name}_network.png`")
        if _exists(re_dir, f"{name}_hub_ranking*"):
            main_figs.append(f"**Fig 2** — Hub gene ranking: `results/{name}_hub_ranking.png`")
        if _exists(re_dir, f"{name}_venn*"):
            main_figs.append(f"**Fig 3** — Venn diagram (herb ∩ disease targets): `results/{name}_venn.png`")
        if _exists(en_dir, f"{name}_KEGG*dotplot*"):
            main_figs.append(f"**Fig 4** — KEGG pathway dotplot: `enrichment/{name}_KEGG_dotplot.png`")
        if _exists(re_dir, f"{name}_docking_heatmap*"):
            main_figs.append(f"**Fig 5** — Molecular docking heatmap: `results/{name}_docking_heatmap.png`")

    # Supplementary Figures
    if _exists(re_dir, f"{name}_admet_network*"):
        supp_figs.append(f"**Supp Fig S1** — ADMET-filtered network: `results/{name}_admet_network.png`")
    if _exists(re_dir, f"{name}_admet_hub_ranking*"):
        supp_figs.append(f"**Supp Fig S2** — ADMET-filtered hub ranking: `results/{name}_admet_hub_ranking.png`")
    for db in ["GO_BP", "GO_MF", "GO_CC", "Reactome"]:
        if _exists(en_dir, f"{name}_{db}*dotplot*"):
            supp_figs.append(f"**Supp Fig** — {db} dotplot: `enrichment/{name}_{db}_dotplot.png`")

    # Main Tables
    main_tables += [
        "**Table 1** — Herb-compound-target summary (Section 2 of this report)",
        "**Table 2** — Hub target genes Top 10~20 (Section 3, Hub Targets)",
        "**Table 3** — Top 10 KEGG pathways (Section 4, Enrichment)",
    ]

    # Supplementary Tables
    if _exists(re_dir, f"{name}_supplementary*"):
        supp_tables.append(f"**Supp Table S1** — Full compound/target/enrichment data: `results/{name}_supplementary.xlsx`")
    if _exists(net_dir, f"{name}_*nodes*"):
        supp_tables.append(f"**Supp Table S2** — Network node list: `networks/{name}_integrated_nodes.csv`")
    if _exists(net_dir, f"{name}_*edges*"):
        supp_tables.append(f"**Supp Table S3** — Network edge list: `networks/{name}_integrated_edges.csv`")

    lines += [
        "## 5. Output Files",
        "",
        f"- Networks: `networks/{name}_*.graphml` / `*_edges.csv` / `*_nodes.csv`",
        f"- Enrichment: `enrichment/{name}_*.csv` / `*_dotplot.png`",
        f"- This report: `results/{name}_report.md`",
        "",
        "---",
        "",
        "## 6. Publication Guide",
        "",
        "> 아래는 저널 투고 시 권장되는 Figure/Table 배치 가이드입니다.",
        "> (Frontiers in Pharmacology, J. Ethnopharmacology 등 기준)",
        "",
        "### Main Figures (본문)",
        "",
    ]
    if main_figs:
        for mf in main_figs:
            lines.append(f"- {mf}")
    else:
        lines.append("- 생성된 Main Figure 없음")

    lines += ["", "### Main Tables (본문)", ""]
    for mt in main_tables:
        lines.append(f"- {mt}")

    lines += ["", "### Supplementary Figures (부록)", ""]
    if supp_figs:
        for sf in supp_figs:
            lines.append(f"- {sf}")
    else:
        lines.append("- 생성된 Supplementary Figure 없음")

    lines += ["", "### Supplementary Tables (부록)", ""]
    if supp_tables:
        for st_ in supp_tables:
            lines.append(f"- {st_}")
    else:
        lines.append("- 생성된 Supplementary Table 없음")

    lines += [
        "",
        "### 참고용 (논문 미포함)",
        "",
        f"- `networks/{name}_*.graphml` — Cytoscape 시각화용",
        f"- `results/{name}_report.md` — 이 보고서 (내부 검토용)",
        f"- `results/{name}_ai_draft.md` / `_manuscript.docx` — AI 초안 (수정 후 사용)",
        "",
        "---",
        "_Generated by netpharm pipeline_",
    ]

    path = RESULTS_DIR / f"{name}_report.md"
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    log.info(f"[보고서] → {path}")
    return path


# ─────────────────────────────────────────────────────────────────────────────
# 메인 파이프라인
# ─────────────────────────────────────────────────────────────────────────────
def run_pipeline(cfg: dict):
    t_start = time.time()
    log.info(f"파이프라인 시작: '{cfg['name']}'")
    log.info(f"한약: {cfg['herbs']}")
    log.info(f"양약: {cfg['drugs']}")
    try:
        from telegram_notify import send as tg_send
        tg_send(f"🔬 <b>[{cfg['name']}]</b> 분석 시작\n• 한약: {len(cfg.get('herbs',{}))}종 / 약물: {', '.join(cfg.get('drugs',[]))}")
    except Exception:
        pass

    collector = NetPharmCollector()

    try:
        # Step 1: 수집
        herb_results, drug_results, all_genes = step_collect(cfg, collector)

        # Step 2: PPI 수집 (타겟 유전자 기반)
        ppi_edges = []
        if all_genes:
            log.info(f"STRING PPI 수집 중 ({len(all_genes)}개 유전자)...")
            ppi_edges = collector.collect_ppi(all_genes, min_score=cfg["ppi_min_score"])

        # Step 3: 네트워크
        net = {}
        if herb_results:
            net = step_network(cfg, herb_results, ppi_edges)
        else:
            log.warning("한약 데이터 없음 — 네트워크 구축 건너뜀")

        # Step 4: 농축분석
        enrich = step_enrichment(cfg, all_genes)

        # Step 5: 보고서
        elapsed = time.time() - t_start
        report_path = generate_report(
            cfg, herb_results, drug_results, all_genes,
            net, enrich, elapsed
        )

        log.info("=" * 60)
        log.info(f"파이프라인 완료 — {elapsed:.1f}초")
        log.info(f"보고서: {report_path}")
        log.info("=" * 60)

        herb_gene_count = len(set(t.get("gene_symbol") for r in herb_results.values()
                                  for c in r.get("compounds", [])
                                  for t in c.get("targets", [])
                                  if t.get("gene_symbol")))
        drug_gene_count = len(set(t.get("gene_symbol","")
                                  for r in drug_results.values()
                                  for t in r.get("targets", [])
                                  if t.get("gene_symbol")))
        try:
            from telegram_notify import notify_analysis_done
            notify_analysis_done(cfg["name"], {
                "herb_genes":  herb_gene_count,
                "drug_genes":  drug_gene_count,
                "intersection": len(set()),  # net에서 추출 가능하나 간단히 0
                "hub_genes":   list(net.get("hub_genes", {}).keys())[:5] if net else [],
            })
        except Exception:
            pass

        return {
            "herb_results": herb_results,
            "drug_results": drug_results,
            "genes":        all_genes,
            "network":      net,
            "enrichment":   enrich,
            "report":       report_path,
        }

    finally:
        collector.close()


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Network Pharmacology Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python src/pipeline.py --herbs "Panax ginseng" "Scutellaria baicalensis" --name test
  python src/pipeline.py --config configs/my_analysis.json
  python src/pipeline.py --template example   # config 템플릿 생성
        """,
    )
    p.add_argument("--config",   type=str, help="JSON config 파일 경로")
    p.add_argument("--name",     type=str, default="analysis", help="분석 이름 (출력 파일 접두사)")
    p.add_argument("--herbs",    nargs="+", default=[], help="한약명 목록 (학명/영문명)")
    p.add_argument("--drugs",    nargs="+", default=[], help="양약명 목록")
    p.add_argument("--ppi-score",type=int,  default=400, help="STRING 최소 신뢰도 점수 (default: 400)")
    p.add_argument("--fdr",      type=float,default=0.05,help="농축분석 FDR 컷오프 (default: 0.05)")
    p.add_argument("--no-plots", action="store_true",    help="시각화 비활성화")
    p.add_argument("--template", type=str, metavar="NAME", help="config 템플릿 파일 생성 후 종료")
    return p


def main():
    parser = build_parser()
    args   = parser.parse_args()

    # 템플릿 생성만 하고 종료
    if args.template:
        save_config_template(args.template)
        sys.exit(0)

    # config 파일 우선, 없으면 CLI 인자 사용
    if args.config:
        cfg = load_config(args.config)
    else:
        if not args.herbs and not args.drugs:
            parser.print_help()
            print("\n오류: --herbs 또는 --drugs 중 하나 이상 필요합니다.")
            sys.exit(1)
        cfg = {
            **DEFAULT_CONFIG,
            "name":          args.name,
            "herbs":         args.herbs,
            "drugs":         args.drugs,
            "ppi_min_score": args.ppi_score,
            "adj_p_cutoff":  args.fdr,
            "make_plots":    not args.no_plots,
        }

    run_pipeline(cfg)


if __name__ == "__main__":
    main()
