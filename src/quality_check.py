"""
분석 결과 학술 품질 진단 모듈
- 항목별 타당성 점검 + 개선 제안 생성
- 결과: Markdown 보고서 파일 저장
"""

from __future__ import annotations
import logging
import math
from datetime import datetime
from pathlib import Path

log = logging.getLogger(__name__)

BASE_DIR = Path(__file__).parent.parent


# ─────────────────────────────────────────────────────────────────────────────
# 개별 진단 함수
# ─────────────────────────────────────────────────────────────────────────────

def _check_target_coverage(pipeline: dict) -> dict:
    """성분당 타깃 수, ChEMBL 미발견률 진단"""
    herb_results = pipeline.get("herb_results", {})
    total_compounds = 0
    zero_target     = 0
    low_target      = 0   # 1~2개
    supplement_used = 0
    per_compound    = []

    for hd in herb_results.values():
        for c in hd.get("compounds", []):
            total_compounds += 1
            tgts = c.get("targets", [])
            n    = len(tgts)
            name = (c.get("molecule") or {}).get("pref_name", "?")
            sources = {t.get("source", "ChEMBL") for t in tgts}
            used_sup = any(s in ("PubChem", "ChEMBL_similarity") for s in sources)
            if used_sup:
                supplement_used += 1
            if n == 0:
                zero_target += 1
            elif n <= 2:
                low_target += 1
            per_compound.append({"name": name, "n_targets": n, "supplement": used_sup})

    if total_compounds == 0:
        return {"status": "skip", "title": "타깃 수집 품질", "issues": [], "suggestions": []}

    zero_pct = zero_target / total_compounds * 100
    low_pct  = (zero_target + low_target) / total_compounds * 100
    avg_n    = sum(p["n_targets"] for p in per_compound) / total_compounds

    issues = []
    suggestions = []
    status = "ok"

    if zero_pct > 30:
        status = "warn"
        issues.append(f"성분의 {zero_pct:.0f}%가 타깃 0개 — 분석 기여도 없음")
        suggestions.append("ChEMBL에 없는 성분은 대체 성분명(동의어)으로 재검증하거나 제거 검토")
    elif zero_pct > 10:
        status = "caution"
        issues.append(f"성분의 {zero_pct:.0f}%가 타깃 0개")
        suggestions.append("DB 검증 패널에서 ❌ 성분 확인 후 대체 후보 적용 권장")

    if avg_n < 5:
        status = max(status, "caution", key=lambda s: {"ok":0,"caution":1,"warn":2}.get(s,0))
        issues.append(f"성분당 평균 타깃 수 {avg_n:.1f}개 — 낮은 편 (권장: ≥10)")
        suggestions.append("PubChem BioAssay 보완이 자동 적용됨; 성분 수 확충 또는 유사도 컷오프 완화(60→50) 고려")

    if supplement_used > 0:
        issues.append(f"{supplement_used}개 성분이 PubChem/유사도 보완 사용 — 논문 Methods에 명시 필요")
        suggestions.append('Methods에 "Targets with no ChEMBL annotation were supplemented via PubChem BioAssay and ChEMBL structural similarity (Tanimoto ≥ 0.65)" 추가')

    zero_list = [p["name"] for p in per_compound if p["n_targets"] == 0]
    if zero_list:
        issues.append(f"타깃 0개 성분: {', '.join(zero_list[:6])}{'...' if len(zero_list)>6 else ''}")

    return {
        "status": status,
        "title": "타깃 수집 품질",
        "stats": {"total": total_compounds, "zero": zero_target,
                  "avg_targets": round(avg_n, 1), "supplement": supplement_used},
        "issues": issues,
        "suggestions": suggestions,
    }


def _check_network(pipeline: dict, admet_filtered: dict | None) -> dict:
    """네트워크 품질 진단 (밀도, 고립 노드, PPI 엣지 수)"""
    active = admet_filtered or pipeline
    net    = active.get("network", {})

    issues      = []
    suggestions = []
    status      = "ok"
    stats       = {}

    try:
        import networkx as nx
        ppi_net = net.get("ppi")
        G = getattr(ppi_net, "G", ppi_net) if ppi_net else None
        if G is None:
            return {"status": "warn", "title": "네트워크 품질",
                    "issues": ["PPI 네트워크 없음"], "suggestions": ["분석 재실행 필요"]}

        n_nodes   = G.number_of_nodes()
        n_edges   = G.number_of_edges()
        density   = nx.density(G)
        n_isolate = len(list(nx.isolates(G)))
        components = list(nx.connected_components(G))
        lcc_frac  = max(len(c) for c in components) / n_nodes if components else 0

        stats = {"nodes": n_nodes, "edges": n_edges,
                 "density": round(density, 4), "isolated": n_isolate,
                 "lcc_fraction": round(lcc_frac, 2)}

        if n_nodes < 20:
            status = "warn"
            issues.append(f"PPI 네트워크 노드 {n_nodes}개 — 너무 작음 (권장: ≥30)")
            suggestions.append("STRING 신뢰도 임계값을 400→300으로 낮추거나 성분/타깃 추가 검토")
        elif n_nodes < 50:
            status = "caution"
            issues.append(f"PPI 네트워크 노드 {n_nodes}개 — 소규모 (권장: ≥50)")

        if n_isolate > 0:
            status = max(status, "caution", key=lambda s: {"ok":0,"caution":1,"warn":2}.get(s,0))
            issues.append(f"고립 노드 {n_isolate}개 — 네트워크와 연결 없음")
            suggestions.append("고립 노드는 실제 상호작용이 없거나 데이터 부재. 분석에서 제외하거나 STRING 점수 완화 고려")

        if lcc_frac < 0.6:
            status = max(status, "caution", key=lambda s: {"ok":0,"caution":1,"warn":2}.get(s,0))
            issues.append(f"최대 연결 컴포넌트가 전체의 {lcc_frac*100:.0f}% — 분절된 네트워크")
            suggestions.append("PPI 신뢰도 완화 또는 타깃 추가로 네트워크 연결성 개선 필요")

        if density < 0.01:
            issues.append(f"네트워크 밀도 {density:.4f} — 희박한 네트워크")
            suggestions.append("희박한 PPI는 STRING 데이터 한계일 수 있음. 논문에서 network sparsity 언급 권장")

    except Exception as e:
        issues.append(f"네트워크 분석 오류: {e}")
        status = "warn"

    return {"status": status, "title": "네트워크 품질",
            "stats": stats, "issues": issues, "suggestions": suggestions}


def _check_disease_intersection(intersection: dict | None, pipeline: dict,
                                admet_filtered: dict | None) -> dict:
    """질환 교차 유의성 — hypergeometric p-value 계산"""
    issues      = []
    suggestions = []
    status      = "ok"
    stats       = {}

    if not intersection:
        return {"status": "skip", "title": "질환 교차분석",
                "issues": ["질환 교차분석 미실행"], "suggestions": ["disease_query 설정 후 재분석"]}

    shared  = intersection.get("shared_count", len(intersection.get("shared", [])))
    h_total = intersection.get("herb_total", 0)
    d_total = intersection.get("disease_total", 0)

    # 전체 인간 단백질 ~ 20,000개 기준으로 hypergeometric test
    N = 20000
    try:
        from scipy.stats import hypergeom
        p_val = hypergeom.sf(shared - 1, N, h_total, d_total)
    except Exception:
        p_val = None

    stats = {"shared": shared, "herb_targets": h_total,
             "disease_targets": d_total, "hypergeom_p": p_val}

    if shared == 0:
        status = "warn"
        issues.append("공유 타깃 0개 — 한약/질환 간 연관성 데이터 없음")
        suggestions.append("disease_query 재설정(영문 질환명 정확히) 또는 OpenTargets 점수 임계값 완화")
    elif shared < 5:
        status = "caution"
        issues.append(f"공유 타깃 {shared}개 — 적은 편 (권장: ≥5)")
        suggestions.append("논문에서 공유 타깃의 생물학적 의의를 개별적으로 상세 기술 필요")

    if p_val is not None:
        if p_val > 0.05:
            # 공유 타깃이 있으면 caution, 없으면 이미 warn 처리됨
            status = max(status, "caution", key=lambda s: {"ok":0,"caution":1,"warn":2}.get(s,0))
            issues.append(f"교차 유의성 p = {p_val:.3f} > 0.05 — 통계적 유의성 미달")
            suggestions.append("질환 키워드를 더 구체적인 영문 질환명으로 변경하거나 타깃 수 확충 검토")
        else:
            issues.append(f"교차 유의성 p = {p_val:.2e} (hypergeometric, N=20,000) — 유의함")
            suggestions.append(f"논문 Methods에 hypergeometric test 결과 (p = {p_val:.2e}) 명시 권장")

    return {"status": status, "title": "질환 교차분석",
            "stats": stats, "issues": issues, "suggestions": suggestions}


def _check_admet(admet: dict | None) -> dict:
    """ADMET 필터링 비율 진단"""
    issues      = []
    suggestions = []
    status      = "ok"

    if not admet:
        return {"status": "skip", "title": "ADMET 스크리닝",
                "issues": ["ADMET 미실행"], "suggestions": []}

    passed = len(admet.get("passed", []))
    failed = len(admet.get("failed", []))
    total  = passed + failed
    if total == 0:
        return {"status": "skip", "title": "ADMET 스크리닝", "issues": [], "suggestions": []}

    pass_rate = passed / total * 100

    if pass_rate == 0:
        status = "warn"
        issues.append("ADMET 통과 성분 0개 — 모든 성분이 탈락")
        suggestions.append("ADMET 기준(Lipinski RO5 등) 완화 또는 성분 재검토 필요")
    elif pass_rate < 30:
        status = "caution"
        issues.append(f"ADMET 통과율 {pass_rate:.0f}% ({passed}/{total}) — 낮은 편")
        suggestions.append("통과율이 낮으면 분석 기반 성분이 제한됨. 기준 완화 고려")
    elif pass_rate > 95:
        issues.append(f"ADMET 통과율 {pass_rate:.0f}% — 필터링 효과 미미")
        suggestions.append("ADMET 기준이 너무 완화되어 있을 수 있음. 기준 재검토 권장")

    if passed > 0 and failed > 0:
        suggestions.append(
            f"탈락 성분 {failed}개는 Methods에서 exclusion criteria로 기술 필요"
        )

    return {"status": status, "title": "ADMET 스크리닝",
            "stats": {"passed": passed, "failed": failed, "rate": round(pass_rate, 1)},
            "issues": issues, "suggestions": suggestions}


def _check_enrichment(pipeline: dict, admet_filtered: dict | None) -> dict:
    """농축분석 품질 — term 수, 경로 관련성"""
    active = admet_filtered or pipeline
    enr    = active.get("enrichment", {})
    issues      = []
    suggestions = []
    status      = "ok"
    stats       = {}

    results = (enr or {}).get("results", {})
    kegg_df = results.get("KEGG")
    gobp_df = results.get("GO_BP")

    try:
        import pandas as pd
        n_kegg = len(kegg_df) if kegg_df is not None else 0
        n_gobp = len(gobp_df) if gobp_df is not None else 0
        stats  = {"kegg_sig_terms": n_kegg, "gobp_sig_terms": n_gobp}

        if n_kegg == 0 and n_gobp == 0:
            status = "warn"
            issues.append("유의한 농축 term 없음 (FDR < 0.05)")
            suggestions.append("FDR 컷오프를 0.05→0.1로 완화하거나 타깃 수를 늘릴 것")
        elif n_kegg < 5:
            status = "caution"
            issues.append(f"KEGG 유의 pathway {n_kegg}개 — 적은 편")
            suggestions.append("타깃 수 부족 또는 경로 데이터베이스 다양화 (Reactome 포함 확인)")

        if kegg_df is not None and n_kegg > 0:
            # 상위 5개 pathway 기록
            top = kegg_df.head(5)
            top_names = top["term"].tolist() if "term" in top.columns else []
            if top_names:
                issues.append(f"상위 KEGG pathway: {', '.join(top_names[:3])}")
                # p값 확인
                if "adj_p_value" in top.columns:
                    min_p = top["adj_p_value"].min()
                    if min_p > 0.05:
                        status = max(status, "caution",
                                     key=lambda s: {"ok":0,"caution":1,"warn":2}.get(s,0))
                        issues.append(f"최소 adj_p = {min_p:.3f} > 0.05 — 유의 경계")

    except Exception as e:
        issues.append(f"농축분석 확인 오류: {e}")

    return {"status": status, "title": "농축분석 (Enrichment)",
            "stats": stats, "issues": issues, "suggestions": suggestions}


def _check_docking(docking) -> dict:
    """도킹 결합 에너지 분포 진단"""
    issues      = []
    suggestions = []
    status      = "ok"

    if docking is None:
        return {"status": "skip", "title": "분자 도킹",
                "issues": ["도킹 미실행"], "suggestions": []}

    try:
        import pandas as pd
        df = docking if isinstance(docking, pd.DataFrame) else pd.DataFrame(docking)
        if df.empty:
            return {"status": "warn", "title": "분자 도킹",
                    "issues": ["도킹 결과 없음 (0쌍 성공)"],
                    "suggestions": ["PDB 구조 다운로드 실패 또는 Vina 오류 확인"]}

        affinities = df["affinity"].dropna()
        n_total  = len(df)
        n_strong = int((affinities < -7.0).sum())   # 강한 결합
        n_mod    = int(((affinities >= -7.0) & (affinities < -5.0)).sum())
        n_weak   = int((affinities >= -5.0).sum())
        best     = float(affinities.min()) if len(affinities) > 0 else 0
        stats    = {"pairs": n_total, "strong(<-7)": n_strong,
                    "moderate(-7~-5)": n_mod, "weak(>-5)": n_weak,
                    "best_affinity": round(best, 2)}

        if n_strong == 0 and n_mod == 0:
            status = "warn"
            issues.append(f"모든 도킹 결과 > -5 kcal/mol — 의미 있는 결합 없음")
            suggestions.append(
                "결합부위(binding pocket) 정의 문제일 수 있음. 공결정 리간드 유무 확인 또는 "
                "exhaustiveness를 8→16으로 높여 재도킹 필요"
            )
        elif n_strong == 0:
            status = "caution"
            issues.append(f"강한 결합(< -7 kcal/mol) 없음. 중간 결합 {n_mod}쌍")
            suggestions.append(
                "중간 결합도 논문에 포함 가능하나 in vitro 실험적 검증 언급 필수"
            )
        else:
            issues.append(f"강한 결합 {n_strong}쌍 (< -7 kcal/mol), 최고 {best} kcal/mol")
            suggestions.append(
                "도킹 결과 신뢰성 향상을 위해 exhaustiveness=16 재도킹 및 "
                "결합 포즈 시각화(PyMOL/UCSF Chimera) 권장"
            )

        issues.append(
            "현재 도킹은 rigid receptor — 수용체 유연성(induced fit) 미고려. "
            "상위 후보는 MD simulation 또는 flexible docking 추가 검증 권장"
        )

    except Exception as e:
        issues.append(f"도킹 결과 확인 오류: {e}")
        status = "warn"

    return {"status": status, "title": "분자 도킹",
            "stats": stats, "issues": issues, "suggestions": suggestions}


def _check_methodology(pipeline: dict, admet: dict | None,
                        intersection: dict | None, docking) -> dict:
    """방법론 완결성 체크리스트"""
    issues      = []
    suggestions = []

    checklist = {
        "타깃 수집 (ChEMBL)":         bool(pipeline.get("herb_results")),
        "PPI 네트워크 구축 (STRING)":  bool(pipeline.get("network")),
        "농축분석 (GO/KEGG)":         bool(pipeline.get("enrichment")),
        "ADMET 스크리닝":             bool(admet),
        "질환 교차분석 (OpenTargets)": bool(intersection),
        "허브 토폴로지 분석":          bool(pipeline.get("network")),
        "분자 도킹 (Vina)":           docking is not None,
    }

    missing = [k for k, v in checklist.items() if not v]
    done    = [k for k, v in checklist.items() if v]

    if missing:
        issues.append(f"미완료 항목: {', '.join(missing)}")
        suggestions.append("미완료 항목은 논문 게재 전 보완 필요 (Reviewer 지적 가능성 높음)")

    # 표준 NP 논문 필수 항목 체크
    required = ["타깃 수집 (ChEMBL)", "PPI 네트워크 구축 (STRING)", "농축분석 (GO/KEGG)"]
    missing_req = [r for r in required if r in missing]
    if missing_req:
        suggestions.append(
            f"필수 항목 미완료: {', '.join(missing_req)} — 네트워크 약리학 논문 기본 요건"
        )
    else:
        issues.append("필수 3개 항목(타깃수집, PPI, 농축분석) 완료")

    suggestions.append(
        "논문 투고 전 체크: (1) 모든 DB 버전 명시 (ChEMBL v33, STRING v12 등) "
        "(2) 분석 재현을 위한 파라미터 기술 (3) 결과 독립적 검증 실험 언급"
    )

    status = "warn" if missing_req else ("caution" if missing else "ok")
    return {"status": status, "title": "방법론 완결성",
            "checklist": checklist, "issues": issues, "suggestions": suggestions}


# ─────────────────────────────────────────────────────────────────────────────
# 종합 보고서 생성
# ─────────────────────────────────────────────────────────────────────────────

_STATUS_ICON = {"ok": "✅", "caution": "⚠️", "warn": "❌", "skip": "⬜"}
_STATUS_LABEL = {"ok": "양호", "caution": "주의", "warn": "개선 필요", "skip": "미실행"}


def run_quality_check(
    pipeline:        dict,
    admet:           dict | None = None,
    admet_filtered:  dict | None = None,
    intersection:    dict | None = None,
    docking          = None,
    cfg:             dict | None = None,
    output_dir:      Path | None = None,
) -> tuple[dict, Path]:
    """
    전체 품질 진단 실행 → dict 결과 + Markdown 파일 반환.

    반환: (report_dict, md_path)
    """
    cfg        = cfg or {}
    output_dir = output_dir or (BASE_DIR / "results")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    checks = [
        _check_target_coverage(pipeline),
        _check_network(pipeline, admet_filtered),
        _check_disease_intersection(intersection, pipeline, admet_filtered),
        _check_admet(admet),
        _check_enrichment(pipeline, admet_filtered),
        _check_docking(docking),
        _check_methodology(pipeline, admet, intersection, docking),
    ]

    # 전체 점수 (ok=2, caution=1, warn=0, skip=-)
    score_map  = {"ok": 2, "caution": 1, "warn": 0, "skip": None}
    scoreable  = [c for c in checks if score_map[c["status"]] is not None]
    total_score = sum(score_map[c["status"]] for c in scoreable)
    max_score   = len(scoreable) * 2
    pct         = total_score / max_score * 100 if max_score else 0
    stars       = "⭐" * round(pct / 20)  # 0~5점

    report = {
        "name":       cfg.get("name", "analysis"),
        "generated":  datetime.now().strftime("%Y-%m-%d %H:%M"),
        "score_pct":  round(pct),
        "stars":      stars,
        "checks":     checks,
    }

    # ── Markdown 생성 ────────────────────────────────────────────────────
    md = _render_markdown(report, cfg)
    name = cfg.get("name", "analysis")
    md_path = output_dir / f"{name}_quality_report.md"
    md_path.write_text(md, encoding="utf-8")
    log.info(f"[QC] 품질 보고서 저장: {md_path}")

    return report, md_path


def _render_markdown(report: dict, cfg: dict) -> str:
    now   = report["generated"]
    name  = report["name"]
    herbs = list((cfg.get("herbs") or {}).keys())
    drugs = cfg.get("drugs", [])
    disease = cfg.get("disease_query", "")
    pct   = report["score_pct"]
    stars = report["stars"]

    lines = [
        f"# NetPharm 분석 품질 진단 보고서",
        f"",
        f"**분석명:** `{name}`  ",
        f"**연구 대상:** {', '.join(herbs)} + {', '.join(drugs)}"
        + (f" | 질환: {disease}" if disease else ""),
        f"**생성일시:** {now}",
        f"",
        f"---",
        f"",
        f"## 종합 판정",
        f"",
        f"| 항목 | 내용 |",
        f"|------|------|",
        f"| 종합 점수 | **{pct}점 / 100점** {stars} |",
    ]

    # 각 항목 상태 요약
    for c in report["checks"]:
        icon  = _STATUS_ICON[c["status"]]
        label = _STATUS_LABEL[c["status"]]
        lines.append(f"| {c['title']} | {icon} {label} |")

    lines += ["", "---", ""]

    # 항목별 상세
    for c in report["checks"]:
        icon   = _STATUS_ICON[c["status"]]
        label  = _STATUS_LABEL[c["status"]]
        lines += [f"## {icon} {c['title']} — {label}", ""]

        # stats
        stats = c.get("stats") or c.get("checklist")
        if stats and isinstance(stats, dict):
            lines.append("**수치 요약**")
            lines.append("")
            lines.append("| 항목 | 값 |")
            lines.append("|------|----|")
            for k, v in stats.items():
                if isinstance(v, bool):
                    v = "✅" if v else "❌"
                elif isinstance(v, float) and v != v:
                    v = "—"
                elif v is None:
                    v = "—"
                lines.append(f"| {k} | {v} |")
            lines.append("")

        if c.get("issues"):
            lines.append("**진단 내용**")
            lines.append("")
            for issue in c["issues"]:
                lines.append(f"- {issue}")
            lines.append("")

        if c.get("suggestions"):
            lines.append("**개선 제안**")
            lines.append("")
            for sug in c["suggestions"]:
                lines.append(f"> {sug}")
            lines.append("")

        lines.append("---")
        lines.append("")

    # 우선순위 action items
    warn_checks = [c for c in report["checks"] if c["status"] == "warn"]
    caution_checks = [c for c in report["checks"] if c["status"] == "caution"]

    lines += ["## 우선 조치 항목", ""]
    if warn_checks:
        lines.append("### 즉시 수정 필요 (❌)")
        for c in warn_checks:
            for sug in c.get("suggestions", [])[:2]:
                lines.append(f"- **[{c['title']}]** {sug}")
        lines.append("")
    if caution_checks:
        lines.append("### 투고 전 검토 권장 (⚠️)")
        for c in caution_checks:
            for sug in c.get("suggestions", [])[:1]:
                lines.append(f"- **[{c['title']}]** {sug}")
        lines.append("")

    lines += [
        "---",
        "",
        "*이 보고서는 NetPharm 자동 생성 — 최종 학술 판단은 연구자가 수행해야 합니다.*",
        f"*생성: {now} | 모델: NetPharm QC v1.0*",
    ]

    return "\n".join(lines)
