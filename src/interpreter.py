"""
Network Pharmacology — AI Interpreter
Gemini API를 이용한 분석 결과 해석 및 논문 초안 생성

사용 모델: gemini-2.0-flash (무료 티어)
  - 15 RPM, 1M tokens/day
"""

import os
import time
import logging
from pathlib import Path

log = logging.getLogger(__name__)

try:
    from google import genai
    from google.genai import types
    from google.genai.errors import ServerError, ClientError
    GENAI_AVAILABLE = True
except ImportError:
    GENAI_AVAILABLE = False
    ServerError = ClientError = Exception
    log.warning("google-genai 미설치. pip install google-genai 실행 필요.")


MODEL = "gemini-3.5-flash"

# 재시도 설정: 1차 30초, 2차 60초, 3차 120초
_RETRY_WAITS = [30, 60, 120]


def _get_client(api_key: str = None):
    key = api_key or os.environ.get("GEMINI_API_KEY", "")
    if not key:
        raise ValueError("GEMINI_API_KEY가 설정되지 않았습니다.")
    return genai.Client(api_key=key)


def _is_retryable(e: Exception) -> bool:
    """429(할당량 초과) / 503(서비스 불가) 일시적 오류 여부"""
    msg = str(e)
    return isinstance(e, (ServerError, ClientError)) and (
        "503" in msg or "429" in msg or "UNAVAILABLE" in msg or "ResourceExhausted" in msg
    )


def _call(client, prompt: str, temperature: float = 0.3,
          on_retry=None) -> str:
    """
    Gemini API 호출. 429/503 일시적 오류 시 최대 3회 자동 재시도.
    on_retry(attempt, wait_sec, error_msg): 재시도 직전 호출되는 콜백 (프론트 상태 표시용)
    """
    last_exc = None
    for attempt in range(1, len(_RETRY_WAITS) + 2):
        try:
            response = client.models.generate_content(
                model=MODEL,
                contents=prompt,
                config=types.GenerateContentConfig(temperature=temperature),
            )
            return response.text.strip()
        except Exception as e:
            last_exc = e
            if attempt > len(_RETRY_WAITS) or not _is_retryable(e):
                raise
            wait = _RETRY_WAITS[attempt - 1]
            log.warning(
                f"[Gemini] 일시적 오류 ({attempt}/3회) — {wait}초 후 재시도: {e}"
            )
            if on_retry:
                on_retry(attempt, wait, str(e))
            time.sleep(wait)
    raise last_exc


# ─────────────────────────────────────────────────────────────────────────────
# 1. 결과 해석
# ─────────────────────────────────────────────────────────────────────────────
def interpret_network(
    herbs:       list,
    drugs:       list,
    hub_genes:   list,
    kegg_terms:  list,
    go_bp_terms: list,
    api_key:     str = None,
    on_retry=None,
) -> str:
    """
    허브 유전자 + 농축분석 결과 → 작용기전 해석

    hub_genes : [{"label": "TP53", "degree": 15, "betweenness": 0.28}, ...]
    kegg_terms: [{"term": "Serotonergic synapse", "gene_count": 4, "adj_p_value": 0.001}, ...]
    go_bp_terms: 동일 형식
    """
    if not GENAI_AVAILABLE:
        return "[오류] google-genai 패키지가 설치되지 않았습니다."

    hub_text  = "\n".join(
        f"  - {g['label']} (degree={g['degree']}, betweenness={g.get('betweenness', 0):.3f})"
        for g in hub_genes[:10]
    )
    kegg_text = "\n".join(
        f"  - {t['term']} (genes={t['gene_count']}, adj.p={t['adj_p_value']:.2e})"
        for t in kegg_terms[:10]
    )
    gobp_text = "\n".join(
        f"  - {t['term']} (genes={t['gene_count']}, adj.p={t['adj_p_value']:.2e})"
        for t in go_bp_terms[:8]
    )

    herb_str = ", ".join(herbs)
    drug_str = ", ".join(drugs) if drugs else "없음"

    prompt = f"""You are an expert in network pharmacology and traditional Chinese medicine research.

Analyze the following network pharmacology results and provide a mechanistic interpretation in Korean.

## Analysis Subject
- Herbs: {herb_str}
- Drugs (co-medication): {drug_str}

## Hub Target Genes (by network degree)
{hub_text}

## Enriched KEGG Pathways (FDR≤0.05)
{kegg_text}

## Enriched GO Biological Processes (FDR≤0.05, top 8)
{gobp_text}

## Task
Write a concise mechanistic interpretation (300-400 words in Korean) covering:
1. 핵심 타겟 유전자의 생물학적 의의
2. 주요 KEGG pathway와의 연관성 및 작용기전 추론
3. 한약-양약 상호작용 가능성 (해당 시)
4. 이 결과의 임상적·약리학적 함의

Use formal academic Korean suitable for a scientific journal."""

    try:
        client = _get_client(api_key)
        return _call(client, prompt, temperature=0.3, on_retry=on_retry)
    except Exception as e:
        log.error(f"[Gemini] 해석 오류: {e}")
        return f"[오류] {e}"


# ─────────────────────────────────────────────────────────────────────────────
# 2. 논문 섹션 초안
# ─────────────────────────────────────────────────────────────────────────────
def draft_methods(
    herbs:          list,
    drugs:          list,
    compounds:      dict,
    ppi_score:      int  = 400,
    adj_p_cutoff:   float = 0.05,
    api_key:        str  = None,
    on_retry=None,
) -> str:
    """
    Methods 섹션 초안 (영문)
    compounds: {herb_name: [compound1, compound2, ...]}
    """
    herb_comp_text = "\n".join(
        f"  - {h}: {', '.join(c)}" for h, c in compounds.items()
    )
    drug_str = ", ".join(drugs) if drugs else "none"

    prompt = f"""Write a Methods section for a network pharmacology study in formal English, suitable for submission to an international journal (e.g., Frontiers in Pharmacology, Journal of Ethnopharmacology).

## Study Parameters
- Herbs analyzed: {', '.join(herbs)}
- Active compounds queried (via ChEMBL database):
{herb_comp_text}
- Western drugs: {drug_str}
- PPI network: STRING database (confidence score ≥ {ppi_score})
- Enrichment analysis: Enrichr API (GO, KEGG; FDR ≤ {adj_p_cutoff})
- Network construction: NetworkX (Python); exported to Cytoscape

## Required Subsections
1. Data Collection and Target Identification
2. Network Construction
3. Protein-Protein Interaction (PPI) Network Analysis
4. Gene Ontology and KEGG Pathway Enrichment Analysis

Write concisely (~400 words total). Use past tense, passive voice where appropriate. Include citations as [Author, Year] placeholders."""

    try:
        client = _get_client(api_key)
        return _call(client, prompt, temperature=0.2, on_retry=on_retry)
    except Exception as e:
        log.error(f"[Gemini] Methods 초안 오류: {e}")
        return f"[오류] {e}"


def draft_results(
    herbs:       list,
    drugs:       list,
    network_stats: dict,
    hub_genes:   list,
    kegg_terms:  list,
    api_key:     str = None,
    on_retry=None,
) -> str:
    """
    Results 섹션 초안 (영문)
    network_stats: {hct_nodes, hct_edges, ppi_nodes, ppi_edges, total_genes}
    """
    hub_text = ", ".join(
        f"{g['label']} (degree={g['degree']})" for g in hub_genes[:5]
    )
    kegg_text = "\n".join(
        f"  - {t['term']} (n={t['gene_count']}, adj.p={t['adj_p_value']:.2e})"
        for t in kegg_terms[:5]
    )

    prompt = f"""Write a Results section for a network pharmacology study in formal English.

## Data to Report
Herbs: {', '.join(herbs)}
Drugs: {', '.join(drugs) if drugs else 'none'}

Network statistics:
- HCT network: {network_stats.get('hct_nodes', '?')} nodes, {network_stats.get('hct_edges', '?')} edges
- PPI network: {network_stats.get('ppi_nodes', '?')} nodes, {network_stats.get('ppi_edges', '?')} edges
- Total unique target genes: {network_stats.get('total_genes', '?')}

Top hub targets: {hub_text}

Top KEGG pathways:
{kegg_text}

## Task
Write a concise Results section (~300 words) reporting these findings objectively.
- Describe the network topology
- Report hub gene identification
- Summarize enrichment results
- No interpretation (save for Discussion)
Use past tense, formal academic English."""

    try:
        client = _get_client(api_key)
        return _call(client, prompt, temperature=0.2, on_retry=on_retry)
    except Exception as e:
        log.error(f"[Gemini] Results 초안 오류: {e}")
        return f"[오류] {e}"


def draft_discussion(
    herbs:       list,
    drugs:       list,
    hub_genes:   list,
    kegg_terms:  list,
    interpretation: str,
    api_key:     str = None,
    on_retry=None,
) -> str:
    """Discussion 섹션 초안 (영문)"""
    hub_text  = ", ".join(g['label'] for g in hub_genes[:5])
    kegg_text = ", ".join(t['term'] for t in kegg_terms[:5])

    prompt = f"""Write a Discussion section for a network pharmacology study in formal English.

## Context
- Herbs: {', '.join(herbs)}
- Key hub genes identified: {hub_text}
- Major KEGG pathways: {kegg_text}
- Korean mechanistic interpretation (for reference):
{interpretation[:600]}

## Task
Write a Discussion (~400 words) that:
1. Interprets the biological significance of hub genes in relation to the herbs
2. Discusses the identified pathways in the context of the herbs' known pharmacological activities
3. Compares with previously published network pharmacology studies on similar herbs
4. Discusses limitations of this computational approach
5. Suggests future experimental validation directions

Use formal academic English. Include [citation needed] placeholders where literature support is expected."""

    try:
        client = _get_client(api_key)
        return _call(client, prompt, temperature=0.4, on_retry=on_retry)
    except Exception as e:
        log.error(f"[Gemini] Discussion 초안 오류: {e}")
        return f"[오류] {e}"


# ─────────────────────────────────────────────────────────────────────────────
# 2-extra. Abstract / Introduction / Conclusion (1 API call)
# ─────────────────────────────────────────────────────────────────────────────
def draft_front_matter(
    herbs:          list,
    drugs:          list,
    hub_genes:      list,
    kegg_terms:     list,
    network_stats:  dict,
    interpretation: str,
    api_key:        str = None,
) -> dict:
    """Abstract, Introduction, Conclusion을 한 번의 API 호출로 생성"""
    herb_str  = ", ".join(herbs)
    drug_str  = ", ".join(drugs) if drugs else "none"
    hub_text  = ", ".join(g["label"] for g in hub_genes[:5])
    kegg_text = ", ".join(t["term"] for t in kegg_terms[:5])

    prompt = f"""Write three sections for a network pharmacology manuscript in formal English.
Return ONLY the three sections separated by the exact markers shown.

## Study Context
- Herbs: {herb_str}
- Co-medications (Western drugs): {drug_str}
- Hub target genes: {hub_text}
- Key KEGG pathways: {kegg_text}
- HCT network: {network_stats.get('hct_nodes','?')} nodes / {network_stats.get('hct_edges','?')} edges
- PPI network: {network_stats.get('ppi_nodes','?')} nodes / {network_stats.get('ppi_edges','?')} edges
- Total target genes: {network_stats.get('total_genes','?')}
- Korean mechanistic summary (reference only): {interpretation[:400]}

===ABSTRACT===
Write a structured abstract (~250 words) with subheadings: Background, Methods, Results, Conclusion.
Use formal academic English. Do NOT include references.

===INTRODUCTION===
Write an Introduction (~400 words) covering:
1. Clinical/pharmacological rationale for studying these herbs
2. Limitations of conventional approaches and rationale for network pharmacology
3. Brief overview of the network pharmacology methodology
4. Aim of this study (last paragraph)
Include [citation needed] placeholders where literature is expected.

===CONCLUSION===
Write a Conclusion (~150 words) summarizing key findings and future directions."""

    try:
        client = _get_client(api_key)
        raw = _call(client, prompt, temperature=0.3)
        sections = {"abstract": "", "introduction": "", "conclusion": ""}
        for key, marker in [("abstract","===ABSTRACT==="),
                             ("introduction","===INTRODUCTION==="),
                             ("conclusion","===CONCLUSION===")]:
            if marker in raw:
                after = raw.split(marker, 1)[1]
                next_markers = [m for m in ["===ABSTRACT===","===INTRODUCTION===","===CONCLUSION==="] if m != marker]
                end = len(after)
                for nm in next_markers:
                    if nm in after:
                        end = min(end, after.index(nm))
                sections[key] = after[:end].strip()
        return sections
    except Exception as e:
        log.error(f"[Gemini] front_matter 오류: {e}")
        return {"abstract": f"[오류] {e}", "introduction": "", "conclusion": ""}


def draft_discovery(
    disease_name: str,
    top_herbs: list,
    disease_genes: list,
    api_key: str = None,
    on_retry=None,
) -> str:
    """
    약재 발굴 결과 AI 해석
    top_herbs: herb_discovery.score_herbs() 반환 리스트
    """
    if not GENAI_AVAILABLE:
        return "[오류] google-genai 패키지가 설치되지 않았습니다."

    client = _get_client(api_key)

    herb_lines = []
    for i, h in enumerate(top_herbs[:10], 1):
        genes = ", ".join(h["overlap_genes"][:8])
        if len(h["overlap_genes"]) > 8:
            genes += f" 외 {len(h['overlap_genes']) - 8}개"
        herb_lines.append(
            f"{i}. {h['pinyin']} ({h.get('cn_name', '')}) — "
            f"오버랩 {h['overlap_count']}개 타깃 | score={h['score']:.3f}\n"
            f"   공유 유전자: {genes}"
        )

    disease_top = ", ".join(list(disease_genes)[:20])

    prompt = (
        f"아래는 '{disease_name}'에 대한 네트워크 약리학적 약재 발굴 결과입니다.\n\n"
        f"분석 방법: TCMSP DB 전체 약재의 활성 성분(OB≥30%, DL≥0.18) 타깃 유전자와 "
        f"OpenTargets '{disease_name}' 관련 유전자를 오버랩 스코어링하여 랭킹.\n\n"
        f"질환 연관 주요 유전자(상위 20개): {disease_top}\n\n"
        f"상위 약재 랭킹:\n" + "\n".join(herb_lines) + "\n\n"
        f"다음을 한국어 학술 문체로 작성하세요:\n"
        f"1. 상위 3~5개 약재의 {disease_name} 관련 잠재 메커니즘 (공유 타깃 유전자 기반)\n"
        f"2. 각 약재가 {disease_name}에 효과적일 수 있는 근거 (전통 사용 + 현대 연구 포함)\n"
        f"3. 후속 연구 방향 제안\n"
        f"총 600~900자 내외."
    )

    return _call(client, prompt, temperature=0.4, on_retry=on_retry)


def _make_figure_legends(
    herbs: list, drugs: list,
    network_stats: dict, hub_genes: list, kegg_terms: list,
) -> str:
    """Figure legends 규칙 기반 자동 생성 (API 호출 없음)"""
    herb_str = ", ".join(herbs)
    drug_str = f" and {', '.join(drugs)}" if drugs else ""
    top_hubs = ", ".join(g["label"] for g in hub_genes[:5])
    top_kegg = "; ".join(t["term"] for t in kegg_terms[:3])

    lines = [
        "**Figure 1. Herb-Compound-Target (HCT) network.**",
        f"The HCT network illustrates the multi-target interactions between active compounds "
        f"from {herb_str}{drug_str} and their protein targets. "
        f"Nodes represent herbs (green hexagons), compounds (blue circles), or target proteins (orange diamonds). "
        f"Edges indicate compound-target binding relationships retrieved from the ChEMBL database. "
        f"The network comprises {network_stats.get('hct_nodes','N/A')} nodes and {network_stats.get('hct_edges','N/A')} edges.",
        "",
        "**Figure 2. Protein-protein interaction (PPI) network and hub target identification.**",
        f"The PPI network was constructed using STRING database (confidence score ≥ 400) "
        f"and comprises {network_stats.get('ppi_nodes','N/A')} nodes and {network_stats.get('ppi_edges','N/A')} edges. "
        f"Node size reflects degree centrality. Hub targets were identified by composite network topology scoring; "
        f"top-ranked hub genes include {top_hubs}.",
        "",
        "**Figure 3. Gene Ontology and KEGG pathway enrichment analysis.**",
        f"Dot plots display significantly enriched GO Biological Process terms and KEGG pathways "
        f"(FDR ≤ 0.05) for the {network_stats.get('total_genes','N/A')} target genes. "
        f"Dot size represents gene count; color intensity indicates adjusted p-value. "
        f"Highlighted pathways include: {top_kegg}.",
        "",
        "**Figure 4. Venn diagram of herb target and disease target overlap.**",
        f"The Venn diagram illustrates the intersection between target genes of {herb_str}{drug_str} "
        f"and known disease-associated genes retrieved from the OpenTargets database. "
        f"Shared targets represent the pharmacological basis of the herb-disease relationship.",
    ]
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# 3. 통합 호출
# ─────────────────────────────────────────────────────────────────────────────
def generate_full_interpretation(pipeline_output: dict, api_key: str = None) -> dict:
    """
    pipeline.run_pipeline() 결과 → 전체 AI 해석 생성

    반환: {interpretation, methods, results, discussion}
    """
    net      = pipeline_output.get("network", {})
    enrich   = pipeline_output.get("enrichment", {})
    cfg_herbs = pipeline_output.get("config", {}).get("herbs", {})
    cfg_drugs = pipeline_output.get("config", {}).get("drugs", [])

    herbs = list(cfg_herbs.keys()) if isinstance(cfg_herbs, dict) else cfg_herbs
    drugs = cfg_drugs if isinstance(cfg_drugs, list) else []

    # 허브 유전자 추출
    hubs_df  = net.get("hubs")
    hub_genes = []
    if hubs_df is not None and not hubs_df.empty:
        hub_genes = hubs_df[["label", "degree", "betweenness"]].to_dict("records")

    # KEGG 결과
    kegg_df   = (enrich.get("results") or {}).get("KEGG")
    kegg_list = []
    if kegg_df is not None and not kegg_df.empty:
        kegg_list = kegg_df.head(10)[["term", "gene_count", "adj_p_value"]].to_dict("records")

    # GO_BP 결과
    gobp_df   = (enrich.get("results") or {}).get("GO_BP")
    gobp_list = []
    if gobp_df is not None and not gobp_df.empty:
        gobp_list = gobp_df.head(8)[["term", "gene_count", "adj_p_value"]].to_dict("records")

    # 네트워크 통계
    hct_g = net.get("hct")
    ppi_g = net.get("ppi")
    network_stats = {
        "hct_nodes":   hct_g.G.number_of_nodes() if hct_g else 0,
        "hct_edges":   hct_g.G.number_of_edges() if hct_g else 0,
        "ppi_nodes":   ppi_g.G.number_of_nodes() if ppi_g else 0,
        "ppi_edges":   ppi_g.G.number_of_edges() if ppi_g else 0,
        "total_genes": len(pipeline_output.get("genes", [])),
    }

    compounds = cfg_herbs if isinstance(cfg_herbs, dict) else {}

    import time as _time
    # 5 RPM 제한: 호출 사이 13초 간격
    interpretation = interpret_network(herbs, drugs, hub_genes, kegg_list, gobp_list, api_key)
    _time.sleep(13)
    methods        = draft_methods(herbs, drugs, compounds, api_key=api_key)
    _time.sleep(13)
    results_text   = draft_results(herbs, drugs, network_stats, hub_genes, kegg_list, api_key)
    _time.sleep(13)
    discussion     = draft_discussion(herbs, drugs, hub_genes, kegg_list, interpretation, api_key)
    _time.sleep(13)
    front_matter   = draft_front_matter(herbs, drugs, hub_genes, kegg_list,
                                        network_stats, interpretation, api_key)

    # Figure legends (API 호출 없이 규칙 기반 생성)
    fig_legends    = _make_figure_legends(herbs, drugs, network_stats, hub_genes, kegg_list)

    return {
        "interpretation": interpretation,
        "methods":        methods,
        "results":        results_text,
        "discussion":     discussion,
        "abstract":       front_matter.get("abstract", ""),
        "introduction":   front_matter.get("introduction", ""),
        "conclusion":     front_matter.get("conclusion", ""),
        "figure_legends": fig_legends,
    }
