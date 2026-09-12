"""
Network Pharmacology Pipeline — Streamlit Frontend v2
포트: 8503  |  실행: streamlit run app.py --server.port 8503
탭 구성: 설정 / 실행 / ADMET / 네트워크+위상 / 질환교차 / 농축분석 / 도킹 / AI해석 / 내보내기
"""

import json, sys, threading, queue, time, os, shutil, pickle
from pathlib import Path
from io import StringIO

import streamlit as st
import pandas as pd

BASE_DIR = Path(__file__).parent
sys.path.insert(0, str(BASE_DIR / "src"))

from settings_manager import load_settings, save_settings

# ── 페이지 설정 ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="NetPharm — 네트워크약리학 파이프라인",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown("""<style>
/* ── 기본 스타일 ── */
.stTabs [data-baseweb="tab"]{font-size:14px;font-weight:600}
div[data-testid="metric-container"]{
    border:1px solid rgba(128,128,128,0.2);
    border-radius:10px;padding:10px 14px;
}
.stButton>button{border-radius:8px}
section[data-testid="stSidebar"]{background:rgba(0,0,0,0.02)}
.stTextInput input,.stTextArea textarea{border-radius:8px}

/* ── 모바일 반응형 (≤768px) ── */
@media (max-width: 768px) {
    /* 탭 — 스크롤 가능하게, 글씨 축소 */
    .stTabs [data-baseweb="tab-list"]{
        overflow-x:auto !important;
        flex-wrap:nowrap !important;
        -webkit-overflow-scrolling:touch;
        scrollbar-width:none;
    }
    .stTabs [data-baseweb="tab-list"]::-webkit-scrollbar{display:none}
    .stTabs [data-baseweb="tab"]{
        font-size:11px !important;
        padding:6px 10px !important;
        white-space:nowrap;
    }

    /* 메인 콘텐츠 패딩 축소 */
    .block-container{
        padding-left:0.75rem !important;
        padding-right:0.75rem !important;
        padding-top:0.5rem !important;
    }

    /* 버튼 — 터치 친화적 높이 */
    .stButton>button{
        min-height:48px !important;
        font-size:15px !important;
    }

    /* 다운로드 버튼 */
    .stDownloadButton>button{
        min-height:48px !important;
        font-size:14px !important;
    }

    /* 2·3열 컬럼 → 수직 스택 */
    [data-testid="column"]{
        width:100% !important;
        flex:1 1 100% !important;
        min-width:100% !important;
    }

    /* 사이드바 — 기본 숨김(Streamlit 기본동작 유지) */
    section[data-testid="stSidebar"] .block-container{
        padding:0.5rem !important;
    }

    /* 슬라이더 트랙 터치 영역 확대 */
    .stSlider [data-baseweb="slider"]{padding:10px 0}

    /* 인풋 폰트 크기 (iOS 자동확대 방지 — 16px 이상) */
    .stTextInput input,
    .stTextArea textarea,
    .stSelectbox select,
    input[type="text"]{font-size:16px !important}

    /* 채팅 인풋 */
    [data-testid="stChatInput"] textarea{font-size:16px !important}

    /* 데이터프레임 — 가로 스크롤 */
    [data-testid="stDataFrame"]{
        overflow-x:auto !important;
        -webkit-overflow-scrolling:touch;
    }

    /* 메트릭 카드 — 더 작게 */
    div[data-testid="metric-container"]{
        padding:6px 10px !important;
    }
    div[data-testid="metric-container"] [data-testid="stMetricValue"]{
        font-size:1.4rem !important;
    }

    /* 헤더 폰트 */
    h1{font-size:1.4rem !important}
    h2{font-size:1.2rem !important}
    h3{font-size:1.05rem !important}

    /* 코드 블록 스크롤 */
    .stCode code,.stCodeBlock{
        font-size:11px !important;
        overflow-x:auto;
    }

    /* 이미지 — 전체 너비 */
    .stImage img{width:100% !important;height:auto !important}

    /* expander 터치 영역 */
    [data-testid="stExpander"] summary{min-height:44px;display:flex;align-items:center}
}

/* ── 초소형 (≤480px, 폴더블/구형 스마트폰) ── */
@media (max-width: 480px) {
    .stTabs [data-baseweb="tab"]{font-size:10px !important;padding:5px 7px !important}
    .block-container{padding-left:0.4rem !important;padding-right:0.4rem !important}
}
</style>""", unsafe_allow_html=True)

# ── 세션 초기화 ───────────────────────────────────────────────────────────────
# ── 저장된 설정 로드 (최초 1회) ──────────────────────────────────────────────
if "settings_loaded" not in st.session_state:
    _saved = load_settings()
    st.session_state.settings_loaded = True
    st.session_state.gemini_key      = _saved.get("gemini_api_key", "")
    st.session_state.setup_model     = _saved.get("setup_model", "gemma-4-31b-it")
    st.session_state.ai_model        = _saved.get("ai_model", "gemini-3.6-flash")
    st.session_state.ppi_score_saved = _saved.get("ppi_min_score", 400)
    st.session_state.fdr_saved       = _saved.get("fdr_cutoff", 0.05)
    st.session_state.top_n_saved     = _saved.get("top_n", 15)
    if st.session_state.gemini_key:
        os.environ["GEMINI_API_KEY"] = st.session_state.gemini_key

DEFAULTS = dict(
    pipeline_result=None, log_lines=[], running=False,
    ai_output={},
    admet_results=None,
    admet_filtered_result=None,
    topology_df=None,
    disease_targets=[], intersection=None,
    docking_df=None,
    interaction_results=[],
    export_paths={},
    current_cfg={},
    cache_stats={},
    chat_messages=[],
    ai_proposed_config=None,
    ai_setup_done=False,
    # 원클릭 전체 분석
    full_run_steps=[],
    full_run_logs=[],
    full_run_done=False,
    full_run_result={},
    full_running=False,
    full_run_state="idle",   # idle | running | done | error
    full_run_error="",
    # DB 검증
    validation_report=None,
    validation_edits={},   # {herb: {comp: new_name}} — 사용자 수정사항
)
for k, v in DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v


# ── 헬퍼: ADMET 필터 적용 여부에 따라 활성 데이터 반환 ────────────────────────
def _active_result() -> dict | None:
    """ADMET 재분석 결과가 있으면 그것을, 없으면 원본 파이프라인 결과를 반환"""
    return st.session_state.admet_filtered_result or st.session_state.pipeline_result


def _active_herb_results() -> dict:
    """ADMET 통과 성분만 담긴 herb_results 반환 (없으면 전체)"""
    admet = st.session_state.admet_results
    base  = (st.session_state.pipeline_result or {}).get("herb_results", {})
    if admet and admet.get("passed"):
        from admet import ADMETScreener, ADMETCriteria
        screener = ADMETScreener()
        return screener.filter_herb_results(base)
    return base


def _active_genes() -> list:
    """현재 활성 herb_results 에서 유전자 심볼 추출"""
    from collect import NetPharmCollector
    hr = _active_herb_results()
    if not hr:
        return (st.session_state.pipeline_result or {}).get("genes", [])
    # NetPharmCollector 없이 직접 추출
    genes = set()
    for hd in hr.values():
        for comp in hd.get("compounds", []):
            for tgt in comp.get("targets", []):
                gs = tgt.get("gene_symbol")
                if gs:
                    genes.add(gs.upper())
    return sorted(genes)

# ── 사이드바 ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🔬 NetPharm")
    st.caption("네트워크약리학 파이프라인 v2")
    st.divider()

    # API 키 상태 표시 (편집은 설정 탭에서)
    if st.session_state.gemini_key:
        st.success("🔑 API 키 설정됨")
    else:
        st.warning("🔑 API 키 미설정\n⚙️ 설정 탭에서 입력")

    st.divider()
    st.subheader("⚙️ 분석 파라미터")
    ppi_score  = st.slider("STRING 신뢰도", 400, 900,
                            int(st.session_state.ppi_score_saved), 100)
    fdr_cutoff = st.select_slider("FDR 컷오프", [0.01, 0.05, 0.1, 0.2],
                                   value=st.session_state.fdr_saved)
    top_n      = st.slider("표시 term 수", 10, 30, int(st.session_state.top_n_saved))
    make_plots = st.toggle("그래프 생성", value=True)

    st.divider()
    if st.session_state.pipeline_result:
        r = st.session_state.pipeline_result
        st.success(f"✅ 분석 완료 ({r.get('elapsed',0):.0f}초)")
        cfg_n = (r.get("config") or {}).get("name","")
        if cfg_n: st.caption(f"프로젝트: {cfg_n}")

# ── 탭 정의 ───────────────────────────────────────────────────────────────────
tabs = st.tabs([
    "💬 AI 설정","▶ 실행","💊 ADMET","🕸 네트워크","🎯 질환교차",
    "📊 농축분석","🧲 도킹","🤖 AI 해석","📦 내보내기",
    "📁 프로젝트","🗄️ DB 관리","⚙️ 설정","🔍 약재발굴","📋 CPG","📝 논문작성"
])
(tab_setup, tab_run, tab_admet, tab_net, tab_disease,
 tab_enrich, tab_dock, tab_ai, tab_export,
 tab_projects, tab_dbmgr, tab_settings, tab_discover, tab_cpg, tab_paper) = tabs


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 1: AI 대화 설정
# ═══════════════════════════════════════════════════════════════════════════════
with tab_setup:
    st.header("💬 AI 연구 설정")
    st.caption("연구 목적을 자유롭게 설명하면 AI가 분석 설정을 자동으로 구성합니다.")

    # ── 채팅 히스토리 표시 (네이티브 chat_message — 다크/라이트 자동 대응) ──
    if not st.session_state.chat_messages:
        st.info(
            "💡 **예시:** \"가미소요산과 스틸녹스(졸피뎀)의 불면증 약물상호작용 분석해줘\"  |  "
            "\"황련해독탕 + 메트포르민 당뇨 분석\"  |  \"보중익기탕 단독 암 표적 분석\""
        )

    for msg in st.session_state.chat_messages:
        role = "user" if msg["role"] == "user" else "assistant"
        with st.chat_message(role):
            st.markdown(msg["content"])

    # ── 입력창 (chat_input: 화면 하단 고정) ──────────────────────────────────
    user_input = st.chat_input("연구 내용을 자유롭게 입력하세요...")

    if user_input and user_input.strip():
        st.session_state.chat_messages.append({"role": "user", "content": user_input})

        with st.chat_message("user"):
            st.markdown(user_input)

        with st.chat_message("assistant"):
            with st.spinner("분석 중..."):
                from ai_setup import ResearchSetupParser
                parser = ResearchSetupParser(api_key=st.session_state.gemini_key or "", model=st.session_state.get("setup_model", "gemma-4-31b-it"))
                ai_result = parser.parse(
                    user_input,
                    chat_history=st.session_state.chat_messages[:-1],
                )

            config  = ai_result.get("config")
            summary = ai_result.get("summary_ko", "")
            clarify = ai_result.get("clarification")

            if config:
                st.session_state.ai_proposed_config = config
                herb_count = len(config.get("herbs", {}))
                drug_count = len(config.get("drugs", []))
                reply = (
                    f"{summary}\n\n"
                    f"**추출된 설정:** 생약 {herb_count}종 · 약물 {drug_count}종 · "
                    f"질환: `{config.get('disease_query','—')}` · 연구명: `{config.get('name','')}`\n\n"
                    f"아래에서 설정을 확인하고 **분석 시작** 버튼을 눌러주세요. ✅"
                )
            elif clarify:
                reply = f"❓ {clarify}"
            else:
                reply = f"⚠️ {summary}\n\n사이드바에 Gemini API 키를 입력하거나 수동 설정을 이용해주세요."

            st.markdown(reply)
            st.session_state.chat_messages.append({"role": "ai", "content": reply})
        st.rerun()

    # ── 제안된 설정 미리보기 + 확정 ─────────────────────────────────────────
    proposed = st.session_state.ai_proposed_config
    if proposed:
        st.divider()
        st.subheader("📋 추출된 설정 미리보기")

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**🌿 생약 구성**")
            for herb, comps in proposed.get("herbs", {}).items():
                st.markdown(f"- **{herb}**: {', '.join(comps[:3])}{'...' if len(comps)>3 else ''}")
        with c2:
            st.markdown("**💊 약물**")
            for d in proposed.get("drugs", []) or ["(없음)"]:
                st.markdown(f"- {d}")
            st.markdown(f"**🏥 질환**: {proposed.get('disease_query','—')}")
            st.markdown(f"**📛 연구명**: `{proposed.get('name','')}`")

        col_ok, col_edit, col_clear = st.columns(3)
        with col_ok:
            if st.button("✅ 설정 확정", type="primary", use_container_width=True, key="btn_apply_config"):
                proposed["ppi_min_score"]    = ppi_score
                proposed["adj_p_cutoff"]     = fdr_cutoff
                proposed["top_n_enrichment"] = top_n
                proposed["make_plots"]       = make_plots
                st.session_state.current_cfg         = proposed
                st.session_state.ai_setup_done       = True
                st.session_state.full_run_done        = False
                st.session_state.full_run_state       = "idle"
                st.session_state.full_run_steps       = []
                st.session_state.full_run_result      = {}
                st.session_state.validation_report    = None  # 자동 검증 트리거
                st.session_state.validation_edits     = {}
                st.rerun()
        with col_edit:
            if st.button("✏️ JSON 직접 수정", use_container_width=True, key="btn_json_edit"):
                st.session_state._show_json_editor = True
        with col_clear:
            if st.button("🗑 초기화", use_container_width=True, key="btn_chat_clear"):
                st.session_state.chat_messages      = []
                st.session_state.ai_proposed_config = None
                st.session_state.ai_setup_done      = False
                st.session_state.validation_report  = None
                st.rerun()

        # JSON 직접 편집
        if st.session_state.get("_show_json_editor"):
            edited = st.text_area(
                "Config JSON 편집", height=300,
                value=json.dumps(proposed, ensure_ascii=False, indent=2),
                key="json_editor_area",
            )
            if st.button("적용", key="btn_json_apply"):
                try:
                    st.session_state.ai_proposed_config = json.loads(edited)
                    st.session_state._show_json_editor = False
                    st.rerun()
                except json.JSONDecodeError as e:
                    st.error(f"JSON 오류: {e}")

    # ── DB 검증 패널 (설정 확정 직후 자동 실행) ───────────────────────────────
    if st.session_state.ai_setup_done and st.session_state.current_cfg:
        cfg_now = st.session_state.current_cfg
        st.divider()
        st.subheader("🔍 DB 검증")

        if not st.session_state.validation_report:
            # 자동 실행
            with st.spinner("DB 검증 중... ChEMBL + OpenTargets (캐시 없을 경우 30~60초)"):
                from validator import validate_config
                from collect   import CacheDB
                _vc = CacheDB()
                _rpt = validate_config(cfg_now, cache=_vc)
                _vc.close()
                st.session_state.validation_report = _rpt
                st.rerun()

        s = st.session_state.validation_report.get("summary", {})
        ok, warn, fail = s.get("ok", 0), s.get("warn", 0), s.get("fail", 0)
        if fail == 0 and warn == 0:
            st.success(f"✅ 전체 {ok}개 정확 매칭 — 바로 분석 시작 가능")
        elif fail == 0:
            st.warning(f"✅ {ok}개 정확 / ⚠️ {warn}개 이름 확인 권장 / ❌ {fail}개 미발견")
        else:
            st.error(f"✅ {ok}개 / ⚠️ {warn}개 / ❌ {fail}개 미발견 — 아래에서 수정하거나 그대로 진행")

        # 재검증 버튼 (수정 후에만 필요)
        if st.button("🔄 재검증", key="btn_revalidate", help="성분명 수정 후 다시 검증"):
            st.session_state.validation_report = None
            st.rerun()

        # ── 검증 결과 표시 + 인라인 수정 ─────────────────────────────────
        report = st.session_state.validation_report
        if report:
            edits = st.session_state.validation_edits

            # 성분 검증 결과
            st.markdown("#### 🌿 생약 성분")
            for herb, comps in report.get("herbs", {}).items():
                with st.expander(f"**{herb}** ({len(comps)}개 성분)", expanded=True):
                    for orig_name, v in comps.items():
                        edit_key = f"edit_{herb}_{orig_name}"
                        status_icon = ("✅" if v["status"] == "exact"
                                       else "⚠️" if v["found"]
                                       else "❌")
                        col_s, col_n, col_m, col_act = st.columns([1, 3, 3, 2])
                        col_s.markdown(f"### {status_icon}")
                        col_n.markdown(f"**입력:** `{orig_name}`")
                        if v["found"]:
                            col_m.markdown(
                                f"**매칭:** `{v['matched_name']}` ({v['chembl_id']})"
                            )
                        else:
                            col_m.markdown(f"*ChEMBL 미발견*")

                        # ── 이름 수정 입력 ──
                        with col_act:
                            new_name = st.text_input(
                                "수정 (선택)", key=edit_key,
                                value=edits.get(herb, {}).get(orig_name, ""),
                                placeholder=v.get("matched_name", "") or "대체 이름",
                                label_visibility="collapsed",
                            )
                            if new_name and new_name != orig_name:
                                if herb not in edits: edits[herb] = {}
                                edits[herb][orig_name] = new_name
                            elif not new_name and herb in edits:
                                edits[herb].pop(orig_name, None)
                        st.session_state.validation_edits = edits

                        # ── 대체 후보 버튼 (not_found 일 때만) ──
                        suggestions = v.get("suggestions", [])
                        if not v["found"] and suggestions:
                            st.caption("💡 대체 후보 (클릭하면 자동 입력):")
                            btn_cols = st.columns(len(suggestions))
                            for i, sug in enumerate(suggestions):
                                label = f"{sug['matched_name']} ({sug['chembl_id']})"
                                if btn_cols[i].button(label, key=f"sug_{herb}_{orig_name}_{i}",
                                                      use_container_width=True):
                                    if herb not in edits: edits[herb] = {}
                                    edits[herb][orig_name] = sug["matched_name"]
                                    st.session_state.validation_edits = edits
                                    st.rerun()
                        elif not v["found"] and not suggestions:
                            st.caption("💡 대체 후보 없음 — 분석 시 PubChem으로 자동 보완")

            # 약물 검증 결과
            if report.get("drugs"):
                st.markdown("#### 💊 약물")
                for drug, v in report["drugs"].items():
                    status_icon = ("✅" if v["status"] == "exact"
                                   else "⚠️" if v["found"] else "❌")
                    col_s, col_n, col_m = st.columns([1, 3, 5])
                    col_s.markdown(f"### {status_icon}")
                    col_n.markdown(f"**입력:** `{drug}`")
                    if v["found"]:
                        col_m.markdown(f"**매칭:** `{v['matched_name']}` ({v['chembl_id']})")
                    else:
                        col_m.markdown(f"*{v['note']}*")

            # 질환 검증 결과
            dv = report.get("disease")
            if dv:
                st.markdown("#### 🏥 질환")
                dicon = "✅" if dv["found"] else "❌"
                if dv["found"]:
                    st.success(
                        f"{dicon} **{dv['disease_name']}** ({dv['disease_id']}) "
                        f"— OpenTargets 관련 타깃 **{dv['n_targets']}개**"
                    )
                else:
                    st.error(f"❌ `{dv['query']}` — {dv['note']}")

            # 수정사항 적용 버튼
            has_edits = any(v for v in edits.values())
            if has_edits:
                st.divider()
                st.info("이름 수정사항이 있습니다. 아래 버튼으로 config에 반영하세요.")
                if st.button("✏️ 수정사항 config에 적용", key="btn_apply_edits",
                             type="secondary", use_container_width=True):
                    updated_cfg = dict(cfg_now)
                    updated_herbs = {}
                    for herb, comps in (cfg_now.get("herbs") or {}).items():
                        new_comps = []
                        for c in comps:
                            replacement = edits.get(herb, {}).get(c, "")
                            new_comps.append(replacement if replacement else c)
                        updated_herbs[herb] = new_comps
                    updated_cfg["herbs"] = updated_herbs
                    st.session_state.current_cfg      = updated_cfg
                    st.session_state.validation_report = None  # 재검증 유도
                    st.session_state.validation_edits  = {}
                    st.success("수정 완료. 'DB 검증 실행'으로 재검증하세요.")
                    st.rerun()

            st.divider()

    # ── 원클릭 전체 분석 ──────────────────────────────────────────────────────
    if st.session_state.ai_setup_done and st.session_state.current_cfg:
        cfg_now = st.session_state.current_cfg
        st.subheader("🚀 논문용 전체 자동 분석")

        if not st.session_state.validation_report:
            st.warning("DB 검증을 먼저 실행하세요. 검증 없이도 분석할 수 있지만 "
                       "성분명 불일치 시 타깃 0개가 될 수 있습니다.")

        # ── 상태별 분기 ──────────────────────────────────────────────────
        _run_state = st.session_state.get("full_run_state", "idle")
        # 하위 호환: 기존 full_run_done/full_running 변수 마이그레이션
        if _run_state == "idle":
            if st.session_state.full_running:
                _run_state = "running"
            elif st.session_state.full_run_done:
                _run_state = "done" if st.session_state.full_run_result.get("pipeline") else "error"

        # ── A: 대기 상태 — 옵션 + 시작 버튼 ─────────────────────────────
        if _run_state == "idle":
            has_disease = bool(cfg_now.get("disease_query"))
            has_key     = bool(st.session_state.gemini_key)

            col_opt1, col_opt2 = st.columns(2)
            with col_opt1:
                opt_dock = st.checkbox(
                    "🧲 분자 도킹 포함",
                    value=True,
                    help="허브 Top3 단백질 × ADMET 통과 상위 5개 성분. 도킹 1쌍당 약 1~3분.",
                    key="full_opt_dock",
                )
                opt_ai = st.checkbox(
                    "🤖 AI 논문 초안 생성",
                    value=has_key,
                    disabled=not has_key,
                    help="Gemini API 키 필요. Methods/Results/Discussion 자동 작성.",
                    key="full_opt_ai",
                )
            with col_opt2:
                total_steps = 5 + int(opt_dock) + int(has_key and opt_ai)
                dock_time   = "~30분" if opt_dock else ""
                st.markdown("**포함 단계:**")
                st.markdown(
                    "✅ 타깃 수집 (ChEMBL + PubChem 보완)  \n"
                    "✅ ADMET 스크리닝  \n"
                    "✅ 네트워크 재구성  \n"
                    f"{'✅' if has_disease else '⬜'} 질환 교차분석  \n"
                    "✅ 허브 토폴로지  \n"
                    f"{'✅' if opt_dock else '⬜'} 분자 도킹 {dock_time}  \n"
                    f"{'✅' if (has_key and opt_ai) else '⬜'} AI 논문 초안  \n"
                    "✅ Excel + Figure 생성"
                )

            if st.button("🚀 전체 분석 시작", type="primary",
                         use_container_width=True, key="btn_full_run"):
                st.session_state.full_run_state  = "running"
                st.session_state.full_running    = True
                st.session_state.full_run_steps  = []
                st.session_state.full_run_logs   = []
                st.session_state.full_run_done   = False
                st.session_state.full_run_result = {}
                st.session_state.full_run_error  = ""
                st.session_state._run_start_ts   = time.time()

                full_q = queue.Queue()
                # 스레드 진입 전에 session_state 값 캡처 (스레드 내부에서 접근 불가)
                _api_key   = st.session_state.get("gemini_key", "") or ""
                _cfg_snap  = cfg_now
                _do_dock   = opt_dock
                _do_ai     = opt_ai

                # ── 텔레그램 봇 연동: 분석 시작/종료 기록 ─────────────
                _proc_file = BASE_DIR / "logs" / "bot_processes.json"
                _analysis_name = _cfg_snap.get("name", "analysis")
                def _bot_procs_write(name, started=True):
                    try:
                        procs = json.loads(_proc_file.read_text()) if _proc_file.exists() else {}
                        if started:
                            procs[name] = {
                                "pid":        os.getpid(),
                                "name":       name,
                                "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                "log":        "",
                                "source":     "streamlit",
                            }
                        else:
                            procs.pop(name, None)
                        _proc_file.write_text(json.dumps(procs, ensure_ascii=False, indent=2))
                    except Exception:
                        pass
                _bot_procs_write(_analysis_name, started=True)

                def _full_run():
                    import logging
                    class FQH(logging.Handler):
                        def emit(self, r): full_q.put({"type": "log", "msg": self.format(r)})
                    root = logging.getLogger()
                    h = FQH()
                    h.setFormatter(logging.Formatter("%(asctime)s  %(message)s", "%H:%M:%S"))
                    root.addHandler(h)
                    try:
                        from full_pipeline import run_full_analysis
                        run_full_analysis(
                            cfg        = _cfg_snap,
                            q          = full_q,
                            api_key    = _api_key,
                            do_docking = _do_dock,
                            do_ai      = _do_ai,
                        )
                    except Exception as e:
                        import traceback
                        full_q.put({"type": "error",
                                    "msg": f"{e}\n{traceback.format_exc()[-600:]}"})
                    finally:
                        root.removeHandler(h)
                        _bot_procs_write(_analysis_name, started=False)

                threading.Thread(target=_full_run, daemon=True).start()

                # ── 실시간 진행 표시 ──────────────────────────────────────
                with st.status("분석 실행 중...", expanded=True) as status_box:
                    step_box = st.empty()
                    log_box  = st.empty()
                    done_flag    = False
                    final_result = {}
                    error_msg    = ""

                    while not done_flag:
                        time.sleep(0.4)
                        while not full_q.empty():
                            msg = full_q.get()
                            mtype = msg.get("type", "")
                            if mtype == "progress":
                                steps = st.session_state.full_run_steps
                                existing = [i for i, s in enumerate(steps)
                                            if s["step"] == msg["step"]]
                                if existing:
                                    steps[existing[0]] = msg
                                else:
                                    steps.append(msg)
                                st.session_state.full_run_steps = steps
                            elif mtype == "log":
                                st.session_state.full_run_logs.append(msg["msg"])
                            elif mtype == "done":
                                done_flag    = True
                                final_result = msg.get("result", {})
                            elif mtype == "error":
                                error_msg = msg["msg"]
                                st.session_state.full_run_logs.append(f"❌ {error_msg}")
                                done_flag = True

                        # 진행 현황 렌더링
                        elapsed = int(time.time() - st.session_state._run_start_ts)
                        step_lines = []
                        for s in st.session_state.full_run_steps:
                            icon = "✅" if s.get("done") and not s.get("warn") else \
                                   "⚠️" if s.get("done") and s.get("warn") else "🔄"
                            step_lines.append(f"{icon} [{s['step']}/{s['total']}] {s['msg']}")
                        if step_lines:
                            step_box.markdown("\n\n".join(step_lines) +
                                              f"\n\n⏱ 경과: {elapsed//60}분 {elapsed%60}초")
                        last_logs = st.session_state.full_run_logs[-6:]
                        if last_logs:
                            log_box.code("\n".join(last_logs), language="text")

                    # 루프 종료 후 status 업데이트
                    if error_msg:
                        status_box.update(label="❌ 분석 중 오류 발생", state="error", expanded=True)
                    else:
                        status_box.update(label="✅ 분석 완료!", state="complete", expanded=False)

                # 결과 세션 저장
                st.session_state.full_run_error = error_msg
                pipe_res = final_result.get("pipeline")
                if pipe_res:
                    st.session_state.pipeline_result = pipe_res
                if final_result.get("admet"):
                    st.session_state.admet_results = final_result["admet"]
                if final_result.get("admet_filtered"):
                    st.session_state.admet_filtered_result = final_result["admet_filtered"]
                if final_result.get("topology") is not None:
                    st.session_state.topology_df = final_result["topology"]
                if final_result.get("intersection"):
                    st.session_state.intersection = final_result["intersection"]
                if final_result.get("disease_targets"):
                    st.session_state.disease_targets = final_result["disease_targets"]
                if final_result.get("docking") is not None:
                    st.session_state.docking_df = final_result["docking"]
                if final_result.get("ai_output"):
                    st.session_state.ai_output = final_result["ai_output"]
                if final_result.get("excel_path"):
                    st.session_state.export_paths["excel"] = final_result["excel_path"]
                if final_result.get("figure_png"):
                    st.session_state.export_paths["multipanel_png"] = final_result["figure_png"]
                if final_result.get("figure_svg"):
                    st.session_state.export_paths["multipanel_svg"] = final_result["figure_svg"]

                _proj_name = cfg_now.get("name", "analysis")
                if cfg_now:
                    _cfg_dir = BASE_DIR / "configs"
                    _cfg_dir.mkdir(exist_ok=True)
                    # config JSON 저장 (오류 여부 무관 — 프로젝트 탭에 표시)
                    try:
                        (_cfg_dir / f"{_proj_name}.json").write_text(
                            json.dumps(cfg_now, ensure_ascii=False, indent=2), encoding="utf-8")
                    except Exception as _e:
                        log.warning(f"[config 저장 실패] {_e}")
                    # pipeline 결과 pkl 캐시
                    if pipe_res:
                        try:
                            with open(BASE_DIR / "results" / f"{_proj_name}_pipeline_cache.pkl", "wb") as _pf:
                                pickle.dump(pipe_res, _pf)
                        except Exception as _e:
                            log.warning(f"[pkl 저장 실패] {_e}")
                    # 원클릭 전체 결과 pkl 캐시 (프로젝트 복원용)
                    try:
                        _safe = {k: v for k, v in final_result.items() if k != "pipeline"}
                        _safe["pipeline"] = pipe_res
                        with open(BASE_DIR / "results" / f"{_proj_name}_fullrun_cache.pkl", "wb") as _pf:
                            pickle.dump(_safe, _pf)
                    except Exception as _e:
                        log.warning(f"[fullrun pkl 저장 실패] {_e}")

                # 품질 진단 자동 실행 (에러가 없고 pipeline 결과 있을 때)
                if not error_msg and final_result.get("pipeline"):
                    try:
                        from quality_check import run_quality_check
                        _qc_report, _qc_md = run_quality_check(
                            pipeline       = final_result["pipeline"],
                            admet          = final_result.get("admet"),
                            admet_filtered = final_result.get("admet_filtered"),
                            intersection   = final_result.get("intersection"),
                            docking        = final_result.get("docking"),
                            cfg            = cfg_now,
                            output_dir     = BASE_DIR / "results",
                        )
                        st.session_state[f"qc_report_{_proj_name}"] = str(_qc_md)
                        st.session_state["_last_qc_report"] = _qc_report
                    except Exception as _e:
                        log.warning(f"[QC 실패] {_e}")

                st.session_state.full_run_state  = "error" if error_msg else "done"
                st.session_state.full_run_done   = True
                st.session_state.full_running    = False
                st.session_state.full_run_result = final_result
                st.rerun()

        # ── B: 실행 중 (페이지 새로고침 등으로 while loop 이탈한 경우) ────
        elif _run_state == "running":
            st.info("🔄 분석이 백그라운드에서 실행 중입니다. 이 페이지를 닫지 마세요.")
            steps = st.session_state.full_run_steps
            if steps:
                lines = [f"{'✅' if s.get('done') else '🔄'} [{s['step']}/{s['total']}] {s['msg']}"
                         for s in steps]
                st.code("\n".join(lines))

        # ── C: 완료 패널 ─────────────────────────────────────────────────
        elif _run_state in ("done", "error"):
            error_msg = st.session_state.get("full_run_error", "")
            steps     = st.session_state.full_run_steps
            r         = st.session_state.full_run_result

            # 단계 요약
            if steps:
                lines = []
                for s in steps:
                    icon = "✅" if s.get("done") and not s.get("warn") else \
                           "⚠️" if s.get("done") and s.get("warn") else "🔄"
                    lines.append(f"{icon} [{s['step']}/{s['total']}] {s['msg']}")
                with st.expander("단계별 결과 보기", expanded=(_run_state == "error")):
                    st.code("\n".join(lines))

            if _run_state == "error":
                st.error(f"❌ 분석 중 오류가 발생했습니다.\n\n```\n{error_msg[:600]}\n```")
            else:
                st.success("✅ 전체 분석 완료! 아래에서 결과를 다운로드하세요.")
                _zero_herbs = r.get("zero_target_herbs", [])
                if _zero_herbs:
                    st.warning(f"⚠ 타깃 0개 약재 감지: **{', '.join(_zero_herbs)}** — DB 검증 탭에서 성분명을 재확인하거나 해당 약재를 제외하세요. 논문 투고 전 반드시 처리 필요.")

                dl_cols = st.columns(4)
                xl = r.get("excel_path") or st.session_state.export_paths.get("excel", "")
                if xl and Path(xl).exists():
                    with dl_cols[0]:
                        st.download_button("📊 Supplementary Excel",
                                           Path(xl).read_bytes(), Path(xl).name,
                                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                                           use_container_width=True, key="dl_full_xl")
                fig_png = r.get("figure_png") or st.session_state.export_paths.get("multipanel_png", "")
                if fig_png and Path(fig_png).exists():
                    with dl_cols[1]:
                        st.download_button("🖼 Figure PNG 300DPI",
                                           Path(fig_png).read_bytes(), Path(fig_png).name,
                                           "image/png", use_container_width=True, key="dl_full_png")
                fig_svg = r.get("figure_svg") or st.session_state.export_paths.get("multipanel_svg", "")
                if fig_svg and Path(fig_svg).exists():
                    with dl_cols[2]:
                        st.download_button("🖼 Figure SVG (편집용)",
                                           Path(fig_svg).read_bytes(), Path(fig_svg).name,
                                           "image/svg+xml", use_container_width=True, key="dl_full_svg")
                ai_md = r.get("ai_draft_path", "")
                if ai_md and Path(ai_md).exists():
                    with dl_cols[3]:
                        st.download_button("📝 AI 논문 초안 (MD)",
                                           Path(ai_md).read_bytes(), Path(ai_md).name,
                                           "text/markdown", use_container_width=True, key="dl_full_ai")
                ai_docx = r.get("manuscript_docx", "")
                if ai_docx and Path(ai_docx).exists():
                    st.download_button(
                        "📄 AI 논문 원고 Word (.docx)",
                        Path(ai_docx).read_bytes(), Path(ai_docx).name,
                        "application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                        use_container_width=True, key="dl_full_docx"
                    )
                # ── 인라인 결과 미리보기 ──────────────────────────────────
                st.divider()
                st.markdown("### 📋 결과 미리보기")

                # 핵심 수치 요약
                _pipe  = r.get("pipeline", {})
                _admet = r.get("admet", {})
                _topo  = r.get("topology")
                _dock  = r.get("docking")
                _inters = r.get("intersection", {})
                _enr   = (r.get("admet_filtered") or _pipe).get("enrichment", {})

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("수집 유전자", f"{len(_pipe.get('genes', []))}개")
                m2.metric("ADMET 통과", f"{len(_admet.get('passed', []))}개" if _admet else "—")
                m3.metric("질환 교차 타깃", f"{_inters.get('shared_count', 0)}개" if _inters else "—")
                n_dock_ok = int((_dock["affinity"] < -5).sum()) if _dock is not None and not _dock.empty else 0
                m4.metric("유효 도킹 결합", f"{n_dock_ok}쌍 (< -5)")

                # 허브 유전자 Top 10
                if _topo is not None and not _topo.empty:
                    with st.expander("🔬 허브 유전자 Top 10", expanded=True):
                        _show = [c for c in ["label","degree","betweenness","closeness"] if c in _topo.columns]
                        st.dataframe(_topo.head(10)[_show].rename(columns={"label":"Gene"}),
                                     use_container_width=True, height=320)

                # KEGG 상위 경로
                _enr_res = (_enr or {}).get("results", {})
                _kegg = _enr_res.get("KEGG") if _enr_res else None
                if _kegg is not None and not _kegg.empty:
                    with st.expander("📊 KEGG Pathway Top 10", expanded=True):
                        _kcols = [c for c in ["term","gene_count","adj_p_value"] if c in _kegg.columns]
                        st.dataframe(_kegg.head(10)[_kcols], use_container_width=True, height=320)

                # 도킹 결과
                if _dock is not None and not _dock.empty:
                    with st.expander("🧲 도킹 결과 (결합에너지 kcal/mol)", expanded=False):
                        st.dataframe(_dock.sort_values("affinity"), use_container_width=True)

                # 통합 Figure
                if fig_png and Path(fig_png).exists():
                    with st.expander("🖼 논문용 통합 Figure", expanded=True):
                        st.image(fig_png, use_container_width=True)

            # ── 학술 품질 진단 (자동 실행 결과 표시) ────────────────────────
            st.divider()
            qc_path_key = f"qc_report_{cfg_now.get('name','')}"
            qc_report   = st.session_state.get("_last_qc_report")
            qc_md_path  = st.session_state.get(qc_path_key, "")

            if qc_report:
                pct   = qc_report["score_pct"]
                stars = qc_report["stars"]
                warns = [c for c in qc_report["checks"] if c["status"] == "warn"]
                caut  = [c for c in qc_report["checks"] if c["status"] == "caution"]
                col_qcs, col_dl_qc = st.columns([3, 1])
                with col_qcs:
                    st.markdown(f"**🔬 학술 품질 진단: {pct}점/100점 {stars}**")
                    if warns:
                        st.error(f"즉시 개선 필요: {', '.join(c['title'] for c in warns)}")
                    elif caut:
                        st.warning(f"투고 전 검토 권장: {', '.join(c['title'] for c in caut)}")
                    else:
                        st.success("모든 항목 양호 — 논문 투고 준비 완료")
                with col_dl_qc:
                    if qc_md_path and Path(qc_md_path).exists():
                        st.download_button(
                            "⬇ 진단 보고서",
                            Path(qc_md_path).read_bytes(),
                            Path(qc_md_path).name,
                            "text/markdown",
                            use_container_width=True,
                            key="dl_qc_md",
                        )
            elif qc_md_path and Path(qc_md_path).exists():
                # 이전 실행 결과만 있는 경우
                st.download_button("⬇ 품질 진단 보고서", Path(qc_md_path).read_bytes(),
                                   Path(qc_md_path).name, "text/markdown", key="dl_qc_md")

            # 다시 분석 버튼 (완료/에러 공통)
            st.divider()
            if st.button("🔄 다시 분석하기", key="btn_full_reset", use_container_width=True):
                for k in ("full_run_state", "full_run_done", "full_running",
                          "full_run_steps", "full_run_logs", "full_run_result", "full_run_error"):
                    st.session_state[k] = DEFAULTS.get(k, False if "run" in k else [] if "steps" in k or "logs" in k else {})
                st.session_state.full_run_state = "idle"
                st.rerun()

    # ── 수동 설정 (접이식 백업) ───────────────────────────────────────────────
    st.divider()
    with st.expander("🔧 수동 설정 (AI 없이 직접 입력)", expanded=not bool(proposed)):
        c1, c2 = st.columns(2)
        with c1:
            st.subheader("🌿 한약")
            analysis_name = st.text_input("분석 이름", "my_analysis")
            n_herbs = st.number_input("한약 수", 1, 6, 2)
            herb_entries = {}
            for i in range(int(n_herbs)):
                with st.expander(f"한약 {i+1}", expanded=(i==0)):
                    hn = st.text_input("학명/영문명", key=f"hn{i}",
                                       placeholder="Panax ginseng")
                    hc = st.text_area("활성 성분 (줄바꿈/쉼표)", key=f"hc{i}", height=70,
                                      placeholder="ginsenoside Rb1\nginsenoside Rg1")
                    if hn and hc:
                        herb_entries[hn] = [x.strip() for x in hc.replace(",","\n").split("\n") if x.strip()]
        with c2:
            st.subheader("💊 양약")
            drugs_inp = st.text_area("약물명", height=70, placeholder="metformin\naspirin")
            drugs_list = [d.strip() for d in drugs_inp.replace(",","\n").split("\n") if d.strip()] if drugs_inp else []
            disease_query   = st.text_input("질환명 (영문)", placeholder="rheumatoid arthritis", key="disease_query_setup")
            disease_min_score = st.slider("OpenTargets 최소 점수", 0.0, 1.0, 0.1, 0.05)

        manual_cfg = dict(
            name=analysis_name, herbs=herb_entries, drugs=drugs_list,
            ppi_min_score=ppi_score, adj_p_cutoff=fdr_cutoff,
            top_n_enrichment=top_n, make_plots=make_plots,
            disease_query=disease_query, disease_min_score=disease_min_score,
            admet=dict(ob_min=30.0, dl_min=0.18, mw_max=500.0,
                       logp_max=5.0, hbd_max=5, hba_max=10),
        )

        col_sv, col_ld = st.columns(2)
        with col_sv:
            if st.button("💾 이 설정 적용", use_container_width=True, key="btn_manual_apply"):
                st.session_state.current_cfg = manual_cfg
                st.success("수동 설정 적용됨 — '▶ 실행' 탭으로 이동하세요.")
        with col_ld:
            cfgs = sorted((BASE_DIR/"configs").glob("*.json"))
            if cfgs:
                sel = st.selectbox("저장된 설정 불러오기", ["—"]+[f.stem for f in cfgs])
                if sel != "—":
                    loaded = json.loads((BASE_DIR/"configs"/f"{sel}.json").read_text())
                    st.session_state.current_cfg    = loaded
                    st.session_state.ai_proposed_config = loaded
                    st.info(f"'{sel}' 로드됨")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 2: 실행
# ═══════════════════════════════════════════════════════════════════════════════
with tab_run:
    st.header("파이프라인 실행")
    cfg = st.session_state.current_cfg
    n_h = len(cfg.get("herbs", {}))

    if n_h:
        src = "AI 설정 ✨" if st.session_state.ai_setup_done else "수동 설정"
        st.info(
            f"**{cfg.get('name','?')}** ({src}) &nbsp;|&nbsp; "
            f"생약 {n_h}종 &nbsp;|&nbsp; 약물 {len(cfg.get('drugs',[]))}종 &nbsp;|&nbsp; "
            f"질환: {cfg.get('disease_query','—') or '—'}"
        )
    else:
        st.warning("← 'AI 설정' 탭에서 연구 내용을 먼저 입력하세요.")

    # ── 실행 옵션 ─────────────────────────────────────────────────────────────
    auto_all = st.toggle(
        "🔄 전체 자동 분석 (ADMET + 질환교차 + 농축분석 순차 실행)",
        value=False,
        help="ON: 기본 파이프라인 완료 후 ADMET·질환교차·농축분석 자동 실행\nOFF: 각 탭에서 수동 실행"
    )

    run_btn = st.button("▶ 분석 시작", type="primary", key="btn_run",
                        disabled=st.session_state.running or n_h == 0,
                        use_container_width=True)
    log_box = st.empty()

    if run_btn and not st.session_state.running:
        st.session_state.running    = True
        st.session_state.log_lines  = []
        st.session_state.pipeline_result = None

        log_q = queue.Queue()

        def _run():
            import logging
            class QH(logging.Handler):
                def emit(self, r): log_q.put(self.format(r))
            root = logging.getLogger()
            h = QH()
            h.setFormatter(logging.Formatter("%(asctime)s  %(message)s", "%H:%M:%S"))
            root.addHandler(h)
            try:
                from pipeline import run_pipeline
                t0  = time.time()
                res = run_pipeline(cfg)
                res["elapsed"] = time.time() - t0
                res["config"]  = cfg

                # 전체 자동 분석
                if auto_all:
                    log_q.put("── ADMET 스크리닝 자동 실행 ──")
                    from admet import ADMETScreener, ADMETCriteria
                    admet_cfg = cfg.get("admet", {})
                    criteria  = ADMETCriteria(**{k: v for k, v in admet_cfg.items()
                                                  if k in ADMETCriteria.__dataclass_fields__})
                    screener  = ADMETScreener(criteria)
                    admet_res = screener.screen_all(res.get("herb_results", {}))
                    res["admet"] = admet_res

                    if cfg.get("disease_query"):
                        log_q.put("── 질환 교차분석 자동 실행 ──")
                        from disease import DiseaseTargetCollector, IntersectionAnalyzer, IntersectionPlotter
                        from collect import CacheDB
                        cache     = CacheDB()
                        collector = DiseaseTargetCollector(cache)
                        did, d_targets = collector.get_targets_by_name(
                            cfg["disease_query"], min_score=cfg.get("disease_min_score", 0.1)
                        )
                        cache.close()
                        herb_genes = set(res.get("genes", []))
                        d_genes    = {t["gene_symbol"] for t in d_targets}
                        analyzer   = IntersectionAnalyzer()
                        inters     = analyzer.intersect(herb_genes, d_genes,
                                                         herb_label="Herbs",
                                                         disease_label=cfg["disease_query"])
                        res["disease_intersection"] = inters
                        res["disease_targets"]      = d_targets
                        IntersectionPlotter(Path("results")).venn_diagram(inters, prefix=cfg.get("name","analysis"))

                log_q.put("__DONE__")
                log_q.put(res)
            except Exception as e:
                import traceback
                log_q.put(f"[오류] {e}\n{traceback.format_exc()[-500:]}")
                log_q.put("__DONE__")
                log_q.put(None)
            finally:
                root.removeHandler(h)

        threading.Thread(target=_run, daemon=True).start()

        done, res_obj = False, None
        while not done:
            time.sleep(0.3)
            while not log_q.empty():
                item = log_q.get()
                if item == "__DONE__":
                    done = True
                elif isinstance(item, dict):
                    res_obj = item
                else:
                    st.session_state.log_lines.append(item)
            log_box.code("\n".join(st.session_state.log_lines[-60:]), language="text")

        st.session_state.pipeline_result = res_obj
        if res_obj and res_obj.get("admet"):
            st.session_state.admet_results = res_obj["admet"]
        if res_obj and res_obj.get("disease_intersection"):
            st.session_state.intersection    = res_obj["disease_intersection"]
            st.session_state.disease_targets = res_obj.get("disease_targets", [])
        # 설정 파일 + 파이프라인 결과 캐시 저장
        if res_obj:
            _proj_name = (st.session_state.current_cfg or {}).get("name", "analysis")
            # configs/ 에 설정 JSON 자동 저장 (프로젝트 관리 탭에 표시)
            _cfg_dir = BASE_DIR / "configs"
            _cfg_dir.mkdir(exist_ok=True)
            try:
                (_cfg_dir / f"{_proj_name}.json").write_text(
                    json.dumps(st.session_state.current_cfg, ensure_ascii=False, indent=2),
                    encoding="utf-8"
                )
            except Exception:
                pass
            # results/ 에 pickle 캐시 저장 (불러오기 복원용)
            _cache_path = BASE_DIR / "results" / f"{_proj_name}_pipeline_cache.pkl"
            try:
                with open(_cache_path, "wb") as _pf:
                    pickle.dump(res_obj, _pf)
            except Exception:
                pass
        st.session_state.running = False
        st.rerun()

    if st.session_state.log_lines:
        log_box.code("\n".join(st.session_state.log_lines[-80:]), language="text")

    if st.session_state.pipeline_result:
        st.success("✅ 기본 파이프라인 완료! 각 탭에서 추가 분석을 실행하세요.")
        rp = st.session_state.pipeline_result.get("report")
        if rp and Path(rp).exists():
            st.download_button("📄 보고서 다운로드",
                               Path(rp).read_text(encoding="utf-8"),
                               Path(rp).name, "text/markdown")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 3: ADMET
# ═══════════════════════════════════════════════════════════════════════════════
with tab_admet:
    st.header("💊 ADMET 약물성 스크리닝")
    result = st.session_state.pipeline_result

    if not result:
        st.info("먼저 '실행' 탭에서 파이프라인을 실행하세요.")
    else:
        cfg      = result.get("config", {})
        admet_cfg = cfg.get("admet", {})

        if st.button("ADMET 스크리닝 실행", type="primary", key="btn_admet_run"):
            from admet import ADMETScreener, ADMETCriteria
            criteria = ADMETCriteria(**{k:v for k,v in admet_cfg.items()
                                        if k in ADMETCriteria.__dataclass_fields__})
            screener = ADMETScreener(criteria)
            with st.spinner("스크리닝 중..."):
                st.session_state.admet_results = screener.screen_all(result["herb_results"])
            st.rerun()

        admet = st.session_state.admet_results
        if admet:
            n_pass = len(admet["passed"])
            n_fail = len(admet["failed"])
            c1,c2,c3 = st.columns(3)
            c1.metric("전체 성분",    n_pass+n_fail)
            c2.metric("✅ 통과",       n_pass, delta=f"{n_pass/(n_pass+n_fail)*100:.0f}%")
            c3.metric("❌ 탈락",       n_fail)

            df = admet["summary"]
            if not df.empty:
                st.subheader("전체 결과")
                show_cols = [c for c in ["compound_name","ob","dl","mw","logp",
                                          "hbd","hba","tpsa","pass_filter","fail_reasons"]
                             if c in df.columns]
                st.dataframe(
                    df[show_cols].style.map(
                        lambda v: "background-color:#c8e6c9" if v is True
                                  else ("background-color:#ffcdd2" if v is False else ""),
                        subset=["pass_filter"] if "pass_filter" in show_cols else []
                    ),
                    use_container_width=True, height=400,
                )

                st.subheader("❌ 탈락 사유")
                fail_df = df[df["pass_filter"]==False][["compound_name","fail_reasons"]] if "pass_filter" in df else pd.DataFrame()
                if not fail_df.empty:
                    st.dataframe(fail_df, use_container_width=True)

                st.download_button("⬇ ADMET 결과 CSV",
                                   df.to_csv(index=False).encode(),
                                   f"{cfg.get('name','analysis')}_admet.csv")

            st.divider()
            st.subheader("🔄 ADMET 필터 적용 → 네트워크·농축분석 재구성")
            st.caption(
                f"통과 성분({n_pass}개)의 타깃 유전자만으로 네트워크와 GO/KEGG를 다시 계산합니다. "
                "이후 네트워크·농축분석·도킹 탭에서 필터된 결과가 자동 사용됩니다."
            )
            if n_pass == 0:
                st.warning("ADMET 통과 성분이 없어 재분석할 수 없습니다.")
            elif st.button("🔄 ADMET 필터 기반 재분석 실행", type="primary", key="btn_admet_rerun",
                           use_container_width=True):
                from admet import ADMETScreener
                from collect import NetPharmCollector
                from network import build_full_network
                from enrichment import run_enrichment

                with st.spinner("ADMET 필터 적용 후 네트워크·농축분석 재구성 중..."):
                    base_hr  = result["herb_results"]
                    screener = ADMETScreener()
                    filt_hr  = screener.filter_herb_results(base_hr)

                    # 유전자 추출
                    genes = set()
                    for hd in filt_hr.values():
                        for comp in hd.get("compounds", []):
                            for tgt in comp.get("targets", []):
                                gs = tgt.get("gene_symbol")
                                if gs: genes.add(gs.upper())
                    genes = sorted(genes)

                    # PPI 재수집
                    nc = NetPharmCollector()
                    ppi_edges = nc.collect_ppi(genes, min_score=cfg.get("ppi_min_score", 400)) if genes else []
                    nc.close()

                    # 네트워크 재구성
                    name = cfg.get("name","analysis") + "_admet"
                    net_f  = build_full_network(filt_hr, ppi_edges, prefix=name,
                                               min_ppi_score=cfg.get("ppi_min_score",400))
                    enr_f  = run_enrichment(genes, prefix=name,
                                            adj_p_cutoff=cfg.get("adj_p_cutoff",0.05),
                                            top_n=cfg.get("top_n_enrichment",20),
                                            make_plots=cfg.get("make_plots",True))

                st.session_state.admet_filtered_result = {
                    **result,
                    "herb_results": filt_hr,
                    "genes":        genes,
                    "network":      net_f,
                    "enrichment":   enr_f,
                    "config":       {**cfg, "name": name},
                }
                st.session_state.topology_df = None   # 재분석 시 토폴로지 초기화
                st.success(f"재분석 완료: 유전자 {len(genes)}개 | 네트워크·농축분석 업데이트됨")
                st.rerun()

            if st.session_state.admet_filtered_result:
                st.info("✅ ADMET 필터 기반 재분석 결과가 네트워크·농축분석·도킹 탭에 적용 중")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 4: 네트워크 + 다중 토폴로지
# ═══════════════════════════════════════════════════════════════════════════════
with tab_net:
    st.header("🕸 네트워크 분석 + 허브 유전자 랭킹")
    result = _active_result()

    if st.session_state.admet_filtered_result:
        st.info("✅ **ADMET 필터 적용 결과** 표시 중 — 통과 성분 타깃만 포함")
    elif st.session_state.admet_results:
        st.warning("⚠️ ADMET 스크리닝 완료됐지만 재분석 미실행 — 아래는 전체 성분 기반 네트워크")

    if not result or not result.get("network"):
        st.info("분석을 먼저 실행하세요.")
    else:
        net = result["network"]
        hct, ppi, itg = net["hct"].G, net["ppi"].G, net["integrated"].G
        c1,c2,c3,c4 = st.columns(4)
        c1.metric("HCT 노드", hct.number_of_nodes())
        c2.metric("HCT 엣지", hct.number_of_edges())
        c3.metric("PPI 노드", ppi.number_of_nodes())
        c4.metric("PPI 엣지", ppi.number_of_edges())

        st.divider()
        st.subheader("허브 유전자 다중 토폴로지 분석")
        top_n_hub = st.slider("허브 표시 수", 5, 20, 10, key="hub_n")

        if st.button("🔍 다중 토폴로지 분석 실행", type="primary", key="btn_topo_run"):
            from topology import TopologyAnalyzerExtended
            with st.spinner("Degree / Betweenness / Closeness / Eigenvector / MCC 계산 중..."):
                analyzer = TopologyAnalyzerExtended(net["integrated"].G)
                topo_df  = analyzer.compute(node_type_filter="target")
                st.session_state.topology_df = topo_df
            st.rerun()

        topo = st.session_state.topology_df
        if topo is not None and not topo.empty:
            top_df = topo.head(top_n_hub)
            _show_cols = [c for c in
                          ["rank","label","degree","betweenness","closeness",
                           "eigenvector","mcc","composite_score"]
                          if c in top_df.columns]
            st.dataframe(
                top_df[_show_cols].rename(columns={"label":"Gene","composite_score":"Score"}),
                use_container_width=True, height=380,
            )

            cfg_name = (result.get("config") or {}).get("name","analysis")
            img_path = BASE_DIR / "results" / f"{cfg_name}_hub_ranking.png"
            if not img_path.exists():
                from topology import TopologyAnalyzerExtended
                analyzer = TopologyAnalyzerExtended(net["integrated"].G)
                img_path = analyzer.plot_hub_ranking(topo, top_n=top_n_hub, prefix=cfg_name)
            if img_path and Path(img_path).exists():
                st.image(str(img_path), use_container_width=True)

            st.download_button("⬇ 토폴로지 CSV",
                               topo.to_csv(index=False).encode(),
                               f"{cfg_name}_topology_extended.csv")

        st.divider()
        st.subheader("📂 Cytoscape 파일")
        name = (result.get("config") or {}).get("name","analysis")
        for fn in [f"{name}_full.graphml",f"{name}_edges.csv",f"{name}_nodes.csv"]:
            fp = BASE_DIR/"networks"/fn
            if fp.exists():
                st.download_button(f"⬇ {fn}", fp.read_bytes(), fn, key=f"dl_{fn}")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 5: 질환 교차 분석
# ═══════════════════════════════════════════════════════════════════════════════
with tab_disease:
    st.header("🎯 질환-타깃 교차 분석")
    base_result = st.session_state.pipeline_result
    act_result  = _active_result()   # ADMET 필터 적용 결과 우선

    if not base_result:
        st.info("분석을 먼저 실행하세요.")
    else:
        if st.session_state.admet_filtered_result:
            st.info("✅ **ADMET 필터 적용 유전자** 기준으로 교차 분석합니다.")

        cfg = (act_result or base_result).get("config", {})

        col_q, col_s = st.columns([3, 1])
        with col_q:
            default_dq = cfg.get("disease_query", "")
            disease_q  = st.text_input("질환명 (영문)", value=default_dq,
                                        placeholder="rheumatoid arthritis")
        with col_s:
            st.markdown("")  # spacing
            max_tgt = st.number_input("최대 타겟 수", 50, 500, 200, 50)

        # OpenTargets 필터 옵션
        st.subheader("📊 OpenTargets 연관 점수 필터")
        st.caption(
            "OpenTargets Overall Association Score (0~1): 유전자-질환 연관 근거 강도를 통합한 점수. "
            "0.1 이상 = 일반 분석, 0.3 이상 = 고신뢰도 타겟 선별 권장"
        )
        col_f1, col_f2 = st.columns(2)
        with col_f1:
            min_score = st.slider(
                "최소 연관 점수 (min_score)",
                min_value=0.0, max_value=1.0,
                value=float(cfg.get("disease_min_score", 0.1)),
                step=0.05,
                help="낮을수록 더 많은 타겟 포함 (TP 증가, 노이즈도 증가)"
            )
        with col_f2:
            st.metric("설정된 컷오프", f"≥ {min_score:.2f}",
                      delta="고신뢰" if min_score >= 0.3 else "표준")

        if st.button("🔍 질환 타겟 조회 + 교차 분석", type="primary", key="btn_disease_run",
                     use_container_width=True) and disease_q:
            from disease import DiseaseTargetCollector, IntersectionAnalyzer, IntersectionPlotter
            from collect import CacheDB
            cache = CacheDB()
            collector = DiseaseTargetCollector(cache)
            with st.spinner(f"OpenTargets에서 '{disease_q}' 타겟 조회 중 (score ≥ {min_score})..."):
                did, d_targets = collector.get_targets_by_name(
                    disease_q, min_score=min_score
                )
                # max_tgt 클라이언트 측 추가 필터
                d_targets = sorted(d_targets, key=lambda x: x.get("score",0), reverse=True)[:max_tgt]
            cache.close()

            st.session_state.disease_targets = d_targets
            # ADMET 필터 적용 유전자 우선 사용
            herb_genes = set(_active_genes())
            d_genes    = {t["gene_symbol"] for t in d_targets}

            analyzer  = IntersectionAnalyzer()
            herbs_str = " + ".join(list((cfg.get("herbs") or {}).keys())[:2])
            inters = analyzer.intersect(herb_genes, d_genes,
                                         herb_label=herbs_str or "Herbs",
                                         disease_label=disease_q)
            st.session_state.intersection = inters

            plotter = IntersectionPlotter(BASE_DIR/"results")
            plotter.venn_diagram(inters, prefix=cfg.get("name","analysis"))
            plotter.shared_gene_barplot(inters, d_targets, prefix=cfg.get("name","analysis"))
            st.rerun()

        inters    = st.session_state.intersection
        d_targets = st.session_state.disease_targets

        if inters:
            c1,c2,c3 = st.columns(3)
            c1.metric("한약 타겟", inters["herb_total"])
            c2.metric("질환 타겟", inters["disease_total"])
            c3.metric("🎯 공유 타겟", inters["shared_count"])

            st.subheader(f"공유 유전자 ({inters['shared_count']}개)")
            if inters["shared"]:
                st.write(", ".join(f"**{g}**" for g in inters["shared"]))

            col_v, col_b = st.columns(2)
            with col_v:
                vp = BASE_DIR/"results"/f"{cfg.get('name','analysis')}_venn.png"
                if vp.exists(): st.image(str(vp), use_container_width=True)
            with col_b:
                bp = BASE_DIR/"results"/f"{cfg.get('name','analysis')}_shared_genes.png"
                if bp.exists(): st.image(str(bp), use_container_width=True)

            if d_targets:
                with st.expander(f"질환 타겟 전체 ({len(d_targets)}개) — score ≥ {min_score}"):
                    df_dt = pd.DataFrame(d_targets).sort_values("score", ascending=False)
                    st.dataframe(df_dt, use_container_width=True, height=300)


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 6: 농축분석
# ═══════════════════════════════════════════════════════════════════════════════
with tab_enrich:
    st.header("📊 GO / KEGG 농축분석")
    result = _active_result()

    if st.session_state.admet_filtered_result:
        st.info("✅ **ADMET 필터 적용 유전자** 기반 농축분석 결과")
    elif st.session_state.admet_results:
        st.warning("⚠️ ADMET 재분석 미실행 — 전체 성분 유전자 기반 결과")

    if not result or not result.get("enrichment"):
        st.info("분석을 먼저 실행하세요.")
    else:
        enrich = result["enrichment"]
        results_e = enrich.get("results", {})
        name = (result.get("config") or {}).get("name","analysis")
        enr_dir = BASE_DIR/"enrichment"

        db_tabs = st.tabs(list(results_e.keys()))
        for dbt, (alias, df) in zip(db_tabs, results_e.items()):
            with dbt:
                if df is None or df.empty:
                    st.warning("유의한 term 없음"); continue
                st.metric("FDR≤0.05 term", len(df))
                sc = [c for c in ["term","gene_count","adj_p_value","combined_score","genes"] if c in df.columns]
                st.dataframe(df[sc].head(30), use_container_width=True, height=300)
                cd, cb = st.columns(2)
                for col, sfx in [(cd,"dotplot"),(cb,"barplot")]:
                    ip = enr_dir/f"{name}_{alias}_{sfx}.png"
                    if ip.exists():
                        with col: st.image(str(ip), use_container_width=True)

        cp = enr_dir/f"{name}_combined_dotplot.png"
        if cp.exists():
            st.divider(); st.subheader("통합 Dot Plot")
            st.image(str(cp), use_container_width=True)

        st.divider()
        sp = enr_dir/f"{name}_enrichment_summary.csv"
        if sp.exists():
            st.download_button("⬇ 통합 요약 CSV", sp.read_bytes(), sp.name)


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 7: 분자 도킹
# ═══════════════════════════════════════════════════════════════════════════════
with tab_dock:
    st.header("🧲 자동 분자 도킹 (AutoDock Vina)")
    result = st.session_state.pipeline_result

    from docking import check_tools
    tools = check_tools()
    # 미설치 시에만 경고 표시
    if not tools["vina"] or not tools["obabel"]:
        missing = []
        if not tools["vina"]:   missing.append("AutoDock Vina")
        if not tools["obabel"]: missing.append("OpenBabel")
        st.error(f"필수 도구 미설치: {', '.join(missing)}")
        if not tools["vina"]:
            st.code("sudo apt-get install -y autodock-vina", language="bash")
        if not tools["obabel"]:
            st.code("sudo apt-get install -y openbabel", language="bash")

    result = _active_result()
    if not result:
        st.info("분석을 먼저 실행하세요.")
    elif tools["vina"] and tools["obabel"]:
        topo = st.session_state.topology_df
        if topo is None or topo.empty:
            st.warning("'네트워크' 탭에서 다중 토폴로지 분석을 먼저 실행하세요.")
        else:
            # ── 리간드 소스 표시 ──────────────────────────────────────────────
            admet_ok = st.session_state.admet_results and st.session_state.admet_results.get("passed")
            if admet_ok:
                n_pass = len(st.session_state.admet_results["passed"])
                st.success(f"✅ ADMET 통과 성분 {n_pass}개를 리간드로 사용합니다.")
            else:
                st.info("ADMET 스크리닝 미실행 — 전체 수집 성분을 리간드로 사용합니다.")

            st.divider()
            c_lig, c_rec, c_exh = st.columns(3)
            max_lig = c_lig.slider("최대 리간드 수", 1, 20, 5)
            max_rec = c_rec.slider("최대 수용체 수", 1, 10, 3)
            exh     = c_exh.slider("Exhaustiveness", 1, 16, 8,
                                    help="높을수록 정확하지만 느림 (N150: 8 권장)")

            # ── Grid Box 설정 ─────────────────────────────────────────────────
            st.subheader("📦 Grid Box 설정")
            grid_mode = st.radio(
                "중심 좌표 설정 방식",
                ["자동 (수용체 Cα 무게중심)", "수동 입력"],
                horizontal=True,
                help="자동: PDB에서 Cα 잔기들의 무게중심을 계산. 수동: 결합 포켓 좌표를 직접 입력"
            )

            box_center = None
            if grid_mode == "수동 입력":
                st.caption("결합 포켓 중심 좌표 (Å) — Chimera/PyMOL에서 확인 후 입력")
                col_cx, col_cy, col_cz = st.columns(3)
                cx = col_cx.number_input("Center X", value=0.0, format="%.2f")
                cy = col_cy.number_input("Center Y", value=0.0, format="%.2f")
                cz = col_cz.number_input("Center Z", value=0.0, format="%.2f")
                box_center = (cx, cy, cz)

            col_sx, col_sy, col_sz = st.columns(3)
            sx = col_sx.number_input("Size X (Å)", value=20.0, min_value=10.0, max_value=60.0, step=5.0)
            sy = col_sy.number_input("Size Y (Å)", value=20.0, min_value=10.0, max_value=60.0, step=5.0)
            sz = col_sz.number_input("Size Z (Å)", value=20.0, min_value=10.0, max_value=60.0, step=5.0)
            box_size = (sx, sy, sz)

            if grid_mode == "자동 (수용체 Cα 무게중심)":
                st.caption(
                    "수용체마다 독립적으로 Cα 무게중심을 계산합니다. "
                    "전체 단백질을 덮으므로 결합 포켓 밖을 포함할 수 있습니다. "
                    "결합 포켓을 알고 있다면 수동 입력을 권장합니다."
                )

            # ── 도킹 실행 ─────────────────────────────────────────────────────
            if st.button("🧲 도킹 실행", type="primary", use_container_width=True, key="btn_dock_run"):
                from docking import BatchDocking

                # ADMET 필터 적용 리간드 우선
                if admet_ok:
                    passed = st.session_state.admet_results["passed"]
                    compounds = [
                        {"chembl_id": r.compound_id, "pref_name": r.compound_name,
                         "smiles": r.smiles}
                        for r in passed if r.smiles
                    ]
                else:
                    compounds = []
                    for hd in result.get("herb_results", {}).values():
                        for c in hd.get("compounds", []):
                            mol = c.get("molecule", {})
                            if mol.get("smiles"):
                                compounds.append(mol)

                _topo_targets = (topo[topo["node_type"] == "target"]
                                 if "node_type" in topo.columns else topo)
                hub_targets = [
                    {"label": row["label"], "uniprot_id": row.get("uniprot_id","") or ""}
                    for _, row in _topo_targets.head(max_rec).iterrows()
                ]

                batcher = BatchDocking(max_ligands=max_lig, max_receptors=max_rec)
                batcher.vina.exhaustiveness = exh

                n_lig = min(len(compounds), max_lig)
                n_rec = min(len(hub_targets), max_rec)
                st.info(f"⏳ 도킹 실행 중: 리간드 {n_lig}개 × 수용체 {n_rec}개\n"
                        f"Vina exhaustiveness={exh} | 성분당 약 {exh}~{exh*2}분 소요")
                dock_log = st.empty()
                dock_q   = queue.Queue()

                def _dock_run():
                    import logging
                    class DockQH(logging.Handler):
                        def emit(self, r): dock_q.put(self.format(r))
                    root = logging.getLogger()
                    h = DockQH()
                    h.setFormatter(logging.Formatter("%(asctime)s  %(message)s", "%H:%M:%S"))
                    root.addHandler(h)
                    try:
                        df = batcher.run(compounds, hub_targets,
                                         box_size=box_size, box_center=box_center)
                        dock_q.put("__DOCK_DONE__")
                        dock_q.put(df)
                    except Exception as e:
                        dock_q.put(f"[오류] {e}")
                        dock_q.put("__DOCK_DONE__")
                        dock_q.put(None)
                    finally:
                        root.removeHandler(h)

                threading.Thread(target=_dock_run, daemon=True).start()

                dock_lines, dock_done, dock_result = [], False, None
                while not dock_done:
                    time.sleep(0.5)
                    while not dock_q.empty():
                        item = dock_q.get()
                        if item == "__DOCK_DONE__": dock_done = True
                        elif isinstance(item, pd.DataFrame): dock_result = item
                        else: dock_lines.append(item)
                    dock_log.code("\n".join(dock_lines[-30:]) or "대기 중...", language="text")

                st.session_state.docking_df = dock_result
                st.rerun()

    dock_df = st.session_state.docking_df
    if dock_df is not None and not dock_df.empty:
        st.divider()
        st.subheader("도킹 결과 (결합 에너지, kcal/mol)")
        st.caption("결합 에너지 < -7 kcal/mol: 강한 결합, < -5: 중간 결합 (일반 기준)")
        st.dataframe(dock_df.sort_values("affinity"), use_container_width=True)
        st.download_button("⬇ 도킹 결과 CSV",
                           dock_df.to_csv(index=False).encode(),
                           "docking_results.csv")

        st.subheader("결합 에너지 히트맵")
        pivot = dock_df.pivot_table(index="compound", columns="target",
                                     values="affinity", aggfunc="min")
        if not pivot.empty:
            import matplotlib.pyplot as plt
            import numpy as np
            fig, ax = plt.subplots(figsize=(min(14, len(pivot.columns)*2+2),
                                             min(10, len(pivot)*1+2)))
            im = ax.imshow(pivot.values.astype(float), cmap="RdYlGn_r", aspect="auto")
            ax.set_xticks(range(len(pivot.columns)))
            ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=9)
            ax.set_yticks(range(len(pivot.index)))
            ax.set_yticklabels(pivot.index, fontsize=9)
            plt.colorbar(im, ax=ax, label="Affinity (kcal/mol)")
            ax.set_title("Docking Affinity Heatmap", fontsize=12, fontweight="bold")
            for i in range(len(pivot.index)):
                for j in range(len(pivot.columns)):
                    v = pivot.iloc[i, j]
                    if pd.notna(v):
                        ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=8)
            plt.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

        # ── 2D 상호작용 분석 ─────────────────────────────────────────────────
        st.divider()
        st.subheader("🔬 단백질-리간드 2D 상호작용 분석")
        st.caption(
            "도킹 포즈에서 H-Bond · 소수성 · Pi-Stacking · Salt Bridge 접촉 잔기를 "
            "자동 검출합니다. 수용체 PDB와 Vina 출력 파일이 필요합니다."
        )

        success_rows = dock_df[dock_df.get("success", pd.Series([False]*len(dock_df))) == True] \
                       if "success" in dock_df.columns else dock_df.dropna(subset=["affinity"])

        if success_rows.empty:
            st.warning("성공한 도킹 결과가 없어 상호작용 분석을 건너뜁니다.")
        else:
            if st.button("🔬 상호작용 분석 실행", type="primary", use_container_width=True, key="btn_interact_run"):
                from interaction import batch_analyze
                from docking import DOCKING_DIR

                # 수용체 PDB 경로 수집
                rec_pdbs = {}
                topo = st.session_state.topology_df
                if topo is not None:
                    for _, row in topo.iterrows():
                        gene = row.get("label","")
                        for suffix in [f"{gene}.pdb", f"AF_{gene}.pdb"]:
                            p = DOCKING_DIR / "receptors" / suffix
                            if p.exists():
                                rec_pdbs[gene] = p; break

                # 리간드 SMILES 수집
                admet_ok = st.session_state.admet_results and st.session_state.admet_results.get("passed")
                lig_smiles = {}
                if admet_ok:
                    for r in st.session_state.admet_results["passed"]:
                        lig_smiles[r.compound_id] = r.smiles
                        lig_smiles[r.compound_name] = r.smiles
                else:
                    for hd in (result or {}).get("herb_results",{}).values():
                        for c in hd.get("compounds",[]):
                            mol = c.get("molecule",{})
                            if mol.get("smiles"):
                                lig_smiles[mol.get("chembl_id","")] = mol["smiles"]
                                lig_smiles[mol.get("pref_name","")] = mol["smiles"]

                cfg_name = (result.get("config") or {}).get("name","analysis")
                with st.spinner("상호작용 분석 중..."):
                    iresults = batch_analyze(
                        success_rows, rec_pdbs, lig_smiles,
                        DOCKING_DIR, prefix=cfg_name
                    )
                st.session_state.interaction_results = iresults
                st.rerun()

        # 결과 표시
        iresults = st.session_state.interaction_results
        if iresults:
            st.success(f"✅ 상호작용 분석 완료: {len(iresults)}개 쌍")
            pair_labels = [f"{r['compound'][:20]} vs {r['target']}" for r in iresults]
            sel = st.selectbox("쌍 선택", pair_labels)
            ir  = next((r for r in iresults
                        if f"{r['compound'][:20]} vs {r['target']}" == sel), None)
            if ir:
                n_cont = ir.get("n_contacts", 0)
                st.metric("접촉 잔기 수", n_cont)
                c_w, c_b = st.columns(2)
                with c_w:
                    wp = ir.get("wheel_png")
                    if wp and Path(wp).exists():
                        st.image(str(wp), caption="Interaction Wheel", use_container_width=True)
                with c_b:
                    bp = ir.get("bar_png")
                    if bp and Path(bp).exists():
                        st.image(str(bp), caption="Interaction Summary", use_container_width=True)

                if ir.get("ligand_svg") and Path(ir["ligand_svg"]).exists():
                    st.download_button(
                        "⬇ 리간드 2D SVG",
                        Path(ir["ligand_svg"]).read_bytes(),
                        Path(ir["ligand_svg"]).name, "image/svg+xml",
                    )

                cdf = ir.get("contacts_df")
                if cdf is not None and not cdf.empty:
                    with st.expander("접촉 잔기 전체 목록"):
                        st.dataframe(cdf, use_container_width=True, height=300)
                    cfg_name = (result.get("config") or {}).get("name","analysis")
                    st.download_button(
                        "⬇ 접촉 잔기 CSV",
                        cdf.to_csv(index=False).encode(),
                        f"{cfg_name}_{ir['compound'][:15]}_vs_{ir['target']}_contacts.csv",
                    )


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 8: AI 해석
# ═══════════════════════════════════════════════════════════════════════════════
with tab_ai:
    st.header("🤖 AI 해석 — Gemini 2.0 Flash")
    result = st.session_state.pipeline_result

    if not result:
        st.info("분석을 먼저 실행하세요.")
    elif not st.session_state.gemini_key:
        st.warning("사이드바에서 Gemini API Key를 입력하세요.")
    else:
        st.caption(f"모델: {st.session_state.get('ai_model', 'gemini-3.6-flash')} | 무료 티어")

        def _get_context():
            net    = result.get("network",{})
            enrich = result.get("enrichment",{})
            cfg_   = result.get("config",{})
            herbs  = list((cfg_.get("herbs") or {}).keys())
            drugs  = cfg_.get("drugs",[])
            hubs_df = net.get("hubs") if net.get("hubs") is not None else st.session_state.topology_df
            hub_genes = hubs_df[["label","degree","betweenness"]].to_dict("records") \
                        if hubs_df is not None and not hubs_df.empty else []
            er = (enrich.get("results") or {})
            kegg_df = er.get("KEGG"); gobp_df = er.get("GO_BP")
            kegg = kegg_df.head(10)[["term","gene_count","adj_p_value"]].to_dict("records") \
                   if kegg_df is not None and not kegg_df.empty else []
            gobp = gobp_df.head(8)[["term","gene_count","adj_p_value"]].to_dict("records") \
                   if gobp_df is not None and not gobp_df.empty else []
            return herbs, drugs, hub_genes, kegg, gobp

        c1,c2,c3,c4 = st.columns(4)
        for col, label, section in [
            (c1,"작용기전\n(한국어)","interpretation"),
            (c2,"Methods\n(영문)","methods"),
            (c3,"Results\n(영문)","results"),
            (c4,"Discussion\n(영문)","discussion"),
        ]:
            with col:
                if st.button(label, use_container_width=True, key=f"ai_{section}"):
                    from interpreter import (interpret_network, draft_methods,
                                              draft_results, draft_discussion)
                    herbs, drugs, hubs, kegg, gobp = _get_context()
                    key = st.session_state.gemini_key
                    with st.status("Gemini 생성 중...", expanded=True) as _st:
                        _msg = st.empty()
                        _msg.write("Gemini API 호출 중...")

                        def _on_retry(attempt, wait_sec, err_msg):
                            _st.update(label=f"Gemini 재시도 중... ({attempt}/3)")
                            _msg.warning(
                                f"일시적 오류 (503/429) — {wait_sec}초 후 자동 재시도 ({attempt}/3)\n\n"
                                f"`{err_msg[:120]}`"
                            )

                        if section == "interpretation":
                            st.session_state.ai_output[section] = interpret_network(
                                herbs, drugs, hubs, kegg, gobp, key, on_retry=_on_retry)
                        elif section == "methods":
                            cfg_ = result.get("config",{})
                            st.session_state.ai_output[section] = draft_methods(
                                herbs, drugs, cfg_.get("herbs",{}), api_key=key, on_retry=_on_retry)
                        elif section == "results":
                            net = result.get("network",{})
                            hg = net.get("hct"); pg = net.get("ppi")
                            stats = dict(
                                hct_nodes=hg.G.number_of_nodes() if hg else 0,
                                hct_edges=hg.G.number_of_edges() if hg else 0,
                                ppi_nodes=pg.G.number_of_nodes() if pg else 0,
                                ppi_edges=pg.G.number_of_edges() if pg else 0,
                                total_genes=len(result.get("genes",[])),
                            )
                            st.session_state.ai_output[section] = draft_results(
                                herbs, drugs, stats, hubs, kegg, key, on_retry=_on_retry)
                        elif section == "discussion":
                            interp = st.session_state.ai_output.get("interpretation","")
                            st.session_state.ai_output[section] = draft_discussion(
                                herbs, drugs, hubs, kegg, interp, key, on_retry=_on_retry)
                        _st.update(label="Gemini 생성 완료", state="complete")
                    st.rerun()

        st.divider()
        full = []
        for sec, lbl in [("interpretation","🔍 작용기전 해석"),("methods","📝 Methods"),
                          ("results","📝 Results"),("discussion","📝 Discussion")]:
            txt = st.session_state.ai_output.get(sec,"")
            if txt:
                st.subheader(lbl); st.markdown(txt); st.divider()
                full.append(f"## {lbl}\n\n{txt}\n")

        if full:
            name = (result.get("config") or {}).get("name","analysis")
            st.download_button("📄 AI 초안 전체 다운로드",
                               "\n".join(full).encode(),
                               f"{name}_ai_draft.md","text/markdown")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 9: 내보내기
# ═══════════════════════════════════════════════════════════════════════════════
with tab_export:
    st.header("📦 투고용 데이터 패키징")
    result = st.session_state.pipeline_result

    if not result:
        st.info("분석을 먼저 실행하세요.")
    else:
        cfg  = result.get("config",{})
        name = cfg.get("name","analysis")

        st.subheader("Excel Supplementary Data (.xlsx)")
        st.caption("모든 분석 결과를 시트별로 패키징합니다.")

        col_ex, col_svg = st.columns(2)
        with col_ex:
            if st.button("📊 Excel 생성", type="primary", use_container_width=True, key="btn_excel_gen"):
                from exporter import SupplementaryExporter, MetadataLogger
                supp = SupplementaryExporter()
                meta = MetadataLogger()
                with st.spinner("Excel 생성 중..."):
                    xl_path = supp.build_excel(
                        prefix          = name,
                        admet_df        = (st.session_state.admet_results or {}).get("summary"),
                        topology_df     = st.session_state.topology_df,
                        enrich_results  = (result.get("enrichment") or {}).get("results"),
                        intersection    = st.session_state.intersection,
                        disease_targets = st.session_state.disease_targets or None,
                        metadata        = meta,
                    )
                st.session_state.export_paths["excel"] = xl_path
                st.rerun()

            xl = st.session_state.export_paths.get("excel")
            if xl and Path(xl).exists():
                st.download_button("⬇ 다운로드 (.xlsx)",
                                   Path(xl).read_bytes(), Path(xl).name,
                                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        with col_svg:
            if st.button("🖼 SVG 벡터 이미지 생성", use_container_width=True, key="btn_svg_gen"):
                from exporter import FigureExporter
                fe = FigureExporter()
                enrich_r = (result.get("enrichment") or {}).get("results", {})
                with st.spinner("SVG 생성 중..."):
                    svgs = fe.export_enrichment_svg(enrich_r, name)
                st.session_state.export_paths["svgs"] = [str(p) for p in svgs]
                st.success(f"SVG {len(svgs)}개 생성 완료")
                st.rerun()

            svgs = st.session_state.export_paths.get("svgs",[])
            for sp in svgs:
                p = Path(sp)
                if p.exists():
                    st.download_button(f"⬇ {p.name}", p.read_bytes(), p.name, "image/svg+xml",
                                       key=f"svg_{p.name}")

        # ── 논문용 멀티패널 Figure ─────────────────────────────────────────────
        st.divider()
        st.subheader("🖼 논문용 통합 멀티패널 Figure (A/B/C/D)")
        st.caption("300 DPI PNG + SVG로 저장. 개별 그래프를 자동 탐색합니다.")

        with st.expander("패널 경로 설정 (자동 탐색 / 수동 지정)"):
            en_dir = BASE_DIR / "enrichment"
            re_dir = BASE_DIR / "results"
            dk_dir = BASE_DIR / "docking"

            def _find(folder, patterns):
                for pat in patterns:
                    hits = sorted(folder.glob(pat))
                    if hits: return str(hits[0])
                return ""

            path_A = st.text_input("Panel A — Venn/교차",
                _find(re_dir, [f"{name}_venn.png","*_venn.png"]))
            path_B = st.text_input("Panel B — 허브 유전자 랭킹",
                _find(re_dir, [f"{name}_hub_ranking.png","*_hub_ranking.png"]))
            path_C = st.text_input("Panel C — GO/KEGG Dot Plot",
                _find(en_dir, [f"{name}_combined_dotplot.png","*_combined_dotplot.png"]))

            # 도킹 상호작용 wheel 자동 탐색
            iresults = st.session_state.interaction_results
            best_wheel = ""
            if iresults:
                best = min(iresults, key=lambda x: x.get("affinity") or 0)
                best_wheel = str(best.get("wheel_png",""))
            path_D = st.text_input("Panel D — 도킹 상호작용 Wheel / 히트맵",
                best_wheel or _find(dk_dir, ["*_interaction_wheel.png","*_wheel.png"]))

            fig_title = st.text_input("Figure 제목 (선택)",
                f"Network Pharmacology Analysis — {name}")
            fig_dpi   = st.slider("DPI", 150, 600, 300, 50)

        if st.button("🖼 통합 Figure 생성", type="primary", use_container_width=True, key="btn_fig_gen"):
            from figure_composer import make_publication_figure
            with st.spinner("멀티패널 Figure 생성 중..."):
                out = make_publication_figure(
                    venn_path    = path_A or None,
                    hub_path     = path_B or None,
                    dotplot_path = path_C or None,
                    docking_path = path_D or None,
                    title        = fig_title,
                    prefix       = name,
                    dpi          = fig_dpi,
                    output_dir   = re_dir,
                )
            st.session_state.export_paths["multipanel_png"] = str(out.get("png",""))
            st.session_state.export_paths["multipanel_svg"] = str(out.get("svg",""))
            st.success(f"생성 완료 ({out.get('panels_used',0)}/4 패널 사용)")
            st.rerun()

        for key_e, label_e, mime_e in [
            ("multipanel_png", "⬇ PNG (논문용 고해상도)", "image/png"),
            ("multipanel_svg", "⬇ SVG (벡터, 편집용)",   "image/svg+xml"),
        ]:
            fp_e = st.session_state.export_paths.get(key_e,"")
            if fp_e and Path(fp_e).exists():
                col1, col2 = st.columns([2,1])
                with col1:
                    if mime_e == "image/png":
                        st.image(fp_e, use_container_width=True)
                with col2:
                    st.download_button(label_e, Path(fp_e).read_bytes(),
                                       Path(fp_e).name, mime_e, key=f"dl_{key_e}")

        # ── 상호작용 데이터 Excel 통합 ────────────────────────────────────────
        st.divider()
        st.subheader("📊 상호작용 분석 Excel 포함 패키징")
        if st.button("📊 상호작용 포함 Excel 재생성", use_container_width=True, key="btn_excel_interact"):
            from exporter import SupplementaryExporter, MetadataLogger
            iresults = st.session_state.interaction_results
            with st.spinner("Excel 생성 중..."):
                # 상호작용 요약 DataFrame 합치기
                contact_dfs = []
                for ir in iresults:
                    cdf = ir.get("contacts_df")
                    if cdf is not None and not cdf.empty:
                        cdf = cdf.copy()
                        cdf.insert(0, "compound", ir["compound"])
                        cdf.insert(1, "target",   ir["target"])
                        cdf.insert(2, "affinity", ir.get("affinity",""))
                        contact_dfs.append(cdf)
                combined_contacts = pd.concat(contact_dfs, ignore_index=True) \
                                    if contact_dfs else None

                supp = SupplementaryExporter()
                meta = MetadataLogger()
                xl_path = supp.build_excel(
                    prefix          = f"{name}_full",
                    admet_df        = (st.session_state.admet_results or {}).get("summary"),
                    topology_df     = st.session_state.topology_df,
                    enrich_results  = (_active_result() or {}).get("enrichment",{}).get("results"),
                    intersection    = st.session_state.intersection,
                    disease_targets = st.session_state.disease_targets or None,
                    metadata        = meta,
                )
                # 접촉 잔기 시트 추가
                if combined_contacts is not None:
                    import openpyxl
                    wb = openpyxl.load_workbook(xl_path)
                    ws = wb.create_sheet("Docking_Interactions")
                    headers = list(combined_contacts.columns)
                    ws.append(headers)
                    for row in combined_contacts.itertuples(index=False):
                        ws.append(list(row))
                    wb.save(xl_path)
            st.session_state.export_paths["excel_full"] = str(xl_path)
            st.rerun()

        fp_xl = st.session_state.export_paths.get("excel_full","")
        if fp_xl and Path(fp_xl).exists():
            st.download_button("⬇ 전체 Supplementary Excel",
                               Path(fp_xl).read_bytes(), Path(fp_xl).name,
                               "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                               key="dl_excel_full")

        # ── 캐시 통계 ─────────────────────────────────────────────────────────
        st.divider()
        st.subheader("📋 API 캐시 & 재현성 메타데이터")
        col_c1, col_c2 = st.columns(2)
        with col_c1:
            if st.button("🗄 캐시 통계 조회", key="btn_cache_stats"):
                from collect import get_cache_stats, DB_PATH
                st.session_state.cache_stats = get_cache_stats(DB_PATH)

            cs = st.session_state.cache_stats
            if cs:
                if "error" in cs:
                    st.error(cs["error"])
                else:
                    st.metric("전체 캐시 항목", cs.get("total_entries",0))
                    st.metric("DB 크기 (KB)", cs.get("db_size_kb",0))
                    st.caption(f"최초 기록: {cs.get('oldest_entry','')}")
                    for src, cnt in (cs.get("by_source") or {}).items():
                        st.write(f"- {src}: {cnt}개")

        with col_c2:
            if st.button("📋 메타데이터 생성", key="btn_meta_gen"):
                from exporter import MetadataLogger
                meta   = MetadataLogger()
                report = meta.to_text_report()
                mp     = BASE_DIR/"results"/f"{name}_metadata.txt"
                mp.write_text(report, encoding="utf-8")
                st.session_state.export_paths["metadata"] = str(mp)
                st.rerun()

            mp = st.session_state.export_paths.get("metadata")
            if mp and Path(mp).exists():
                txt = Path(mp).read_text(encoding="utf-8")
                st.code(txt[:2000], language="text")
                st.download_button("⬇ 메타데이터", txt.encode(), Path(mp).name, "text/plain")

        st.divider()
        st.subheader("📂 전체 결과 파일 목록")
        for folder, label in [
            (BASE_DIR/"networks",   "Networks"),
            (BASE_DIR/"enrichment", "Enrichment"),
            (BASE_DIR/"results",    "Results"),
            (BASE_DIR/"docking",    "Docking"),
        ]:
            files = sorted(folder.glob(f"{name}_*")) if folder.exists() else []
            if files:
                with st.expander(f"{label} ({len(files)}개 파일)"):
                    for fp in files:
                        st.download_button(f"⬇ {fp.name}", fp.read_bytes(),
                                           fp.name, key=f"all_{fp.name}")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 10: 프로젝트 관리
# ═══════════════════════════════════════════════════════════════════════════════
with tab_projects:
    st.header("📁 프로젝트 관리")
    st.caption("분석 프로젝트별 파일 현황 확인 및 삭제로 디스크 용량을 관리합니다.")

    PROJECT_DIRS = {
        "configs":    BASE_DIR / "configs",
        "results":    BASE_DIR / "results",
        "enrichment": BASE_DIR / "enrichment",
        "networks":   BASE_DIR / "networks",
        "docking":    BASE_DIR / "docking",
    }

    def _get_projects() -> list[dict]:
        cfg_dir = BASE_DIR / "configs"
        if not cfg_dir.exists():
            return []
        projects = []
        for cfg_file in sorted(cfg_dir.glob("*.json"), key=lambda p: -p.stat().st_mtime):
            name = cfg_file.stem
            total_bytes = cfg_file.stat().st_size
            file_count  = 1
            file_list   = [cfg_file]
            for folder in [BASE_DIR/"results", BASE_DIR/"enrichment",
                           BASE_DIR/"networks", BASE_DIR/"docking"]:
                if folder.exists():
                    for f in folder.glob(f"{name}_*"):
                        total_bytes += f.stat().st_size
                        file_count  += 1
                        file_list.append(f)
            try:
                cfg_data = json.loads(cfg_file.read_text())
            except Exception:
                cfg_data = {}
            projects.append({
                "name":        name,
                "cfg_file":    cfg_file,
                "size_mb":     round(total_bytes / 1024 / 1024, 2),
                "file_count":  file_count,
                "file_list":   file_list,
                "herbs":       len(cfg_data.get("herbs", {})),
                "drugs":       len(cfg_data.get("drugs", [])),
                "disease":     cfg_data.get("disease_query", "—") or "—",
                "mtime":       cfg_file.stat().st_mtime,
            })
        return projects

    if st.button("🔄 새로고침", key="proj_refresh"):
        st.rerun()

    projects = _get_projects()

    if not projects:
        st.info("저장된 프로젝트가 없습니다. AI 설정 탭에서 분석을 실행하면 프로젝트가 생성됩니다.")
    else:
        # 전체 통계
        total_size = sum(p["size_mb"] for p in projects)
        c1, c2, c3 = st.columns(3)
        c1.metric("총 프로젝트", len(projects))
        c2.metric("총 사용 용량", f"{total_size:.1f} MB")
        c3.metric("평균 프로젝트 크기", f"{total_size/len(projects):.1f} MB")

        st.divider()

        for proj in projects:
            import datetime
            mtime_str = datetime.datetime.fromtimestamp(proj["mtime"]).strftime("%Y-%m-%d %H:%M")
            with st.expander(
                f"**{proj['name']}** — {proj['size_mb']} MB · {proj['file_count']}개 파일 · {mtime_str}",
                expanded=False
            ):
                col_info, col_act = st.columns([3, 1])
                with col_info:
                    st.markdown(
                        f"- 생약 **{proj['herbs']}**종 · 약물 **{proj['drugs']}**종 · 질환: **{proj['disease']}**\n"
                        f"- 파일 수: **{proj['file_count']}**개 · 용량: **{proj['size_mb']} MB**"
                    )
                    with st.expander("파일 목록 및 다운로드"):
                        _MIME = {
                            ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            ".png":  "image/png",
                            ".svg":  "image/svg+xml",
                            ".md":   "text/markdown",
                            ".csv":  "text/csv",
                            ".json": "application/json",
                            ".pkl":  "application/octet-stream",
                        }
                        for f in sorted(proj["file_list"]):
                            size_kb = f.stat().st_size // 1024
                            ext = f.suffix.lower()
                            fc1, fc2, fc3 = st.columns([4, 1, 1])
                            fc1.caption(f"`{f.parent.name}/{f.name}` — {size_kb} KB")
                            with fc2:
                                st.download_button(
                                    "⬇", f.read_bytes(), f.name,
                                    _MIME.get(ext, "application/octet-stream"),
                                    key=f"dl_{proj['name']}_{f.name}",
                                    use_container_width=True,
                                )
                            with fc3:
                                if ext == ".png":
                                    if st.button("👁", key=f"prev_{proj['name']}_{f.name}",
                                                 use_container_width=True):
                                        st.session_state[f"_preview_{proj['name']}_{f.name}"] = \
                                            not st.session_state.get(f"_preview_{proj['name']}_{f.name}", False)
                                elif ext == ".md":
                                    if st.button("👁", key=f"prev_{proj['name']}_{f.name}",
                                                 use_container_width=True):
                                        st.session_state[f"_preview_{proj['name']}_{f.name}"] = \
                                            not st.session_state.get(f"_preview_{proj['name']}_{f.name}", False)
                                else:
                                    st.write("")
                            # 미리보기 토글
                            _pkey = f"_preview_{proj['name']}_{f.name}"
                            if st.session_state.get(_pkey):
                                if ext == ".png":
                                    st.image(str(f), use_container_width=True)
                                elif ext == ".md":
                                    st.markdown(f.read_text(encoding="utf-8"))

                with col_act:
                    # 현재 활성 프로젝트로 로드
                    if st.button("📂 불러오기", key=f"load_{proj['name']}",
                                 use_container_width=True):
                        loaded = json.loads(proj["cfg_file"].read_text())
                        st.session_state.current_cfg         = loaded
                        st.session_state.ai_proposed_config  = loaded
                        st.session_state.ai_setup_done       = True
                        _restored = False
                        # 원클릭 전체 결과 복원 (우선)
                        _fullcache = BASE_DIR / "results" / f"{proj['name']}_fullrun_cache.pkl"
                        if _fullcache.exists():
                            try:
                                _fr = pickle.loads(_fullcache.read_bytes())
                                if _fr.get("pipeline"):
                                    st.session_state.pipeline_result = _fr["pipeline"]
                                if _fr.get("admet"):
                                    st.session_state.admet_results = _fr["admet"]
                                if _fr.get("admet_filtered"):
                                    st.session_state.admet_filtered_result = _fr["admet_filtered"]
                                if _fr.get("topology") is not None:
                                    st.session_state.topology_df = _fr["topology"]
                                if _fr.get("intersection"):
                                    st.session_state.intersection    = _fr["intersection"]
                                    st.session_state.disease_targets = _fr.get("disease_targets", [])
                                if _fr.get("docking") is not None:
                                    st.session_state.docking_df = _fr["docking"]
                                if _fr.get("ai_output"):
                                    st.session_state.ai_output = _fr["ai_output"]
                                st.session_state.full_run_result = _fr
                                st.session_state.full_run_state  = "done"
                                _restored = True
                                st.success(f"'{proj['name']}' 로드됨 (전체 분석 결과 복원)")
                            except Exception:
                                pass
                        # fallback: 기본 pipeline 캐시
                        if not _restored:
                            _cache = BASE_DIR / "results" / f"{proj['name']}_pipeline_cache.pkl"
                            if _cache.exists():
                                try:
                                    _res = pickle.loads(_cache.read_bytes())
                                    st.session_state.pipeline_result = _res
                                    st.success(f"'{proj['name']}' 로드됨 (기본 분석 결과 복원)")
                                except Exception:
                                    st.success(f"'{proj['name']}' 설정만 로드됨 — 재실행 필요")
                            else:
                                st.info(f"'{proj['name']}' 설정 로드됨 — 재실행 필요")
                        st.rerun()

                    # 삭제 확인 상태 관리
                    confirm_key = f"confirm_del_{proj['name']}"
                    if st.session_state.get(confirm_key):
                        st.error("⚠️ 정말 삭제하시겠습니까?")
                        col_y, col_n = st.columns(2)
                        with col_y:
                            if st.button("✅ 확인", key=f"yes_{proj['name']}",
                                         use_container_width=True):
                                for f in proj["file_list"]:
                                    try: f.unlink()
                                    except Exception: pass
                                st.session_state[confirm_key] = False
                                st.success(f"'{proj['name']}' 삭제 완료")
                                st.rerun()
                        with col_n:
                            if st.button("❌ 취소", key=f"no_{proj['name']}",
                                         use_container_width=True):
                                st.session_state[confirm_key] = False
                                st.rerun()
                    else:
                        if st.button("🗑 삭제", key=f"del_{proj['name']}",
                                     use_container_width=True, type="secondary"):
                            st.session_state[confirm_key] = True
                            st.rerun()


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 11: 설정
# ═══════════════════════════════════════════════════════════════════════════════
with tab_settings:
    st.header("⚙️ 설정")

    st.subheader("🔑 Google AI Studio API 키")
    st.caption(
        "Google AI Studio (aistudio.google.com) 에서 발급. "
        "Gemma 4 31B (무료) 및 Gemini 2.5 Flash 사용에 필요합니다."
    )

    api_key_input = st.text_input(
        "API Key",
        value=st.session_state.gemini_key,
        type="password",
        placeholder="AIzaSy...",
        label_visibility="collapsed",
    )
    col_save, col_clear = st.columns([2, 1])
    with col_save:
        if st.button("💾 저장", type="primary", use_container_width=True, key="btn_settings_save"):
            st.session_state.gemini_key = api_key_input
            os.environ["GEMINI_API_KEY"] = api_key_input
            save_settings({
                "gemini_api_key": api_key_input,
                "setup_model":    st.session_state.setup_model,
                "ai_model":       st.session_state.ai_model,
                "ppi_min_score":  ppi_score,
                "fdr_cutoff":     fdr_cutoff,
                "top_n":          top_n,
            })
            st.success("✅ 저장되었습니다. 앱을 재시작해도 유지됩니다.")
    with col_clear:
        if st.button("🗑 초기화", use_container_width=True, key="btn_apikey_clear"):
            st.session_state.gemini_key = ""
            os.environ.pop("GEMINI_API_KEY", None)
            save_settings({"gemini_api_key": ""})
            st.info("API 키가 삭제되었습니다.")
            st.rerun()

    st.divider()
    st.subheader("🤖 AI 모델 설정")

    col_m1, col_m2 = st.columns(2)
    with col_m1:
        st.markdown("**💬 AI 설정 파서 (연구 내용 인식)**")
        st.info(f"모델: `{st.session_state.setup_model}`")
        st.caption("Gemma 4 31B — 무료, 자연어 연구 설명 → 설정 자동 추출")
    with col_m2:
        st.markdown("**📝 AI 해석기 (논문 초안 작성)**")
        st.info(f"모델: `{st.session_state.ai_model}`")
        st.caption("Gemini 2.5 Flash — Methods/Results/Discussion 자동 생성")

    st.divider()
    st.subheader("📊 기본 분석 파라미터")
    st.caption("사이드바 슬라이더와 동일. 여기서 저장하면 다음 실행 시 기본값으로 적용됩니다.")
    col_p1, col_p2, col_p3 = st.columns(3)
    with col_p1:
        st.metric("STRING 신뢰도 임계값", ppi_score)
    with col_p2:
        st.metric("FDR 컷오프", fdr_cutoff)
    with col_p3:
        st.metric("표시 term 수", top_n)

    if st.button("📊 현재 파라미터를 기본값으로 저장", use_container_width=True, key="btn_param_save"):
        save_settings({
            "gemini_api_key": st.session_state.gemini_key,
            "setup_model":    st.session_state.setup_model,
            "ai_model":       st.session_state.ai_model,
            "ppi_min_score":  ppi_score,
            "fdr_cutoff":     fdr_cutoff,
            "top_n":          top_n,
        })
        st.session_state.ppi_score_saved = ppi_score
        st.session_state.fdr_saved       = fdr_cutoff
        st.session_state.top_n_saved     = top_n
        st.success("저장 완료")

    st.divider()
    st.subheader("🗄 캐시 관리")
    from collect import get_cache_stats
    stats = get_cache_stats()
    if "error" not in stats:
        col_s1, col_s2, col_s3 = st.columns(3)
        col_s1.metric("캐시 항목", stats.get("total_entries", 0))
        col_s2.metric("DB 크기", f"{stats.get('db_size_kb', 0):.0f} KB")
        col_s3.metric("최초 수집일", str(stats.get("oldest_entry", "—"))[:10])
        with st.expander("소스별 상세"):
            for src, n in (stats.get("by_source") or {}).items():
                st.write(f"- **{src}**: {n}개")
    if st.button("🗑 캐시 전체 삭제", type="secondary", key="btn_cache_clear"):
        db_path = BASE_DIR / "db" / "cache.db"
        if db_path.exists():
            import sqlite3
            conn = sqlite3.connect(db_path)
            conn.execute("DELETE FROM cache")
            conn.commit()
            conn.close()
            st.success("캐시 삭제 완료 (다음 실행 시 API를 새로 호출합니다)")


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 11: DB 관리 (TCMSP 캐시 빌드 모니터)
# ═══════════════════════════════════════════════════════════════════════════════
with tab_dbmgr:
    st.header("🗄️ TCMSP 로컬 DB 관리")
    st.caption("TCMSP 500종 약재 데이터를 로컬에 캐시합니다. 한 번 구축하면 분석이 수초 내 완료됩니다.")

    STATUS_FILE = BASE_DIR / "logs" / "build_tcmsp_cache_status.json"
    ERRORS_FILE = BASE_DIR / "logs" / "build_tcmsp_cache_errors.json"
    CACHE_DIR   = BASE_DIR / "data" / "tcmsp_cache"
    HERB_LIST_F = BASE_DIR / "data" / "tcmsp_herb_list.json"

    # ── 현재 캐시 현황 ───────────────────────────────────────────────────────
    cached_files = list(CACHE_DIR.glob("*.json")) if CACHE_DIR.exists() else []
    n_cached = len(cached_files)
    n_total  = 500

    col1, col2, col3 = st.columns(3)
    col1.metric("전체 약재", f"{n_total}종")
    col2.metric("캐시 완료", f"{n_cached}종", delta=f"+{n_cached}" if n_cached else None)
    col3.metric("남은 약재", f"{n_total - n_cached}종")

    st.progress(n_cached / n_total, text=f"{n_cached}/{n_total} ({n_cached/n_total*100:.1f}%)")

    st.divider()

    # ── 빌드 상태 ────────────────────────────────────────────────────────────
    st.subheader("빌드 진행 상황")

    status = {}
    if STATUS_FILE.exists():
        try:
            status = json.loads(STATUS_FILE.read_text())
        except Exception:
            pass

    build_status = status.get("status", "idle")

    if build_status == "running":
        st.info(f"🔄 빌드 진행 중 (PID: {status.get('pid', '?')})")
        proc = status.get("processed", 0)
        rem  = status.get("remaining", n_total)
        if rem > 0:
            st.progress(proc / rem, text=f"처리: {proc}/{rem}")
        cur = status.get("current_herb", "")
        if cur:
            st.caption(f"현재 처리 중: **{cur}**")
        eta = status.get("eta_minutes")
        if eta is not None:
            h, m = divmod(int(eta), 60)
            st.caption(f"예상 남은 시간: {h}시간 {m}분")
        cols = st.columns(3)
        cols[0].metric("성공", status.get("success", 0))
        cols[1].metric("실패", status.get("fail", 0))
        cols[2].metric("마지막 업데이트", status.get("last_updated", "")[-8:] if status.get("last_updated") else "-")
        if st.button("🔄 새로고침", key="btn_db_refresh"):
            st.rerun()

    elif build_status == "completed":
        st.success(f"✅ 빌드 완료! 성공 {status.get('success',0)}개 / 실패 {status.get('fail',0)}개")

    else:
        st.info("대기 중 — 아직 빌드가 시작되지 않았습니다.")

    # ── 빌드 시작 버튼 ───────────────────────────────────────────────────────
    st.divider()
    st.subheader("빌드 실행")

    if build_status == "running":
        st.warning("이미 빌드가 실행 중입니다. 터미널에서 `tail -f logs/build_tcmsp_cache.log` 로 확인하세요.")
    else:
        st.markdown("""
터미널에서 아래 명령으로 백그라운드 실행하세요. 앱/세션을 닫아도 계속 진행됩니다.
```bash
cd /home/beanalogue/netpharm
nohup python3 build_tcmsp_cache.py > logs/build_tcmsp_cache.log 2>&1 &
echo $!
```
        """)

    # ── 실패 목록 ────────────────────────────────────────────────────────────
    if ERRORS_FILE.exists():
        st.divider()
        st.subheader("실패 약재 목록")
        try:
            errors = json.loads(ERRORS_FILE.read_text())
            if errors:
                err_df = pd.DataFrame(errors)[["pinyin", "cn_name", "error", "timestamp"]]
                err_df.columns = ["병음", "한자", "오류", "시각"]
                st.dataframe(err_df, use_container_width=True, hide_index=True)
            else:
                st.success("실패 없음")
        except Exception:
            st.warning("실패 로그를 읽을 수 없습니다.")

    # ── 캐시된 약재 목록 ─────────────────────────────────────────────────────
    st.divider()
    st.subheader(f"캐시된 약재 ({n_cached}종)")
    if HERB_LIST_F.exists():
        herb_db_data = json.loads(HERB_LIST_F.read_text())
        cached_names = {f.stem.replace("_", " ") for f in cached_files}
        rows = []
        for en_name, hdata in herb_db_data.items():
            py  = hdata.get("herb_pinyin", "")
            cn  = hdata.get("herb_cn_name", "")
            is_cached = py in cached_names or py.replace(" ", "_") in {f.stem for f in cached_files}
            rows.append({"병음": py, "한자": cn, "약전명": en_name, "캐시": "✅" if is_cached else "⬜"})
        herb_table = pd.DataFrame(rows)
        col_filter = st.columns([1, 3])
        show_all   = col_filter[0].checkbox("전체 보기", value=False, key="db_show_all")
        search_kw  = col_filter[1].text_input("검색", placeholder="병음 또는 한자 이름", key="db_search", label_visibility="collapsed")
        if not show_all:
            herb_table = herb_table[herb_table["캐시"] == "✅"]
        if search_kw:
            mask = herb_table["병음"].str.contains(search_kw, case=False) | herb_table["한자"].str.contains(search_kw, case=False)
            herb_table = herb_table[mask]
        st.dataframe(herb_table, use_container_width=True, hide_index=True, height=400)


# ═══════════════════════════════════════════════════════════════════════════════
# TAB 13: 약재 발굴 모드 (Reverse Network Pharmacology)
# ═══════════════════════════════════════════════════════════════════════════════
with tab_discover:
    st.header("🔍 약재 발굴 모드")
    st.caption(
        "질환/건강상태를 입력하면 TCMSP DB 전체 약재를 스코어링하여 "
        "잠재적으로 효과적인 약재를 랭킹합니다."
    )

    # ── 캐시 현황 알림 ───────────────────────────────────────────────────────
    _disc_cache_dir = BASE_DIR / "data" / "tcmsp_cache"
    _disc_n_cached  = len(list(_disc_cache_dir.glob("*.json"))) if _disc_cache_dir.exists() else 0
    if _disc_n_cached < 400:
        st.warning(
            f"TCMSP 캐시가 {_disc_n_cached}/500종만 완성되어 있습니다. "
            "캐시가 충분히 쌓일수록 발굴 결과의 커버리지가 높아집니다. "
            "지금도 사용 가능합니다."
        )
    else:
        st.success(f"TCMSP 캐시 {_disc_n_cached}/500종 준비됨 ✅")

    st.divider()

    # ── 입력 폼 ─────────────────────────────────────────────────────────────
    with st.form("discover_form"):
        _disc_disease = st.text_input(
            "질환 / 건강 상태",
            placeholder="예: fatigue, hypertension, type 2 diabetes, depression",
            help="영어로 입력하면 OpenTargets 검색 정확도가 높습니다.",
        )
        _col1, _col2, _col3, _col4 = st.columns(4)
        _disc_top_n      = _col1.number_input("상위 약재 수",   min_value=5,  max_value=50, value=20, step=5)
        _disc_min_ol     = _col2.number_input("최소 오버랩",    min_value=1,  max_value=10, value=2,  step=1,
                                              help="질환 타깃과 공유 유전자가 이 수 이상인 약재만 표시")
        _disc_ob         = _col3.number_input("OB 기준 (%)",    min_value=10, max_value=60, value=30, step=5,
                                              help="경구 생체이용률 필터")
        _disc_dl         = _col4.number_input("DL 기준",        min_value=0.10, max_value=0.50, value=0.18, step=0.01,
                                              help="약물 유사성 필터")
        _col_ppi1, _col_ppi2 = st.columns([1, 3])
        _disc_use_ppi = _col_ppi1.toggle(
            "🔗 PPI 시너지 분석",
            value=False,
            help="STRING DB를 이용해 약재 조합 간 단백질 상호작용 시너지를 계산합니다. 조합 분석 필요, 추가 시간 소요."
        )
        if _disc_use_ppi:
            _col_ppi2.caption("STRING DB (string-db.org) 연동 — 약재 타깃 간 PPI 강도로 조합 우선순위 보정. 2-herb 조합에 적용됩니다.")
        _disc_run = st.form_submit_button("🔍 발굴 시작", use_container_width=True, type="primary")

    # ── 실행 ────────────────────────────────────────────────────────────────
    if _disc_run:
        if not _disc_disease.strip():
            st.error("질환명을 입력해주세요.")
        else:
            _disc_result = None
            with st.status("약재 발굴 분석 중...", expanded=True) as _disc_st:
                _disc_msg = st.empty()

                def _disc_progress(msg: str):
                    _disc_msg.write(f"⏳ {msg}")
                    _disc_st.update(label=msg)

                try:
                    from herb_discovery import run_discovery
                    _disc_result = run_discovery(
                        disease_name          = _disc_disease.strip(),
                        top_n                 = int(_disc_top_n),
                        ob_cutoff             = float(_disc_ob),
                        dl_cutoff             = float(_disc_dl),
                        min_overlap           = int(_disc_min_ol),
                        combo_sizes           = [2] if _disc_use_ppi else None,
                        top_n_single_for_combo= 20,
                        use_ppi               = _disc_use_ppi,
                        ppi_top_n             = 10,
                        on_progress           = _disc_progress,
                    )
                    if "error" in _disc_result:
                        _disc_st.update(label="오류 발생", state="error")
                        st.error(_disc_result["error"])
                        _disc_result = None
                    else:
                        _disc_st.update(label="발굴 완료 ✅", state="complete")
                except Exception as _disc_e:
                    _disc_st.update(label="오류 발생", state="error")
                    st.error(f"오류: {_disc_e}")

            # ── 결과 표시 ─────────────────────────────────────────────────
            if _disc_result:
                herbs_found = _disc_result["herb_results"]
                st.markdown(
                    f"**{_disc_result['disease_name']}** | "
                    f"질환 타깃: {_disc_result['disease_count']}개 | "
                    f"분석 약재: {_disc_result['cache_count']}개 | "
                    f"조건 충족: **{len(herbs_found)}개**"
                )

                if not herbs_found:
                    st.warning("조건을 충족하는 약재가 없습니다. 최소 오버랩 값을 낮춰보세요.")
                else:
                    # 요약 테이블
                    _disc_rows = []
                    for rank, h in enumerate(herbs_found, 1):
                        _disc_rows.append({
                            "순위":       rank,
                            "약재(병음)": h["pinyin"],
                            "한자":       h.get("cn_name", ""),
                            "활성성분":   h["active_count"],
                            "오버랩":     h["overlap_count"],
                            "Coverage":   f"{h['coverage']*100:.1f}%",
                            "Score":      h["score"],
                            "주요 공유유전자": ", ".join(h["overlap_genes"][:5])
                                              + (f" 외 {len(h['overlap_genes'])-5}개"
                                                 if len(h["overlap_genes"]) > 5 else ""),
                        })
                    _disc_df = pd.DataFrame(_disc_rows)
                    st.dataframe(_disc_df, use_container_width=True, hide_index=True)

                    # 바 차트
                    import matplotlib
                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt

                    _fig, _ax = plt.subplots(figsize=(8, max(4, len(herbs_found) * 0.4)))
                    _names  = [f"{h['pinyin']}\n({h.get('cn_name','')})" for h in herbs_found]
                    _scores = [h["overlap_count"] for h in herbs_found]
                    _bars   = _ax.barh(_names[::-1], _scores[::-1],
                                       color="#4C72B0", alpha=0.8, edgecolor="white")
                    _ax.set_xlabel("Disease Target Overlap Count", fontsize=10)
                    _ax.set_title(
                        f"Herb Ranking for '{_disc_result['disease_name']}'\n"
                        f"(OB≥{int(_disc_ob)}%, DL≥{_disc_dl}, min overlap≥{int(_disc_min_ol)})",
                        fontsize=11, fontweight="bold"
                    )
                    _ax.tick_params(axis="y", labelsize=8)
                    plt.tight_layout()
                    st.pyplot(_fig)
                    plt.close(_fig)

                    # PPI 시너지 조합 결과
                    _disc_combos = _disc_result.get("combinations", {}).get(2, [])
                    if _disc_combos:
                        st.subheader("🔗 2-herb 조합 우선순위 (PPI 시너지 포함)" if _disc_use_ppi else "🔗 2-herb 조합 우선순위")
                        _combo_rows = []
                        for c in _disc_combos[:10]:
                            row = {
                                "조합": " + ".join(c["herbs"]),
                                "한자": " + ".join(c.get("cn_names", c.get("herbs", []))),
                                "타깃오버랩": c["overlap_count"],
                                "시너지(+타깃)": c["synergy"],
                            }
                            if _disc_use_ppi:
                                row["PPI교차수"]  = c.get("ppi_cross_n", 0)
                                row["PPI강도"]   = round(c.get("ppi_cross_score", 0), 2)
                                row["최종점수"]   = round(c.get("combined_score", 0), 4)
                                top_p = c.get("top_pairs", [])
                                row["주요 PPI"] = f"{top_p[0][0]}↔{top_p[0][1]}({top_p[0][2]:.2f})" if top_p else "-"
                            _combo_rows.append(row)
                        st.dataframe(pd.DataFrame(_combo_rows), use_container_width=True, hide_index=True)
                        if _disc_use_ppi:
                            st.caption("PPI 강도: STRING DB 단백질 상호작용 점수 합계 (높을수록 약재 간 타깃이 생물학적으로 연결됨)")

                    # 약재별 상세 (접이식)
                    st.subheader("약재별 공유 유전자 상세")
                    for h in herbs_found[:10]:
                        with st.expander(
                            f"**{h['pinyin']}** ({h.get('cn_name','')}) — "
                            f"오버랩 {h['overlap_count']}개 | score {h['score']:.4f}"
                        ):
                            st.write(f"**활성 성분 ({h['active_count']}개)**:")
                            st.write(", ".join(h["active_compounds"][:10])
                                     + (f" 외 {h['active_count']-10}개" if h["active_count"] > 10 else ""))
                            st.write(f"**공유 유전자 ({h['overlap_count']}개)**:")
                            st.write(", ".join(h["overlap_genes"]))

                    # AI 해석
                    st.divider()
                    _settings_d = load_settings()
                    _api_key_d  = _settings_d.get("gemini_api_key", "")
                    if not _api_key_d:
                        st.info("AI 메커니즘 해석을 보려면 사이드바에서 Gemini API 키를 입력하세요.")
                    else:
                        if st.button("🤖 AI 메커니즘 해석 생성", key="disc_ai_btn"):
                            with st.status("Gemini 생성 중...", expanded=True) as _disc_ai_st:
                                _disc_ai_msg = st.empty()
                                _disc_ai_msg.write("Gemini API 호출 중...")

                                def _disc_on_retry(attempt, wait_sec, err_msg):
                                    _disc_ai_st.update(label=f"Gemini 재시도 중... ({attempt}/3)")
                                    _disc_ai_msg.warning(
                                        f"일시적 오류 (503/429) — {wait_sec}초 후 자동 재시도 ({attempt}/3)\n\n"
                                        f"`{err_msg[:120]}`"
                                    )

                                try:
                                    from interpreter import draft_discovery as _draft_disc
                                    _disc_ai_text = _draft_disc(
                                        disease_name  = _disc_result["disease_name"],
                                        top_herbs     = herbs_found,
                                        disease_genes = _disc_result["disease_genes"],
                                        api_key       = _api_key_d,
                                        on_retry      = _disc_on_retry,
                                    )
                                    _disc_ai_st.update(label="Gemini 생성 완료 ✅", state="complete")
                                    st.markdown(_disc_ai_text)
                                except Exception as _disc_ai_e:
                                    _disc_ai_st.update(label="오류", state="error")
                                    st.error(f"AI 해석 오류: {_disc_ai_e}")

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 14: 논문 작성
# ═══════════════════════════════════════════════════════════════════════════════
with tab_paper:
    import os
    from pathlib import Path

    PAPER_DIR = Path("paper")
    PAPER_DIR.mkdir(exist_ok=True)
    DRAFT_PATH = PAPER_DIR / "draft_v1.md"

    st.header("📝 논문 작성 / Paper Writing")
    st.caption("Claude AI와 대화하며 논문을 작성하고, 편집기에서 직접 다듬을 수 있습니다.")

    # ── 서브탭 ───────────────────────────────────────────────────
    _ptab_chat, _ptab_edit = st.tabs(["💬 AI 채팅", "📄 편집 / 미리보기"])

    # ════════════════════════════════════════════════════════════
    # 서브탭 1: Claude 채팅
    # ════════════════════════════════════════════════════════════
    with _ptab_chat:

        # ── API 키 ───────────────────────────────────────────────
        _p_settings = load_settings()
        _p_ant_key  = _p_settings.get("anthropic_api_key", "")

        if not _p_ant_key:
            _p_ant_key = st.text_input(
                "Anthropic API 키",
                type="password",
                placeholder="sk-ant-...",
                help="Claude API 키를 입력하세요.",
                key="paper_ant_key_input",
            )
            if _p_ant_key:
                save_settings({"anthropic_api_key": _p_ant_key})
                st.rerun()
        else:
            st.caption(f"✅ Anthropic API 키 설정됨  |  모델: claude-sonnet-4-6")

        if not _p_ant_key:
            st.stop()

        # ── 논문 연구 컨텍스트 시스템 프롬프트 ──────────────────
        _PAPER_SYSTEM = """You are a scientific writing assistant specializing in computational pharmacology and traditional medicine (TM) research. You are helping write an academic paper with the following details:

## Paper Concept
Title: "A Reverse Network Pharmacology Pipeline for Molecular Evidence Enrichment in East Asian Traditional Medicine Clinical Practice Guidelines: A Multi-Disease Proof-of-Concept Study"

Core claim: In silico reverse network pharmacology can serve as a supplementary molecular evidence layer for TM Clinical Practice Guidelines (CPGs) — demonstrated with Korean TM CPGs across three diseases.

## Methodology
- Disease gene sets: OpenTargets (MONDO ontology), top 499–500 genes per disease
- Herb-compound-target: TCMSP database, 500 herbs, OB≥30% / DL≥0.18 filter
- Scoring: hypergeometric p-value + normalized overlap score
- Validation metric: Pool AUROC (Wilcoxon rank-sum + permutation test, n=1,000)
- Multiple testing correction: Bonferroni (α = 0.05/3 = 0.0167), 3 diseases
- Per-formula analysis: exploratory/descriptive only (no significance stars)

## Key Results
PRIMARY OUTCOME (Disease Pool AUROC, Bonferroni-corrected):
| Disease | MONDO | Genes | Pool AUROC | Permutation p | Significant |
|---------|-------|-------|-----------|--------------|-------------|
| Essential Hypertension | MONDO_0001134 | 500 | 0.6245 | 0.010 | YES (p<0.0167) |
| Insomnia | MONDO_0013600 | 500 | 0.6070 | 0.048 | borderline |
| Dementia | MONDO_0001627 | 499 | 0.5566 | 0.188 | NO |

EXPLORATORY (per-formula, descriptive only):
Insomnia top 3: Suanzaoren-tang AUROC=0.762, Guipi-tang 0.722, Wendan-tang 0.720
Hypertension: Qiju-Dihuang-wan 0.689, Tianma-Gouteng-yin 0.657 (Tianma excluded)
Dementia: all AUROC ≈ 0.48–0.57 (null pattern)

## Sensitivity Analysis (10 scenarios)
- Robust finding: HTN significant in 8/10 scenarios, dementia null in 9/10
- OB≥40 anomaly: removes low-ranking formula herbs → artificially inflates AUROC (documented)
- OT≥0.30 weakens HTN signal (only 202 genes) → OT≥0.10 pre-specified as baseline

## Key Limitations
1. TCMSP coverage: Tianma (Gastrodia) absent from TCMSP's 500 herbs; its main compound gastrodin also fails OB/DL filter → systematic exclusion of phenolic glycosides
2. Herb name mapping: Longdancao resolved to Longdan (confirmed in cache)
3. CPG composition: Korean CPGs list formula names only, not herb compositions → standard ChP 2020 used as proxy
4. Circular reasoning risk: TCMSP compound-target data may be biased toward traditional use indications
5. Dementia MONDO term: pan-dementia gene set is heterogeneous → reduces specificity

## Target Journal
Chinese Medicine (BMC), IF=5.8 — or Frontiers in Pharmacology (Network Pharmacology section)

## Writing Style
- Academic English
- Journal: Chinese Medicine (BMC) style
- Word count target: ~4,000–5,000 words (excl. references)
- Sections: Abstract, Introduction, Methods, Results, Discussion, Conclusion

When generating sections, embed actual numbers from the results above. Be scientifically accurate and honest about limitations. Write in a way that clearly defends the methodology against anticipated reviewer critiques."""

        # ── 세션 상태 초기화 ─────────────────────────────────────
        if "paper_chat_history" not in st.session_state:
            st.session_state["paper_chat_history"] = []

        # ── 상단 컨트롤 ──────────────────────────────────────────
        _pc_col1, _pc_col2, _pc_col3 = st.columns([2, 1, 1])
        with _pc_col1:
            _pc_section = st.selectbox(
                "빠른 섹션 생성",
                ["직접 입력", "Abstract", "Introduction", "Methods", "Results",
                 "Discussion", "Conclusion", "Limitations", "Figure legends"],
                key="paper_chat_section_sel",
                label_visibility="collapsed",
            )
        with _pc_col2:
            if st.button("🚀 섹션 요청", use_container_width=True, key="paper_chat_quick_btn"):
                if _pc_section != "직접 입력":
                    _quick_prompt = f"Please write the **{_pc_section}** section for this paper."
                    st.session_state["paper_chat_history"].append(
                        {"role": "user", "content": _quick_prompt}
                    )
                    st.session_state["paper_chat_trigger"] = True
                    st.rerun()
        with _pc_col3:
            if st.button("🗑️ 대화 초기화", use_container_width=True, key="paper_chat_clear_btn"):
                st.session_state["paper_chat_history"] = []
                st.rerun()

        # ── 채팅 히스토리 표시 ────────────────────────────────────
        _chat_container = st.container(height=520)
        with _chat_container:
            if not st.session_state["paper_chat_history"]:
                st.info(
                    "👋 논문 작성을 시작하세요.\n\n"
                    "예시:\n"
                    "- *'Abstract 써줘'*\n"
                    "- *'Methods 섹션을 한국어로 설명해줘'*\n"
                    "- *'Discussion에서 치매 null result를 어떻게 설명하면 좋을까?'*\n"
                    "- *'Reviewer가 다중검정 문제를 지적할 때 반박 논리 써줘'*"
                )
            for _msg in st.session_state["paper_chat_history"]:
                with st.chat_message(_msg["role"]):
                    st.markdown(_msg["content"])

        # ── 자동 응답 (빠른 섹션 생성 후) ─────────────────────────
        if st.session_state.get("paper_chat_trigger"):
            st.session_state.pop("paper_chat_trigger")
            try:
                from anthropic import Anthropic
                _ant_client = Anthropic(api_key=_p_ant_key)
                with _chat_container:
                    with st.chat_message("assistant"):
                        _resp_placeholder = st.empty()
                        _full_resp = ""
                        with _ant_client.messages.stream(
                            model="claude-sonnet-4-6",
                            max_tokens=4096,
                            system=_PAPER_SYSTEM,
                            messages=st.session_state["paper_chat_history"],
                        ) as _stream:
                            for _chunk in _stream.text_stream:
                                _full_resp += _chunk
                                _resp_placeholder.markdown(_full_resp + "▌")
                        _resp_placeholder.markdown(_full_resp)
                st.session_state["paper_chat_history"].append(
                    {"role": "assistant", "content": _full_resp}
                )
                st.session_state["paper_last_response"] = _full_resp
            except Exception as _e:
                st.error(f"Claude API 오류: {_e}")

        # ── 채팅 입력 ─────────────────────────────────────────────
        if _user_input := st.chat_input("논문에 대해 질문하거나 섹션 작성을 요청하세요...", key="paper_chat_input"):
            st.session_state["paper_chat_history"].append(
                {"role": "user", "content": _user_input}
            )
            try:
                from anthropic import Anthropic
                _ant_client = Anthropic(api_key=_p_ant_key)
                with _chat_container:
                    with st.chat_message("user"):
                        st.markdown(_user_input)
                    with st.chat_message("assistant"):
                        _resp_placeholder = st.empty()
                        _full_resp = ""
                        with _ant_client.messages.stream(
                            model="claude-sonnet-4-6",
                            max_tokens=4096,
                            system=_PAPER_SYSTEM,
                            messages=st.session_state["paper_chat_history"],
                        ) as _stream:
                            for _chunk in _stream.text_stream:
                                _full_resp += _chunk
                                _resp_placeholder.markdown(_full_resp + "▌")
                        _resp_placeholder.markdown(_full_resp)
                st.session_state["paper_chat_history"].append(
                    {"role": "assistant", "content": _full_resp}
                )
                st.session_state["paper_last_response"] = _full_resp
            except Exception as _e:
                st.error(f"Claude API 오류: {_e}")
            st.rerun()

        # ── 마지막 응답 저장 버튼 ────────────────────────────────
        if st.session_state.get("paper_last_response"):
            st.divider()
            _save_col1, _save_col2 = st.columns([3, 1])
            with _save_col1:
                _save_fname = st.text_input(
                    "저장할 파일명",
                    value="draft_v1.md",
                    key="paper_chat_save_fname",
                    label_visibility="collapsed",
                )
            with _save_col2:
                if st.button("💾 초안에 추가", use_container_width=True, key="paper_chat_save_btn"):
                    _sp = PAPER_DIR / _save_fname
                    _existing = _sp.read_text(encoding="utf-8") if _sp.exists() else ""
                    _sp.write_text(
                        _existing + "\n\n" + st.session_state["paper_last_response"],
                        encoding="utf-8",
                    )
                    st.success(f"저장됨: {_sp.name}")

    # ════════════════════════════════════════════════════════════
    # 서브탭 2: 편집 / 미리보기 (기존 기능 유지)
    # ════════════════════════════════════════════════════════════
    with _ptab_edit:
        col_left, col_right = st.columns([1, 2])

        with col_left:
            st.subheader("파일 관리")
            md_files = sorted(PAPER_DIR.glob("*.md"))
            if md_files:
                sel_file = st.selectbox(
                    "초안 선택",
                    options=[f.name for f in md_files],
                    index=0,
                    key="paper_sel_file",
                )
                current_path = PAPER_DIR / sel_file
            else:
                current_path = DRAFT_PATH
                st.info("저장된 초안이 없습니다.")

            if current_path.exists():
                raw_md = current_path.read_text(encoding="utf-8")
                st.download_button(
                    "⬇️ Markdown",
                    data=raw_md,
                    file_name=current_path.name,
                    mime="text/markdown",
                    use_container_width=True,
                )
                if st.button("⬇️ Word (.docx)", use_container_width=True, key="paper_docx_btn"):
                    try:
                        from docx import Document
                        import io, re
                        doc = Document()
                        for line in raw_md.split("\n"):
                            line = line.rstrip()
                            if line.startswith("# "):   doc.add_heading(line[2:], 1)
                            elif line.startswith("## "): doc.add_heading(line[3:], 2)
                            elif line.startswith("### "):doc.add_heading(line[4:], 3)
                            elif line == "":             doc.add_paragraph("")
                            else:
                                clean = re.sub(r'\*\*(.+?)\*\*', r'\1',
                                        re.sub(r'\*(.+?)\*',   r'\1', line))
                                doc.add_paragraph(clean)
                        buf = io.BytesIO()
                        doc.save(buf); buf.seek(0)
                        st.download_button("📄 Word 저장", data=buf,
                            file_name=current_path.stem+".docx",
                            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                            key="paper_docx_dl")
                    except Exception as _e:
                        st.error(f"Word 오류: {_e}")

        with col_right:
            st.subheader("편집 / 미리보기")
            view_mode = st.radio("모드", ["미리보기", "편집"], horizontal=True, key="paper_view_mode")
            current_content = current_path.read_text(encoding="utf-8") if current_path.exists() else "# 새 논문 초안\n"

            if view_mode == "미리보기":
                with st.container(height=750):
                    st.markdown(current_content)
            else:
                edited = st.text_area("내용 편집 (Markdown)", value=current_content,
                                      height=650, key="paper_editor",
                                      label_visibility="collapsed")
                c1, c2 = st.columns(2)
                with c1:
                    if st.button("💾 저장", use_container_width=True, key="paper_save_btn"):
                        current_path.write_text(edited, encoding="utf-8")
                        st.success("저장 완료"); st.rerun()
                with c2:
                    new_fname = st.text_input("버전 파일명", value="draft_v2.md",
                                              key="paper_new_fname", label_visibility="collapsed")
                    if st.button("📋 다른 이름 저장", use_container_width=True, key="paper_saveas_btn"):
                        (PAPER_DIR / new_fname).write_text(edited, encoding="utf-8")
                        st.success(f"저장: {new_fname}"); st.rerun()

            if current_path.exists():
                txt = current_path.read_text(encoding="utf-8")
                st.caption(f"단어: **{len(txt.split()):,}** | 섹션: **{txt.count(chr(10)+'## ')}** | `{current_path.name}`")


with tab_cpg:
    st.subheader("📋 한의 CPG 관리")
    st.caption("임상진료지침(CPG) 문서를 업로드하여 약재발굴 결과 검증에 활용합니다.")

    from pathlib import Path as _PPath
    _cpg_settings = load_settings()
    _cpg_api_key  = _cpg_settings.get("gemini_api_key", "")

    _cpg_col1, _cpg_col2 = st.columns([1, 1])

    with _cpg_col1:
        st.markdown("**CPG 문서 입력**")
        _cpg_input_mode = st.radio(
            "입력 방식",
            ["📄 파일 업로드", "📋 텍스트 붙여넣기"],
            horizontal=True,
            key="cpg_input_mode",
        )
        _cpg_file = None
        _cpg_paste_text = ""
        if _cpg_input_mode == "📄 파일 업로드":
            _cpg_file = st.file_uploader(
                "CPG 파일 선택 (PDF / DOCX / TXT)",
                type=["pdf", "docx", "doc", "txt"],
                key="cpg_uploader",
            )
        else:
            _cpg_paste_text = st.text_area(
                "PDF에서 복사한 텍스트를 붙여넣으세요",
                height=200,
                key="cpg_paste",
                placeholder="Ctrl+A → Ctrl+C 후 여기에 붙여넣기",
            )

        _cpg_disease_key = st.text_input(
            "질환 키 (영문 소문자, 예: insomnia)",
            value="insomnia",
            key="cpg_disease_key",
        )

        _cpg_has_input = bool(_cpg_file or _cpg_paste_text.strip())
        if _cpg_has_input and _cpg_api_key:
            if st.button("🤖 AI로 CPG 구조화 추출", use_container_width=True, key="cpg_extract_btn"):
                with st.status("CPG 분석 중...", expanded=True) as _cpg_st:
                    try:
                        from cpg_parser import extract_text, extract_cpg_with_gemini, map_herbs_to_pinyin, save_cpg_entry

                        _cpg_st.update(label="텍스트 준비 중...")
                        if _cpg_paste_text.strip():
                            _cpg_text = _cpg_paste_text.strip()
                            st.write(f"붙여넣기 텍스트: {len(_cpg_text):,}자")
                        else:
                            _cpg_bytes = _cpg_file.getvalue()
                            _cpg_text  = extract_text(_cpg_bytes, _cpg_file.name)
                            st.write(f"텍스트 추출 완료: {len(_cpg_text):,}자")

                        _cpg_st.update(label="Gemini 구조화 추출 중...")
                        _cpg_result = extract_cpg_with_gemini(
                            _cpg_text, _cpg_api_key,
                            on_progress=lambda m: _cpg_st.update(label=m)
                        )

                        if not _cpg_result["success"]:
                            _cpg_st.update(label="추출 실패", state="error")
                            st.error(f"오류: {_cpg_result['error']}")
                        else:
                            data = _cpg_result["data"]

                            # 한글명 → Pinyin 자동 매핑
                            _cpg_st.update(label="약재명 Pinyin 매핑 중...")
                            for formula in data.get("formulas", []):
                                ko_herbs = formula.get("herbs_ko", [])
                                mapping  = map_herbs_to_pinyin(ko_herbs)
                                formula["herbs_pinyin"] = [mapping.get(h) for h in ko_herbs]
                                formula["pinyin_map"]   = mapping

                            for herb in data.get("single_herbs", []):
                                ko = herb.get("name_ko", "")
                                herb["pinyin"] = map_herbs_to_pinyin([ko]).get(ko)

                            # 저장
                            save_cpg_entry(_cpg_disease_key, data)
                            _cpg_st.update(label="CPG 저장 완료 ✅", state="complete")
                            st.session_state["cpg_extracted"] = data
                            st.success(f"CPG 저장 완료 → data/korean_cpg.json [{_cpg_disease_key}]")

                    except Exception as _cpg_e:
                        _cpg_st.update(label="오류", state="error")
                        st.error(f"오류: {_cpg_e}")

        elif _cpg_has_input and not _cpg_api_key:
            st.warning("Gemini API 키가 필요합니다 (AI 설정 탭에서 입력).")

    with _cpg_col2:
        st.markdown("**저장된 CPG 확인**")
        try:
            from cpg_parser import load_cpg_data, GRADE_COLORS
            _saved_cpg = load_cpg_data()

            if not _saved_cpg:
                st.info("저장된 CPG 없음")
            else:
                for _dk, _dv in _saved_cpg.items():
                    with st.expander(f"**{_dk}** — {_dv.get('cpg_source','')}", expanded=(_dk == "hypertension")):
                        st.caption(f"발행: {_dv.get('publisher','')} {_dv.get('year','')}")

                        # 처방 목록
                        _formulas = _dv.get("formulas", [])
                        if _formulas:
                            st.markdown("**권고 처방:**")
                            for _f in _formulas:
                                _gc = GRADE_COLORS.get(_f.get("grade","?"), "⚫")
                                _herbs_str = ", ".join(
                                    f"{k}({v})" if v else k
                                    for k, v in (_f.get("pinyin_map") or {}).items()
                                )
                                st.markdown(
                                    f"{_gc} **{_f.get('name_ko','')}** "
                                    f"[Grade {_f.get('grade','?')}] — {_f.get('evidence_type','')}\n\n"
                                    f"&nbsp;&nbsp;{_herbs_str}"
                                )

                        # 단미 약재
                        _singles = _dv.get("single_herbs", [])
                        if _singles:
                            st.markdown("**단미 약재 권고:**")
                            _s_rows = []
                            for _s in _singles:
                                _gc = GRADE_COLORS.get(_s.get("grade","?"), "⚫")
                                _s_rows.append({
                                    "등급": f"{_gc} {_s.get('grade','?')}",
                                    "약재(한글)": _s.get("name_ko",""),
                                    "Pinyin": _s.get("pinyin",""),
                                    "근거": _s.get("evidence_type",""),
                                    "적응증": _s.get("indication",""),
                                })
                            st.dataframe(pd.DataFrame(_s_rows), use_container_width=True, hide_index=True)

                        # 수동 편집 버튼
                        if st.button(f"🗑 {_dk} 삭제", key=f"cpg_del_{_dk}"):
                            _saved_cpg.pop(_dk)
                            from cpg_parser import CPG_PATH
                            CPG_PATH.write_text(json.dumps(_saved_cpg, ensure_ascii=False, indent=2))
                            st.rerun()

        except Exception as _cpg_load_e:
            st.error(f"CPG 로드 오류: {_cpg_load_e}")

    # ── 추출 결과 미리보기 (추출 직후) ──────────────────────────────────────
    if "cpg_extracted" in st.session_state:
        st.divider()
        st.markdown("**추출 결과 미리보기** (저장 완료)")
        _prev = st.session_state["cpg_extracted"]

        # 매핑 미확인 약재 경고
        _unmapped = []
        for _f in _prev.get("formulas", []):
            for _ko, _pin in (_f.get("pinyin_map") or {}).items():
                if not _pin:
                    _unmapped.append(_ko)
        for _s in _prev.get("single_herbs", []):
            if not _s.get("pinyin"):
                _unmapped.append(_s.get("name_ko",""))

        if _unmapped:
            st.warning(
                f"Pinyin 매핑 실패 약재 {len(_unmapped)}개: {', '.join(set(_unmapped))}\n\n"
                "→ `src/cpg_parser.py`의 `KO_TO_PINYIN` 테이블에 추가하거나 직접 JSON 편집하세요."
            )

        st.json(_prev, expanded=False)
