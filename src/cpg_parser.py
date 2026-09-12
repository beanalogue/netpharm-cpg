"""
cpg_parser.py
한의 임상진료지침(CPG) 문서 파싱 및 구조화 모듈

PDF/DOCX/TXT → Gemini AI → 구조화된 JSON (약재 + 권고등급)
"""

import json
import logging
import re
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

CPG_PATH = Path("data/korean_cpg.json")

# 권고등급 한글/영문 매핑
GRADE_MAP = {
    "강력권고": "A", "권고": "B", "약한권고": "C", "전문가합의": "D",
    "적극권고": "A", "조건부권고": "B",
    "grade a": "A", "grade b": "B", "grade c": "C", "grade d": "D",
    "strong": "A", "moderate": "B", "weak": "C", "consensus": "D",
    "1a": "A", "1b": "A", "2a": "B", "2b": "B", "3": "C", "4": "D",
    "a": "A", "b": "B", "c": "C", "d": "D",
}

GRADE_COLORS = {"A": "🟢", "B": "🔵", "C": "🟡", "D": "⚪", "?": "⚫"}


# ─────────────────────────────────────────────────────────────────────────────
# 텍스트 추출
# ─────────────────────────────────────────────────────────────────────────────
def extract_text_from_pdf(file_bytes: bytes) -> str:
    import io

    # 1차: pymupdf
    try:
        import pymupdf
        doc = pymupdf.open(stream=file_bytes, filetype="pdf")
        pages = []
        for page in doc:
            t = page.get_text()
            if not t.strip():
                blocks = page.get_text("blocks")
                if isinstance(blocks, list):
                    t = "\n".join(b[4] for b in blocks if len(b) > 4 and isinstance(b[4], str))
            pages.append(t)
        result = "\n".join(pages).strip()
        if result:
            log.info(f"pymupdf 추출: {len(result)}자")
            return result
    except Exception as e:
        log.warning(f"pymupdf 실패: {e}")

    # 2차: pdfminer.six (한글 CID 폰트 대응)
    try:
        from pdfminer.high_level import extract_text as pdfminer_extract
        result = pdfminer_extract(io.BytesIO(file_bytes)).strip()
        if result:
            log.info(f"pdfminer 추출: {len(result)}자")
            return result
    except Exception as e:
        log.warning(f"pdfminer 실패: {e}")

    log.error("PDF 텍스트 추출 전부 실패")
    return ""


def extract_text_from_docx(file_bytes: bytes) -> str:
    try:
        import io
        from docx import Document
        doc = Document(io.BytesIO(file_bytes))
        return "\n".join(p.text for p in doc.paragraphs if p.text.strip())
    except Exception as e:
        log.error(f"DOCX 파싱 오류: {e}")
        return ""


def extract_text(file_bytes: bytes, filename: str) -> str:
    ext = Path(filename).suffix.lower()
    if ext == ".pdf":
        return extract_text_from_pdf(file_bytes)
    elif ext in (".docx", ".doc"):
        return extract_text_from_docx(file_bytes)
    else:  # txt, md 등
        return file_bytes.decode("utf-8", errors="replace")


# ─────────────────────────────────────────────────────────────────────────────
# Gemini로 구조화 추출
# ─────────────────────────────────────────────────────────────────────────────
EXTRACTION_PROMPT = """당신은 한의학 임상진료지침(CPG) 분석 전문가입니다.
아래는 한의 CPG 문서의 텍스트입니다.

이 문서에서 다음 정보를 추출하여 **반드시 JSON 형식으로만** 응답하세요.
다른 설명 없이 JSON 블록만 출력하세요.

추출 목표:
1. 질환명 (한글, 영문)
2. CPG 출처 (발행기관, 연도)
3. 권고 처방 목록 — 처방명, 권고등급, **구성 약재 한글명 목록** (반드시 추출), 근거 유형, 적응증
4. 단미 약재 권고 목록 — 약재명, 권고등급, 근거 유형, 적응증

중요: herbs_ko 필드에 구성 약재를 반드시 채우세요. 문서 내 처방 구성, 약재 목록, 조성 등의 표나 본문에서 찾아 채우세요.

권고등급 기준:
- A: 강력권고 / Strong recommendation (고품질 RCT 근거)
- B: 권고 / Moderate recommendation (중등도 근거)
- C: 약한권고 / Weak recommendation (낮은 근거)
- D: 전문가합의 / Expert consensus (근거 불충분)

출력 JSON 형식 (중괄호 포함, 이 형식 그대로):
{{
  "disease_ko": "고혈압",
  "disease_en": "Hypertension",
  "cpg_source": "대한한의학회 고혈압 한의 임상진료지침 2023",
  "publisher": "대한한의학회",
  "year": 2023,
  "formulas": [
    {{
      "name_ko": "천마구등음",
      "name_en": "Tianma Gouteng Yin",
      "grade": "A",
      "evidence_type": "RCT",
      "herbs_ko": ["천마", "구등", "황금", "두충", "우슬", "익모초"],
      "indication": "간양상항형 고혈압"
    }}
  ],
  "single_herbs": [
    {{
      "name_ko": "두충",
      "name_en": "Eucommia",
      "grade": "B",
      "evidence_type": "SR",
      "indication": "혈압강하"
    }}
  ]
}}

처방이나 약재가 없으면 빈 배열 []로 표시하세요.
권고등급이 명확하지 않으면 "?" 로 표시하세요.

--- CPG 문서 텍스트 ---
{text}
--- 끝 ---
"""


def _extract_json_from_text(raw: str) -> dict:
    """Gemini 응답에서 JSON 객체를 최대한 robust하게 추출."""
    # 1) ```json ... ``` 블록
    m = re.search(r"```(?:json)?\s*([\s\S]+?)```", raw)
    if m:
        return json.loads(m.group(1).strip())

    # 2) 첫 번째 { 부터 마지막 } 까지
    start = raw.find("{")
    end   = raw.rfind("}")
    if start != -1 and end != -1 and end > start:
        return json.loads(raw[start:end + 1])

    raise ValueError(f"JSON을 찾을 수 없습니다. 응답 앞부분: {raw[:200]}")


def _sample_text(text: str, max_chars: int = 300000) -> str:
    """Gemini 장문 컨텍스트 지원 — 가능한 전체 전달."""
    return text[:max_chars]


def extract_cpg_with_gemini(
    text: str,
    api_key: str,
    on_progress=None,
) -> dict:
    """Gemini로 CPG 텍스트 → 구조화 JSON 추출 (텍스트 모드 전용)."""
    raw = ""
    try:
        if not text or len(text.strip()) < 50:
            return {"success": False, "error": "텍스트 추출 실패: PDF에서 텍스트를 읽을 수 없습니다."}

        from google import genai as genai_new
        from google.genai import types as genai_types

        client = genai_new.Client(api_key=api_key)

        if on_progress:
            on_progress(f"Gemini로 CPG 구조화 추출 중... ({len(text):,}자)")

        sampled = _sample_text(text, max_chars=20000)
        prompt  = EXTRACTION_PROMPT.format(text=sampled)

        _models = ["gemini-2.5-flash", "gemini-3.6-flash", "gemini-2.5-pro"]
        response = None
        last_err = None
        for _m in _models:
            try:
                if on_progress:
                    on_progress(f"Gemini 구조화 추출 중... ({_m})")
                response = client.models.generate_content(
                    model=_m,
                    contents=prompt,
                    config=genai_types.GenerateContentConfig(
                        response_mime_type="application/json",
                    ),
                )
                break
            except Exception as _me:
                last_err = _me
                log.warning(f"모델 {_m} 실패: {_me}")
                continue

        if response is None:
            raise last_err

        raw    = response.text.strip()
        result = json.loads(raw)

        if not isinstance(result, dict):
            return {"success": False, "error": f"예상치 못한 응답 형식: {type(result).__name__}\n원문: {raw[:400]}"}

        return {"success": True, "data": result}

    except json.JSONDecodeError as e:
        return {"success": False, "error": f"JSON 파싱 실패: {e}\n\n응답원문:\n{raw[:800]}"}
    except Exception as e:
        return {"success": False, "error": f"{type(e).__name__}: {e}"}


# ─────────────────────────────────────────────────────────────────────────────
# TCMSP 병음 자동 매핑
# ─────────────────────────────────────────────────────────────────────────────
# 한글 약재명 → TCMSP Pinyin 매핑 테이블 (주요 약재)
KO_TO_PINYIN = {
    "구등": "Gouteng", "조구등": "Gouteng",
    "황금": "Huangqin",
    "두충": "Duzhong",
    "우슬": "Niuxi", "천우슬": "Niuxi",
    "익모초": "Yimucao",
    "상기생": "Sangjisheng",
    "복령": "Fuling",
    "치자": "Zhizi",
    "시호": "Chaihu",
    "택사": "Zexie",
    "차전자": "Cheqianzi",
    "당귀": "Danggui",
    "감초": "Gancao",
    "용담": "Longdan", "용담초": "Longdan",
    "백작약": "Baishao",
    "현삼": "Xuanshen",
    "천련자": "Chuanlianzi",
    "맥아": "Maiya",
    "인진": "Yinchen", "인진호": "Yinchen",
    "구기자": "Gouqizi",
    "국화": "Juhua",
    "숙지황": "Shudihuang",
    "산수유": "Shanzhuyu",
    "산약": "Shanyao",
    "목단피": "Mudanpi",
    "반하": "Banxia",
    "백출": "Baizhu",
    "진피": "Chenpi",
    "생강": "Shengjiang",
    "대조": "Dazao",
    "황기": "Huangqi",
    "인삼": "Renshen",
    "단삼": "Danshen",
    "천마": "Tianma",
    "지모": "Zhimu",
    "갈근": "Gegen",
    "오미자": "Beiwuweizi",
    "황련": "Huanglian",
    "천화분": "Tianhuafen",
    "지골피": "Digupi",
    "희첨초": "Xixiancao",
    "독활": "Duhuo",
    "진교": "Qinjiao",
    "방풍": "Fangfeng",
    "계지": "Guizhi",
    "강활": "Qianghuo",
    "천궁": "Chuanxiong",
    "속단": "Xuduan",
    "육계": "Rougui",
    "마황": "Mahuang",
    "부자": "Fuzi",
    "향부자": "Xiangfuzi",
    "오약": "Wuyao",
    "현호색": "Yanhusuo",
    "삼릉": "Sanleng",
    "아출": "Ezhu",
    "홍화": "Honghua",
    "도인": "Taoren",
    "적작약": "Chishao",
    "천산갑": "Chuanshanjia",
    "수질": "Shuizhi",
    "망충": "Mengchong",
    "지룡": "Dilong",
    "산사": "Shanzha",
    "신곡": "Shenqu",
    "래복자": "Laifuzi",
    "지실": "Zhishi",
    "후박": "Houpo",
    "대황": "Dahuang",
    "망초": "Mangxiao",
    "빈랑": "Binlang",
    "목향": "Muxiang",
    "침향": "Chenxiang",
    "오수유": "Wuzhuyu",
    "소회향": "Xiaohuixiang",
    "천오": "Chuanwu",
    "초오": "Caowu",
    "세신": "Xixin",
    "고량강": "Gaoliangjiang",
    "필발": "Bibo",
    "화초": "Huajiao",
    "정향": "Dingxiang",
    "용안육": "Longyanrou",
    "원지": "Yuanzhi",
    "산조인": "Suanzaoren",
    "백자인": "Baiziren",
    "합환피": "Hehuanpi",
    "하수오": "Heshouwu",
    "구판": "Guiban",
    # ── 불면증 CPG 추가 약재 ──
    "건강": "Ganjiang",
    "소엽": "Zisuye",
    "울금": "Yujin",
    "청피": "Qingpi",
    "택란": "Zelan",
    "별갑": "Biejia",
    "전갈": "Quanxie",
    "자충": "Zhechong",
    "승마": "Shengma",
    "지유": "Diyu",
    "형방": "Jingfang",       # 형개+방풍 합방 약재
    "변향부자": "Fuzi",        # 포부자 계열 변형
    # ── 기타 자주 쓰이는 약재 ──
    "포부자": "Fuzi",
    "백복령": "Fuling",
    "맥문동": "Maimendong",
    "황정": "Huangjing",
    "석창포": "Shichangpu",
    "원육": "Longyanrou",
    "우황": "Niuhuang",
    "침향": "Chenxiang",
    "자소엽": "Zisuye",
    "형개": "Jingjie",
    "패모": "Beimu",          # 천패모
    "천패모": "Beimu",
    "절패모": "Beimu",
    "모려": "Muli",
    "용골": "Longgu",
    "주사": "Zhusha",
    "야교등": "Yejiaoteng",
    "합환화": "Hehuanhua",
    "죽여": "Zhuru",
    "박하": "Bohe",
    "지황": "Dihuang", "생지황": "Dihuang", "건지황": "Dihuang",
    "지각": "Zhike",
}


def map_herbs_to_pinyin(herbs_ko: list[str]) -> dict:
    """한글 약재명 → TCMSP Pinyin 매핑. 매핑 실패 시 None 반환"""
    result = {}
    for herb in herbs_ko:
        herb_clean = herb.strip()
        pinyin = KO_TO_PINYIN.get(herb_clean)
        result[herb_clean] = pinyin
    return result


# ─────────────────────────────────────────────────────────────────────────────
# CPG 데이터 저장/로드
# ─────────────────────────────────────────────────────────────────────────────
def load_cpg_data() -> dict:
    if CPG_PATH.exists():
        try:
            return json.loads(CPG_PATH.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def save_cpg_entry(disease_key: str, cpg_data: dict) -> None:
    """CPG 데이터를 질환 키로 저장"""
    existing = load_cpg_data()
    existing[disease_key] = cpg_data
    CPG_PATH.parent.mkdir(parents=True, exist_ok=True)
    CPG_PATH.write_text(json.dumps(existing, ensure_ascii=False, indent=2))


def build_cpg_herb_pool(disease_key: str, min_grade: str = "B") -> dict:
    """
    CPG에서 검증용 약재 풀 추출.
    min_grade: 이 등급 이상만 포함 (A > B > C > D)
    반환: {"herbs": set[pinyin], "grade_map": {pinyin: grade}, "source": str}
    """
    grade_order = {"A": 0, "B": 1, "C": 2, "D": 3, "?": 4}
    min_order = grade_order.get(min_grade, 4)

    cpg = load_cpg_data().get(disease_key, {})
    if not cpg:
        return {"herbs": set(), "grade_map": {}, "source": ""}

    herbs = set()
    grade_map = {}

    # 처방에서 약재 추출
    for formula in cpg.get("formulas", []):
        fgrade = formula.get("grade", "?")
        if grade_order.get(fgrade, 4) <= min_order:
            for pinyin in formula.get("herbs_pinyin", []):
                if pinyin:
                    herbs.add(pinyin)
                    # 더 높은 등급으로 업데이트
                    if pinyin not in grade_map or \
                       grade_order.get(fgrade, 4) < grade_order.get(grade_map[pinyin], 4):
                        grade_map[pinyin] = fgrade

    # 단미 약재
    for herb in cpg.get("single_herbs", []):
        hgrade = herb.get("grade", "?")
        pinyin = herb.get("pinyin")
        if pinyin and grade_order.get(hgrade, 4) <= min_order:
            herbs.add(pinyin)
            if pinyin not in grade_map or \
               grade_order.get(hgrade, 4) < grade_order.get(grade_map.get(pinyin, "?"), 4):
                grade_map[pinyin] = hgrade

    return {
        "herbs":     herbs,
        "grade_map": grade_map,
        "source":    cpg.get("cpg_source", "Korean CPG"),
        "total":     len(herbs),
    }
