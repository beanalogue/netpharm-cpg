"""
AI 기반 연구 설정 파서
자연어 연구 설명 → 파이프라인 설정 JSON 자동 추출

지원 모델 (Google AI Studio 무료 tier):
  - gemini-2.0-flash  (기본값, 빠름·정확)
  - gemma-3-27b-it    (오픈소스, 무료)

사용:
    parser = ResearchSetupParser(api_key="AIza...")
    result = parser.parse("가미소요산과 스틸녹스의 불면증 관련 상호작용 분석")
    # result: {"config": {...}, "summary": "...", "clarification": None}
"""

import json
import logging
import re
import time
from typing import Optional

_RETRY_WAITS = [30, 60, 120]   # 503/429 재시도 대기(초)

log = logging.getLogger(__name__)

# ── 알려진 한약 처방 데이터베이스 ─────────────────────────────────────────────
KNOWN_FORMULAS = {
    "가미소요산": {
        "en": "Jiawei Xiaoyao San / Kamishoyosan",
        "herbs": {
            "Bupleurum chinense":        ["saikosaponin A", "saikosaponin D", "rutin", "quercetin"],
            "Paeonia lactiflora":        ["paeoniflorin", "albiflorin", "oxypaeoniflorin"],
            "Angelica sinensis":         ["ligustilide", "ferulic acid", "Z-butylidenephthalide"],
            "Glycyrrhiza uralensis":     ["liquiritigenin", "isoliquiritigenin", "licochalcone A"],
            "Atractylodes macrocephala": ["atractylenolide I", "atractylenolide III"],
            "Poria cocos":               ["pachymic acid", "poricoic acid A"],
            "Gardenia jasminoides":      ["genipin", "geniposide", "chlorogenic acid"],
            "Mentha haplocalyx":         ["menthol", "menthone", "rosmarinic acid"],
            "Zingiber officinale":       ["6-gingerol", "zingerone"],
            "Paeonia suffruticosa":      ["paeonol", "mudanpioside C"],
        }
    },
    "소요산": {
        "en": "Xiaoyao San",
        "herbs": {
            "Bupleurum chinense":        ["saikosaponin A", "rutin", "quercetin"],
            "Paeonia lactiflora":        ["paeoniflorin", "albiflorin"],
            "Angelica sinensis":         ["ligustilide", "ferulic acid"],
            "Glycyrrhiza uralensis":     ["liquiritigenin", "isoliquiritigenin"],
            "Atractylodes macrocephala": ["atractylenolide I", "atractylenolide III"],
            "Poria cocos":               ["pachymic acid", "poricoic acid A"],
            "Zingiber officinale":       ["6-gingerol"],
            "Mentha haplocalyx":         ["menthol", "menthone"],
        }
    },
    "보중익기탕": {
        "en": "Hochuekkito / Bu Zhong Yi Qi Tang",
        "herbs": {
            "Astragalus membranaceus":   ["astragaloside IV", "calycosin", "formononetin"],
            "Codonopsis pilosula":       ["lobetyolin", "atractylodin"],
            "Atractylodes macrocephala": ["atractylenolide I", "atractylenolide III"],
            "Glycyrrhiza uralensis":     ["liquiritigenin", "glycyrrhizin"],
            "Angelica sinensis":         ["ligustilide", "ferulic acid"],
            "Cimicifuga foetida":        ["cimicifugoside", "actein"],
            "Bupleurum chinense":        ["saikosaponin A", "rutin"],
            "Zingiber officinale":       ["6-gingerol"],
        }
    },
    "황련해독탕": {
        "en": "Orengedokuto / Huang Lian Jie Du Tang",
        "herbs": {
            "Coptis chinensis":          ["berberine", "coptisine", "palmatine"],
            "Scutellaria baicalensis":   ["baicalein", "baicalin", "wogonin"],
            "Phellodendron amurense":    ["berberine", "phellodendrine"],
            "Gardenia jasminoides":      ["genipin", "geniposide", "chlorogenic acid"],
        }
    },
    "육미지황탕": {
        "en": "Yukmijihwang-tang / Liu Wei Di Huang Tang",
        "herbs": {
            "Rehmannia glutinosa":       ["catalpol", "acteoside", "aucubin"],
            "Cornus officinalis":        ["loganin", "morroniside", "cornuside"],
            "Dioscorea opposita":        ["diosgenin", "allantoin"],
            "Alisma orientale":          ["alisol A", "alisol B", "alismoxide"],
            "Poria cocos":               ["pachymic acid", "poricoic acid A"],
            "Paeonia suffruticosa":      ["paeonol", "paeoniflorin"],
        }
    },
    "반하사심탕": {
        "en": "Banhasasim-tang / Ban Xia Xie Xin Tang",
        "herbs": {
            "Pinellia ternata":          ["homogentisic acid", "ephedrine"],
            "Scutellaria baicalensis":   ["baicalein", "baicalin", "wogonin"],
            "Coptis chinensis":          ["berberine", "coptisine"],
            "Zingiber officinale":       ["6-gingerol", "zingerone"],
            "Panax ginseng":             ["ginsenoside Rb1", "ginsenoside Rg1", "ginsenoside Re"],
            "Glycyrrhiza uralensis":     ["liquiritigenin", "glycyrrhizin"],
            "Ziziphus jujuba":           ["jujuboside A", "betulinic acid"],
        }
    },
    "사물탕": {
        "en": "Samul-tang / Si Wu Tang",
        "herbs": {
            "Rehmannia glutinosa":       ["catalpol", "acteoside"],
            "Paeonia lactiflora":        ["paeoniflorin", "albiflorin"],
            "Angelica sinensis":         ["ligustilide", "ferulic acid"],
            "Ligusticum chuanxiong":     ["senkyunolide A", "tetramethylpyrazine", "ligustilide"],
        }
    },
}

# ── 약물 이름 매핑 (상품명/한글명 → 성분명) ────────────────────────────────────
DRUG_ALIASES = {
    "스틸녹스": "zolpidem", "졸피뎀": "zolpidem", "zolpidem": "zolpidem",
    "아스피린": "aspirin", "aspirin": "aspirin",
    "메트포르민": "metformin", "metformin": "metformin", "글루코파지": "metformin",
    "세로퀄": "quetiapine", "쿠에티아핀": "quetiapine", "quetiapine": "quetiapine",
    "리스페달": "risperidone", "리스페리돈": "risperidone",
    "렉사프로": "escitalopram", "에스시탈로프람": "escitalopram",
    "자낙스": "alprazolam", "알프라졸람": "alprazolam",
    "리보트릴": "clonazepam", "클로나제팜": "clonazepam",
    "데파킨": "valproic acid", "발프로산": "valproic acid",
    "리피토": "atorvastatin", "아토르바스타틴": "atorvastatin",
    "크레스토": "rosuvastatin", "로수바스타틴": "rosuvastatin",
    "타미플루": "oseltamivir", "오셀타미비르": "oseltamivir",
    "글리벡": "imatinib", "이마티닙": "imatinib",
    "허셉틴": "trastuzumab", "트라스투주맙": "trastuzumab",
    "독소루비신": "doxorubicin", "아드리아마이신": "doxorubicin",
    "시스플라틴": "cisplatin", "cisplatin": "cisplatin",
    "메토트렉세이트": "methotrexate", "methotrexate": "methotrexate",
    "프레드니손": "prednisone", "prednisone": "prednisone",
    "덱사메타손": "dexamethasone", "dexamethasone": "dexamethasone",
    "이부프로펜": "ibuprofen", "부루펜": "ibuprofen",
    "나프록센": "naproxen", "탁센": "naproxen",
    "셀레콕시브": "celecoxib", "쎄레브렉스": "celecoxib",
    "암로디핀": "amlodipine", "노바스크": "amlodipine",
    "리시노프릴": "lisinopril", "로트릴": "lisinopril",
    "옴프라졸": "omeprazole", "오메프라졸": "omeprazole",
    "판토프라졸": "pantoprazole", "판토록": "pantoprazole",
}

# ── 질환 영어 매핑 ────────────────────────────────────────────────────────────
DISEASE_ALIASES = {
    "불면증": "insomnia", "수면장애": "sleep disorder",
    "당뇨": "diabetes", "당뇨병": "diabetes mellitus",
    "고혈압": "hypertension", "혈압": "hypertension",
    "암": "cancer", "종양": "tumor",
    "대장암": "colorectal cancer", "간암": "hepatocellular carcinoma",
    "유방암": "breast cancer", "폐암": "lung cancer",
    "위암": "gastric cancer", "췌장암": "pancreatic cancer",
    "우울증": "depression", "우울": "depression",
    "불안": "anxiety", "불안장애": "anxiety disorder",
    "치매": "Alzheimer's disease", "알츠하이머": "Alzheimer's disease",
    "파킨슨": "Parkinson's disease",
    "간염": "hepatitis", "간경화": "liver cirrhosis", "지방간": "fatty liver",
    "관절염": "arthritis", "류마티스": "rheumatoid arthritis",
    "천식": "asthma", "폐질환": "pulmonary disease", "COPD": "COPD",
    "신부전": "renal failure", "신장병": "kidney disease",
    "심부전": "heart failure", "심근경색": "myocardial infarction",
    "뇌졸중": "stroke", "뇌경색": "cerebral infarction",
    "비만": "obesity", "대사증후군": "metabolic syndrome",
    "골다공증": "osteoporosis",
    "아토피": "atopic dermatitis", "건선": "psoriasis",
    "크론": "Crohn's disease", "궤양성대장염": "ulcerative colitis",
    "루푸스": "lupus", "전신홍반루푸스": "systemic lupus erythematosus",
}

SYSTEM_PROMPT = """You are an expert in traditional Korean/Chinese medicine (TKM/TCM) and network pharmacology.
Your task: extract a structured research configuration from the user's natural language description.

Known herb formulas (use these exact herb/compound lists when the formula name matches):
""" + json.dumps({k: v["herbs"] for k, v in KNOWN_FORMULAS.items()}, ensure_ascii=False, indent=2) + """

Known drug aliases (map to generic name):
""" + json.dumps(DRUG_ALIASES, ensure_ascii=False) + """

Rules:
1. If the user mentions a known formula (가미소요산, 보중익기탕, etc.), use its full herb list.
2. If individual herbs are mentioned without a known formula, use top 3-5 known active compounds per herb.
3. Map drug brand names / Korean names to generic English names.
4. Extract disease keyword in English.
5. Generate a study name in snake_case (e.g., kamishoyosan_zolpidem).
6. If any information is ambiguous or missing, set clarification to a specific Korean question.
7. Respond ONLY with a JSON object in this exact format:

{
  "name": "study_name_snake_case",
  "herbs": {
    "Latin herb name": ["compound1", "compound2", "compound3"],
    ...
  },
  "drugs": ["generic_drug_name"],
  "disease_query": "english disease keyword",
  "summary_ko": "한국어로 이해한 연구 내용 2-3줄 요약",
  "clarification": null
}

If you cannot extract enough information, set clarification to a Korean question asking for the missing info,
and set other fields to reasonable defaults or empty values."""


class ResearchSetupParser:
    """Gemini/Gemma 기반 연구 설정 자동 추출"""

    def __init__(self, api_key: str, model: str = "gemini-2.0-flash"):
        self.api_key = api_key
        self.model   = model
        self._client = None

    def _get_client(self):
        if self._client is None:
            try:
                from google import genai
                self._client = genai.Client(api_key=self.api_key)
            except Exception as e:
                raise RuntimeError(f"google-genai SDK 오류: {e}")
        return self._client

    def parse(self, user_message: str, chat_history: list = None) -> dict:
        """
        자연어 → 파이프라인 설정 추출

        반환:
        {
            "config":        dict | None,   # 파이프라인 config
            "summary_ko":    str,           # 한국어 요약
            "clarification": str | None,    # 추가 질문 (None이면 추출 성공)
            "raw_response":  str,           # LLM 원본 응답
        }
        """
        # 로컬 처방 DB에서 먼저 매칭 시도
        quick = self._quick_match(user_message)
        if quick:
            return quick

        # LLM 호출 (503/429 시 최대 3회 자동 재시도)
        client  = self._get_client()
        history = self._build_history(chat_history or [])
        prompt  = SYSTEM_PROMPT + "\n\nUser: " + user_message
        last_exc = None

        for attempt in range(1, len(_RETRY_WAITS) + 2):
            try:
                response = client.models.generate_content(
                    model    = self.model,
                    contents = history + [{"role": "user", "parts": [{"text": prompt}]}],
                )
                return self._parse_response(response.text.strip())
            except Exception as e:
                last_exc = e
                msg = str(e)
                is_temp = ("503" in msg or "429" in msg or
                           "UNAVAILABLE" in msg or "ResourceExhausted" in msg)
                if attempt > len(_RETRY_WAITS) or not is_temp:
                    break
                wait = _RETRY_WAITS[attempt - 1]
                log.warning(f"[AI Setup] 일시적 오류 ({attempt}/3) — {wait}초 후 재시도: {e}")
                time.sleep(wait)

        log.error(f"[AI Setup] LLM 호출 실패: {last_exc}")
        return {
            "config": None,
            "summary_ko": f"LLM 호출 오류: {last_exc}",
            "clarification": "API 키를 확인해주세요. 아니면 수동 설정 탭을 이용해주세요.",
            "raw_response": str(last_exc),
        }

    def _quick_match(self, text: str) -> dict | None:
        """로컬 처방 DB 빠른 매칭 (API 호출 없음)"""
        found_formula = None
        for name, data in KNOWN_FORMULAS.items():
            if name in text:
                found_formula = (name, data)
                break
        if not found_formula:
            return None

        formula_name, formula_data = found_formula

        # 약물 찾기
        drugs = []
        for alias, generic in DRUG_ALIASES.items():
            if alias in text and generic not in drugs:
                drugs.append(generic)

        # 질환 찾기
        disease = "not specified"
        for kor, eng in DISEASE_ALIASES.items():
            if kor in text:
                disease = eng
                break

        # 영문에서도 질환 검색
        for eng in ["insomnia", "depression", "anxiety", "cancer", "diabetes",
                    "hypertension", "arthritis", "alzheimer", "parkinson"]:
            if eng.lower() in text.lower() and disease == "not specified":
                disease = eng

        study_name = formula_name.replace(" ", "_")
        if drugs:
            study_name += "_" + drugs[0].replace(" ", "_")

        config = {
            "name":            study_name,
            "herbs":           formula_data["herbs"],
            "drugs":           drugs,
            "disease_query":   disease,
            "disease_min_score": 0.1,
            "ppi_min_score":   400,
            "adj_p_cutoff":    0.05,
            "top_n_enrichment": 20,
            "make_plots":      True,
            "enrichment_dbs":  None,
            "admet": {
                "ob_min": 30.0, "dl_min": 0.18, "mw_max": 500.0,
                "logp_max": 5.0, "hbd_max": 5, "hba_max": 10,
                "tpsa_max": 140.0, "rotbonds_max": 10,
            },
        }

        drug_str = f"및 {', '.join(drugs)}" if drugs else ""
        dis_str  = f"({disease})" if disease != "not specified" else ""

        return {
            "config":        config,
            "summary_ko":    f"{formula_name} {drug_str} 상호작용 네트워크 약리학 분석 {dis_str}\n"
                             f"구성 생약 {len(formula_data['herbs'])}종, "
                             f"약물 {len(drugs)}종으로 분석을 구성했습니다.",
            "clarification": None,
            "raw_response":  "[로컬 DB 빠른 매칭]",
        }

    def _parse_response(self, raw: str) -> dict:
        """LLM 응답에서 JSON 추출"""
        # 마크다운 코드블록 제거
        cleaned = re.sub(r"```(?:json)?\s*", "", raw).strip().rstrip("```").strip()
        # 첫 { ... } 블록 추출
        m = re.search(r"\{[\s\S]+\}", cleaned)
        if not m:
            return {
                "config": None,
                "summary_ko": "응답 파싱 실패",
                "clarification": "다시 시도하거나 더 구체적으로 설명해 주세요.",
                "raw_response": raw,
            }
        try:
            data = json.loads(m.group())
        except json.JSONDecodeError:
            return {
                "config": None,
                "summary_ko": "JSON 파싱 오류",
                "clarification": "다시 시도해 주세요.",
                "raw_response": raw,
            }

        # 약물 alias 보정
        drugs = [DRUG_ALIASES.get(d.lower(), d) for d in data.get("drugs", [])]

        config = {
            "name":            data.get("name", "analysis"),
            "herbs":           data.get("herbs", {}),
            "drugs":           drugs,
            "disease_query":   data.get("disease_query", ""),
            "disease_min_score": 0.1,
            "ppi_min_score":   400,
            "adj_p_cutoff":    0.05,
            "top_n_enrichment": 20,
            "make_plots":      True,
            "enrichment_dbs":  None,
            "admet": {
                "ob_min": 30.0, "dl_min": 0.18, "mw_max": 500.0,
                "logp_max": 5.0, "hbd_max": 5, "hba_max": 10,
                "tpsa_max": 140.0, "rotbonds_max": 10,
            },
        }

        return {
            "config":        config if config["herbs"] else None,
            "summary_ko":    data.get("summary_ko", ""),
            "clarification": data.get("clarification"),
            "raw_response":  raw,
        }

    @staticmethod
    def _build_history(chat_history: list) -> list:
        """Gemini API 형식 히스토리 변환"""
        result = []
        for msg in chat_history[-6:]:  # 최근 6개만
            role = "user" if msg["role"] == "user" else "model"
            result.append({"role": role, "parts": [{"text": msg["content"]}]})
        return result

    def chat_reply(self, user_message: str, chat_history: list = None) -> str:
        """
        일반 대화 응답 (설정 추출이 아닌 추가 질문/수정 대화용)
        """
        try:
            client = self._get_client()
            system = (
                "당신은 네트워크 약리학 연구 보조 AI입니다. "
                "한국어로 짧고 명확하게 대답하세요. "
                "연구 설정 수정 요청이 있으면 정확히 무엇을 바꾸는지 확인하고, "
                "변경된 설정을 JSON으로 제시하세요."
            )
            history = self._build_history(chat_history or [])
            response = client.models.generate_content(
                model    = self.model,
                contents = history + [{"role": "user", "parts": [{"text": system + "\n\n" + user_message}]}],
            )
            return response.text.strip()
        except Exception as e:
            return f"오류: {e}"
