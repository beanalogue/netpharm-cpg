"""
NetPharm 설정 영구 저장 모듈
~/.netpharm_settings.json 에 저장 (API 키, 모델 등)
"""
import json
from pathlib import Path

SETTINGS_PATH = Path.home() / ".netpharm_settings.json"

DEFAULTS = {
    "gemini_api_key": "",
    "setup_model":    "gemma-4-31b-it",      # AI 설정 파서 모델
    "ai_model":       "gemini-3.5-flash",      # AI 해석 모델
    "ppi_min_score":  400,
    "fdr_cutoff":     0.05,
    "anthropic_api_key": "",
    "top_n":          15,
}


def load_settings() -> dict:
    if SETTINGS_PATH.exists():
        try:
            saved = json.loads(SETTINGS_PATH.read_text())
            return {**DEFAULTS, **saved}
        except Exception:
            pass
    return dict(DEFAULTS)


def save_settings(settings: dict):
    existing = load_settings()          # 기존 저장값 먼저 로드
    merged = {**existing, **settings}   # 새 값으로 덮어쓰되 나머지는 보존
    SETTINGS_PATH.write_text(json.dumps(merged, ensure_ascii=False, indent=2))
