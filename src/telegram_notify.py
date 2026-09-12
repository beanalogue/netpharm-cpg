"""
telegram_notify.py — 텔레그램 알림 모듈
분석 진행 상황 및 완료/에러 알림 전송

설정:
    환경변수 또는 .env 파일에 아래 항목 설정
    TELEGRAM_TOKEN  = "봇 토큰 (BotFather에서 발급)"
    TELEGRAM_CHAT_ID = "채팅 ID (봇에게 메시지 후 /getUpdates로 확인)"
"""

import logging
import os
import time
from pathlib import Path

import requests

log = logging.getLogger(__name__)

_BASE_URL = "https://api.telegram.org/bot{token}/sendMessage"
_LAST_SENT: float = 0
_MIN_INTERVAL = 2.0  # 최소 전송 간격 (초) — rate limit 방지


def _load_config() -> tuple[str, str]:
    """토큰·chat_id 로드: 환경변수 → .env 파일 순서"""
    token   = os.environ.get("TELEGRAM_TOKEN", "")
    chat_id = os.environ.get("TELEGRAM_CHAT_ID", "")

    if not token or not chat_id:
        env_path = Path(__file__).parent.parent / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                line = line.strip()
                if line.startswith("TELEGRAM_TOKEN="):
                    token = line.split("=", 1)[1].strip().strip('"')
                elif line.startswith("TELEGRAM_CHAT_ID="):
                    chat_id = line.split("=", 1)[1].strip().strip('"')

    return token, chat_id


def send(message: str, silent: bool = False) -> bool:
    """
    텔레그램 메시지 전송.
    설정이 없으면 조용히 무시(False 반환).
    silent=True 면 알림음 없이 전송.
    """
    global _LAST_SENT

    token, chat_id = _load_config()
    if not token or not chat_id:
        return False

    # rate limit 방지
    gap = time.time() - _LAST_SENT
    if gap < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - gap)

    try:
        url = _BASE_URL.format(token=token)
        resp = requests.post(url, json={
            "chat_id":              chat_id,
            "text":                 message,
            "parse_mode":           "HTML",
            "disable_notification": silent,
        }, timeout=10)
        resp.raise_for_status()
        _LAST_SENT = time.time()
        return True
    except Exception as e:
        log.debug(f"[Telegram] 전송 실패: {e}")
        return False


def notify_start(name: str, total: int, remaining: int):
    send(
        f"🔬 <b>분석 시작</b>\n"
        f"• 작업: {name}\n"
        f"• 전체: {total}개 / 남은 작업: {remaining}개"
    )


def notify_progress(name: str, processed: int, total: int,
                    success: int, fail: int, eta_min: int | None, current: str = ""):
    eta_str = f"{eta_min // 60}h {eta_min % 60}m" if eta_min else "계산 중"
    cur_str = f"\n• 현재: {current}" if current else ""
    send(
        f"📊 <b>{name} 진행 중</b>\n"
        f"• 진행: {processed}/{total} ({processed/total*100:.0f}%)\n"
        f"• 성공 {success} / 실패 {fail}\n"
        f"• 예상 남은 시간: {eta_str}"
        f"{cur_str}",
        silent=True,
    )


def notify_done(name: str, success: int, fail: int, elapsed_h: float):
    fail_str = f"\n⚠️ 실패 {fail}개 — logs/build_tcmsp_cache_errors.json 확인" if fail else ""
    send(
        f"✅ <b>{name} 완료!</b>\n"
        f"• 성공: {success}개\n"
        f"• 소요 시간: {elapsed_h:.1f}시간"
        f"{fail_str}"
    )


def notify_error(name: str, error: str):
    send(f"❌ <b>{name} 오류 발생</b>\n<code>{error[:300]}</code>")


def notify_analysis_step(project: str, step: str, detail: str = ""):
    detail_str = f"\n{detail}" if detail else ""
    send(
        f"▶ <b>[{project}]</b> {step}{detail_str}",
        silent=True,
    )


def notify_analysis_done(project: str, summary: dict):
    """분석 완료 알림 — summary dict에서 핵심 지표 추출"""
    herb_genes  = summary.get("herb_genes", "?")
    drug_genes  = summary.get("drug_genes", "?")
    intersect   = summary.get("intersection", "?")
    top_genes   = summary.get("hub_genes", [])
    top_str     = ", ".join(top_genes[:5]) if top_genes else "-"

    send(
        f"🎉 <b>[{project}] 분석 완료!</b>\n"
        f"• 한약 타깃: {herb_genes}개\n"
        f"• 약물 타깃: {drug_genes}개\n"
        f"• 교차 타깃: {intersect}개\n"
        f"• 허브 유전자: {top_str}"
    )
