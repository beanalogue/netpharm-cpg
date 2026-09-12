"""
TCMSP 전체 500종 약재 캐시 빌더
- 이미 캐시된 약재는 자동 스킵 (재실행하면 이어서 진행)
- 진행 상태를 logs/build_tcmsp_cache_status.json 에 실시간 기록 (프론트엔드 표시용)
- 실패 내역을 logs/build_tcmsp_cache_errors.json 에 별도 저장

실행 방법:
    cd /home/beanalogue/netpharm
    nohup python3 build_tcmsp_cache.py > logs/build_tcmsp_cache.log 2>&1 &
    echo $!                              # PID 기록
    tail -f logs/build_tcmsp_cache.log   # 진행 모니터링
"""

import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

BASE = Path(__file__).parent
sys.path.insert(0, str(BASE / "src"))

from tcmsp_scraper import scrape_herb, _CACHE_DIR
from telegram_notify import notify_start, notify_progress, notify_done, notify_error

# ── 로깅 ─────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.StreamHandler()],
)
log = logging.getLogger(__name__)

# ── 경로 ─────────────────────────────────────────────────────────────────────
LOGS_DIR   = BASE / "logs"
LOGS_DIR.mkdir(exist_ok=True)

STATUS_PATH = LOGS_DIR / "build_tcmsp_cache_status.json"
ERRORS_PATH = LOGS_DIR / "build_tcmsp_cache_errors.json"
HERB_LIST   = BASE / "data" / "tcmsp_herb_list.json"


# ── 유틸 ─────────────────────────────────────────────────────────────────────
def already_cached(pinyin: str) -> bool:
    safe = pinyin.strip().replace(" ", "_")
    p = _CACHE_DIR / f"{safe}.json"
    if not p.exists():
        return False
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        # OB 필터 적용된 구버전 캐시는 유효하지 않음
        return data.get("ob_cutoff", 0) == 0
    except Exception:
        return False


def count_cached(herbs: list) -> int:
    return sum(1 for _, v in herbs if already_cached(v.get("herb_pinyin", "")))


def write_status(status: dict):
    status["last_updated"] = datetime.now().isoformat()
    STATUS_PATH.write_text(json.dumps(status, ensure_ascii=False, indent=2))


def write_errors(errors: list):
    ERRORS_PATH.write_text(json.dumps(errors, ensure_ascii=False, indent=2))


# ── 메인 ─────────────────────────────────────────────────────────────────────
def main():
    herb_db = json.loads(HERB_LIST.read_text(encoding="utf-8"))
    herbs   = list(herb_db.items())   # [(herb_en_name, {pinyin, cn_name}), ...]
    total   = len(herbs)

    cached_count = count_cached(herbs)
    remaining    = total - cached_count

    log.info("=== TCMSP 전체 캐시 빌드 시작 ===")
    log.info(f"전체: {total}개 / 이미 캐시됨: {cached_count}개 / 남은 작업: {remaining}개")
    notify_start("TCMSP 캐시 빌드", total, remaining)

    status = {
        "status":       "running",
        "pid":          os.getpid(),
        "started_at":   datetime.now().isoformat(),
        "total":        total,
        "cached":       cached_count,
        "remaining":    remaining,
        "processed":    0,
        "success":      0,
        "fail":         0,
        "current_herb": "",
        "current_idx":  0,
        "eta_minutes":  None,
        "last_updated": "",
        "failed_herbs": [],
    }
    write_status(status)

    errors  = []
    t_start = time.time()

    import tcmsp_scraper as _scraper
    _scrape_count = 0   # 실제 스크래핑 횟수 (토큰 갱신 주기용)
    _consec_fail  = 0   # 연속 빠른 실패 횟수 (서버 다운 감지)
    _CONSEC_LIMIT = 5   # 이 이상 연속 실패 시 서버 다운으로 간주
    _WAIT_MINUTES = 10  # 서버 다운 시 대기 시간(분)

    for idx, (herb_en_name, herb_data) in enumerate(herbs, 1):
        pinyin = herb_data.get("herb_pinyin", "")
        cn     = herb_data.get("herb_cn_name", "")

        if already_cached(pinyin):
            log.info(f"[{idx}/{total}] 스킵: {pinyin} ({cn})")
            continue

        # 50개마다 토큰 강제 갱신 (세션 만료 방지)
        if _scrape_count > 0 and _scrape_count % 50 == 0:
            log.info("[토큰 갱신] 50개 처리 완료 — 세션 토큰 재발급")
            _scraper._shared_token = ""
            _scraper._shared_sess  = None

        log.info(f"[{idx}/{total}] 스크래핑: {pinyin} / {cn}")
        status["current_herb"] = f"{pinyin} ({cn})"
        status["current_idx"]  = idx
        write_status(status)

        t_req = time.time()
        try:
            result   = scrape_herb(pinyin, use_cache=True)
            req_time = time.time() - t_req
            n_active = result.get("active_compounds", 0)
            n_total  = result.get("total_compounds", 0)

            if result.get("error"):
                raise RuntimeError(result["error"])

            log.info(f"  -> 성공: 활성 {n_active}개 / 전체 {n_total}개")
            status["success"] += 1
            _scrape_count += 1
            _consec_fail   = 0  # 성공 시 연속 실패 카운터 초기화

        except Exception as e:
            req_time  = time.time() - t_req
            fast_fail = req_time < 1.0 and ("herb_not_found" in str(e) or "검색 결과 없음" in str(e))
            _consec_fail += 1
            log.error(f"  -> 실패: {e} (응답 {req_time:.2f}s, 연속실패 {_consec_fail})")

            # 연속 빠른 실패 → 서버 다운으로 판단, 대기 후 재시도
            if fast_fail and _consec_fail >= _CONSEC_LIMIT:
                log.warning(
                    f"[서버 다운 감지] {_consec_fail}회 연속 빠른 실패 — "
                    f"토큰 초기화 후 {_WAIT_MINUTES}분 대기"
                )
                _scraper._shared_token = ""
                _scraper._shared_sess  = None
                time.sleep(_WAIT_MINUTES * 60)
                _consec_fail = 0
                log.info("[재시도] 대기 완료 — 현재 약재 재시도")
                try:
                    result   = scrape_herb(pinyin, use_cache=True)
                    if result.get("error"):
                        raise RuntimeError(result["error"])
                    n_active = result.get("active_compounds", 0)
                    n_total  = result.get("total_compounds", 0)
                    log.info(f"  -> 재시도 성공: 활성 {n_active}개 / 전체 {n_total}개")
                    status["success"] += 1
                    _scrape_count += 1
                    status["processed"] += 1
                    status["cached"] = count_cached(herbs)
                    elapsed = time.time() - t_start
                    done    = status["processed"]
                    if done > 0:
                        eta_sec = (elapsed / done) * (remaining - done)
                        status["eta_minutes"] = round(eta_sec / 60)
                    write_status(status)
                    time.sleep(0.3)
                    continue
                except Exception as e2:
                    log.error(f"  -> 재시도도 실패: {e2}")
                    e = e2

            status["fail"] += 1
            err_entry = {
                "herb_en_name": herb_en_name,
                "pinyin":       pinyin,
                "cn_name":      cn,
                "error":        str(e),
                "timestamp":    datetime.now().isoformat(),
            }
            errors.append(err_entry)
            status["failed_herbs"].append(pinyin)
            write_errors(errors)

        status["processed"] += 1
        status["cached"]     = count_cached(herbs)

        # ETA 계산
        elapsed = time.time() - t_start
        done    = status["processed"]
        if done > 0:
            eta_sec = (elapsed / done) * (remaining - done)
            status["eta_minutes"] = round(eta_sec / 60)

        write_status(status)

        # 50개마다 요약 로그 + 텔레그램 알림
        if status["processed"] % 50 == 0:
            log.info(
                f"--- 진행: {done}/{remaining} | "
                f"성공 {status['success']} / 실패 {status['fail']} | "
                f"예상 {status['eta_minutes']}분 남음 ---"
            )
            notify_progress(
                "TCMSP 캐시 빌드",
                done, remaining,
                status["success"], status["fail"],
                status["eta_minutes"],
                status["current_herb"],
            )
        elif status["processed"] % 10 == 0:
            log.info(
                f"--- 진행: {done}/{remaining} | "
                f"성공 {status['success']} / 실패 {status['fail']} | "
                f"예상 {status['eta_minutes']}분 남음 ---"
            )

        time.sleep(0.3)

    # 완료
    elapsed_total = time.time() - t_start
    log.info(
        f"=== 완료: 성공 {status['success']}개 / "
        f"실패 {status['fail']}개 / 총 {elapsed_total/3600:.1f}시간 ==="
    )
    notify_done("TCMSP 캐시 빌드", status["success"], status["fail"], elapsed_total / 3600)
    status["status"]    = "completed"
    status["cached"]    = count_cached(herbs)
    status["eta_minutes"] = 0
    write_status(status)


if __name__ == "__main__":
    main()
