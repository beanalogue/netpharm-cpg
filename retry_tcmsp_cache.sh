#!/bin/bash
# TCMSP 캐시 재시도 스크립트
# - 서버 접속 가능 시 빌드 자동 실행
# - 195개 모두 완료 시 텔레그램 알림 후 크론 자체 삭제

BASE="/home/beanalogue/netpharm"
LOG="$BASE/logs/retry_tcmsp_cron.log"
STATUS="$BASE/logs/build_tcmsp_cache_status.json"
PY="/usr/bin/python3"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 크론 실행" >> "$LOG"

# 1. 이미 빌드 실행 중이면 스킵
if pgrep -f "build_tcmsp_cache.py" > /dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 이미 실행 중 — 스킵" >> "$LOG"
    exit 0
fi

# 2. 실패 약재 남아있는지 확인
FAILED=$(python3 -c "
import json
try:
    d = json.load(open('$STATUS'))
    print(len(d.get('failed_herbs', [])))
except:
    print(999)
" 2>/dev/null)

if [ "$FAILED" = "0" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 실패 약재 없음 — 크론 삭제" >> "$LOG"
    crontab -l 2>/dev/null | grep -v "retry_tcmsp_cache.sh" | crontab -
    exit 0
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 실패 약재 ${FAILED}개 남음" >> "$LOG"

# 3. TCMSP 서버 접속 확인
HTTP=$(curl -s --max-time 10 "https://old.tcmsp-e.com/" -o /dev/null -w "%{http_code}" 2>/dev/null)
if [ "$HTTP" != "200" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 서버 불가 (HTTP $HTTP) — 다음 시간 재시도" >> "$LOG"
    exit 0
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 서버 정상 (HTTP $HTTP) — 빌드 시작" >> "$LOG"

# 4. 빌드 실행
BUILD_NUM=$(ls $BASE/logs/build_tcmsp_cache*.log 2>/dev/null | wc -l)
BUILD_NUM=$((BUILD_NUM + 1))
BUILDLOG="$BASE/logs/build_tcmsp_cache${BUILD_NUM}.log"

cd "$BASE"
nohup $PY build_tcmsp_cache.py > "$BUILDLOG" 2>&1 &
BPID=$!
echo "[$(date '+%Y-%m-%d %H:%M:%S')] 빌드 시작 PID=$BPID 로그=$BUILDLOG" >> "$LOG"

# 5. 빌드 완료 대기 후 결과 확인
wait $BPID

SUCCESS=$(python3 -c "
import json
try:
    d = json.load(open('$STATUS'))
    print(d.get('success', 0))
except:
    print(0)
" 2>/dev/null)

FAIL=$(python3 -c "
import json
try:
    d = json.load(open('$STATUS'))
    print(len(d.get('failed_herbs', [])))
except:
    print(999)
" 2>/dev/null)

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 빌드 완료 — 성공 ${SUCCESS}개, 실패 ${FAIL}개" >> "$LOG"

# 6. 텔레그램 알림
cd "$BASE"
$PY -c "
import sys; sys.path.insert(0, 'src')
from telegram_notify import send_message
success = $SUCCESS
fail = $FAIL
if fail == 0:
    send_message(f'✅ TCMSP 캐시 빌드 완전 완료!\n성공: {success}개 / 실패: 0개\n500개 전체 캐시 완성!')
else:
    send_message(f'🔄 TCMSP 캐시 재시도 완료\n성공: {success}개 / 실패: {fail}개 남음\n다음 시간에 자동 재시도')
" 2>/dev/null

# 7. 모두 성공했으면 크론 삭제
if [ "$FAIL" = "0" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 전체 완료 — 크론 삭제" >> "$LOG"
    crontab -l 2>/dev/null | grep -v "retry_tcmsp_cache.sh" | crontab -
fi
