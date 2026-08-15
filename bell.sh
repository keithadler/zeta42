#!/bin/zsh
# The 6-hour bell for the T6 aggregate, then the T7 sequence.
#
#   1. wait until the T6 kissat hits 6h total OR delivers a verdict
#   2. record/kill accordingly
#   3. run kissat on the T7 directed-pair instance (the B-question; no cap)
#   4. then kissat on the T7 aggregate instance
#
# All output flushed to bell.log; every verdict also appended there.

cd /Users/keith/cnp
log() { echo "[$(date '+%H:%M:%S')] $*" >> bell.log }

BELL_SECONDS=21600

log "bell armed: T6 kissat cap ${BELL_SECONDS}s"
while true; do
  PID=$(pgrep -x kissat | head -1)
  if [ -z "$PID" ]; then
    V=$(grep -m1 -E 's (UN)?SATISFIABLE' t6core_kissat.log 2>/dev/null)
    log "T6 kissat finished on its own: ${V:-no verdict line}"
    break
  fi
  ET=$(ps -o etimes= -p "$PID" | tr -d ' ')
  if [ "$ET" -ge "$BELL_SECONDS" ]; then
    kill "$PID"
    log "BELL: T6 aggregate killed at ${ET}s without verdict (open, expensive)"
    break
  fi
  sleep 60
done

# wait for the pair instance export if it is still being written
while [ ! -s t7pair_dial.cnf ]; do sleep 30; done
sleep 5

log "launching kissat on t7pair_dial.cnf (directed B-question, no cap)"
kissat -q t7pair_dial.cnf > t7pair_kissat.log 2>&1
V=$(grep -m1 -E 's (UN)?SATISFIABLE' t7pair_kissat.log)
log "T7 PAIR VERDICT: ${V:-none}"

log "launching kissat on t7core_dial.cnf (aggregate)"
kissat -q t7core_dial.cnf > t7core_kissat.log 2>&1
V=$(grep -m1 -E 's (UN)?SATISFIABLE' t7core_kissat.log)
log "T7 AGGREGATE VERDICT: ${V:-none}"
log "bell sequence complete"
