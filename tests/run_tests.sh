#!/bin/bash
# =====================================================================
# run_tests.sh -- tests with the real program (two processes, TCP)
# One honest run with --check, then one run per cheating point: the honest
# server must print "BLAME S<cheater>".
# =====================================================================
set -u
cd "$(dirname "$0")/.."
T=test_run
rm -rf $T && mkdir -p $T/base && cp params.txt $T/base/
PORT=9500
pass=0; fail=0

run2() {   # run2 <dir> <cmd> <cheating server: 0|1|-> <point> [args...]
    local dir=$1 cmd=$2 cs=$3 pt=$4; shift 4
    local a0=() a1=()
    [ "$cs" = 0 ] && a0=(--cheat "$pt")
    [ "$cs" = 1 ] && a1=(--cheat "$pt")
    PORT=$((PORT + 1))
    ./mpcprf p0 "$cmd" --dir "$dir" --port $PORT "${a0[@]}" "$@" > "$dir/out0.txt" 2>&1 &
    local pid=$!
    ./mpcprf p1 "$cmd" --dir "$dir" --port $PORT "${a1[@]}" "$@" > "$dir/out1.txt" 2>&1
    wait $pid
}
ok()  { echo "  PASS  $1"; pass=$((pass + 1)); }
bad() { echo "  FAIL  $1"; fail=$((fail + 1)); tail -4 "$2/out0.txt" "$2/out1.txt" | sed 's/^/        /'; }

echo "== init and one honest run =="
run2 $T/base init - - || { echo "init failed"; cat $T/base/out*.txt; exit 1; }
run2 $T/base run - - --n 3 --check
if grep -q "CHECK OK" $T/base/out0.txt && grep -q "CHECK OK" $T/base/out1.txt; then ok "honest run: shared and delivery mode give (k + x)^d"; else bad "honest run" $T/base; fi

echo "== a cheating server is blamed by the honest one =="
for case in "0 triple.psi" "1 mta.E" "0 mta.big" "1 mta.nu" "0 mta.Ep" "1 mta.dec" "0 mask.pow" "1 omega" "0 delivery.mul"; do
    set -- $case
    d=$T/cheat_$2_$1
    mkdir -p $d && cp $T/base/params.txt $T/base/s*.state $d/
    run2 $d run "$1" "$2" --n 1
    honest=$((1 - $1))
    if grep -q "BLAME S$1" "$d/out$honest.txt"; then ok "$2 by S$1 -> S$honest blames S$1"; else bad "$2 by S$1" $d; fi
done
d=$T/cheat_deliver
mkdir -p $d && cp $T/base/params.txt $T/base/s*.state $d/
run2 $d run 1 setup.deliver --n 1
if grep -q "BLAME S1" "$d/out0.txt"; then ok "wrong share to a writer by S1 -> the writer's complaint blames S1"; else bad "setup.deliver" $d; fi

d=$T/cheat_silent
mkdir -p $d && cp $T/base/params.txt $T/base/s*.state $d/
run2 $d run 1 silent --n 1 --deadline 3
if grep -q "round deadline.*BLAME S1" "$d/out0.txt"; then ok "S1 stops sending -> S0 blames S1 at the round deadline"; else bad "silent" $d; fi

echo
echo "run tests: $pass passed, $fail failed"
[ $fail = 0 ]
