#!/usr/bin/env bash
# Timed multi-force loop for RgpotPot (requires potserv + engine).
# SocketNWChem path is only timed if NWCHEM_SOCKET_BENCH=1 and nwchem is available
# (full i-PI harness is heavier; default documents that path as untimed here).
set -euo pipefail
OUT="${1:-bench_rgpot_vs_socket_forces.txt}"
N_FORCE="${N_FORCE:-20}"
TEST_BIN="${TEST_BIN:-}"
HOST="${RGPOT_POTSERV_HOST:-127.0.0.1}"
PORT="${RGPOT_POTSERV_PORT:-19111}"

{
  echo "=== Timed force-loop comparison scaffold ==="
  echo "Date: $(date -Is)"
  echo "Host: $(hostname)"
  echo "N_FORCE=$N_FORCE"
  echo
  echo "--- RgpotPot (NWChem backend via potserv) ---"
  if [[ -z "$TEST_BIN" || ! -x "$TEST_BIN" ]]; then
    echo "SKIP: set TEST_BIN to test_rgpot_pot (built with -Dwith_rgpot=true)"
  else
    export RGPOT_POTSERV_HOST="$HOST" RGPOT_POTSERV_PORT="$PORT" RGPOT_BACKEND=NWChem
    # Warm-up
    "$TEST_BIN" "[rgpot]" >/dev/null 2>&1 || true
    start=$(date +%s.%N)
    for ((i=0; i<N_FORCE; i++)); do
      "$TEST_BIN" "[rgpot]" >/dev/null 2>&1
    done
    end=$(date +%s.%N)
    python3 - <<PY
start=float("$start"); end=float("$end"); n=int("$N_FORCE")
print(f"rgpot_force_loops={n}")
print(f"rgpot_wall_s_total={end-start:.6f}")
print(f"rgpot_wall_s_per_force_call_incl_process={((end-start)/n):.6f}")
print("note: each iteration re-executes test binary (process overhead dominates vs in-process Opt)")
PY
  fi
  echo
  echo "--- SocketNWChem ---"
  if [[ "${NWCHEM_SOCKET_BENCH:-0}" == "1" ]] && command -v nwchem >/dev/null; then
    echo "NWCHEM_SOCKET_BENCH requested; run client/unit_tests/run_nwchem_test.sh externally and paste times."
  else
    echo "NOT TIMED: requires nwchem + i-PI socket setup (SocketNWChemPotTest / run_nwchem_test.sh)."
    echo "Without that, cannot claim RGPOT faster than SocketNWChem for real NWChem optimize."
  fi
  echo
  echo "--- Optimize (minimization) ---"
  echo "NOT TIMED: no eonclient minimization job timed for both potentials on the same deck in this harness."
} | tee "$OUT"
