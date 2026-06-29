#!/usr/bin/env bash
# Compare setup notes: SocketNWChem vs RgpotPot(NWChem). Times are wall-clock
# for a single force call harness when servers are up.
set -euo pipefail
OUT="${1:-rgpot_nwchem_vs_socket.txt}"
{
  echo "=== Rgpot NWChem vs SocketNWChem (setup comparison) ==="
  echo "Date: $(date -Is)"
  echo
  echo "SocketNWChem: eOn listens; NWChem connects (i-PI). See"
  echo "  client/potentials/SocketNWChem/ and unit_tests/run_nwchem_test.sh"
  echo "  Requires nwchem on PATH and a .nwi with driver socket."
  echo
  echo "RgpotPot NWChem: potserv hosts NWChemPot (dlopen libnwchemc)."
  echo "  Start: NWCHEMC_LIBRARY=... potserv <port> NWChem"
  echo "  eOn: potential=RGPOT backend=NWChem host/port -> Cap'n Proto configure+calculate"
  echo "  Fake CI path: libnwchemc_fake_engine.so + tests/test_rpc_e2e style"
  echo
  if [[ -n "${RGPOT_POTSERV_PORT:-}" ]]; then
    echo "RGPOT_POTSERV_PORT=${RGPOT_POTSERV_PORT} (run test_rgpot_pot for timings)"
  else
    echo "Set RGPOT_POTSERV_HOST/PORT and run: meson test -C <build> test_rgpot_pot --suite rgpot"
  fi
} | tee "$OUT"
