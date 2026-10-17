#!/usr/bin/env bash
# Two endpoint structures, one shared cell, and one CPMD oracle for the band.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
reactant=$(readlink -f "${1:?usage: run_neb.sh REACTANT.xyz PRODUCT.xyz}")
product=$(readlink -f "${2:?usage: run_neb.sh REACTANT.xyz PRODUCT.xyz}")
set -a
source "${EON_ENV_FILE:-$here/.env}"
set +a
client=$(command -v "${EON_CLIENT:-eonclient}" || true)
test -x "$client" || { echo 'Activate eOn or set EON_CLIENT in .env' >&2; exit 2; }
work=$here/runs/neb-$(date +%Y%m%dT%H%M%S)-$$
mkdir -p "$work"
uv run --script "$here/prepare_neb.py" "$reactant" "$product" "$work"
python3 - "$here" "$work" <<'CONFIG'
from pathlib import Path
import sys
here, work = map(Path, sys.argv[1:])
text = (here / "config-neb.ini").read_text()
(work / "config.ini").write_text(text.replace("./ext_pot", str(here / "ext_pot")))
CONFIG
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EON_CLIENT_STANDALONE=1
printf 'NEB record: %s\n' "$work"
cd "$work"
"$client" > eon.log 2>&1 || { tail -80 eon.log >&2; exit 1; }
python3 - <<'CHECK'
from pathlib import Path
rows = [line.split(maxsplit=1) for line in Path("results.dat").read_text().splitlines()]
values = {p[1]:p[0] for p in rows if len(p)==2}
print(Path("results.dat").read_text())
if values.get("termination_reason_text") not in {"GOOD", "Success"}:
    raise SystemExit("NEB has not reported convergence; inspect eon.log and neb.dat")
print("Record:", Path.cwd())
CHECK
