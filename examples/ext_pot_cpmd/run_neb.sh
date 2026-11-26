#!/usr/bin/env bash
# Two endpoint structures, one shared cell, and one CPMD oracle for the band.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
test $# -eq 2 || { echo 'usage: run_neb.sh REACTANT.xyz PRODUCT.xyz' >&2; exit 2; }
# BSD readlink has no -f, and GNU readlink -f accepts a missing file.
endpoint() {
  test -f "$1" || { echo "No endpoint file: $1" >&2; exit 2; }
  printf '%s/%s\n' "$(cd "$(dirname "$1")" && pwd -P)" "$(basename "$1")"
}
reactant=$(endpoint "$1")
product=$(endpoint "$2")
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
