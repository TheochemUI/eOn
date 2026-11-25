#!/usr/bin/env bash
# One Born-Oppenheimer point through the external CPMD executable.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
input=$(readlink -f "${1:-$here/pos.con}")
export EON_ENV_FILE=${EON_ENV_FILE:-$here/.env}
set -a
source "$EON_ENV_FILE"
set +a
client=$(command -v "${EON_CLIENT:-eonclient}" || true)
test -x "$client" || { echo 'eonclient missing; run pixi install or set EON_CLIENT in .env' >&2; exit 2; }
: "${CPMD_BIN:?Set CPMD_BIN in .env}"
: "${PP_LIBRARY_PATH:?Set PP_LIBRARY_PATH in .env}"
: "${CPMD_PP:?Set CPMD_PP species mapping in .env}"
: "${CPMD_MPI_LAUNCH:?Set CPMD_MPI_LAUNCH in .env}"
test -x "$CPMD_BIN" || { echo "CPMD not executable: $CPMD_BIN" >&2; exit 2; }
test -f "$input" || { echo "Input missing: $input" >&2; exit 2; }
mkdir -p "$here/runs"
work=$(mktemp -d "$here/runs/point-$(date +%Y%m%dT%H%M%S)-XXXXXX")
if [[ ${input,,} == *.xyz ]]; then
    uv run --script "$here/xyz_to_con.py" "$input" "$work/pos.con"
else
    uv run --script "$here/con_to_input.py" "$input" "$work/pos.con"
fi
cat > "$work/config.ini" <<CONFIG
[Main]
job = point

[Potential]
potential = ext_pot
ext_pot_path = $here/ext_pot
CONFIG
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EON_CLIENT_STANDALONE=1
printf 'Record: %s\nCPMD: %s\n' "$work" "$CPMD_BIN"
printf 'E[eV] = E[Hartree] * 27.211386245988\nF[eV/A] = F[Hartree/bohr] * 51.4220674763\n'
cd "$work"
"$client" > eon.log 2>&1 || { tail -60 eon.log >&2; exit 1; }
python3 - <<'CHECK'
import math
from pathlib import Path
rows = [line.split(maxsplit=1) for line in Path("results.dat").read_text().splitlines()]
values = {p[1]: p[0] for p in rows if len(p) == 2}
energy = float(values["potential_energy"])
if not math.isfinite(energy):
    raise SystemExit("Non-finite point energy; inspect eon.log")
print(f"Energy: {energy:.10f} eV")
print("Record:", Path.cwd())
CHECK
