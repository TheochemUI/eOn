#!/bin/bash
#SBATCH --job-name=eon-cpmd-point
#SBATCH --partition=grant
#SBATCH --account=g2025a360c
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=eon-cpmd-point-%j.out
# Submit from the eOn checkout with its Pixi environment active:
#   sbatch /path/to/run-cpmd-point.sbatch /absolute/path/to/pos.con
# The default structure is the public seven-silicon vacancy example.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
: "${SLURM_JOB_ID:?Submit this file with sbatch; do not run on the login node}"
EON_REPO=${EON_REPO:-$PWD}
EON_CLIENT=${EON_CLIENT:-$(command -v eonclient || true)}
CPMD_BIN=${CPMD_BIN:-$HOME/myCPMD/bin/cpmd4.3-oneAPI-2022-AMD.x}
PP_LIBRARY_PATH=${PP_LIBRARY_PATH:-$HOME/pseudo}
CPMD_RANKS=${CPMD_RANKS:-${SLURM_NTASKS:-4}}
source_structure=${1:-$EON_REPO/examples/akmc-cpmd-slurm/pos.con}
wrapper=$EON_REPO/examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py
if [ ! -x "$EON_CLIENT" ]; then
  echo 'eonclient missing: activate the installed eOn Pixi environment before sbatch, or set EON_CLIENT to its absolute path' >&2
  exit 2
fi
for file in "$CPMD_BIN" "$source_structure" "$wrapper"; do
  test -f "$file" || { echo "Required file missing: $file" >&2; exit 2; }
done
test -x "$CPMD_BIN" || { echo "CPMD is not executable: $CPMD_BIN" >&2; exit 2; }
test -d "$PP_LIBRARY_PATH" || { echo "Missing pseudopotential directory: $PP_LIBRARY_PATH" >&2; exit 2; }
CPMD_BIN=$(readlink -f "$CPMD_BIN")
PP_LIBRARY_PATH=$(readlink -f "$PP_LIBRARY_PATH")
EON_CLIENT=$(readlink -f "$EON_CLIENT")
work=$PWD/eon-cpmd-point-$SLURM_JOB_ID
mkdir "$work"
cp "$source_structure" "$work/pos.con"
cp "$wrapper" "$work/cpmd_extpot.py"
chmod +x "$work/cpmd_extpot.py"

# Only the CPMD child loads Intel MPI; eOn keeps its own Pixi libraries.
cat > "$work/cpmd-launch" <<'LAUNCH'
#!/bin/bash -l
set -euo pipefail
unset LD_LIBRARY_PATH LD_PRELOAD LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH
module purge -s
module load intel/oneAPI-2022
module load intelmpi/oneAPI-2022
mpirun --version >&2
exec mpirun -np "${CPMD_RANKS:?}" "$@"
LAUNCH
chmod +x "$work/cpmd-launch"
cat > "$work/config.ini" <<CONFIG
[Main]
job = point

[Potential]
potential = ext_pot
ext_pot_path = $work/cpmd_extpot.py
CONFIG
export CPMD_BIN PP_LIBRARY_PATH CPMD_RANKS
export CPMD_LAUNCH="$work/cpmd-launch"
# Public vacancy example defaults: silicon, BLYP, 20 Ry.
# Supply a comma-separated Z:file:options map for any other structure.
export CPMD_PP=${CPMD_PP:-14:Si_MT_BLYP.psp:LMAX=P}
export CPMD_CUTOFF=${CPMD_CUTOFF:-20}
export CPMD_FUNCTIONAL=${CPMD_FUNCTIONAL:-BLYP}
export CPMD_CONVERGENCE=${CPMD_CONVERGENCE:-1e-6}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export EON_CLIENT_STANDALONE=1
printf 'eOn: %s\nCPMD: %s\nPseudopotentials: %s\nRecord: %s\n' \
  "$EON_CLIENT" "$CPMD_BIN" "$PP_LIBRARY_PATH" "$work"
printf 'Units: E[eV] = E[Ha] * 27.211386245988\n'
printf 'Units: F[eV/A] = F[Ha/bohr] * 51.4220674763\n'
cd "$work"
"$EON_CLIENT" > eon.log 2>&1
python3 - "$work/results.dat" <<'VERIFY'
import math
from pathlib import Path
import sys
rows = [line.split(maxsplit=1) for line in Path(sys.argv[1]).read_text().splitlines()]
values = {p[1]: p[0] for p in rows if len(p) == 2}
energy = float(values["potential_energy"])
if not math.isfinite(energy):
    raise SystemExit("CPMD point returned a non-finite energy; inspect eon.log")
print(f"CPMD point energy: {energy:.9f} eV")
print("Input and output remain in", Path(sys.argv[1]).parent)
VERIFY
