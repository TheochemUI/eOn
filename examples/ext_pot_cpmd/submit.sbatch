#!/bin/bash
#SBATCH --job-name=eon-cpmd
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=eon-cpmd-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
pixi_bin=${PIXI_BIN:-}
if [ -z "$pixi_bin" ]; then
    pixi_bin=$(command -v pixi || true)
fi
pixi_bin=${pixi_bin:-$HOME/.pixi/bin/pixi}
test -x "$pixi_bin" || { echo 'Pixi missing; install it or set PIXI_BIN to its absolute path' >&2; exit 2; }
exec "$pixi_bin" run --manifest-path "$PWD/pixi.toml" bash ./run.sh "$@"
