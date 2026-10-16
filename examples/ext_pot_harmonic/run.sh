#!/usr/bin/env bash
set -euo pipefail
example=$(cd "$(dirname "$0")" && pwd)
work=${1:-"$example/runs/point-$(date +%Y%m%dT%H%M%S)-$$"}
mkdir -p "$work"
work=$(cd "$work" && pwd)
if [ -f "$work/results.dat" ]; then
  echo "Result exists in $work; select a different work directory" >&2
  exit 2
fi
cp "$example/pos.con" "$example/config.ini" "$example/harmonic_extpot.py" "$work/"
chmod +x "$work/harmonic_extpot.py"
cd "$work"
pixi run --manifest-path "$example/pixi.toml" eonclient > eon.log 2>&1
pixi run --manifest-path "$example/pixi.toml" python "$example/check_result.py" "$work/results.dat"
printf 'Record: %s\n' "$work"
