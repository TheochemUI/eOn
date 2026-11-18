#!/usr/bin/env bash
# Serial force benchmark on the 4,000-atom Cu V/SIA cell used by amsel.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BENCH_BIN="${BENCH_BIN:-$ROOT/build-inprocess/client/cu_force_benchmark}"
EONCLIENT="${EONCLIENT:-$ROOT/build-inprocess/client/eonclient}"
TABLE="${TABLE:-$ROOT/examples/CECAM_examples/advance_topic/FCC_GB/potfiles/Cu01.eam.alloy}"
WORK="${WORK:-$ROOT/benchmark-cu}"
REPEATS="${REPEATS:-8}"

if ((REPEATS < 6)); then
  echo "REPEATS must be at least 6" >&2
  exit 2
fi
if [[ ! -x "$BENCH_BIN" || ! -f "$TABLE" ]]; then
  echo "BENCH_BIN must be executable and TABLE must name Cu01.eam.alloy" >&2
  exit 2
fi

mkdir -p "$WORK"
cp "$TABLE" "$WORK/Cu01.eam.alloy"
cat >"$WORK/in.lammps" <<'EOF'
units metal
pair_style eam/alloy
pair_coeff * * Cu01.eam.alloy Cu
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
EOF

# Same construction as amsel scripts/cu-vsia/build_catalog.py:
# a 10^3 FCC cell, one vacancy, and one split [100] dumbbell. Removing two
# sites and appending the two dumbbell atoms leaves exactly 4,000 atoms.
python3 - "$WORK/pos.con" <<'PY'
import sys

nc = 10
a0 = 3.615
d = 0.75
center = (10, 10, 10)
vacancy = (12, 10, 10)
sites = [
    (i, j, k)
    for i in range(2 * nc)
    for j in range(2 * nc)
    for k in range(2 * nc)
    if (i + j + k) % 2 == 0
]
positions = [
    tuple(a0 * q / 2 for q in site)
    for site in sites
    if site not in (center, vacancy)
]
site = tuple(a0 * q / 2 for q in center)
positions.extend(
    [(site[0] - d, site[1], site[2]), (site[0] + d, site[1], site[2])]
)
assert len(positions) == 4000
with open(sys.argv[1], "w", newline="\n") as out:
    out.write("Cu V/SIA benchmark\n\n")
    out.write(f"{nc*a0:.15f} {nc*a0:.15f} {nc*a0:.15f}\n")
    out.write("90.0 90.0 90.0\n\n\n1\n4000\n63.546\nCu\n")
    out.write("Coordinates of Component 1\n")
    for index, (x, y, z) in enumerate(positions):
        out.write(f"{x:.15f} {y:.15f} {z:.15f} 0 {index}\n")
PY

(
  cd "$WORK"
  for mode in worker library rgpot; do
    OMP_NUM_THREADS=1 "$BENCH_BIN" "$mode" pos.con "$REPEATS" \
      >"$mode.raw"
    awk -v mode="$mode" \
      '$0 ~ /^mode,/ || index($0, mode ",") == 1' \
      "$mode.raw" >"$mode.csv"
  done
  OMP_NUM_THREADS=1 "$BENCH_BIN" compare pos.con "$REPEATS" >parity.raw
  awk '/^(energy_difference|force_linf_difference)/' parity.raw >parity.csv
)

python3 - "$WORK" <<'PY'
import csv
import math
import pathlib
import statistics
import sys

work = pathlib.Path(sys.argv[1])
for mode in ("worker", "library", "rgpot"):
    rows = list(csv.DictReader((work / f"{mode}.csv").open()))
    values = [float(row["latency_ms"]) for row in rows]
    mean = statistics.mean(values)
    se = statistics.stdev(values) / math.sqrt(len(values))
    print(f"{mode}: n={len(values)} mean_ms={mean:.9f} se_ms={se:.9f}")
print((work / "parity.csv").read_text(), end="")
PY

cat >"$WORK/saddle-config.template" <<'EOF'
[Main]
job = process_search
temperature = 600
random_seed = SEED
write_log = false

[Potential]
potential = lammps

[Process Search]
minimize_first = true

[Optimizer]
converged_force = 0.01
max_iterations = 50

[Saddle Search]
min_mode_method = dimer
client_displace_type = random
displace_magnitude = 0.1
max_iterations = 10
max_energy = 20.0
EOF

echo "mode,seed,wall_s,force_calls" >"$WORK/saddle.csv"
for mode in worker rgpot; do
  for seed in 8003 8004 8005 8006 8007 8008; do
    run="$WORK/saddle-$mode-$seed"
    mkdir -p "$run"
    cp "$WORK/pos.con" "$WORK/in.lammps" "$WORK/Cu01.eam.alloy" "$run/"
    sed "s/SEED/$seed/" "$WORK/saddle-config.template" >"$run/config.ini"
    start=$(date +%s.%N)
    if [[ "$mode" == worker ]]; then
      (cd "$run" && OMP_NUM_THREADS=1 EON_LAMMPS_WORKER=1 \
        "$EONCLIENT" >run.log 2>&1)
    else
      (cd "$run" && OMP_NUM_THREADS=1 "$EONCLIENT" >run.log 2>&1)
    fi
    end=$(date +%s.%N)
    python3 - "$mode" "$seed" "$start" "$end" "$run/results.dat" \
      >>"$WORK/saddle.csv" <<'PY'
import pathlib
import sys

mode, seed, start, end, result = sys.argv[1:]
force_calls = None
for line in pathlib.Path(result).read_text().splitlines():
    fields = line.split()
    if len(fields) == 2 and fields[1] == "total_force_calls":
        force_calls = int(fields[0])
        break
if force_calls is None:
    raise RuntimeError("results.dat has no total_force_calls")
print(f"{mode},{seed},{float(end)-float(start):.12f},{force_calls}")
PY
  done
done

python3 - "$WORK/saddle.csv" <<'PY'
import csv
import math
import statistics
import sys

rows = list(csv.DictReader(open(sys.argv[1])))
for mode in ("worker", "rgpot"):
    values = [float(row["wall_s"]) for row in rows if row["mode"] == mode]
    mean = statistics.mean(values)
    se = statistics.stdev(values) / math.sqrt(len(values))
    print(f"{mode} saddle: n={len(values)} mean_s={mean:.9f} se_s={se:.9f}")
PY
