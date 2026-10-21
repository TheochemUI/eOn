#!/usr/bin/env bash
# Build libcpmdc.so against a pinned OpenCPMD commit and run a point job.
# The default inputs are the Si3N4 isomer 1 kit: the energy line must match
# -1396.269526 eV within 1e-6 eV. Set BOND_LENGTHS to rebuild STRUCTURE at
# those H-H distances, check ENERGY_AT against ENERGY_EV, and check the
# quadratic bond-length minimum against BOND_REF_ANG.
set -euo pipefail
OPENCPMD_COMMIT=${OPENCPMD_COMMIT:-062582b7cfd832d36f88f504cd08e4ead42eb404}
CPMDC_COMMIT=${CPMDC_COMMIT:-95d779a839a2787e022148216baa331c2a6f412b}
OPENCPMD_REPO=${OPENCPMD_REPO:-https://github.com/OpenCPMD/CPMD.git}
CPMDC_REPO=${CPMDC_REPO:-https://github.com/OmniPotentRPC/cpmdc.git}
PREFIX=${PREFIX:-/usr}
STACK=${STACK:?set STACK to a build directory}
KIT=${KIT:?set KIT to the Si3N4 inputs}
EONCLIENT=${EONCLIENT:?set EONCLIENT to eonclient}
JOBS=${JOBS:-$(nproc)}
PATCHES=${PATCHES:-"embed_rinitwf embed_geometry embed_rwfopt converged_state kpoints_inputfile stopgm_return tistopgm c_mem_addrs embed_teardown"}
here=$(cd "$(dirname "$0")" && pwd)
src=$STACK/src/opencpmd
dest=$STACK/opencpmd
cpmdc=$STACK/cpmdc-src

if [ ! -d "$src/.git" ]; then
  mkdir -p "$STACK/src"
  git clone "$OPENCPMD_REPO" "$src"
fi
git -C "$src" fetch --depth 1 origin "$OPENCPMD_COMMIT"
git -C "$src" checkout --force --detach FETCH_HEAD
git -C "$src" clean -fd
test "$(git -C "$src" rev-parse HEAD)" = "$OPENCPMD_COMMIT"

if [ -n "${CPMDC_SRC:-}" ]; then
  cpmdc=$CPMDC_SRC
else
  if [ ! -d "$cpmdc/.git" ]; then
    git clone "$CPMDC_REPO" "$cpmdc"
  fi
  git -C "$cpmdc" fetch --depth 1 origin "$CPMDC_COMMIT"
  git -C "$cpmdc" checkout --detach FETCH_HEAD
fi
test "$(git -C "$cpmdc" rev-parse HEAD)" = "$CPMDC_COMMIT"

for p in $PATCHES; do
  patch -d "$src" -p1 --forward < "$cpmdc/tools/opencpmd_$p.patch"
done
# GCC 16 rejects KIND() of an imported BIND(C) enumerator in this file.
python3 - "$src/src/cuda_interfaces.mod.F90" << 'PY'
import pathlib, sys
path = pathlib.Path(sys.argv[1])
text = path.read_text()
text = text.replace(
    "INTEGER( KIND( cudaMemcpyHostToHost ) )",
    "INTEGER( C_INT )",
)
lines = []
for line in text.splitlines(True):
    if "IMPORT ::" in line and "cudaMemcpyHostToHost" in line and "C_INT" not in line:
        line = line.replace("cudaMemcpyHostToHost", "cudaMemcpyHostToHost, C_INT", 1)
    lines.append(line)
path.write_text("".join(lines))
PY
sed "s|@PREFIX@|$PREFIX|g" "$here/LINUX-CONDA-PIXI" > "$src/configure/LINUX-CONDA-PIXI"
(cd "$src" && ./configure.sh -DEST="$dest" LINUX-CONDA-PIXI)
make -C "$dest/obj" -f "$dest/Makefile" -j "$JOBS" "$dest/lib/libcpmd.a" timetag.o "$dest/bin/cpmd.x"
test -s "$dest/lib/libcpmd.a"
cpmd_bin=$(find "$dest" -name cpmd.x -type f | head -1)
test -x "$cpmd_bin"

build=$STACK/cpmdc-build
meson setup "$build" -Dwith_cpmd=true -Dcpmd_root="$dest" -Dwith_tests=false "$cpmdc"
meson compile -C "$build" -j "$JOBS"
so=$(find "$build" -name 'libcpmdc.so*' -type f | head -1)
test -n "$so"
cp -f "$so" "$STACK/libcpmdc.so"
test -s "$STACK/libcpmdc.so"

EXT_POT=${EXT_POT:-$KIT/bin/ext_pot}
CONFIG=${CONFIG:-$KIT/configs/point/config.ini}
STRUCTURE=${STRUCTURE:-$KIT/structures/si3n4_isomer1.con}
# Si3N4's ext_pot reads CPMD_PP as a directory. The in-repo H2 wrapper
# reads CPMD_PP as "Z:file[:options]" and PP_LIBRARY_PATH as the directory.
export CPMD_PP=${CPMD_PP:-$KIT/pseudo}
export PP_LIBRARY_PATH=${PP_LIBRARY_PATH:-$KIT/pseudo}
export CPMD_BIN=$cpmd_bin
export CPMD_LAUNCH=${CPMD_LAUNCH:-"mpirun -np 4 --bind-to none --map-by :OVERSUBSCRIBE"}
export CPMDC_LIBRARY=$STACK/libcpmdc.so
export OMP_NUM_THREADS=1
chmod +x "$EXT_POT"
test -x "$EXT_POT"
test -f "$CONFIG"
test -f "$STRUCTURE"

run_point() {
  local structure=$1
  local work=$2
  mkdir -p "$work"
  sed "s|@EXT_POT@|$EXT_POT|" "$CONFIG" > "$work/config.ini"
  cp "$structure" "$work/pos.con"
  (cd "$work" && "$EONCLIENT" > eon.log 2>&1) || {
    echo "point job failed ($work)" >&2
    tail -40 "$work/eon.log" >&2 || true
    exit 1
  }
}

if [ -z "${BOND_LENGTHS:-}" ]; then
  ENERGY_EV=${ENERGY_EV:--1396.269526}
  ENERGY_TOL_EV=${ENERGY_TOL_EV:-1e-6}
  work=$STACK/point
  run_point "$STRUCTURE" "$work"
  export POINT_RESULTS=$work/results.dat
  export ENERGY_EV ENERGY_TOL_EV
  python3 - << 'PY'
import os
import sys
from pathlib import Path
text = Path(os.environ["POINT_RESULTS"]).read_text()
found = None
for line in text.splitlines():
    parts = line.split()
    if len(parts) == 2 and parts[1] == "potential_energy":
        found = float(parts[0])
if found is None:
    sys.exit("results.dat has no potential_energy line")
target = float(os.environ["ENERGY_EV"])
tol = float(os.environ["ENERGY_TOL_EV"])
if abs(found - target) > tol:
    sys.exit(f"potential_energy {found} is not {target} within {tol}")
print(f"{found:.9f} potential_energy")
PY
else
  : "${ENERGY_EV:?set ENERGY_EV to the reference potential energy in eV}"
  : "${ENERGY_AT:?set ENERGY_AT to the bond length whose energy is checked}"
  : "${BOND_REF_ANG:?set BOND_REF_ANG to the reference H-H distance in Angstrom}"
  : "${BOND_TOL_ANG:?set BOND_TOL_ANG to the bond-length tolerance in Angstrom}"
  ENERGY_TOL_EV=${ENERGY_TOL_EV:-1e-3}
  geom=$STACK/geoms
  rm -rf "$geom" "$STACK"/point_r_*
  mkdir -p "$geom"
  export STRUCTURE GEOM_DIR=$geom BOND_LENGTHS
  python3 - << 'PY'
import math
import os
from pathlib import Path
template = Path(os.environ["STRUCTURE"]).read_text()
lines = template.splitlines()
start = None
for i, line in enumerate(lines):
    if line.startswith("Coordinates of Component"):
        start = i
        break
if start is None:
    raise SystemExit("pos.con has no coordinates")
atoms = []
j = start + 1
while j < len(lines) and lines[j].strip():
    parts = lines[j].split()
    if len(parts) < 3:
        break
    atoms.append((j, [float(parts[0]), float(parts[1]), float(parts[2])], parts[3:]))
    j += 1
if len(atoms) != 2:
    raise SystemExit(f"bond scan needs 2 atoms, found {len(atoms)}")
p1 = atoms[0][1]
p2 = atoms[1][1]
center = [(a + b) / 2 for a, b in zip(p1, p2)]
direction = [b - a for a, b in zip(p1, p2)]
norm = math.sqrt(sum(v * v for v in direction))
if norm < 1e-8:
    raise SystemExit("template bond length is zero")
direction = [v / norm for v in direction]
out_dir = Path(os.environ["GEOM_DIR"])
for raw in os.environ["BOND_LENGTHS"].split():
    distance = float(raw)
    if distance <= 0:
        raise SystemExit(f"bond length {raw} is not positive")
    out = list(lines)
    for sign, (line_no, _, tail) in zip((-0.5, 0.5), atoms):
        xyz = [c + sign * distance * d for c, d in zip(center, direction)]
        fields = [f"{v:.10f}" for v in xyz] + tail
        out[line_no] = " ".join(fields)
    path = out_dir / f"r_{distance:.4f}.con"
    path.write_text("\n".join(out) + "\n")
    print(f"wrote {path.name} bond {distance:.4f} A")
PY
  for con in "$geom"/r_*.con; do
    tag=$(basename "$con" .con)
    run_point "$con" "$STACK/point_${tag}"
  done
  export POINT_DIR=$STACK ENERGY_EV ENERGY_TOL_EV ENERGY_AT BOND_REF_ANG BOND_TOL_ANG
  python3 - << 'PY'
import os
import sys
from pathlib import Path

def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )

def value(text, key):
    found = None
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1] == key:
            found = float(parts[0])
    return found

root = Path(os.environ["POINT_DIR"])
points = []
for work in sorted(p for p in root.glob("point_r_*") if p.is_dir()):
    tag = work.name[len("point_"):]
    if not tag.startswith("r_"):
        continue
    distance = float(tag[2:])
    text = (work / "results.dat").read_text()
    energy = value(text, "potential_energy")
    force = value(text, "Max_Force")
    if energy is None:
        sys.exit(f"{work}/results.dat has no potential_energy line")
    points.append((distance, energy, force))
    print(f"{distance:.4f} {energy:.9f} potential_energy max_force {force}")
if len(points) < 3:
    sys.exit("bond scan needs at least 3 lengths")
# e = a (r-r0)^2 + b (r-r0) + c, with r0 the mean length so the
# normal matrix is not a Vandermonde in r ~ 1. The minimum is r0 - b/(2a).
r0 = sum(r for r, _e, _force in points) / len(points)
s = [0.0] * 5
t = [0.0, 0.0, 0.0]
for r, e, _force in points:
    x = r - r0
    p = 1.0
    for k in range(5):
        s[k] += p
        p *= x
    t[0] += e
    t[1] += e * x
    t[2] += e * x * x
matrix = [
    [s[0], s[1], s[2]],
    [s[1], s[2], s[3]],
    [s[2], s[3], s[4]],
]
denom = det3(matrix)
if abs(denom) < 1e-30:
    sys.exit("bond lengths are degenerate")

def column(index):
    copied = [row[:] for row in matrix]
    for i in range(3):
        copied[i][index] = t[i]
    return det3(copied) / denom

a = column(2)
b = column(1)
if a <= 0:
    sys.exit(f"bond scan curvature is not positive (a={a})")
minimum = r0 - b / (2.0 * a)
print(f"{minimum:.6f} bond_minimum_angstrom")
energy_at = float(os.environ["ENERGY_AT"])
energy_tol = float(os.environ["ENERGY_TOL_EV"])
energy_ref = float(os.environ["ENERGY_EV"])
bond_ref = float(os.environ["BOND_REF_ANG"])
bond_tol = float(os.environ["BOND_TOL_ANG"])
matched = [row for row in points if abs(row[0] - energy_at) <= 1e-6]
if len(matched) != 1:
    sys.exit(f"ENERGY_AT {energy_at} is not one of the scanned bond lengths")
found = matched[0][1]
if abs(found - energy_ref) > energy_tol:
    sys.exit(
        f"potential_energy {found} at {energy_at} A is not {energy_ref} within {energy_tol}"
    )
if abs(minimum - bond_ref) > bond_tol:
    sys.exit(
        f"bond minimum {minimum} A is not {bond_ref} within {bond_tol}"
    )
print(f"{found:.9f} potential_energy")
print(f"{minimum:.6f} bond_length_angstrom")
PY
fi
