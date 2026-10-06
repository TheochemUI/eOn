#!/usr/bin/env bash
# Build libcpmdc.so against a pinned OpenCPMD commit and run the Si3N4
# isomer 1 point job. The energy line must match -1396.269526 eV.
set -euo pipefail
OPENCPMD_COMMIT=${OPENCPMD_COMMIT:-062582b7cfd832d36f88f504cd08e4ead42eb404}
CPMDC_COMMIT=${CPMDC_COMMIT:-8439c25cf0bd5f4caefd6ed07bebf142e2b181f5}
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

work=$STACK/point
mkdir -p "$work"
sed "s|@EXT_POT@|$KIT/bin/ext_pot|" "$KIT/configs/point/config.ini" > "$work/config.ini"
cp "$KIT/structures/si3n4_isomer1.con" "$work/pos.con"
export CPMD_BIN=$cpmd_bin
export CPMD_PP=$KIT/pseudo
export CPMD_LAUNCH=${CPMD_LAUNCH:-"mpirun -np 4 --bind-to none --map-by :OVERSUBSCRIBE"}
export CPMDC_LIBRARY=$STACK/libcpmdc.so
export OMP_NUM_THREADS=1
(cd "$work" && "$EONCLIENT" > eon.log 2>&1) || {
  echo "point job failed" >&2
  tail -40 "$work/eon.log" >&2 || true
  exit 1
}
python3 - "$work/results.dat" << 'PY'
import sys
from pathlib import Path
text = Path(sys.argv[1]).read_text()
found = None
for line in text.splitlines():
    parts = line.split()
    if len(parts) == 2 and parts[1] == "potential_energy":
        found = float(parts[0])
if found is None:
    sys.exit("results.dat has no potential_energy line")
target = -1396.269526
if abs(found - target) > 1e-6:
    sys.exit(f"potential_energy {found} is not {target}")
print(f"{found:.9f} potential_energy")
PY
