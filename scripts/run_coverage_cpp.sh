#!/usr/bin/env bash
# Meson b_coverage + Catch2 unit tests → LCOV for client/ C++ sources.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="${1:-cpp_lcov.info}"
BUILD="${EON_COV_BUILDDIR:-bbdir-cov}"

if ! command -v meson >/dev/null 2>&1; then
  echo "ERROR: meson not on PATH (pixi shell / activate)" >&2
  exit 1
fi
if ! command -v lcov >/dev/null 2>&1; then
  echo "ERROR: lcov not on PATH" >&2
  exit 1
fi

echo "==> meson setup ${BUILD} (coverage, tests, fortran pots)"
rm -rf "$BUILD"
# Lean-ish default feature set used by ci_build_akmc, plus tests + coverage.
meson setup "$BUILD" \
  -Dwith_tests=true \
  -Dwith_fortran=true \
  -Dwith_cuh2=true \
  -Dwith_mpi=disabled \
  -Dwith_xtb=false \
  -Dwith_metatomic=false \
  -Dwith_gprd=disabled \
  -Dwith_serve=false \
  -Db_coverage=true \
  --buildtype=debug

meson compile -C "$BUILD"
# Soft-continue individual test failures; still collect gcov.
set +e
meson test -C "$BUILD" --print-errorlogs --timeout-multiplier 2
set -e

echo "==> lcov capture (C++ client)"
RAW=$(mktemp)
lcov --directory "$BUILD" --capture --output-file "$RAW" \
  --rc lcov_branch_coverage=1 2>/dev/null \
  || lcov --directory "$BUILD" --capture --output-file "$RAW"

TMP=$(mktemp)
if lcov --extract "$RAW" \
  '*/client/*.cpp' '*/client/*.cc' '*/client/*.hpp' '*/client/*.h' \
  '*/client/*/*.cpp' '*/client/*/*.cc' \
  --output-file "$TMP" 2>/dev/null; then
  # Prefer extract of client/ (exclude unit_tests via remove)
  lcov --remove "$TMP" '*/unit_tests/*' '*/thirdparty/*' '*/tests/*' \
    --output-file "$OUT" 2>/dev/null || mv "$TMP" "$OUT"
else
  lcov --remove "$RAW" '/usr/*' '*/.pixi/*' '*/unit_tests/*' '*/thirdparty/*' \
    --output-file "$OUT" || cp "$RAW" "$OUT"
fi
rm -f "$RAW" "$TMP"

test -s "$OUT"
python3 - "$OUT" <<'PY'
import sys
hits = tot = 0
for line in open(sys.argv[1], encoding="utf-8", errors="replace"):
    if line.startswith("DA:"):
        h = int(line.split(":")[1].split(",")[1])
        tot += 1
        hits += h > 0
print(f"cpp lcov {100*hits/tot:.1f}% {hits}/{tot}" if tot else "empty")
if tot == 0:
    raise SystemExit(1)
PY
echo "OK wrote $OUT"
