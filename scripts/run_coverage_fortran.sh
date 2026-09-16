#!/usr/bin/env bash
# LCOV for Fortran potentials (client/potentials/* and fortcuh2) from a
# Meson coverage build tree, or standalone fortcuh2 instrumentation.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="${1:-fortran_lcov.info}"
BUILD="${EON_COV_BUILDDIR:-bbdir-cov}"

if ! command -v lcov >/dev/null 2>&1; then
  echo "ERROR: lcov not on PATH" >&2
  exit 1
fi

if [[ -d "$BUILD" ]]; then
  echo "==> extract Fortran from $BUILD"
  RAW=$(mktemp)
  # inih wrap builds unittest.c twice; lcov 2.x treats that as fatal
  # inconsistent unless ignored. Those TUs are not Fortran coverage.
  lcov --directory "$BUILD" --capture --output-file "$RAW" \
    --ignore-errors inconsistent,inconsistent \
    --exclude '*/subprojects/inih*' \
    2>/dev/null \
    || lcov --directory "$BUILD" --capture --output-file "$RAW" \
      --ignore-errors inconsistent,inconsistent \
      --exclude '*/subprojects/inih*'
  if lcov --extract "$RAW" \
    '*/client/potentials/*' '*/fortcuh2/*' '*/subprojects/fortcuh2/*' \
    --output-file "$OUT" 2>/dev/null; then
    :
  else
    # Fall back to .f90 SF records only
    python3 - "$RAW" "$OUT" <<'PY'
import sys
inp, outp = sys.argv[1], sys.argv[2]
keep = False
buf, out = [], []
for line in open(inp, encoding="utf-8", errors="replace"):
    if line.startswith("SF:"):
        if buf and keep:
            out.extend(buf)
        buf = [line]
        p = line[3:].strip().lower()
        keep = p.endswith((".f90", ".f", ".f95", ".for")) or "fortcuh2" in p or "/potentials/" in p
    elif line.startswith("end_of_record"):
        buf.append(line)
        if keep:
            out.extend(buf)
        buf, keep = [], False
    else:
        buf.append(line)
if buf and keep:
    out.extend(buf)
text = "".join(out)
if "SF:" not in text:
    raise SystemExit("no Fortran SF records in LCOV")
open(outp, "w", encoding="utf-8").write(text)
PY
  fi
  rm -f "$RAW"
else
  echo "==> standalone fortcuh2 instrumented build"
  command -v gfortran >/dev/null || { echo "need gfortran"; exit 1; }
  STG=$(mktemp -d)
  trap 'rm -rf "$STG"' EXIT
  SRC="$ROOT/subprojects/fortcuh2"
  if [[ ! -d "$SRC" ]]; then
    echo "    cloning fortcuh2 wrap source"
    git clone --depth 1 --branch v0.0.2       https://github.com/TheochemUI/fortran_cuh2_src.git "$SRC" || true
  fi
  [[ -d "$SRC" ]] || { echo "missing $SRC"; exit 1; }
  gfortran -c -O0 -g --coverage -cpp -J"$STG" -I"$STG" \
    "$SRC"/eam_dat.f90 -o "$STG"/eam_dat.o
  gfortran -c -O0 -g --coverage -cpp -J"$STG" -I"$STG" \
    "$SRC"/eam_isoc.f90 -o "$STG"/eam_isoc.o
  gfortran -c -O0 -g --coverage -cpp -J"$STG" -I"$STG" \
    "$SRC"/eamroutines.f90 -o "$STG"/eamroutines.o
  gfortran -O0 -g --coverage -cpp -J"$STG" -I"$STG" \
    "$SRC"/eam_dat.f90 "$SRC"/eam_isoc.f90 "$SRC"/eamroutines.f90 "$SRC"/main.f90 \
    -o "$STG"/fortcuh2_main 2>/dev/null || true
  if [[ -x "$STG/fortcuh2_main" ]]; then
    (cd "$STG" && ./fortcuh2_main) || true
  fi
  lcov --directory "$STG" --capture --output-file "$OUT"
  python3 - "$OUT" "$SRC" <<'PY'
import os, sys
outp, src = sys.argv[1], sys.argv[2]
lines = open(outp, encoding="utf-8", errors="replace").read().splitlines(True)
fixed = []
for line in lines:
    if line.startswith("SF:"):
        p = line[3:].strip()
        cand = os.path.join(src, os.path.basename(p))
        if os.path.isfile(cand):
            fixed.append(f"SF:{os.path.abspath(cand)}\n")
        else:
            fixed.append(line)
    else:
        fixed.append(line)
open(outp, "w", encoding="utf-8").writelines(fixed)
PY
fi

test -s "$OUT"
python3 - "$OUT" <<'PY'
import sys
hits = tot = 0
for line in open(sys.argv[1], encoding="utf-8", errors="replace"):
    if line.startswith("DA:"):
        h = int(line.split(":")[1].split(",")[1])
        tot += 1
        hits += h > 0
print(f"fortran lcov {100*hits/tot:.1f}% {hits}/{tot}" if tot else "empty")
if tot == 0:
    raise SystemExit(1)
PY
echo "OK wrote $OUT"
