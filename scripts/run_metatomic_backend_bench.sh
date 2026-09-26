#!/usr/bin/env bash
# Reproducible metatomic backend bench (fat · RGPOT · ASE) → docs JSON + figure.
#
# Intended entry point (from repo root)::
#
#   pixi run -e mta-bench mta-backend-bench
#
# What this does:
#   1. meson setup/compile fat pyeonclient (with_metatomic; RGPOT is always linked)
#   2. meson setup/compile ASE-safe pyeonclient (no C++ metatomic; TORCH_LIBRARY)
#   3. locate PET-MAD model + geometry under subprojects/gpr_optim/bench_data
#   4. run scripts/compare_metatomic_backends.py → docs JSON SSoT
#   5. run scripts/plot_metatomic_backend_bench.py → fig/generated/
#
# Env overrides:
#   EON_PET_MAD_MODEL, EON_PET_MAD_POS, RGPOT_METATOMIC_ENGINE
#   EON_MTA_BENCH_BUILD_FAT   (default: bbdir-mta-bench)
#   EON_MTA_BENCH_BUILD_ASE   (default: bbdir-mta-bench-ase)
#   EON_MTA_BENCH_SKIP_BUILD=1  skip meson if builds already present
#   EON_MTA_BENCH_JSON        output JSON path
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BUILD_FAT="${EON_MTA_BENCH_BUILD_FAT:-bbdir-mta-bench}"
BUILD_ASE="${EON_MTA_BENCH_BUILD_ASE:-bbdir-mta-bench-ase}"
JSON_OUT="${EON_MTA_BENCH_JSON:-$ROOT/docs/source/fig/data/metatomic_backend_bench.json}"
TORCH_VER="${EON_MTA_BENCH_TORCH:-2.10}"
PREFIX="${CONDA_PREFIX:-${PIXI_ENVIRONMENT_HOME:-}}"
if [[ -z "${PREFIX}" ]]; then
  echo "error: CONDA_PREFIX / PIXI_ENVIRONMENT_HOME unset; run via pixi run -e mta-bench" >&2
  exit 2
fi

log() { printf '==> %s\n' "$*"; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "error: missing command: $1" >&2
    exit 2
  }
}

need_cmd meson
need_cmd ninja
need_cmd python

ensure_gpr_optim() {
  if [[ -d "$ROOT/subprojects/gpr_optim/bench_data/petmad" ]]; then
    return 0
  fi
  if [[ -d "$ROOT/../gpr_optim" ]]; then
    log "rsync gpr_optim → subprojects/ (bench data + headers)"
    rsync -a --delete --exclude .git "$ROOT/../gpr_optim/" "$ROOT/subprojects/gpr_optim/"
    return 0
  fi
  echo "error: need subprojects/gpr_optim/bench_data/petmad (rsync gpr_optim or set EON_PET_MAD_MODEL)" >&2
  exit 2
}

meson_setup_compile() {
  local builddir="$1"
  shift
  local -a opts=("$@")
  if [[ "${EON_MTA_BENCH_SKIP_BUILD:-0}" == "1" && -f "$builddir/client/python/_core.abi3.so" ]]; then
    log "skip build (exists): $builddir"
    return 0
  fi
  if [[ ! -f "$builddir/build.ninja" ]]; then
    log "meson setup $builddir"
    # shellcheck disable=SC2086
    meson setup "$builddir" \
      --prefix="$PREFIX" \
      --libdir=lib \
      --buildtype=release \
      ${NATIVE_FILES:+$NATIVE_FILES} \
      "${opts[@]}"
  else
    log "meson reconfigure $builddir"
    meson setup --reconfigure "$builddir" \
      --prefix="$PREFIX" \
      --libdir=lib \
      --buildtype=release \
      ${NATIVE_FILES:+$NATIVE_FILES} \
      "${opts[@]}" || true
  fi
  log "meson compile -C $builddir"
  meson compile -C "$builddir"
}

# Optional mold/ccache native files when present (same as feature.metatomic.setupeon).
NATIVE_FILES=""
if [[ -f "$ROOT/nativeFiles/mold.ini" ]]; then
  NATIVE_FILES="--native-file $ROOT/nativeFiles/mold.ini"
fi
if [[ -f "$ROOT/nativeFiles/ccache_gnu.ini" ]]; then
  NATIVE_FILES="$NATIVE_FILES --native-file $ROOT/nativeFiles/ccache_gnu.ini"
fi

ensure_gpr_optim

log "build fat pyeonclient (metatomic + rgpot)"
meson_setup_compile "$BUILD_FAT" \
  -Dwith_metatomic=true \
  -Dpip_metatomic=true \
  -Dtorch_version="$TORCH_VER" \
  -Dwith_pyeonclient=true \
  -Dwith_tests=false \
  -Dwith_xtb=false

log "build ASE-safe pyeonclient (no C++ metatomic)"
meson_setup_compile "$BUILD_ASE" \
  -Dwith_metatomic=false \
  -Dpip_metatomic=false \
  -Dwith_pyeonclient=true \
  -Dwith_tests=false \
  -Dwith_xtb=false

# Resolve assets (repo-relative first).
MODEL="${EON_PET_MAD_MODEL:-}"
POS="${EON_PET_MAD_POS:-}"
ENGINE="${RGPOT_METATOMIC_ENGINE:-}"

if [[ -z "$MODEL" || ! -f "$MODEL" ]]; then
  for cand in \
    "$ROOT/subprojects/gpr_optim/bench_data/petmad/pet-mad-s-v1.5.0.pt" \
    "$ROOT/../gpr_optim/bench_data/petmad/pet-mad-s-v1.5.0.pt"; do
    if [[ -f "$cand" ]]; then MODEL="$cand"; break; fi
  done
fi
if [[ -z "$POS" || ! -f "$POS" ]]; then
  for cand in \
    "$ROOT/subprojects/gpr_optim/bench_data/petmad/d016_pos.con" \
    "${MODEL:+$(dirname "$MODEL")/d016_pos.con}"; do
    if [[ -n "$cand" && -f "$cand" ]]; then POS="$cand"; break; fi
  done
fi
if [[ -z "$ENGINE" || ! -f "$ENGINE" ]]; then
  for cand in \
    "$ROOT/$BUILD_FAT/client/libmetatomic_engine.so" \
    "$ROOT/bbdir/client/libmetatomic_engine.so"; do
    if [[ -f "$cand" ]]; then ENGINE="$cand"; break; fi
  done
fi

if [[ -z "$MODEL" || ! -f "$MODEL" ]]; then
  echo "error: PET-MAD model not found; set EON_PET_MAD_MODEL" >&2
  exit 2
fi
if [[ -z "$POS" || ! -f "$POS" ]]; then
  echo "error: geometry not found; set EON_PET_MAD_POS" >&2
  exit 2
fi
if [[ -z "$ENGINE" || ! -f "$ENGINE" ]]; then
  echo "error: libmetatomic_engine.so not found under $BUILD_FAT" >&2
  exit 2
fi

CORE_FAT="$ROOT/$BUILD_FAT/client/python/_core.abi3.so"
CORE_ASE="$ROOT/$BUILD_ASE/client/python/_core.abi3.so"
if [[ ! -f "$CORE_FAT" ]]; then
  echo "error: missing $CORE_FAT" >&2
  exit 2
fi
if [[ ! -f "$CORE_ASE" ]]; then
  echo "error: missing $CORE_ASE" >&2
  exit 2
fi

# Runtime library path: fat client tree + pixi env + torch/metatomic wheels.
export LD_LIBRARY_PATH="\
$ROOT/$BUILD_FAT/client:\
$ROOT/$BUILD_FAT/client/python:\
$ROOT/$BUILD_ASE/client:\
$ROOT/$BUILD_ASE/client/python:\
$PREFIX/lib:\
${LD_LIBRARY_PATH:-}"

# Prefer source package + inject fat _core for the parent process import.
export PYTHONPATH="$ROOT/client/python:$ROOT/$BUILD_FAT/client/python:${PYTHONPATH:-}"
export EON_ASE_PYEON_CORE="$CORE_ASE"
export EON_PET_MAD_MODEL="$MODEL"
export EON_PET_MAD_POS="$POS"
export RGPOT_METATOMIC_ENGINE="$ENGINE"

# Point pyeonclient at the fat extension for parent process (list_backends etc.).
# Workers pick ASE core via EON_ASE_PYEON_CORE.
if [[ ! -e "$ROOT/client/python/pyeonclient/_core.abi3.so" ]] || \
   [[ "$(readlink -f "$ROOT/client/python/pyeonclient/_core.abi3.so" 2>/dev/null || true)" != "$(readlink -f "$CORE_FAT")" ]]; then
  # Install a symlink for import without clobbering a real file permanently:
  # only if missing or already a symlink.
  if [[ -L "$ROOT/client/python/pyeonclient/_core.abi3.so" || ! -e "$ROOT/client/python/pyeonclient/_core.abi3.so" ]]; then
    ln -sfn "$CORE_FAT" "$ROOT/client/python/pyeonclient/_core.abi3.so"
  fi
fi

mkdir -p "$(dirname "$JSON_OUT")"
log "compare backends → $JSON_OUT"
python "$ROOT/scripts/compare_metatomic_backends.py" \
  --model "$MODEL" \
  --pos "$POS" \
  --engine "$ENGINE" \
  --json "$JSON_OUT" \
  --repeat "${EON_MTA_BENCH_REPEAT:-3}"

log "plot figure from JSON"
python "$ROOT/scripts/plot_metatomic_backend_bench.py" \
  --json "$JSON_OUT"

log "done"
log "  JSON:  $JSON_OUT"
log "  figure: $ROOT/docs/source/fig/generated/metatomic_backend_bench.svg"
log "  model:  $MODEL"
log "  engine: $ENGINE"
