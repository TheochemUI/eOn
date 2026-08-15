#!/usr/bin/env bash
# Build a pyeonclient wheel for the *current* interpreter (PyPI-style).
#
# Variants (env):
#   PYEONCLIENT_VARIANT=base        # default: with_rgpot, no torch link → PyPI
#   PYEONCLIENT_VARIANT=metatomic   # with_metatomic + with_rgpot; wheel tagged
#                                   # 0.X.Y+metatomic (GitHub artifacts / local;
#                                   # not flattened onto PyPI as 0.X.Y)
#
# Examples:
#   pip install -U build nanobind numpy meson ninja meson-python
#   PYEONCLIENT_VARIANT=base ./scripts/pyeonclient_build_wheel.sh
#   pip install torch metatomic-torch metatensor-torch vesin
#   PYEONCLIENT_VARIANT=metatomic ./scripts/pyeonclient_build_wheel.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
VARIANT="${PYEONCLIENT_VARIANT:-base}"

bash scripts/pyeonclient_prepare_pyproject.sh

cleanup() {
  # restore version rewrite if any
  if [[ -f pyproject-pyeonclient.toml.bak-version ]]; then
    mv -f pyproject-pyeonclient.toml.bak-version pyproject-pyeonclient.toml
    # re-activate after restore if still in prepare mode
    if [[ -f pyproject.eon-akmc.toml.bak ]]; then
      cp -a pyproject-pyeonclient.toml pyproject.toml
    fi
  fi
  bash scripts/pyeonclient_prepare_pyproject.sh restore || true
}
trap cleanup EXIT

python -m pip install -U pip build nanobind 'numpy>=1.26.4' meson ninja meson-python

# A conda or pixi environment ships its .pc files under $CONDA_PREFIX but
# leaves PKG_CONFIG_PATH empty, so pkg-config finds none of them. The base
# wheel builds with -Dwith_rgpot=true, which turns on rgpot's RPC arm and
# needs capnp-rpc; without this the build stops at
# subprojects/rgpot/CppCore/rgpot/rpc/meson.build.
if [ -n "${CONDA_PREFIX:-}" ] && [ -d "${CONDA_PREFIX}/lib/pkgconfig" ]; then
  export PKG_CONFIG_PATH="${CONDA_PREFIX}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
fi

# The readcon-core wrap builds Rust with whatever rustflags the host
# ~/.cargo/config.toml carries. Those are chosen for the host toolchain: a
# developer pointing rustc at an absolute linker ('-fuse-ld=/usr/bin/mold',
# say) makes cargo pass it to the conda cc that meson selects here, which
# rejects it, and the wheel build stops at libreadcon_core.a. A wheel is a
# redistributable artifact, so it is built without host-tuned flags either
# way. RUSTFLAGS present but empty is what tells cargo to ignore the config
# ones; unset would let them back in.
export RUSTFLAGS="${PYEONCLIENT_RUSTFLAGS:-}"

# readcon-core wrap needs cbindgen (same as manylinux before-all)
export PATH="${HOME}/.cargo/bin:/root/.cargo/bin:${PATH:-}"
if ! command -v cbindgen >/dev/null 2>&1; then
  if command -v cargo >/dev/null 2>&1; then
    cargo install cbindgen
  else
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    # shellcheck disable=SC1091
    source "${HOME}/.cargo/env"
    cargo install cbindgen
  fi
fi
command -v cbindgen
command -v rustc || true

# Metatomic builds: PEP 440 local version so filenames never clobber base 0.X.Y
# wheels in a shared wheelhouse. Local versions are NOT uploaded to PyPI.
if [[ "$VARIANT" == "metatomic" ]]; then
  cp -a pyproject-pyeonclient.toml pyproject-pyeonclient.toml.bak-version
  python3 - <<'PY'
from pathlib import Path
p = Path("pyproject-pyeonclient.toml")
t = p.read_text()
# version = "0.3.0" -> version = "0.3.0+metatomic"
import re
t2, n = re.subn(
    r'(?m)^version = "([^"+]+)"\s*$',
    r'version = "\1+metatomic"',
    t,
    count=1,
)
if n != 1:
    raise SystemExit(f"failed to rewrite version (n={n})")
p.write_text(t2)
print("version ->", re.search(r'^version = "(.*)"', t2, re.M).group(1))
PY
  # re-copy into active pyproject
  cp -a pyproject-pyeonclient.toml pyproject.toml
fi

SETUP_ARGS=(
  "-Dwith_pyeonclient=true"
  "-Dinstall_eon_server=false"
  "-Dwith_tests=false"
  "-Dwith_ase=false"
  "-Dwith_xtb=false"
  "-Dwith_serve=false"
  "-Dwith_fortran=true"
  "-Dwith_cuh2=true"
  "-Dwith_rgpot=true"
  "-Dwrap_mode=forcefallback"
  # Explicit: manylinux/PyPI wheels never need LTO; avoids GCC format LTO traps
  "-Db_lto=false"
)

case "$VARIANT" in
  base)
    SETUP_ARGS+=("-Dwith_metatomic=false")
    ;;
  metatomic)
    SETUP_ARGS+=("-Dwith_metatomic=true")
    if [[ -n "${TORCH_PATH:-}" ]]; then
      SETUP_ARGS+=("-Dtorch_path=${TORCH_PATH}" "-Dpip_metatomic=false")
    else
      SETUP_ARGS+=("-Dpip_metatomic=true")
    fi
    python -c "import torch; print('torch', torch.__version__)" || {
      echo "metatomic variant needs torch on PYTHONPATH (pip install torch)" >&2
      exit 2
    }
    ;;
  *)
    echo "unknown PYEONCLIENT_VARIANT=$VARIANT (base|metatomic)" >&2
    exit 2
    ;;
esac

if [[ -n "${PYEONCLIENT_SETUP_ARGS:-}" ]]; then
  # shellcheck disable=SC2206
  SETUP_ARGS+=($PYEONCLIENT_SETUP_ARGS)
fi

CFG=()
for a in "${SETUP_ARGS[@]}"; do
  CFG+=(--config-setting="setup-args=${a}")
done

rm -rf dist build
echo "building pyeonclient variant=${VARIANT} args=${SETUP_ARGS[*]}"
python -m build --wheel "${CFG[@]}" 2>&1
ls -la dist/

# Vendor non-system NEEDED libs (capnp/kj for RGPOT) and strip host RPATH.
# Same idea as auditwheel / torch: portable wheel, no build-host RUNPATH.
if [[ "${PYEONCLIENT_SKIP_REPAIR:-0}" != "1" ]]; then
  python -m pip install -U auditwheel patchelf 2>/dev/null || true
  bash scripts/pyeonclient_repair_wheel.sh dist/pyeonclient-*.whl
fi

python - <<'PY'
import glob, zipfile, sys, os, re
whls = sorted(glob.glob("dist/*.whl"))
if not whls:
    sys.exit("no wheel produced")
whl = whls[-1]
print("wheel", whl)
variant = os.environ.get("PYEONCLIENT_VARIANT", "base")
if variant == "metatomic" and "+metatomic" not in whl:
    sys.exit(f"metatomic wheel filename missing +metatomic local version: {whl}")
if variant == "base" and "+metatomic" in whl:
    sys.exit(f"base wheel must not use +metatomic local version: {whl}")
with zipfile.ZipFile(whl) as z:
    names = z.namelist()
    sample = [n for n in names if n.startswith("pyeonclient/")][:15]
    print("sample members", sample)
    caps = [n for n in names if re.search(r"libcapnp|libkj", n)]
    # Base+RGPOT wheels NEEDED capnp — after repair they must be vendored.
    if variant == "base" and not caps:
        # tolerate only if extension truly has no capnp NEEDED (read later)
        print("NOTE: no libcapnp/libkj members in wheel (repair may have been skipped)")
    else:
        print("vendored_capnp_kj", caps)
print("WHEEL_BUILD_OK", whl)
PY
