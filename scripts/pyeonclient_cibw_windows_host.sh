#!/usr/bin/env bash
# Host toolchain for the Windows pyeonclient wheel, printed as GITHUB_ENV
# lines. Git bash on the runner is MSYS, so this script keeps a POSIX PATH
# while it runs and the Python helper writes the Windows form.
set -euo pipefail

prefix=${1:?pixi prefix}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

if [ ! -d "${prefix}/Library" ]; then
  echo "wheel host prefix is missing: ${prefix}" >&2
  find "${root}/scripts/ci/wheel-host/.pixi" -maxdepth 3 -type d >&2 || true
  exit 1
fi

# MinGW and Strawberry break the CMake probe (ld, windres, gmake rewriting
# /nologo, Pure-Perl pkg-config). Git usr/bin stays so sh.exe remains;
# ilammy already put MSVC link.exe ahead of GNU link.
PATH=$(printf '%s\n' "$PATH" | tr ':' '\n' \
  | grep -vE '[Mm]ingw|[Ss]trawberry' \
  | paste -sd:)
export PATH="${prefix}/Library/bin:${prefix}/bin:${PATH}"

if ! command -v flang-new >/dev/null 2>&1; then
  echo "flang-new not on PATH after adding ${prefix}" >&2
  exit 1
fi
export FC=flang-new
if ! command -v rustc >/dev/null 2>&1 || ! command -v cargo >/dev/null 2>&1; then
  echo "rustc/cargo missing under ${prefix}" >&2
  exit 1
fi
# cbindgen is not packaged. Install it into the host prefix so the
# readcon-core wrap can generate headers. A cached prefix may already
# have it.
if ! command -v cbindgen >/dev/null 2>&1; then
  cargo install --root "${prefix}/Library" cbindgen >&2
fi

# shellcheck source=scripts/ci/win_flang_rt.sh
source "${root}/scripts/ci/win_flang_rt.sh"
# Status lines stay off stdout. Only the Python helper may write GITHUB_ENV.
eon_export_flang_rt "$prefix" --required >&2

py=(python)
if ! command -v python >/dev/null 2>&1; then
  py=(py -3)
fi
"${py[@]}" "${root}/scripts/pyeonclient_cibw_windows_host.py" \
  --emit --prefix "$prefix"
