#!/usr/bin/env bash
# One procedure for eonclient and libcpmdc.so.
# The run path of each binary must not contain a /home/ directory.
#
#   build-client.sh --check BIN [BIN...]
#   build-client.sh --sanitize BIN [BIN...]
#   build-client.sh STACK EON_URL EON_COMMIT CPMDC_SRC
#
# STACK is the build root. It is not read from a home scratch directory.
# CPMDC_SRC is a cpmdc checkout. STACK/opencpmd must already hold the
# configured OpenCPMD tree. The last form calls build_cpmdc.sh and
# build_eon.sh beside this file, then rewrites both run paths.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)

runpath_of() {
  readelf -d "$1" | sed -n 's/.*Library \(runpath\|rpath\): \[\(.*\)\]/\2/p' | head -1
}

check_one() {
  local rp
  rp=$(runpath_of "$1")
  if printf '%s\n' "$rp" | grep -F '/home/' >/dev/null; then
    echo "RUNPATH contains /home/: $1" >&2
    echo "$rp" >&2
    return 1
  fi
}

sanitize_one() {
  local rp kept=() part joined
  rp=$(runpath_of "$1")
  if [ -n "$rp" ]; then
    IFS=':' read -ra parts <<< "$rp"
    for part in "${parts[@]}"; do
      case "$part" in
        *'/home/'*) ;;
        *) kept+=("$part") ;;
      esac
    done
  fi
  if [ "${#kept[@]}" -eq 0 ]; then
    joined='$ORIGIN'
  else
    joined=$(IFS=:; printf '%s' "${kept[*]}")
  fi
  patchelf --set-rpath "$joined" "$1"
  check_one "$1"
}

case "${1:-}" in
  --check)
    shift
    [ "$#" -ge 1 ]
    status=0
    for bin in "$@"; do
      check_one "$bin" || status=1
    done
    exit "$status"
    ;;
  --sanitize)
    shift
    [ "$#" -ge 1 ]
    for bin in "$@"; do
      sanitize_one "$bin"
    done
    ;;
  *)
    [ "$#" -eq 4 ]
    stack=$1
    eon_url=$2
    eon_commit=$3
    cpmdc_src=$4
    jobs=${JOBS:-4}
    bash "$here/build_cpmdc.sh" "$cpmdc_src" "$stack" "$jobs"
    bash "$here/build_eon.sh" "$eon_url" "$eon_commit" "$stack" "$jobs"
    sanitize_one "$stack/eon/prefix/bin/eonclient"
    so=$(find "$stack/cpmdc" -name 'libcpmdc.so*' -type f | head -1)
    sanitize_one "$so"
    ;;
esac
