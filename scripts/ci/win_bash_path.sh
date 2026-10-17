# shellcheck shell=bash
# Git bash stores PATH with ':' separators. A Windows C:\ entry is split
# on the drive colon, so every prepend has to be a POSIX path (/c/...).

eon_posix_path() {
  local p="$1"
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -u -- "$p"
    return
  fi
  printf '%s\n' "$p"
}

eon_drop_cmake_poison() {
  # MinGW (ld, windres, gfortran), Git's ucrt64 toolchain, and Strawberry
  # (Pure-Perl pkg-config, gmake) break the MSVC CMake probe. Git usr/bin
  # stays so sh.exe remains, but Git bash prepends it and its link.exe is
  # GNU. Place that directory immediately after the MSVC linker.
  local -a kept=() deferred=()
  local p low msvc_at=-1 i out=""
  local IFS=':'
  for p in $1; do
    [ -n "$p" ] || continue
    low=$(printf '%s' "$p" | tr '[:upper:]' '[:lower:]')
    case "$low" in
      *mingw*|*strawberry*|*ucrt64*) continue ;;
    esac
    case "$low" in
      */git/usr/bin|/usr/bin)
        # Git bash mounts that directory as /usr/bin and prepends it.
        deferred+=("$p")
        continue
        ;;
    esac
    case "$low" in
      *vc/tools/msvc*hostx64/x64) msvc_at=${#kept[@]} ;;
    esac
    kept+=("$p")
  done
  if [ "$msvc_at" -ge 0 ] && [ "${#deferred[@]}" -gt 0 ]; then
    local -a merged=()
    for ((i = 0; i <= msvc_at; i++)); do
      merged+=("${kept[i]}")
    done
    for p in "${deferred[@]}"; do
      merged+=("$p")
    done
    for ((i = msvc_at + 1; i < ${#kept[@]}; i++)); do
      merged+=("${kept[i]}")
    done
    kept=("${merged[@]}")
  elif [ "${#deferred[@]}" -gt 0 ]; then
    kept+=("${deferred[@]}")
  fi
  for p in "${kept[@]}"; do
    if [ -z "$out" ]; then
      out="$p"
    else
      out="${out}:${p}"
    fi
  done
  printf '%s' "$out"
}
