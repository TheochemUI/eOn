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
  # MinGW (ld, windres, gfortran) and Strawberry (Pure-Perl pkg-config,
  # gmake) break the MSVC CMake probe. Git usr/bin stays: sh.exe is
  # there, and MSVC link.exe is already earlier on PATH.
  local out="" p low
  local IFS=':'
  for p in $1; do
    [ -n "$p" ] || continue
    low=$(printf '%s' "$p" | tr '[:upper:]' '[:lower:]')
    case "$low" in
      *mingw*|*strawberry*) continue ;;
    esac
    if [ -z "$out" ]; then
      out="$p"
    else
      out="${out}:${p}"
    fi
  done
  printf '%s' "$out"
}
