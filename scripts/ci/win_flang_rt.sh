# shellcheck shell=bash
# Locate conda-forge flang_rt for:
#   - MSVC link: import libs (FortranRuntime.dynamic.lib / flang_rt*.lib) via LIB
#   - LoadLibrary of pot plugins: runtime DLLs via PATH
#
# Source: source scripts/ci/win_flang_rt.sh
# Use:    eon_export_flang_rt "${CONDA_PREFIX:-}" --required

eon_win_to_unix_path() {
  local p="${1//\\//}"
  if [[ "$p" =~ ^([A-Za-z]):(.*)$ ]]; then
    local drive rest
    drive="$(echo "${BASH_REMATCH[1]}" | tr 'A-Z' 'a-z')"
    rest="${BASH_REMATCH[2]}"
    echo "/${drive}${rest}"
  else
    echo "$p"
  fi
}

eon_unix_to_win_path() {
  # MSVC link.exe reads LIB with Windows paths. MSYS converts PATH for native
  # children but does NOT rewrite LIB — Unix /d/... entries are ignored → LNK1104.
  local p="$1"
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -w "$p"
    return
  fi
  if [[ "$p" =~ ^/([a-zA-Z])/(.*)$ ]]; then
    local drive rest
    drive="$(echo "${BASH_REMATCH[1]}" | tr 'a-z' 'A-Z')"
    rest="${BASH_REMATCH[2]//\//\\}"
    echo "${drive}:\\${rest}"
  else
    echo "${p//\//\\}"
  fi
}

eon_dir_has() {
  # $1=dir $2=glob-ish name prefix (FortranRuntime or flang_rt) $3=ext (lib|dll)
  local d="$1" base="$2" ext="$3" f
  [ -d "$d" ] || return 1
  for f in "$d"/"${base}".dynamic."${ext}" \
           "$d"/"${base}".dynamic_dbg."${ext}" \
           "$d"/"${base}"."${ext}" \
           "$d"/"${base}"*."${ext}"
  do
    case "$f" in *\**) continue ;; esac
    [ -f "$f" ] && return 0
  done
  ls "$d"/${base}*."${ext}" >/dev/null 2>&1
}

eon_has_import_lib() {
  eon_dir_has "$1" FortranRuntime lib || eon_dir_has "$1" flang_rt lib
}

eon_has_runtime_dll() {
  eon_dir_has "$1" FortranRuntime dll || eon_dir_has "$1" flang_rt dll \
    || eon_dir_has "$1" FortranDecimal dll
}

eon_candidate_list() {
  local root res candidate
  local roots=()
  for root in "$@"; do
    [ -n "$root" ] || continue
    roots+=("$(eon_win_to_unix_path "$root")")
  done

  res=""
  if command -v flang-new >/dev/null 2>&1; then
    res="$(flang-new -print-resource-dir 2>/dev/null || true)"
  fi
  if [ -z "$res" ] && command -v flang >/dev/null 2>&1; then
    res="$(flang -print-resource-dir 2>/dev/null || true)"
  fi
  if [ -n "$res" ]; then
    res="$(eon_win_to_unix_path "$res")"
    echo "flang -print-resource-dir -> $res" >&2
    # Prefer the MSVC triple layout first (feedstock / clang resource dir).
    echo "${res}/lib/x86_64-pc-windows-msvc"
    echo "${res}/lib"
    echo "${res}"
  fi

  for root in "${roots[@]}"; do
    # shellcheck disable=SC2086
    for candidate in \
      "${root}/Library/lib/clang"/*/lib/x86_64-pc-windows-msvc \
      "${root}/lib/clang"/*/lib/x86_64-pc-windows-msvc \
      "${root}/Library/bin" \
      "${root}/Library/lib" \
      "${root}/bin" \
      "${root}/lib"
    do
      case "$candidate" in *\**) continue ;; esac
      echo "$candidate"
    done
  done
}

eon_find_first() {
  # $1=predicate function name; remaining args = roots
  local pred="$1"
  shift
  local c
  while IFS= read -r c; do
    [ -n "$c" ] || continue
    if "$pred" "$c"; then
      echo "$c"
      return 0
    fi
  done < <(eon_candidate_list "$@")
  # find(1) fallback under roots
  local root found name
  if [ "$pred" = eon_has_import_lib ]; then
    name='FortranRuntime.dynamic.lib'
  else
    name='FortranRuntime.dynamic.dll'
  fi
  for root in "$@"; do
    [ -n "$root" ] || continue
    root="$(eon_win_to_unix_path "$root")"
    [ -d "$root" ] || continue
    found="$(find "$root" -name "$name" 2>/dev/null | head -1 || true)"
    if [ -n "$found" ]; then
      dirname "$found"
      return 0
    fi
    # alternate names
    found="$(find "$root" \( -name 'flang_rt*.dll' -o -name 'flang_rt*.lib' \
      -o -name 'FortranRuntime*.dll' -o -name 'FortranRuntime*.lib' \) \
      2>/dev/null | head -1 || true)"
    if [ -n "$found" ]; then
      if [ "$pred" = eon_has_import_lib ] && [[ "$found" == *.lib ]]; then
        dirname "$found"
        return 0
      fi
      if [ "$pred" = eon_has_runtime_dll ] && [[ "$found" == *.dll ]]; then
        dirname "$found"
        return 0
      fi
    fi
  done
  return 1
}

eon_export_flang_rt() {
  # --required: must find import libs (MSVC link). Runtime DLLs are optional:
  # conda-forge flang often ships only .lib; pot plugins should link static
  # FortranRuntime (see client/meson.build _fortran_pot_link_args).
  local required=0
  local require_dll=0
  local roots=()
  local arg
  for arg in "$@"; do
    if [ "$arg" = "--required" ]; then
      required=1
    elif [ "$arg" = "--require-dll" ]; then
      require_dll=1
    elif [ -n "$arg" ]; then
      roots+=("$arg")
    fi
  done

  local lib_dir dll_dir
  lib_dir="$(eon_find_first eon_has_import_lib "${roots[@]}")" || lib_dir=""
  dll_dir="$(eon_find_first eon_has_runtime_dll "${roots[@]}")" || dll_dir=""

  # Library/bin often holds the runtime DLLs even when import libs live in Library/lib
  if [ -z "$dll_dir" ]; then
    local root
    for root in "${roots[@]}"; do
      root="$(eon_win_to_unix_path "$root")"
      if eon_has_runtime_dll "${root}/Library/bin"; then
        dll_dir="${root}/Library/bin"
        break
      fi
    done
  fi

  if [ -z "$lib_dir" ] && [ -z "$dll_dir" ]; then
    echo "WARNING: flang_rt import libs and runtime DLLs not found (roots: ${roots[*]:-none})" >&2
    if command -v flang-new >/dev/null 2>&1; then
      echo "flang-new=$(command -v flang-new)" >&2
      flang-new -print-resource-dir 2>&1 || true
    fi
    for root in "${roots[@]}"; do
      root="$(eon_win_to_unix_path "$root")"
      echo "find under $root:" >&2
      find "$root" \( -name 'FortranRuntime*' -o -name 'flang_rt*' \) 2>/dev/null | head -30 >&2 || true
    done
    if [ "$required" -eq 1 ]; then
      return 1
    fi
    return 0
  fi

  if [ -n "$lib_dir" ]; then
    export FLANG_RT_LIB_DIR="$lib_dir"
    local lib_win
    lib_win="$(eon_unix_to_win_path "$lib_dir")"
    export LIB="${lib_win};${LIB:-}"
    echo "Using flang_rt LIB: $lib_dir (win=$lib_win)"
    ls "$lib_dir"/FortranRuntime* "$lib_dir"/flang_rt* 2>/dev/null | head -20 || true
  else
    echo "WARNING: flang_rt import lib dir not found (link may LNK1104)" >&2
    if [ "$required" -eq 1 ]; then
      return 1
    fi
  fi

  if [ -n "$dll_dir" ]; then
    export FLANG_RT_DLL_DIR="$dll_dir"
    export PATH="${dll_dir}:${PATH}"
    echo "Using flang_rt PATH (runtime DLLs): $dll_dir"
    ls "$dll_dir"/FortranRuntime*.dll "$dll_dir"/flang_rt*.dll \
       "$dll_dir"/FortranDecimal*.dll 2>/dev/null | head -20 || true
  else
    echo "NOTE: no flang_rt runtime DLLs in env (pot plugins should use static flang_rt)" >&2
    if [ "$require_dll" -eq 1 ]; then
      return 1
    fi
  fi

  # Keep prefix bin on PATH (conda puts many runtime DLLs there).
  local root
  for root in "${roots[@]}"; do
    root="$(eon_win_to_unix_path "$root")"
    if [ -d "${root}/Library/bin" ]; then
      export PATH="${root}/Library/bin:${PATH}"
    fi
    if [ -d "${root}/bin" ]; then
      export PATH="${root}/bin:${PATH}"
    fi
  done

  export FLANG_RT_DIR="${dll_dir:-$lib_dir}"
  return 0
}
