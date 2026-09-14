# shellcheck shell=bash
# Windows GHA meson/test env: conda first for torch DLLs; MSVC via absolute tools.

_path_without_usr_bin() {
  local _out="" _p
  local IFS=':'
  for _p in $1; do
    [ "$_p" = "/usr/bin" ] && continue
    [ -z "$_p" ] && continue
    if [ -z "$_out" ]; then _out="$_p"; else _out="${_out}:${_p}"; fi
  done
  printf '%s' "$_out"
}

PATH="$(_path_without_usr_bin "$PATH")"
export PATH

if [ -n "${CONDA_PREFIX:-}" ]; then
  export PATH="${CONDA_PREFIX}/Scripts:${CONDA_PREFIX}:${CONDA_PREFIX}/Library/bin:${PATH}"
  _torch_lib="${CONDA_PREFIX}/Lib/site-packages/torch/lib"
  if [ -d "$_torch_lib" ]; then
    export PATH="${_torch_lib}:${PATH}"
  fi
fi

# Absolute MSVC tools from activate_msvc_gha.ps1
if [ -n "${EON_MSVC_CL:-}" ] && [ -x "$EON_MSVC_CL" ]; then
  export CC="$EON_MSVC_CL"
  export CXX="$EON_MSVC_CL"
  # Append (not prepend) MSVC bin so LoadLibrary still prefers conda
  export PATH="${PATH}:${EON_MSVC_BIN:-$(dirname "$EON_MSVC_CL")}"
else
  export CC=cl
  export CXX=cl
fi
if [ -n "${EON_MSVC_LINK:-}" ] && [ -x "$EON_MSVC_LINK" ]; then
  export LD="$EON_MSVC_LINK"
else
  export LD=link
fi

unset LDFLAGS || true
if [ -n "${FFLAGS:-}" ]; then
  FFLAGS="$(printf '%s' "$FFLAGS" | sed 's/-fuse-ld=lld//g')"
  export FFLAGS
fi
