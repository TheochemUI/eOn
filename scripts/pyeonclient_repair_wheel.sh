#!/usr/bin/env bash
# Vendor non-system shared libs into a pyeonclient wheel and strip host RPATH/RUNPATH.
# Torch-style portable wheel (same role as auditwheel for non-manylinux deps like capnp).
#
# Usage:
#   ./scripts/pyeonclient_repair_wheel.sh dist/pyeonclient-*.whl
# Env:
#   PYEONCLIENT_VENDOR_SEARCH_PATHS  colon-separated dirs to resolve missing libs
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <wheel.whl> [more.whl...]" >&2
  exit 2
fi

LIBS_DIR_NAME="${PYEONCLIENT_LIBS_DIR:-.pyeonclient.mesonpy.libs}"

is_system_lib() {
  local base="$1"
  case "$base" in
    linux-vdso.so*|ld-linux*.so*|ld-linux-x86-64.so*) return 0 ;;
    libc.so*|libm.so*|libmvec.so*|libdl.so*|librt.so*|libpthread.so*|libresolv.so*|libutil.so*) return 0 ;;
    libgcc_s.so*|libstdc++.so*|libgomp.so*) return 0 ;;
    libgfortran.so*|libquadmath.so*) return 0 ;;
    libz.so*|libbz2.so*|liblzma.so*|libzstd.so*) return 0 ;;
    libpython*.so*) return 0 ;;
    *) return 1 ;;
  esac
}

needed_libs() {
  readelf -d "$1" 2>/dev/null | sed -n 's/.*Shared library: \[\(.*\)\]/\1/p' || true
}

repair_one() {
  local whl_in="$1"
  [[ -f "$whl_in" ]] || { echo "missing wheel: $whl_in" >&2; exit 2; }
  whl_in="$(cd "$(dirname "$whl_in")" && pwd)/$(basename "$whl_in")"
  local work
  work="$(mktemp -d "${TMPDIR:-/tmp}/pyeon-repair.XXXXXX")"
  # shellcheck disable=SC2064
  trap "rm -rf '$work'" RETURN

  echo "repair: extracting $whl_in"
  python3 -m zipfile -e "$whl_in" "$work"
  local libs_dir="$work/$LIBS_DIR_NAME"
  mkdir -p "$libs_dir"

  local search="${PYEONCLIENT_VENDOR_SEARCH_PATHS:-}"
  if [[ -z "$search" ]]; then
    search="${LD_LIBRARY_PATH:-}:${LIBRARY_PATH:-}:/usr/local/lib:/usr/local/lib64:/usr/lib:/usr/lib64"
  fi
  # meson-python leaves the wrap cdylib under .mesonpy-*/ or cargo-target/.
  local found
  while IFS= read -r found; do
    [[ -z "$found" ]] && continue
    search="${search}:$(dirname "$found")"
  done < <(find "${PYEONCLIENT_BUILD_ROOT:-$PWD}" /project \
    \( -name 'libreadcon_core.so' -o -name 'libreadcon_core.so.*' \) \
    2>/dev/null | head -20)
  export LD_LIBRARY_PATH="$libs_dir:${search}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

  # Pull auditwheel / mesonpy lib packs into one dir so a single $ORIGIN rpath
  # covers hash-suffixed capnp/kj/gfortran that auditwheel already staged.
  local pack
  while IFS= read -r -d '' pack; do
    [[ "$(cd "$pack" && pwd)" == "$(cd "$libs_dir" && pwd)" ]] && continue
    local f
    while IFS= read -r -d '' f; do
      local base
      base="$(basename "$f")"
      if [[ ! -f "$libs_dir/$base" ]]; then
        echo "consolidate: $base <- $pack"
        cp -aL "$f" "$libs_dir/$base"
        chmod u+w "$libs_dir/$base" 2>/dev/null || true
      fi
    done < <(find "$pack" -maxdepth 1 -type f \( -name '*.so' -o -name '*.so.*' \) -print0)
  done < <(find "$work" -type d \( -name '*.libs' -o -name '.libs' -o -name '.pyeonclient.mesonpy.libs' \) -print0)

  local changed=1 round=0
  while [[ $changed -eq 1 && $round -lt 30 ]]; do
    changed=0
    round=$((round + 1))
    local so
    while IFS= read -r -d '' so; do
      local needed
      while IFS= read -r needed; do
        [[ -z "$needed" ]] && continue
        is_system_lib "$needed" && continue
        if [[ -f "$libs_dir/$needed" ]]; then
          continue
        fi
        local resolved=""
        resolved="$(ldd "$so" 2>/dev/null | awk -v n="$needed" '$1 == n { print $3; exit }')" || true
        if [[ -z "$resolved" || "$resolved" == "not" || ! -f "${resolved:-}" ]]; then
          local d
          local IFS=':'
          for d in $search; do
            [[ -z "$d" || ! -d "$d" ]] && continue
            if [[ -f "$d/$needed" ]]; then
              resolved="$d/$needed"
              break
            fi
          done
          unset IFS
        fi
        if [[ -z "${resolved:-}" || ! -f "$resolved" ]]; then
          # still missing — try find under search roots
          resolved="$(find ${search//:/ } -maxdepth 2 -name "$needed" 2>/dev/null | head -1 || true)"
        fi
        if [[ -z "${resolved:-}" || ! -f "$resolved" ]]; then
          echo "WARNING: cannot resolve NEEDED $needed (from $(basename "$so"))" >&2
          continue
        fi
        case "$resolved" in
          "$work"/*)
            # already inside wheel (e.g. other pack dir) — copy into libs_dir
            if [[ ! -f "$libs_dir/$needed" ]]; then
              echo "vendor-in-wheel: $needed <- $resolved"
              cp -aL "$resolved" "$libs_dir/$needed"
              chmod u+w "$libs_dir/$needed" 2>/dev/null || true
              changed=1
            fi
            continue
            ;;
        esac
        echo "vendor: $needed <- $resolved"
        cp -aL "$resolved" "$libs_dir/$(basename "$resolved")"
        # keep the NEEDED basename even if source had a different soname file
        if [[ "$(basename "$resolved")" != "$needed" ]]; then
          cp -aL "$resolved" "$libs_dir/$needed"
        fi
        chmod u+w "$libs_dir/$needed" 2>/dev/null || true
        changed=1
      done < <(needed_libs "$so")
    done < <(find "$work" -type f \( -name '*.so' -o -name '*.so.*' \) -print0 | sort -z)
  done

  # Strip host RPATH/RUNPATH → $ORIGIN-relative only (point everything at libs_dir)
  local so so_dir libs_abs rel new_rpath
  libs_abs="$(cd "$libs_dir" && pwd)"
  while IFS= read -r -d '' so; do
    so_dir="$(cd "$(dirname "$so")" && pwd)"
    rel="$(python3 -c "import os; print(os.path.relpath('$libs_abs', '$so_dir'))")"
    if [[ "$rel" == "." ]]; then
      new_rpath="\$ORIGIN"
    else
      new_rpath="\$ORIGIN/$rel:\$ORIGIN"
    fi
    echo "patchelf: $(basename "$so") -> rpath=$new_rpath"
    patchelf --set-rpath "$new_rpath" "$so"
  done < <(find "$work" -type f \( -name '*.so' -o -name '*.so.*' \) -print0)

  # Reject host absolute paths in RPATH/RUNPATH
  local bad=0
  while IFS= read -r -d '' so; do
    local rp
    rp="$(readelf -d "$so" 2>/dev/null | sed -n 's/.*Library \(runpath\|rpath\): \[\(.*\)\]/\2/p' || true)"
    if echo "$rp" | grep -E '/home/|/Users/|/opt/conda|/\.pixi/envs' >/dev/null 2>&1; then
      echo "ERROR: host path remains in RPATH of $so: $rp" >&2
      bad=1
    fi
  done < <(find "$work" -type f \( -name '*.so' -o -name '*.so.*' \) -print0)
  [[ $bad -eq 0 ]] || exit 1

  # Every non-system NEEDED of the extension must live in $libs_dir after consolidate.
  local core
  core="$(find "$work/pyeonclient" -name '_core*.so' | head -1 || true)"
  if [[ -n "$core" ]]; then
    local needed
    while IFS= read -r needed; do
      [[ -z "$needed" ]] && continue
      is_system_lib "$needed" && continue
      if [[ ! -f "$libs_dir/$needed" ]]; then
        echo "ERROR: non-system NEEDED not vendored: $needed" >&2
        ls -la "$libs_dir" | head -40 >&2 || true
        bad=1
      fi
    done < <(needed_libs "$core")
  fi
  [[ $bad -eq 0 ]] || exit 1

  if [[ -n "$core" ]] && needed_libs "$core" | grep -q 'libcapnp'; then
    local ncap
    ncap="$(find "$libs_dir" -name 'libcapnp*.so*' | wc -l)"
    if [[ "$ncap" -lt 1 ]]; then
      echo "ERROR: libcapnp NEEDED but not vendored" >&2
      exit 1
    fi
    echo "vendor-ok: $ncap capnp libraries in wheel"
  fi

  # Avoid SONAME collisions with pip rgpot (both ship a liblennard_jones.so).
  # Rename eOn pot libs to libeon_* and rewrite NEEDED across the wheel.
  local old_soname new_soname so
  declare -A soname_map=(
    [liblennard_jones.so]=libeon_lennard_jones.so
    [liblennard_jones_cluster.so]=libeon_lennard_jones_cluster.so
  )
  for old_soname in "${!soname_map[@]}"; do
    new_soname="${soname_map[$old_soname]}"
    if [[ -f "$libs_dir/$old_soname" ]]; then
      mv -f "$libs_dir/$old_soname" "$libs_dir/$new_soname"
      patchelf --set-soname "$new_soname" "$libs_dir/$new_soname"
      echo "soname: $old_soname -> $new_soname"
    fi
  done
  while IFS= read -r -d '' so; do
    for old_soname in "${!soname_map[@]}"; do
      new_soname="${soname_map[$old_soname]}"
      if needed_libs "$so" | grep -qx "$old_soname"; then
        patchelf --replace-needed "$old_soname" "$new_soname" "$so"
        echo "replace-needed: $(basename "$so") $old_soname -> $new_soname"
      fi
    done
  done < <(find "$work" -type f \( -name '*.so' -o -name '*.so.*' \) -print0)

  # Repack with zip (always portable; no wheel.cli dependency)
  python3 - <<PY
import zipfile
from pathlib import Path
root = Path("$work")
out = Path("$whl_in")
# write to temp then replace
tmp = out.with_suffix(".whl.tmp")
with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(root.rglob("*")):
        if p.is_file():
            z.write(p, p.relative_to(root).as_posix())
tmp.replace(out)
print("repacked", out, "bytes", out.stat().st_size)
with zipfile.ZipFile(out) as z:
    caps = [n for n in z.namelist() if "capnp" in n or "libkj" in n]
    print("vendored_members", caps)
    print("n_members", len(z.namelist()))
PY
  echo "REPAIR_OK $whl_in"
  trap - RETURN
  rm -rf "$work"
}

for w in "$@"; do
  repair_one "$w"
done
