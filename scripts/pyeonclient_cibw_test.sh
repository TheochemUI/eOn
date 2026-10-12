#!/usr/bin/env bash
# CIBW_TEST_COMMAND for pyeonclient wheels. Base extra is rgpot-only;
# metatomic extra needs torch/metatensor library dirs on the loader path.
# Linux uses LD_LIBRARY_PATH. Windows searches PATH, and torch DLLs live
# in torch/bin as well as torch/lib.
set -euo pipefail
EXTRA_LIBS="$(
  python -c '
import pathlib, sys
try:
    import torch
except ImportError:
    raise SystemExit(0)
mm = ".".join(torch.__version__.split("+")[0].split(".")[:2])
root = pathlib.Path(sys.prefix)
torch_root = pathlib.Path(torch.__file__).resolve().parent
dirs = []
for sub in ("lib", "bin"):
    d = torch_root / sub
    if d.is_dir():
        dirs.append(str(d))
for pat in (
    "libmetatensor*.so*",
    "libmetatomic*.so*",
    "metatensor*.dll",
    "metatomic*.dll",
    "libmetatensor*.dll",
    "libmetatomic*.dll",
):
    for so in root.rglob(pat):
        parent = str(so.resolve().parent)
        if f"torch-{mm}" in parent or so.suffix.lower() == ".dll":
            dirs.append(parent)
print("\n".join(dict.fromkeys(dirs)))
'
)"
if [ -n "${EXTRA_LIBS}" ]; then
  ld=""
  win=""
  while IFS= read -r dir; do
    [ -z "$dir" ] && continue
    ld="${ld:+$ld:}$dir"
    win="${win:+$win;}$dir"
  done <<< "${EXTRA_LIBS}"
  export LD_LIBRARY_PATH="${ld}:${LD_LIBRARY_PATH:-}"
  echo "cibw-test EXTRA_LIBS=${ld}"
  if [ "${OS:-}" = "Windows_NT" ]; then
    export PATH="${win};${PATH:-}"
    echo "cibw-test PATH_PREFIX=${win}"
  fi
fi
# Static Fortran archives still import the flang runtime. The wheel test
# venv does not ship those DLLs; the pixi prefix used to compile does.
if [ "${OS:-}" = "Windows_NT" ]; then
  _flang_bin="C:/eon-wheel-deps/.pixi/envs/default/Library/bin"
  if [ -d "${_flang_bin}" ]; then
    export PATH="${_flang_bin};${PATH:-}"
  fi
fi
python - <<'PY'
import os
import pathlib
import platform
import site
import subprocess

# Git bash mangles PATH, and Windows does not search LD_LIBRARY_PATH.
# Register the directories with the loader before importing the extension.
if os.name == "nt":
    dll_dirs = []
    flang = pathlib.Path(r"C:\eon-wheel-deps\.pixi\envs\default\Library\bin")
    if flang.is_dir():
        dll_dirs.append(flang)
    roots = [pathlib.Path(p) for p in site.getsitepackages()]
    try:
        import torch
    except ImportError:
        torch = None
    if torch is not None:
        roots.append(pathlib.Path(torch.__file__).resolve().parent)
    mm = None
    if torch is not None:
        mm = ".".join(torch.__version__.split("+")[0].split(".")[:2])
        print("torch", torch.__version__, flush=True)
    # metatensor-torch ships one DLL tree per torch minor. The first
    # add_dll_directory wins, so an older tree hides the one this
    # interpreter was built against and the import fails with a missing
    # procedure.
    needles = ("metatensor", "metatomic", "vesin", "torch", "c10", "libiomp", "libomp")
    for root in roots:
        if not root.is_dir():
            continue
        for dll in root.rglob("*.dll"):
            low = str(dll).lower().replace("\\", "/")
            if not any(n in low for n in needles):
                continue
            if "/torch-" in low and (mm is None or f"/torch-{mm}/" not in low):
                continue
            dll_dirs.append(dll.parent)
    seen = []
    for directory in dll_dirs:
        text = str(directory)
        if text in seen:
            continue
        seen.append(text)
        print("add_dll_directory", text, flush=True)
        os.add_dll_directory(text)

try:
    import pyeonclient as p
except ImportError:
    if os.name == "nt":
        hits = []
        for sp in site.getsitepackages():
            hits.extend(pathlib.Path(sp).rglob("_core*.pyd"))
        for hit in hits:
            print("core", hit, flush=True)
            probe = subprocess.run(
                ["dumpbin", "/dependents", str(hit)],
                text=True,
                capture_output=True,
            )
            print(probe.stdout)
            print(probe.stderr)
    raise

print(p.__version__, "mta", p.built_with_metatomic(), "rgpot", p.built_with_rgpot())
assert p.built_with_rgpot()
parent = pathlib.Path(p.__file__).parent
cores = list(parent.glob("_core*.so")) + list(parent.glob("_core*.pyd"))
assert cores, sorted(parent.iterdir())
if platform.system() != "Windows":
    out = subprocess.check_output(["readelf", "-d", str(cores[0])], text=True)
    assert "libtorch" not in out or p.built_with_metatomic()
    print("NEEDED_OK")
else:
    print("NEEDED_SKIP_WIN", cores[0].name)
assert hasattr(p.NudgedElasticBand, "path_frames")
PY
