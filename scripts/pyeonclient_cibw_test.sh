#!/usr/bin/env bash
# CIBW_TEST_COMMAND for pyeonclient wheels. Base extra is rgpot-only;
# metatomic extra needs torch/metatensor .so dirs on LD_LIBRARY_PATH.
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
torch_lib = pathlib.Path(torch.__file__).resolve().parent / "lib"
dirs = [str(torch_lib)] if torch_lib.is_dir() else []
# metatensor-torch ships one .so per torch minor. Only the matching
# torch-X.Y tree is ABI-compatible.
for pat in ("libmetatensor*.so*", "libmetatomic*.so*"):
    for so in root.rglob(pat):
        p = str(so.resolve().parent).replace("\\\\", "/")
        if f"torch-{mm}" in p:
            dirs.append(str(so.resolve().parent))
print(":".join(dict.fromkeys(dirs)))
'
)"
if [ -n "${EXTRA_LIBS}" ]; then
  export LD_LIBRARY_PATH="${EXTRA_LIBS}:${LD_LIBRARY_PATH:-}"
  echo "cibw-test EXTRA_LIBS=${EXTRA_LIBS}"
fi
python -c "import pyeonclient as p, pathlib, subprocess; print(p.__version__, 'mta', p.built_with_metatomic(), 'rgpot', p.built_with_rgpot()); assert p.built_with_rgpot(); core=next(pathlib.Path(p.__file__).parent.glob('_core*.so')); out=subprocess.check_output(['readelf','-d',str(core)],text=True); assert 'libtorch' not in out or p.built_with_metatomic(); print('NEEDED_OK'); assert hasattr(p.NudgedElasticBand, 'path_frames')"
