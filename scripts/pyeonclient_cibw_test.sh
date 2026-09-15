#!/usr/bin/env bash
# CIBW_TEST_COMMAND for pyeonclient wheels. Base extra is rgpot-only;
# metatomic extra needs torch/metatensor .so dirs on LD_LIBRARY_PATH.
set -euo pipefail
EXTRA_LIBS="$(
  python -c '
import pathlib, sys
# Only the active env. Walking sys.path also hits the manylinux image
# site-packages (torch-2.3 .. 2.14) and dlopens the wrong ABI.
root = pathlib.Path(sys.prefix)
dirs = []
for pat in (
    "libtorch.so*",
    "libmetatensor*.so*",
    "libmetatomic*.so*",
    "libc10.so*",
):
    for so in root.rglob(pat):
        dirs.append(str(so.resolve().parent))
print(":".join(dict.fromkeys(dirs)))
'
)"
if [ -n "${EXTRA_LIBS}" ]; then
  export LD_LIBRARY_PATH="${EXTRA_LIBS}:${LD_LIBRARY_PATH:-}"
  echo "cibw-test EXTRA_LIBS=${EXTRA_LIBS}"
fi
python -c "import pyeonclient as p, pathlib, subprocess; print(p.__version__, 'mta', p.built_with_metatomic(), 'rgpot', p.built_with_rgpot()); assert p.built_with_rgpot(); core=next(pathlib.Path(p.__file__).parent.glob('_core*.so')); out=subprocess.check_output(['readelf','-d',str(core)],text=True); assert 'libtorch' not in out or p.built_with_metatomic(); print('NEEDED_OK'); assert hasattr(p.NudgedElasticBand, 'path_frames')"
