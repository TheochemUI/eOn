#!/usr/bin/env bash
# CIBW_TEST_COMMAND for pyeonclient wheels. Base extra is rgpot-only;
# metatomic extra needs torch/metatensor .so dirs on LD_LIBRARY_PATH.
set -euo pipefail
EXTRA_LIBS="$(
  python -c '
import importlib, pathlib
dirs = []
for n in ("torch", "metatensor", "metatensor.torch", "metatomic.torch"):
    try:
        m = importlib.import_module(n)
        p = pathlib.Path(m.__file__).resolve().parent
        dirs += [str(p), str(p / "lib"), str(p.parent)]
    except Exception:
        pass
print(":".join(dict.fromkeys(d for d in dirs if pathlib.Path(d).is_dir())))
' 2>/dev/null || true
)"
if [ -n "${EXTRA_LIBS}" ]; then
  export LD_LIBRARY_PATH="${EXTRA_LIBS}:${LD_LIBRARY_PATH:-}"
fi
python -c "import pyeonclient as p, pathlib, subprocess; print(p.__version__, 'mta', p.built_with_metatomic(), 'rgpot', p.built_with_rgpot()); assert p.built_with_rgpot(); core=next(pathlib.Path(p.__file__).parent.glob('_core*.so')); out=subprocess.check_output(['readelf','-d',str(core)],text=True); assert 'libtorch' not in out or p.built_with_metatomic(); print('NEEDED_OK'); assert hasattr(p.NudgedElasticBand, 'path_frames')"
