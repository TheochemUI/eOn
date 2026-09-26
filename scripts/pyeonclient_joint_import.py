#!/usr/bin/env python3
"""Import pip rgpot and a bundled pyeonclient wheel in one process.

pyeonclient vendors librgpot as libeon_rgpot.so.3. The rgpot wheel keeps
the public librgpot SONAME. Both imports must succeed, and neither
extension may satisfy the other's NEEDED entry.
"""

from __future__ import annotations

import pathlib
import re
import subprocess
import sys


def needed(path: pathlib.Path) -> str:
    return subprocess.check_output(["readelf", "-d", str(path)], text=True)


def main() -> int:
    import pyeonclient
    import rgpot

    py_core = next(pathlib.Path(pyeonclient.__file__).parent.rglob("_core*.so"))
    rg_core = next(pathlib.Path(rgpot.__file__).parent.rglob("_core*.so"))
    py_need = needed(py_core)
    rg_need = needed(rg_core)
    print(py_core)
    print(py_need)
    print(rg_core)
    print(rg_need)
    if "libeon_rgpot.so.3" not in py_need:
        print("pyeonclient did not bundle libeon_rgpot.so.3", file=sys.stderr)
        return 1
    if re.search(r"\[librgpot\.so", py_need):
        print("pyeonclient still NEEDs the pip rgpot SONAME", file=sys.stderr)
        return 1
    if not re.search(r"\[librgpot\.so", rg_need):
        print("rgpot wheel lost the librgpot SONAME", file=sys.stderr)
        return 1
    print("JOINT_IMPORT_OK", pyeonclient.__version__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
