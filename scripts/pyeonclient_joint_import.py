#!/usr/bin/env python3
"""Import pip rgpot and a bundled pyeonclient wheel in one process.

pyeonclient vendors librgpot as libeon_rgpot.so.N. N is the SONAME
rgpot was linked with. rgpot 3.7.0 is libeon_rgpot.so.6. The rgpot
wheel keeps the public librgpot SONAME. Both imports must succeed,
and neither package may satisfy the other's NEEDED entry.

With rgpot 3.7.0 _core.abi3.so does not NEED librgpot itself;
libeonclib.so does, and the repair step rewrites that entry and leaves
_core alone. So the check walks the NEEDED closure of the pyeonclient
extension modules through the libraries the wheel ships, including
.pyeonclient.mesonpy.libs, and every shared object in the wheel is
checked for a public librgpot entry.
"""

from __future__ import annotations

import importlib.metadata
import pathlib
import re
import subprocess
import sys

PRIVATE_RGPOT = re.compile(r"^libeon_rgpot\.so\.[0-9]+$")
PUBLIC_RGPOT = re.compile(r"^librgpot\.so(\.|$)")


def is_shared_object(path: pathlib.Path) -> bool:
    return path.is_file() and re.search(r"\.so(\.[0-9]+)*$", path.name) is not None


def needed(path: pathlib.Path) -> list[str]:
    out = subprocess.check_output(["readelf", "-d", str(path)], text=True)
    return re.findall(r"\(NEEDED\)\s+Shared library: \[([^\]]+)\]", out)


def wheel_shared_objects(dist: str) -> list[pathlib.Path]:
    files = importlib.metadata.files(dist) or []
    found = {pathlib.Path(f.locate()).absolute() for f in files}
    return sorted(p for p in found if is_shared_object(p))


def closure(
    roots: list[pathlib.Path], shipped: dict[str, pathlib.Path]
) -> dict[pathlib.Path, list[str]]:
    """NEEDED entries of roots and of every shipped library they reach."""
    seen: dict[pathlib.Path, list[str]] = {}
    todo = list(roots)
    while todo:
        so = todo.pop()
        if so in seen:
            continue
        seen[so] = needed(so)
        todo.extend(shipped[lib] for lib in seen[so] if lib in shipped)
    return seen


def mapped_libraries() -> list[str]:
    maps = pathlib.Path("/proc/self/maps")
    if not maps.exists():
        return []
    return sorted(
        {line.split()[-1] for line in maps.read_text().splitlines() if "/" in line}
    )


def main() -> int:
    import pyeonclient
    import rgpot

    print("pyeonclient", pyeonclient.__version__, pyeonclient.__file__)
    print("rgpot", rgpot.__file__)

    py_sos = wheel_shared_objects("pyeonclient")
    rg_sos = wheel_shared_objects("rgpot")
    extensions = [so for so in py_sos if so.name.startswith("_core")]
    if not extensions:
        print("pyeonclient installed no _core extension", file=sys.stderr)
        return 1
    if not rg_sos:
        print("rgpot installed no shared objects", file=sys.stderr)
        return 1

    shipped = {so.name: so for so in py_sos}
    reached = closure(extensions, shipped)
    for so, libs in reached.items():
        print(so)
        for lib in libs:
            print("  NEEDED", lib)

    ok = True
    private = sorted(
        {lib for libs in reached.values() for lib in libs if PRIVATE_RGPOT.match(lib)}
    )
    if not private:
        print(
            "pyeonclient did not bundle libeon_rgpot: no libeon_rgpot.so.N in "
            "the NEEDED closure of " + ", ".join(so.name for so in extensions),
            file=sys.stderr,
        )
        ok = False
    for lib in private:
        if lib not in shipped:
            print(f"pyeonclient NEEDs {lib} but does not ship it", file=sys.stderr)
            ok = False
    for so in py_sos:
        for lib in reached.get(so) or needed(so):
            if PUBLIC_RGPOT.match(lib):
                print(
                    f"{so.name} still NEEDs the pip rgpot SONAME {lib}",
                    file=sys.stderr,
                )
                ok = False
    if not any(PUBLIC_RGPOT.match(lib) for so in rg_sos for lib in needed(so)):
        print("rgpot wheel lost the librgpot SONAME", file=sys.stderr)
        ok = False

    maps = mapped_libraries()
    if maps and private:
        libs_dirs = {shipped[lib].parent for lib in private if lib in shipped}
        if not any(pathlib.Path(m).parent in libs_dirs for m in maps):
            print("the bundled libeon_rgpot is not loaded", file=sys.stderr)
            ok = False
    if not ok:
        return 1
    print("JOINT_IMPORT_OK", pyeonclient.__version__)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
