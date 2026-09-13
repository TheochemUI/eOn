#!/usr/bin/env python3
"""Copy wrap-built librgpot into the install prefix.

rgpot sets install: not meson.is_subproject(). eOn CI uses
meson install --skip-subprojects. eonclient then looks for
librgpot.3.dylib (soversion 3 in rgpot 3.2.0) next to the prefix.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def main() -> int:
    src = Path(sys.argv[1])
    libdir = Path(os.environ["MESON_INSTALL_DESTDIR_PREFIX"]) / sys.argv[2]
    libdir.mkdir(parents=True, exist_ok=True)
    parent = src.parent
    copied = 0
    for path in parent.iterdir():
        name = path.name
        if not (
            name.startswith("librgpot")
            or name.startswith("rgpot.")
            or name == "rgpot.dll"
            or name == "rgpot.lib"
        ):
            continue
        dest = libdir / name
        if path.is_symlink():
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            dest.symlink_to(os.readlink(path))
        else:
            shutil.copy2(path, dest)
        copied += 1
    if copied == 0:
        print(f"install_wrap_rgpot: no librgpot next to {src}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
