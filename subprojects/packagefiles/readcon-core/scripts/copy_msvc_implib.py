#!/usr/bin/env python3
"""Copy the MSVC cdylib import library out of the cargo target directory.

The static archive is a different file. Linking that archive into a DLL
that also imports torch binds a Rust runtime symbol to a torch DLL that
does not export it, and the process exits 0xC0000139 before main.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(
            "usage: copy_msvc_implib.py CARGO_TARGET PROFILE DEST",
            file=sys.stderr,
        )
        return 2
    root = Path(argv[1]) / argv[2]
    dest = Path(argv[3])
    names = ("readcon_core.dll.lib", "libreadcon_core.dll.a")
    for name in names:
        for cand in (root / name, root / "deps" / name):
            if cand.is_file():
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(cand, dest)
                return 0
    listing = sorted(p.name for p in root.iterdir()) if root.is_dir() else []
    print(f"import library not in {root}: {listing}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
