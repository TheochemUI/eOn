#!/usr/bin/env python3
"""Copy the MSVC import library rustc leaves beside a cdylib.

Meson's readcon-db target publishes readcon_db.dll and the static
readcon_db.lib. link.exe cannot take the DLL (LNK1181 / LNK1107). The
import library is readcon_db.dll.lib in the cargo target directory.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) not in (3, 4):
        print(
            "usage: copy_msvc_implib.py DLL OUT_IMPLIB [OUT_STAMP]",
            file=sys.stderr,
        )
        return 2
    dll = Path(argv[1]).resolve()
    out = Path(argv[2])
    root = dll.parent
    cands = [p for p in root.rglob("*.dll.lib") if p.name.lower() == "readcon_db.dll.lib"]
    if not cands:
        print(f"no readcon_db.dll.lib under {root}", file=sys.stderr)
        for p in sorted(root.rglob("*")):
            if p.suffix.lower() in {".lib", ".dll", ".a"} or p.name.endswith(".dll.lib"):
                print(f"  {p}", file=sys.stderr)
        return 1
    cands.sort(key=lambda p: (0 if "cargo-target" in p.parts else 1, len(p.parts)))
    shutil.copy2(cands[0], out)
    # A python -c stamp in meson is rewritten by cmd.exe on Windows, so the
    # C file starts with a stray slash and cl reports C2059.
    if len(argv) == 4:
        Path(argv[3]).write_text("typedef int readcon_db_order;\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
