#!/usr/bin/env python3
"""Build an MSVC import library from a MinGW xtb DLL.

conda-forge win-64 xtb ships libxtb-6.dll plus libxtb.dll.a. MSVC
link.exe cannot open the GNU import library (it asks for xtb.lib).
dumpbin /exports plus lib /DEF writes that import library.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


def _tool(name: str) -> str:
    found = shutil.which(name)
    if found is None:
        sys.exit(f"{name} not on PATH; activate the MSVC developer environment")
    return found


def export_names(text: str) -> list[str]:
    names: list[str] = []
    in_table = False
    for line in text.splitlines():
        lower = line.lower()
        if "ordinal" in lower and "hint" in lower and "name" in lower:
            in_table = True
            continue
        if not in_table:
            continue
        stripped = line.strip()
        if not stripped or stripped.lower().startswith("summary"):
            if stripped.lower().startswith("summary"):
                break
            continue
        parts = stripped.split()
        if len(parts) < 4 or not parts[0].isdigit() or not parts[1].isdigit():
            continue
        name = parts[3]
        if name.lower() in {"name", "rva"}:
            continue
        names.append(name)
    return names


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} DLL OUT.lib MACHINE")
    dll = Path(sys.argv[1])
    out_lib = Path(sys.argv[2])
    machine = sys.argv[3]
    if not dll.is_file():
        sys.exit(f"xtb DLL not found: {dll}")
    dumpbin = _tool("dumpbin.exe")
    lib_exe = _tool("lib.exe")
    text = subprocess.check_output(
        [dumpbin, "/nologo", "/exports", str(dll)],
        text=True,
        errors="replace",
    )
    names = export_names(text)
    if not names:
        sys.exit(f"no exports in {dll}\n{text[:2000]}")
    out_lib.parent.mkdir(parents=True, exist_ok=True)
    out_def = out_lib.with_suffix(".def")
    out_def.write_text("EXPORTS\n" + "\n".join(names) + "\n", encoding="ascii")
    subprocess.check_call(
        [
            lib_exe,
            "/nologo",
            f"/machine:{machine}",
            f"/def:{out_def}",
            f"/out:{out_lib}",
        ]
    )
    if not out_lib.is_file():
        sys.exit(f"lib.exe did not write {out_lib}")


if __name__ == "__main__":
    main()
