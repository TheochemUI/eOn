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


# The C header spells these with internal capitals. MinGW bind(C) exports
# are lowercase, often with a trailing underscore. link.exe does not
# case-fold, so the import library has to name the header spelling.
CANONICAL_XTB = {
    "xtb_newenvironment": "xtb_newEnvironment",
    "xtb_delenvironment": "xtb_delEnvironment",
    "xtb_releaseoutput": "xtb_releaseOutput",
    "xtb_setverbosity": "xtb_setVerbosity",
    "xtb_delmolecule": "xtb_delMolecule",
    "xtb_newcalculator": "xtb_newCalculator",
    "xtb_delcalculator": "xtb_delCalculator",
    "xtb_newresults": "xtb_newResults",
    "xtb_delresults": "xtb_delResults",
    "xtb_checkenvironment": "xtb_checkEnvironment",
    "xtb_geterror": "xtb_getError",
    "xtb_newmolecule": "xtb_newMolecule",
    "xtb_updatemolecule": "xtb_updateMolecule",
    "xtb_loadgfn0xtb": "xtb_loadGFN0xTB",
    "xtb_loadgfn1xtb": "xtb_loadGFN1xTB",
    "xtb_loadgfn2xtb": "xtb_loadGFN2xTB",
    "xtb_loadgfnff": "xtb_loadGFNFF",
    "xtb_setaccuracy": "xtb_setAccuracy",
    "xtb_setmaxiter": "xtb_setMaxIter",
    "xtb_setelectronictemp": "xtb_setElectronicTemp",
    "xtb_singlepoint": "xtb_singlepoint",
    "xtb_getenergy": "xtb_getEnergy",
    "xtb_getgradient": "xtb_getGradient",
}


def alias_lines(names: list[str]) -> list[str]:
    """DEF export lines, including aliases for the C header spellings."""
    export_set = set(names)
    lines = list(names)
    aliased: set[str] = set()

    def add(plain: str, name: str) -> None:
        if (
            plain.startswith("xtb_")
            and plain not in export_set
            and plain not in aliased
            and plain != name
        ):
            lines.append(f"{plain}={name}")
            aliased.add(plain)

    for name in names:
        plain = name[1:] if name.startswith("_") else name
        if plain.endswith("_"):
            plain = plain[:-1]
        add(plain, name)
        canon = CANONICAL_XTB.get(plain.lower())
        if canon is not None:
            add(canon, name)
    return lines


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
        if len(parts) < 3 or not parts[0].isdigit() or not parts[1].isdigit():
            continue
        name = ""
        for token in parts:
            bare = token.split("=", 1)[0]
            if bare.lower().startswith("xtb") or bare.lower().startswith("_xtb"):
                name = bare
        if not name:
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
    # MinGW gfortran often exports bind(C) names in lowercase with a
    # trailing underscore. The xtb header calls the camel-case C name.
    for name in names:
        if "environment" in name.lower():
            print(f"xtb export: {name}", file=sys.stderr)
    lines = alias_lines(names)
    exported = {line.split("=", 1)[0] for line in lines}
    if "xtb_newEnvironment" not in exported:
        shown = "\n".join(names[:80])
        sys.exit(
            "xtb import library has no xtb_newEnvironment. "
            f"DLL exports:\n{shown}"
        )
    out_lib.parent.mkdir(parents=True, exist_ok=True)
    out_def = out_lib.with_suffix(".def")
    out_def.write_text("EXPORTS\n" + "\n".join(lines) + "\n", encoding="ascii")
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
