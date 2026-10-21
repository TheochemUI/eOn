"""Copy a Windows DLL's non-system imports beside it.

The loader searches the executable directory before PATH. A same-named
DLL earlier on the search path can lack the entry point the import lib
recorded, and the process exits 0xC0000139. The copy kept is the first
one that exports every imported symbol.
"""

from __future__ import annotations

import os
import shutil
import struct
import sys
from pathlib import Path

_SYSTEM = {
    "kernel32.dll",
    "kernelbase.dll",
    "ntdll.dll",
    "user32.dll",
    "advapi32.dll",
    "ws2_32.dll",
    "bcrypt.dll",
    "crypt32.dll",
    "ole32.dll",
    "oleaut32.dll",
    "shell32.dll",
    "shlwapi.dll",
    "gdi32.dll",
    "msvcrt.dll",
    "ucrtbase.dll",
    "ntdll.dll",
    "bcryptprimitives.dll",
}


def _system(name: str) -> bool:
    low = name.lower()
    if low in _SYSTEM:
        return True
    return low.startswith(
        ("api-ms-win-", "ext-ms-", "vcruntime", "msvcp", "concrt", "ucrtbase")
    )


def _rva_to_off(data: bytes, sections: list[tuple[int, int, int]], rva: int) -> int:
    for va, vsize, raw in sections:
        if va <= rva < va + max(vsize, 1):
            return raw + (rva - va)
    raise ValueError(f"RVA 0x{rva:x} is not in a section")


def _pe_view(path: Path) -> tuple[bytes, list[tuple[int, int, int]], int, int] | None:
    data = path.read_bytes()
    if data[:2] != b"MZ" or len(data) < 0x40:
        return None
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    if pe + 24 > len(data) or data[pe : pe + 4] != b"PE\0\0":
        return None
    nsect = struct.unpack_from("<H", data, pe + 6)[0]
    opt_size = struct.unpack_from("<H", data, pe + 20)[0]
    opt = pe + 24
    magic = struct.unpack_from("<H", data, opt)[0]
    dd = opt + (112 if magic == 0x20B else 96)
    sect_off = opt + opt_size
    sections: list[tuple[int, int, int]] = []
    for i in range(nsect):
        off = sect_off + i * 40
        if off + 24 > len(data):
            break
        vsize, va, rawsize, raw = struct.unpack_from("<IIII", data, off + 8)
        sections.append((va, max(vsize, rawsize), raw))
    return data, sections, dd, magic


def _cstring(data: bytes, off: int) -> str:
    end = data.index(b"\0", off)
    return data[off:end].decode("ascii", "replace")


def imported_dlls(path: Path) -> list[str]:
    return [name for name, _syms in import_entries(path)]


def import_entries(path: Path) -> list[tuple[str, list[str]]]:
    view = _pe_view(path)
    if view is None:
        return []
    data, sections, dd, magic = view
    if dd + 16 > len(data):
        return []
    import_rva = struct.unpack_from("<I", data, dd + 8)[0]
    if import_rva == 0:
        return []
    pe32plus = magic == 0x20B
    thunk_size = 8 if pe32plus else 4
    ordinal_flag = 1 << (63 if pe32plus else 31)
    entries: list[tuple[str, list[str]]] = []
    try:
        desc = _rva_to_off(data, sections, import_rva)
    except ValueError:
        return []
    while desc + 20 <= len(data):
        name_rva = struct.unpack_from("<I", data, desc + 12)[0]
        if name_rva == 0:
            break
        try:
            name = _cstring(data, _rva_to_off(data, sections, name_rva))
        except (ValueError, IndexError):
            break
        oft = struct.unpack_from("<I", data, desc)[0]
        if oft == 0:
            oft = struct.unpack_from("<I", data, desc + 16)[0]
        syms: list[str] = []
        if oft != 0:
            try:
                toff = _rva_to_off(data, sections, oft)
            except ValueError:
                toff = -1
            while toff >= 0 and toff + thunk_size <= len(data):
                if pe32plus:
                    val = struct.unpack_from("<Q", data, toff)[0]
                else:
                    val = struct.unpack_from("<I", data, toff)[0]
                if val == 0:
                    break
                if val & ordinal_flag:
                    syms.append(f"#{val & 0xFFFF}")
                else:
                    try:
                        no = _rva_to_off(data, sections, val & 0x7FFFFFFF)
                        syms.append(_cstring(data, no + 2))
                    except (ValueError, IndexError):
                        break
                toff += thunk_size
        entries.append((name, syms))
        desc += 20
    return entries


def export_names(path: Path) -> set[str]:
    view = _pe_view(path)
    if view is None:
        return set()
    data, sections, dd, _magic = view
    if dd + 8 > len(data):
        return set()
    export_rva, export_size = struct.unpack_from("<II", data, dd)
    if export_rva == 0 or export_size == 0:
        return set()
    try:
        off = _rva_to_off(data, sections, export_rva)
        _num_funcs, num_names, _addr_funcs, addr_names, _addr_ords = struct.unpack_from(
            "<IIIII", data, off + 20
        )
        name_off = _rva_to_off(data, sections, addr_names)
    except (ValueError, struct.error):
        return set()
    names: set[str] = set()
    for i in range(num_names):
        slot = name_off + 4 * i
        if slot + 4 > len(data):
            break
        rva = struct.unpack_from("<I", data, slot)[0]
        try:
            names.add(_cstring(data, _rva_to_off(data, sections, rva)))
        except (ValueError, IndexError):
            break
    return names


def _missing(path: Path, needed: list[str]) -> list[str]:
    if not needed:
        return []
    have = export_names(path)
    return [sym for sym in needed if not sym.startswith("#") and sym not in have]


def pick_source(name: str, needed: list[str], search: list[Path]) -> Path | None:
    """First search hit that exports every named import.

    A same-named DLL earlier on the path can lack the entry point the
    import lib was linked against. Loading that copy exits 0xC0000139.
    """
    found = [d / name for d in search if (d / name).is_file()]
    if not found:
        return None
    for cand in found:
        if not _missing(cand, needed):
            return cand
    return found[0]


def stage(dll: Path, search: list[Path], stamp: Path) -> None:
    dest = dll.parent
    seen: set[str] = set()
    queue = [dll]
    copied: list[str] = []
    while queue:
        current = queue.pop(0)
        key = str(current.resolve()).lower()
        if key in seen or not current.is_file():
            continue
        seen.add(key)
        for name, syms in import_entries(current):
            if _system(name):
                continue
            target = dest / name
            # An import that is already beside the engine still has its own
            # imports. Walking only newly copied files leaves those behind.
            # A copy that lacks an imported symbol is replaced by a later
            # search hit that exports it.
            if target.is_file() and not _missing(target, syms):
                if str(target.resolve()).lower() not in seen:
                    queue.append(target)
                continue
            source = pick_source(name, syms, search)
            if source is None:
                print(
                    f"win_stage_imports: {name} not in search path",
                    file=sys.stderr,
                )
                continue
            lack = _missing(source, syms)
            if lack:
                shown = ", ".join(lack[:8])
                print(
                    f"win_stage_imports: {name} from {source} lacks {shown}",
                    file=sys.stderr,
                )
            elif target.is_file():
                print(
                    f"win_stage_imports: replaced {name} with {source}",
                    file=sys.stderr,
                )
            replaced = target.is_file()
            shutil.copy2(source, target)
            if replaced:
                seen.discard(str(target.resolve()).lower())
            copied.append(name)
            queue.append(target)
    stamp.write_text("\n".join(copied) + "\n", encoding="utf-8")


def _path_dirs() -> list[Path]:
    raw = os.environ.get("PATH", "")
    sep = ";" if ";" in raw else ":"
    return [Path(p) for p in raw.split(sep) if p]


def report(roots: list[Path]) -> int:
    """Print imports whose resolved DLL does not export the symbol."""
    if not roots:
        return 2
    exe_dir = roots[0].parent
    search = [exe_dir] + _path_dirs()
    seen: set[str] = set()
    queue = [p for p in roots if p.is_file()]
    while queue:
        current = queue.pop(0)
        key = str(current.resolve()).lower()
        if key in seen:
            continue
        seen.add(key)
        for name, syms in import_entries(current):
            if _system(name):
                continue
            found = next((d / name for d in search if (d / name).is_file()), None)
            if found is None:
                print(
                    f"win_stage_imports: MISSING {name} (not on PATH) needed by {current.name}",
                    file=sys.stderr,
                )
                continue
            lack = _missing(found, syms)
            if lack:
                print(
                    f"win_stage_imports: MISSING {name} from {found} "
                    f"needed by {current.name}: {', '.join(lack[:12])}",
                    file=sys.stderr,
                )
            fkey = str(found.resolve()).lower()
            if fkey not in seen:
                queue.append(found)
    return 0


def main(argv: list[str]) -> int:
    if len(argv) >= 2 and argv[1] == "--report":
        return report([Path(p) for p in argv[2:]])
    if len(argv) < 4:
        print(
            "usage: win_stage_imports.py DLL STAMP SEARCH_DIR...",
            file=sys.stderr,
        )
        return 2
    stage(Path(argv[1]), [Path(p) for p in argv[3:]], Path(argv[2]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
