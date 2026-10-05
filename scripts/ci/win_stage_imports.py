"""Copy a Windows DLL's non-system imports beside it.

The loader searches the executable directory before PATH. A same-named
DLL from another prefix is missing the entry point the engine linked,
and the process exits 0xC0000139. The directory that supplied a DLL is
searched first for that DLL's own imports, then the caller search list.
An existing copy is kept and still scanned, so a DLL placed by an
earlier target is not replaced and its imports are not skipped.
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
    "bcryptprimitives.dll",
    "crypt32.dll",
    "ole32.dll",
    "oleaut32.dll",
    "shell32.dll",
    "shlwapi.dll",
    "gdi32.dll",
    "msvcrt.dll",
    "ucrtbase.dll",
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


def _pe_sections(data: bytes) -> tuple[int, list[tuple[int, int, int]], int] | None:
    if len(data) < 0x40 or data[:2] != b"MZ":
        return None
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    if pe + 24 > len(data) or data[pe : pe + 4] != b"PE\0\0":
        return None
    nsect = struct.unpack_from("<H", data, pe + 6)[0]
    opt_size = struct.unpack_from("<H", data, pe + 20)[0]
    opt = pe + 24
    if opt + opt_size > len(data):
        return None
    magic = struct.unpack_from("<H", data, opt)[0]
    dd = opt + (112 if magic == 0x20B else 96)
    sect_off = opt + opt_size
    sections: list[tuple[int, int, int]] = []
    for i in range(nsect):
        off = sect_off + i * 40
        if off + 24 > len(data):
            return None
        vsize, va, rawsize, raw = struct.unpack_from("<IIII", data, off + 8)
        sections.append((va, max(vsize, rawsize), raw))
    return dd, sections, magic


def _cstring(data: bytes, off: int) -> str:
    end = data.index(b"\0", off)
    return data[off:end].decode("ascii")


def imported_dlls(path: Path) -> list[str]:
    return [name for name, _syms in imported_symbols(path)]


def imported_symbols(path: Path) -> list[tuple[str, list[str]]]:
    data = path.read_bytes()
    parsed = _pe_sections(data)
    if parsed is None:
        return []
    dd, sections, magic = parsed
    if dd + 16 > len(data):
        return []
    import_rva = struct.unpack_from("<I", data, dd + 8)[0]
    if import_rva == 0:
        return []
    names: list[tuple[str, list[str]]] = []
    desc = _rva_to_off(data, sections, import_rva)
    wide = magic == 0x20B
    while desc + 20 <= len(data):
        ilt_rva = struct.unpack_from("<I", data, desc)[0]
        name_rva = struct.unpack_from("<I", data, desc + 12)[0]
        iat_rva = struct.unpack_from("<I", data, desc + 16)[0]
        if name_rva == 0:
            break
        name_off = _rva_to_off(data, sections, name_rva)
        dll = _cstring(data, name_off)
        thunk_rva = ilt_rva or iat_rva
        symbols: list[str] = []
        if thunk_rva:
            thunk = _rva_to_off(data, sections, thunk_rva)
            step = 8 if wide else 4
            top = 1 << (63 if wide else 31)
            while thunk + step <= len(data):
                word = struct.unpack_from("<Q" if wide else "<I", data, thunk)[0]
                if word == 0:
                    break
                if word & top:
                    symbols.append(f"#{word & 0xFFFF}")
                else:
                    hint_off = _rva_to_off(data, sections, word & 0x7FFFFFFF)
                    symbols.append(_cstring(data, hint_off + 2))
                thunk += step
        names.append((dll, symbols))
        desc += 20
    return names


def _find(name: str, dirs: list[Path]) -> Path | None:
    for directory in dirs:
        candidate = directory / name
        if candidate.is_file():
            return candidate
    return None


def stage(dll: Path, search: list[Path], stamp: Path) -> None:
    dest = dll.parent
    seen: set[str] = set()
    queue: list[tuple[Path, list[Path]]] = [(dll, [])]
    copied: list[str] = []
    while queue:
        current, siblings = queue.pop(0)
        try:
            key = str(current.resolve()).lower()
        except OSError:
            continue
        if key in seen or not current.is_file():
            continue
        seen.add(key)
        dirs = siblings + search
        try:
            imports = imported_symbols(current)
        except (ValueError, UnicodeDecodeError, IndexError):
            print(f"win_stage_imports: {current.name} has no import table", file=sys.stderr)
            continue
        for name, _symbols in imports:
            if _system(name):
                continue
            target = dest / name
            if target.is_file():
                queue.append((target, [target.parent]))
                continue
            source = _find(name, dirs)
            if source is None:
                print(f"win_stage_imports: {name} not in search path", file=sys.stderr)
                continue
            shutil.copy2(source, target)
            copied.append(f"{name} <- {source}")
            queue.append((target, [source.parent]))
    stamp.write_text("\n".join(copied) + "\n", encoding="utf-8")


def exported_names(path: Path) -> set[str]:
    data = path.read_bytes()
    parsed = _pe_sections(data)
    if parsed is None:
        return set()
    dd, sections, _magic = parsed
    if dd + 8 > len(data):
        return set()
    export_rva, export_size = struct.unpack_from("<II", data, dd)
    if export_rva == 0 or export_size == 0:
        return set()
    exp = _rva_to_off(data, sections, export_rva)
    if exp + 40 > len(data):
        return set()
    nnames = struct.unpack_from("<I", data, exp + 24)[0]
    names_rva = struct.unpack_from("<I", data, exp + 32)[0]
    if nnames == 0 or names_rva == 0:
        return set()
    table = _rva_to_off(data, sections, names_rva)
    found: set[str] = set()
    for i in range(nnames):
        slot = table + i * 4
        if slot + 4 > len(data):
            break
        name_rva = struct.unpack_from("<I", data, slot)[0]
        found.add(_cstring(data, _rva_to_off(data, sections, name_rva)))
    return found


def _system_dirs() -> list[Path]:
    root = os.environ.get("SystemRoot", r"C:\Windows")
    return [Path(root) / "System32", Path(root) / "SysWOW64", Path(root)]


def audit_missing(root: Path, search: list[Path]) -> list[str]:
    """Return lines for imports whose DLL or named symbol is absent."""
    problems: list[str] = []
    seen: set[str] = set()
    queue = [root]
    system_dirs = _system_dirs()
    while queue:
        current = queue.pop(0)
        try:
            key = str(current.resolve()).lower()
        except OSError:
            continue
        if key in seen or not current.is_file():
            continue
        seen.add(key)
        try:
            imports = imported_symbols(current)
        except (ValueError, UnicodeDecodeError, IndexError):
            problems.append(f"{current.name}: import table is unreadable")
            continue
        for name, symbols in imports:
            if _system(name):
                continue
            source = _find(name, [current.parent] + search)
            if source is None:
                if _find(name, system_dirs) is not None:
                    continue
                problems.append(f"{current.name}: {name} not in search path")
                continue
            if symbols:
                try:
                    have = exported_names(source)
                except (ValueError, UnicodeDecodeError, IndexError):
                    problems.append(f"{name}: export table is unreadable ({source})")
                    continue
                if have:
                    missing = [
                        symbol
                        for symbol in symbols
                        if not symbol.startswith("#") and symbol not in have
                    ]
                    if missing:
                        shown = ", ".join(missing[:8])
                        problems.append(
                            f"{current.name}: {name} from {source} missing {shown}"
                        )
            queue.append(source)
    return problems


def main(argv: list[str]) -> int:
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
