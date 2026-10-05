"""Copy a Windows DLL's non-system imports beside it.

The loader searches the executable directory before PATH. Torch's
directory is on PATH ahead of the pixi prefix, so a same-named DLL
from torch is missing the entry point the engine was linked against
and the process exits 0xC0000139. Pot plugins already keep their
runtime DLLs beside the loader. This does that for the engine DLL.
"""

from __future__ import annotations

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


def imported_dlls(path: Path) -> list[str]:
    data = path.read_bytes()
    if data[:2] != b"MZ":
        return []
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    if data[pe : pe + 4] != b"PE\0\0":
        return []
    nsect = struct.unpack_from("<H", data, pe + 6)[0]
    opt_size = struct.unpack_from("<H", data, pe + 20)[0]
    opt = pe + 24
    magic = struct.unpack_from("<H", data, opt)[0]
    dd = opt + (112 if magic == 0x20B else 96)
    import_rva = struct.unpack_from("<I", data, dd + 8)[0]
    if import_rva == 0:
        return []
    sect_off = opt + opt_size
    sections: list[tuple[int, int, int]] = []
    for i in range(nsect):
        off = sect_off + i * 40
        vsize, va, rawsize, raw = struct.unpack_from("<IIII", data, off + 8)
        sections.append((va, max(vsize, rawsize), raw))
    names: list[str] = []
    desc = _rva_to_off(data, sections, import_rva)
    while True:
        name_rva = struct.unpack_from("<I", data, desc + 12)[0]
        if name_rva == 0:
            break
        name_off = _rva_to_off(data, sections, name_rva)
        end = data.index(b"\0", name_off)
        names.append(data[name_off:end].decode("ascii"))
        desc += 20
    return names


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
        for name in imported_dlls(current):
            if _system(name):
                continue
            target = dest / name
            if target.is_file():
                continue
            source = next(
                (d / name for d in search if (d / name).is_file()),
                None,
            )
            if source is None:
                print(f"win_stage_imports: {name} not in search path", file=sys.stderr)
                continue
            shutil.copy2(source, target)
            copied.append(name)
            queue.append(target)
    stamp.write_text("\n".join(copied) + "\n", encoding="utf-8")


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
