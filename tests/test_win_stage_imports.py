"""The Windows import stager reads a PE import table and skips system DLLs."""

import importlib.util
import struct
from pathlib import Path


def _mod():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "ci"
        / "win_stage_imports.py"
    )
    spec = importlib.util.spec_from_file_location("win_stage_imports", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _pe(import_name: bytes) -> bytes:
    # PE32+, one section, one import descriptor.
    dos = bytearray(128)
    dos[:2] = b"MZ"
    struct.pack_into("<I", dos, 0x3C, 128)
    pe = bytearray()
    pe += b"PE\0\0"
    # COFF: machine, 1 section, timestamp, sym, nsym, opt size, characteristics
    pe += struct.pack("<HHIIIHH", 0x8664, 1, 0, 0, 0, 240, 0x22)
    opt = bytearray(240)
    struct.pack_into("<H", opt, 0, 0x20B)
    # Import directory is data directory index 1, at optional+112.
    struct.pack_into("<I", opt, 112 + 8, 0x1000)
    struct.pack_into("<I", opt, 112 + 12, 40)
    pe += opt
    # Section .rdata va 0x1000, raw at 0x200
    sec = bytearray(40)
    sec[:6] = b".rdata"
    struct.pack_into("<IIII", sec, 8, 0x200, 0x1000, 0x200, 0x200)
    pe += sec
    blob = bytearray(dos) + pe
    raw = 0x200
    blob += b"\0" * (raw - len(blob))
    desc = bytearray(20)
    struct.pack_into("<I", desc, 12, 0x1000 + 20)
    blob += desc + import_name + b"\0" + bytearray(20)
    return bytes(blob)


def test_import_table_lists_the_dll_and_skips_kernel32(tmp_path: Path):
    mod = _mod()
    pe = tmp_path / "engine.dll"
    pe.write_bytes(_pe(b"capnp.dll"))
    assert mod.imported_dlls(pe) == ["capnp.dll"]
    assert mod._system("KERNEL32.dll")
    assert not mod._system("capnp.dll")


def test_libomp_comes_from_the_first_search_directory(tmp_path: Path):
    mod = _mod()
    pe = tmp_path / "engine.dll"
    pe.write_bytes(_pe(b"KERNEL32.dll"))
    first = tmp_path / "prefix"
    second = tmp_path / "plugin"
    first.mkdir()
    second.mkdir()
    (first / "libomp.dll").write_bytes(b"llvm")
    (second / "libomp.dll").write_bytes(b"other")
    mod.stage(pe, [first, second], tmp_path / "stamp.txt")
    assert (tmp_path / "libomp.dll").read_bytes() == b"llvm"
