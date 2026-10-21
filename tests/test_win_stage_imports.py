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


def test_pick_source_skips_a_dll_that_lacks_the_symbol(tmp_path: Path):
    mod = _mod()
    bad_dir = tmp_path / "pixi"
    good_dir = tmp_path / "torch"
    bad_dir.mkdir()
    good_dir.mkdir()
    bad = bad_dir / "c10.dll"
    good = good_dir / "c10.dll"
    bad.write_bytes(b"MZ-incomplete")
    good.write_bytes(b"MZ-complete")

    def fake_missing(path: Path, needed: list[str]) -> list[str]:
        if path == bad:
            return ["MissingSym"]
        return []

    mod._missing = fake_missing
    assert mod.pick_source("c10.dll", ["MissingSym"], [bad_dir, good_dir]) == good
    assert mod.pick_source("c10.dll", ["MissingSym"], [bad_dir]) == bad
    # The importing DLL's directory is first. It wins when it exports the
    # symbol, and a later search hit is kept when it does not.
    assert mod.pick_source("c10.dll", ["MissingSym"], [bad_dir], good_dir) == good
    assert mod.pick_source("c10.dll", ["MissingSym"], [good_dir], bad_dir) == good
