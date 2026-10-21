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
    struct.pack_into("<IIII", sec, 8, 0x400, 0x1000, 0x400, 0x200)
    pe += sec
    blob = bytearray(dos) + pe
    raw = 0x200
    blob += b"\0" * (raw - len(blob))
    # Two descriptors so the name cannot be read as the next one.
    desc = bytearray(40)
    struct.pack_into("<I", desc, 12, 0x1000 + 40)
    blob += desc + import_name + b"\0"
    return bytes(blob)


def test_import_table_lists_the_dll_and_skips_kernel32(tmp_path: Path):
    mod = _mod()
    pe = tmp_path / "engine.dll"
    pe.write_bytes(_pe(b"capnp.dll"))
    assert mod.imported_dlls(pe) == ["capnp.dll"]
    assert mod._system("KERNEL32.dll")
    assert mod._system("bcryptprimitives.dll")
    assert not mod._system("capnp.dll")


def test_torch_uv_is_taken_from_the_torch_directory(tmp_path: Path):
    mod = _mod()
    engine = tmp_path / "engine.dll"
    engine.write_bytes(_pe(b"torch_cpu.dll"))
    torch_dir = tmp_path / "torch"
    prefix = tmp_path / "prefix"
    torch_dir.mkdir()
    prefix.mkdir()
    (torch_dir / "torch_cpu.dll").write_bytes(_pe(b"uv.dll"))
    (torch_dir / "uv.dll").write_bytes(b"torch-uv")
    (prefix / "uv.dll").write_bytes(b"prefix-uv")
    mod.stage(engine, [torch_dir, prefix], tmp_path / "stamp.txt")
    assert (tmp_path / "torch_cpu.dll").is_file()
    assert (tmp_path / "uv.dll").read_bytes() == b"torch-uv"


def test_an_import_already_beside_the_engine_is_still_walked(tmp_path: Path):
    mod = _mod()
    engine = tmp_path / "engine.dll"
    engine.write_bytes(_pe(b"torch_cpu.dll"))
    (tmp_path / "torch_cpu.dll").write_bytes(_pe(b"uv.dll"))
    torch_dir = tmp_path / "torch"
    torch_dir.mkdir()
    (torch_dir / "uv.dll").write_bytes(b"torch-uv")
    mod.stage(engine, [torch_dir], tmp_path / "stamp.txt")
    assert (tmp_path / "uv.dll").read_bytes() == b"torch-uv"
