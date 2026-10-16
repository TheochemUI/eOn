"""The MSVC import library keeps the rkrdb_ ABI and drops readcon-core."""

import importlib.util
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "subprojects"
    / "packagefiles"
    / "readcon-db"
    / "scripts"
    / "meson_cargo_build.py"
)

DUMPBIN = """
    ordinal hint RVA      name

          1    0 00001000 create_writer_from_path_c
          2    1 00001100 rkr_frame_new
          3    2 00001200 free_rkr_frame
          4    3 00001300 rkrdb_open
          5    4 00001400 rkrdb_append_trajectory_str
          6    5 00001500 _rkrdb_close
          7    6 00001600 con_frame_iterator_next

  Summary
"""


def _mod():
    spec = importlib.util.spec_from_file_location("meson_cargo_build", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_import_library_drops_readcon_core_symbols():
    mod = _mod()
    names = mod.export_names(DUMPBIN)
    assert "create_writer_from_path_c" in names
    kept = mod.abi_exports(names)
    assert kept == ["rkrdb_open", "rkrdb_append_trajectory_str", "_rkrdb_close"]
    assert all(mod.is_readcon_db_export(name) for name in kept)


def test_def_names_the_readcon_db_dll():
    mod = _mod()
    text = mod.abi_exports(mod.export_names(DUMPBIN))
    # write_implib writes this body; lib.exe is the MSVC half.
    body = f"LIBRARY readcon_db.dll\nEXPORTS\n" + "\n".join(text) + "\n"
    assert body.startswith("LIBRARY readcon_db.dll\nEXPORTS\n")
    assert "create_writer_from_path_c" not in body
    assert "rkrdb_open" in body


def test_built_shared_library_exports_both_abis():
    """The Linux cdylib exports the same two ABIs the Windows DLL does."""
    library = ROOT / "bbdir" / "subprojects" / "readcon-db" / "libreadcon_db.so"
    if not library.is_file():
        pytest.skip("libreadcon_db.so is not built")
    text = subprocess.check_output(
        ["nm", "-D", "--defined-only", str(library)],
        text=True,
    )
    names = [line.split()[-1] for line in text.splitlines() if line.strip()]
    mod = _mod()
    kept = mod.abi_exports(names)
    assert "create_writer_from_path_c" in names
    assert "create_writer_from_path_c" not in kept
    assert "rkrdb_open" in kept
    assert "rkrdb_append_trajectory_str" in kept
