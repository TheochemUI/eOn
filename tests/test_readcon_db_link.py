"""The client links readcon-db the way it links readcon-core."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WRAP = ROOT / "subprojects" / "readcon-db.wrap"
MESON = ROOT / "client" / "meson.build"
HEADER = ROOT / "client" / "ReadconDbMirror.h"
CPP = ROOT / "client" / "ReadconDbMirror.cpp"
CASE = ROOT / "client" / "unit_tests" / "ReadconDbMirrorTest.cpp"
PYPROJECT = ROOT / "pyproject.toml"


def _python_pin() -> str:
    text = PYPROJECT.read_text(encoding="utf-8")
    match = re.search(r"readcon-db>=([0-9.]+)", text)
    assert match is not None
    return match.group(1)


def test_wrap_revision_matches_the_python_package():
    pin = _python_pin()
    text = WRAP.read_text(encoding="utf-8")
    assert f"revision = v{pin}" in text
    assert "https://github.com/lode-org/readcon-db.git" in text
    meson = MESON.read_text(encoding="utf-8")
    # readcon-core is resolved with minimage and linkcell: all three
    # pkg-config deps, or all three wraps. The C ABI is one umbrella
    # staticlib, so the client never links libreadcon_core.a.
    assert "readcon_dep = _rkr_pc" in meson
    assert "readcon_dep = subproject(" in meson
    # The wheel sets wrap_mode=forcefallback. Without this, the
    # headers-only overlays satisfy the probe and the umbrella is skipped.
    mi = meson.split("_mi_pc = dependency(", 1)[1].split("_lc_pc = dependency(", 1)[0]
    lc = meson.split("_lc_pc = dependency(", 1)[1].split("_rkr_pc = dependency(", 1)[0]
    rkr = meson.split("_rkr_pc = dependency(", 1)[1].split(
        "if _mi_pc.found()", 1
    )[0]
    assert "allow_fallback: false" in mi
    assert "allow_fallback: false" in lc
    assert "allow_fallback: false" in rkr
    assert "'readcon-core'" in meson
    assert "version: '>=0.16.0'" in meson
    assert "eon_rust_abi" in meson
    assert "libreadcon_core" not in meson
    db = meson.split("readcon_db_dep = dependency(", 1)[1].split(
        "if not readcon_db_dep.found()", 1
    )[0]
    assert "'readcon-db'" in db
    assert "allow_fallback: false" in db
    assert f"version: '>={pin}'" in db
    assert "subproject(\n        'readcon-db'" in meson
    assert f"'readcon-db',\n        version: '>={pin}'" in meson
    assert "'-lreadcon_db'" not in meson
    # The wrap fallback keeps the cdylib. MSVC links the import library
    # rustc writes beside the DLL, not the DLL itself.
    assert "_readcon_db_link_with = [_readcon_db_so]" in meson
    assert "link_with: _readcon_db_link_with" in meson
    assert "readcon_db.dll.lib" in meson


def test_missing_library_leaves_iostatus_and_fails_the_check():
    header = HEADER.read_text(encoding="utf-8")
    cpp = CPP.read_text(encoding="utf-8")
    case = CASE.read_text(encoding="utf-8")
    assert "bool readcon_db_mirror_ok();" in header
    assert "Failure to load the library" in header
    assert "IoStatus" in header
    assert "EON_READCON_DB_LIBRARY" in cpp
    assert "IoStatus" not in cpp
    assert "rkrdb_open" in cpp
    assert "EON_READCON_DB_VERSION" in cpp
    assert "IoStatus::Ok" in case
    assert "readcon_db_mirror_ok()" in case
    assert "/no/such/libreadcon_db.so" in case
    assert _python_pin() in case
