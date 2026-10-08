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
    core = meson.split("readcon_dep = dependency(", 1)[1].split("if not readcon_dep.found()", 1)[0]
    db = meson.split("readcon_db_dep = dependency(", 1)[1].split(
        "if not readcon_db_dep.found()", 1
    )[0]
    assert "'readcon-core'" in core
    assert "'readcon-db'" in db
    assert "subproject(\n        'readcon-db'" in meson or "subproject(\n    'readcon-db'" in meson
    assert f"'readcon-db',\n    version: '>={pin}'" in meson
    assert "'-lreadcon_db'" not in meson
    fallback = meson.split("if not readcon_db_dep.found()", 1)[1].split("else", 1)[0]
    assert "link_with: _readcon_db_order" in fallback
    assert "_readcon_db_order = [_readcon_db_so]" in fallback


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
