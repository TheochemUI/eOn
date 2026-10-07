"""The INI loader does not mutate the reader."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_ini_loader_takes_a_const_reader():
    header = (ROOT / "include" / "eon" / "ParametersINI.h").read_text(encoding="utf-8")
    generated = (
        ROOT / "include" / "eon" / "generated" / "ParametersSSOTIni.inc"
    ).read_text(encoding="utf-8")
    assert "load_ini(const INIReader &ini, Parameters &params)" in header
    assert "project_ssot_ini(const INIReader &ini, Parameters &params)" in generated
    assert not (ROOT / "client" / "INIFile.cpp").exists()
