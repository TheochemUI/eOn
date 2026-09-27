"""MSVC xtb import library aliases."""

import importlib.util
from pathlib import Path


def _implib():
    path = Path(__file__).resolve().parents[1] / "scripts" / "ci" / "msvc_xtb_implib.py"
    spec = importlib.util.spec_from_file_location("msvc_xtb_implib", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_lowercase_export_aliases_the_header_spelling():
    mod = _implib()
    lines = mod.alias_lines(["xtb_newenvironment_", "xtb_getgradient_"])
    assert "xtb_newEnvironment=xtb_newenvironment_" in lines
    assert "xtb_getGradient=xtb_getgradient_" in lines
    assert "xtb_newenvironment=xtb_newenvironment_" in lines


def test_existing_header_spelling_is_not_aliased():
    mod = _implib()
    lines = mod.alias_lines(["xtb_newEnvironment"])
    assert lines == ["xtb_newEnvironment"]
