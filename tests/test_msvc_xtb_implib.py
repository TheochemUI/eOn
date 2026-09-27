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


def test_hex_hint_does_not_drop_the_c_api_name():
    mod = _implib()
    text = """
    ordinal hint RVA      name

       2870  0B36 0005A580 xtb_newEnvironment
"""
    assert mod.export_names(text) == ["xtb_newEnvironment"]


def test_export_table_uses_the_xtb_token_not_the_rva_column():
    mod = _implib()
    text = """
    ordinal hint RVA      name

          1    0 00001000 00001000 xtb_newenvironment_
"""
    assert mod.export_names(text) == ["xtb_newenvironment_"]


def test_def_names_the_conda_dll_not_the_lib():
    mod = _implib()
    text = mod.def_text("libxtb-6.dll", ["xtb_newEnvironment"])
    assert text.startswith("LIBRARY libxtb-6.dll\nEXPORTS\n")


def test_existing_header_spelling_is_not_aliased():
    mod = _implib()
    lines = mod.alias_lines(["xtb_newEnvironment"])
    assert lines == ["xtb_newEnvironment"]
