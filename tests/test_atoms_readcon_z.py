"""Z/symbol lookups go through readcon's Python helpers (0.14.9+)."""

import pytest

readcon = pytest.importorskip("readcon")
from eon import atoms


def test_readcon_exports_symbol_z_helpers():
    assert hasattr(readcon, "symbol_to_atomic_number")
    assert hasattr(readcon, "atomic_number_to_symbol")
    assert readcon.symbol_to_atomic_number("H") == 1
    assert readcon.atomic_number_to_symbol(6) == "C"


def test_atomic_number_hydrogen():
    assert atoms.atomic_number("H") == 1
    assert atoms.atomic_number(1) == 1


def test_symbol_for_z_carbon():
    assert atoms.symbol_for_z(6) == "C"


def test_no_ctypes_hole():
    assert not hasattr(atoms, "_rkr_lib")


def test_table_still_has_radius_and_color():
    assert "radius" in atoms.elements["Si"]
    assert "color" in atoms.elements["Si"]
    assert atoms.elements[14]["number"] == 14
    assert atoms.elements[14]["symbol"] == "Si"
    assert atoms.elements[0]["symbol"] == "Xx"
