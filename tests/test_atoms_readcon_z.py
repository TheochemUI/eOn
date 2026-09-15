"""Z/symbol lookups prefer readcon and fall back to the local table."""

import pytest

pytest.importorskip("readcon")
from eon import atoms


def test_atomic_number_hydrogen():
    assert atoms.atomic_number("H") == 1
    assert atoms.atomic_number(1) == 1


def test_symbol_for_z_carbon():
    assert atoms.symbol_for_z(6) == "C"


def test_table_still_has_radius_and_color():
    assert "radius" in atoms.elements["Si"]
    assert "color" in atoms.elements["Si"]
    assert atoms.elements[14]["number"] == 14
    assert atoms.elements[14]["symbol"] == "Si"
    assert atoms.elements[0]["symbol"] == "Xx"
