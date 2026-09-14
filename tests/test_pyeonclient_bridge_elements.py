"""pyeonclient.bridge symbol/Z helpers: full table, no server import, no Z79."""

from __future__ import annotations

import importlib.util
import inspect
from pathlib import Path

import pytest

_BRIDGE = (
    Path(__file__).resolve().parents[1]
    / "client"
    / "python"
    / "pyeonclient"
    / "bridge.py"
)


def _load_bridge():
    spec = importlib.util.spec_from_file_location(
        "pyeonclient_bridge_helpers", _BRIDGE
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bridge = _load_bridge()
_symbol_to_z = bridge._symbol_to_z
_z_to_symbol = bridge._z_to_symbol


def test_bridge_does_not_import_eon_atoms():
    src = inspect.getsource(bridge)
    assert "from eon.atoms import" not in src
    assert "from eon import atoms" not in src


@pytest.mark.parametrize(
    "symbol,z",
    [
        ("H", 1),
        ("Ag", 47),
        ("Ni", 28),
        ("Pd", 46),
        ("Ir", 77),
        ("Mo", 42),
        ("W", 74),
        ("Pt", 78),
        ("Au", 79),
        ("Og", 118),
    ],
)
def test_symbol_z_roundtrip_catalysis_and_edges(symbol: str, z: int):
    assert _symbol_to_z(symbol) == z
    assert _z_to_symbol(z) == symbol


def test_symbol_to_z_accepts_odd_case():
    assert _symbol_to_z("ag") == 47
    assert _symbol_to_z("AG") == 47
    assert _symbol_to_z(" w ") == 74


def test_unknown_symbol_raises():
    with pytest.raises(KeyError, match="no atomic number"):
        _symbol_to_z("Xy")
    with pytest.raises(KeyError, match="empty chemical symbol"):
        _symbol_to_z("")
    with pytest.raises(KeyError, match="empty chemical symbol"):
        _symbol_to_z("   ")


def test_unknown_z_raises_instead_of_fabricating():
    with pytest.raises(KeyError, match="no chemical symbol for Z=0"):
        _z_to_symbol(0)
    with pytest.raises(KeyError, match="no chemical symbol for Z=999"):
        _z_to_symbol(999)
    with pytest.raises(KeyError, match="no chemical symbol"):
        _z_to_symbol(-1)
