"""Structure (server) ↔ Matter (client) without requiring ASE.

ASE converters live in :mod:`pyeonclient.ase_bridge` (optional dependency).
"""

from __future__ import annotations

from typing import Any

import numpy as np

# IUPAC symbols Z=1..118. pyeonclient ships as its own distribution and must
# not import eon.atoms. Prefer readcon.symbol_to_z / readcon.z_to_symbol when
# the pymodule exports them (same table the C++ client uses via
# readcon::symbol_to_z / readcon::z_to_symbol); otherwise use this table.
_SYMBOLS: tuple[str, ...] = (
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)

_SYMBOL_Z: dict[str, int] = {sym: z for z, sym in enumerate(_SYMBOLS, start=1)}
_Z_SYMBOL: dict[int, str] = {z: sym for sym, z in _SYMBOL_Z.items()}


def _canonical_symbol(sym: str) -> str:
    s = (sym or "").strip()
    if not s:
        return s
    return s[:1].upper() + s[1:].lower()


def _readcon_symbol_to_z(sym: str) -> int | None:
    try:
        import readcon
    except Exception:
        return None
    fn = getattr(readcon, "symbol_to_z", None) or getattr(
        readcon, "symbol_to_atomic_number", None
    )
    if fn is None:
        return None
    try:
        z = int(fn(sym))
    except Exception:
        return None
    return z if z > 0 else None


def _readcon_z_to_symbol(z: int) -> str | None:
    try:
        import readcon
    except Exception:
        return None
    fn = getattr(readcon, "z_to_symbol", None) or getattr(
        readcon, "atomic_number_to_symbol", None
    )
    if fn is None:
        return None
    try:
        symbol = fn(int(z))
    except Exception:
        return None
    if not symbol or symbol in {"X", "Xx"}:
        return None
    return str(symbol)


def _symbol_to_z(sym: str) -> int:
    s = (sym or "").strip()
    if not s:
        raise KeyError("empty chemical symbol")
    z = _readcon_symbol_to_z(s)
    if z is not None:
        return z
    t = _canonical_symbol(s)
    if t != s:
        z = _readcon_symbol_to_z(t)
        if z is not None:
            return z
    if s in _SYMBOL_Z:
        return _SYMBOL_Z[s]
    if t in _SYMBOL_Z:
        return _SYMBOL_Z[t]
    raise KeyError(f"no atomic number for symbol {sym!r}")


def _z_to_symbol(z: int) -> str:
    zi = int(z)
    symbol = _readcon_z_to_symbol(zi)
    if symbol is not None:
        return symbol
    if zi in _Z_SYMBOL:
        return _Z_SYMBOL[zi]
    raise KeyError(f"no chemical symbol for Z={zi}")


def structure_to_matter(
    structure: Any,
    potential: Any,
    parameters: Any,
    periodic: bool | None = None,
) -> Any:
    """Build a live :class:`Matter` from a server Structure/Atoms-like object.

    *periodic* defaults to the ``periodic`` attribute of *structure* when it
    has one, and to True otherwise (the .con convention). Pass False for a
    cluster, so later ``setPositions`` calls do not wrap through
    ``applyPeriodicBoundary``.
    """
    from pyeonclient import Matter

    if periodic is None:
        periodic = bool(getattr(structure, "periodic", True))
    n = len(structure)
    m = Matter(potential, parameters)
    m.resize(n)
    # Cell before positions: default Matter cell is Zero; with PBC on, setting
    # positions first would wrap through a singular cell.
    m.cell = np.ascontiguousarray(structure.box, dtype=np.float64)
    # Periodicity before positions, for the same reason.
    m.periodic = periodic
    m.positions = np.ascontiguousarray(structure.r, dtype=np.float64)
    m.masses = np.ascontiguousarray(structure.mass, dtype=np.float64)
    free = np.asarray(structure.free, dtype=np.float64).reshape(-1)
    m.fixed = (free < 0.5).astype(np.int64)
    z = np.empty(n, dtype=np.int64)
    for i, name in enumerate(structure.names):
        z[i] = _symbol_to_z(name)
    m.atomic_numbers = z
    return m


def matter_to_structure(matter: Any) -> Any:
    """Convert a live :class:`Matter` to :class:`eon.structure.Structure`."""
    from eon.structure import Structure

    n = int(matter.n_atoms)
    s = Structure(n)
    s.r = np.asarray(matter.positions, dtype=float).copy()
    s.box = np.asarray(matter.cell, dtype=float).copy()
    s.mass = np.asarray(matter.masses, dtype=float).copy()
    fixed = np.asarray(matter.fixed, dtype=np.int64).reshape(-1)
    s.free = np.where(fixed != 0, 0.0, 1.0)
    z = np.asarray(matter.atomic_numbers, dtype=np.int64).reshape(-1)
    s.names = [_z_to_symbol(int(z[i])) for i in range(n)]
    return s


# PEP-style aliases used by tools / notebooks
to_structure = matter_to_structure
from_structure = structure_to_matter
