"""Structure (server) ↔ Matter (client) without requiring ASE.

ASE converters live in :mod:`pyeonclient.ase_bridge` (optional dependency).
"""

from __future__ import annotations

from typing import Any

import numpy as np

# Minimal Z table for common eOn fixtures; extended via eon.atoms.elements when
# the server package is installed alongside.
#
# readcon-core carries the full table (helpers.rs symbol_to_atomic_number /
# atomic_number_to_symbol) and exposes it to C++ as rkr_symbol_to_z /
# rkr_z_to_symbol, which is what client/ConFileIO.cpp calls. The readcon
# pymodule registers neither, so a pyeonclient-only install (readcon is a
# dependency, eon is not) is limited to the symbols listed here. Delete this
# table once readcon exports the two helpers.
_SYMBOL_Z = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Ti": 22,
    "Fe": 26,
    "Cu": 29,
    "Zn": 30,
    "Pt": 78,
    "Au": 79,
}


def _symbol_to_z(sym: str) -> int:
    s = (sym or "").strip()
    if not s:
        raise KeyError("empty chemical symbol")
    if s in _SYMBOL_Z:
        return _SYMBOL_Z[s]
    t = s[:1].upper() + s[1:].lower()
    if t in _SYMBOL_Z:
        return _SYMBOL_Z[t]
    try:
        from eon.atoms import elements

        info = elements.get(s) or elements.get(t)
        if info is not None:
            return int(info["number"])
    except Exception:
        pass
    raise KeyError(
        f"no atomic number for symbol {sym!r}; extend pyeonclient.bridge._SYMBOL_Z"
    )


def _z_to_symbol(z: int) -> str:
    try:
        from eon.atoms import elements

        info = elements.get(int(z))
        if info is not None:
            return str(info["symbol"])
    except Exception:
        pass
    for sym, zi in _SYMBOL_Z.items():
        if zi == int(z):
            return sym
    raise KeyError(
        f"no chemical symbol for Z={int(z)}; extend pyeonclient.bridge._SYMBOL_Z"
    )


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
