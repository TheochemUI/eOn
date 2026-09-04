"""Periodic boundary conditions and distance helpers."""

from __future__ import annotations

from typing import Optional

import numpy as np


def _pbc_numpy(r, box, ibox: Optional[np.ndarray]) -> np.ndarray:
    if ibox is None:
        ibox = np.linalg.inv(box)
    vdir = np.dot(r, ibox)
    vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    return np.dot(vdir, box)


def pbc(r, box, ibox: Optional[np.ndarray] = None) -> np.ndarray:
    """Minimum-image convention for displacement(s).

    Uses :mod:`minimage` when that package is importable so the wrap
    matches the C/Rust cell used by linkcell. Falls back to the numpy
    wrap if minimage is not installed.

    Parameters
    ----------
    r : (3,) or (N, 3)
        Displacement vector(s).
    box : (3, 3)
        Cell matrix (rows = lattice vectors).
    ibox : (3, 3), optional
        Inverse of *box*; computed if omitted.
    """
    r = np.asarray(r, dtype=float)
    box = np.asarray(box, dtype=float)
    try:
        import minimage

        cell = minimage.Cell.from_vesin(box.tolist())
        if r.ndim == 1:
            return np.asarray(cell.displacement([0.0, 0.0, 0.0], r.tolist()), dtype=float)
        out = np.empty_like(r, dtype=float)
        zero = [0.0, 0.0, 0.0]
        for i, row in enumerate(np.atleast_2d(r)):
            out[i] = cell.displacement(zero, row.tolist())
        return out
    except Exception:
        return _pbc_numpy(r, box, ibox)


def per_atom_norm(v, box, ibox: Optional[np.ndarray] = None) -> np.ndarray:
    """Per-row Euclidean norm after PBC (shape ``(N,)``)."""
    diff = pbc(v, box, ibox)
    return np.sqrt(np.sum(diff**2.0, axis=1))


def per_atom_norm_gen(v, box, ibox: Optional[np.ndarray] = None):
    """Yield per-row norms after PBC (legacy generator API)."""
    diff = pbc(v, box, ibox)
    for d in diff:
        yield np.linalg.norm(d)
