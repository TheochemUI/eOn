"""Periodic boundary conditions and distance helpers."""

from __future__ import annotations

from typing import Sequence, Union

import numpy as np

PeriodicSpec = Union[bool, Sequence[bool], np.ndarray]


def _periodic_mask(periodic: PeriodicSpec | None) -> np.ndarray | None:
    """None means wrap all axes. Otherwise a length-3 bool mask."""
    if periodic is None or periodic is True:
        return None
    if periodic is False:
        return np.zeros(3, dtype=bool)
    flags = np.asarray(periodic, dtype=bool)
    if flags.shape == ():
        return None if bool(flags) else np.zeros(3, dtype=bool)
    if flags.shape != (3,):
        raise ValueError(f"periodic must be bool or length-3, got {flags.shape}")
    if bool(flags.all()):
        return None
    return flags


def _pbc_numpy(
    r, box, ibox: np.ndarray | None, periodic: PeriodicSpec | None = None
) -> np.ndarray:
    if ibox is None:
        ibox = np.linalg.inv(box)
    vdir = np.dot(r, ibox)
    wrapped = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    flags = _periodic_mask(periodic)
    if flags is not None:
        wrapped = np.where(flags, wrapped, vdir)
    return np.dot(wrapped, box)


def pbc(
    r,
    box,
    ibox: np.ndarray | None = None,
    periodic: PeriodicSpec | None = None,
) -> np.ndarray:
    """Minimum-image convention for displacement(s).

    The kernel is :mod:`minimage` (same wrap linkcell uses). Packed
    ``(N, 3)`` rows go through ``Cell.wrap_many`` when that method
    exists. The numpy path remains as :func:`pbc_eon_legacy` for the
    eOn/GROMACS/LAMMPS agreement tests.

    Parameters
    ----------
    r : (3,) or (N, 3)
        Displacement vector(s).
    box : (3, 3)
        Cell matrix (rows = lattice vectors).
    ibox : (3, 3), optional
        Inverse of *box*; computed if omitted.
    periodic : bool or (3,) bool, optional
        Axes to wrap. True / omitted wraps all three. False wraps none.
        A length-3 mask wraps only the True axes. Partial masks skip
        minimage and use the numpy wrap so a free axis is left alone.
    """
    r = np.asarray(r, dtype=float)
    box = np.asarray(box, dtype=float)
    flags = _periodic_mask(periodic)
    if flags is not None:
        return _pbc_numpy(r, box, ibox, periodic)
    try:
        import minimage

        cell = minimage.Cell.from_vesin(box.tolist())
        if r.ndim == 1:
            wrap = getattr(cell, "wrap", None)
            if wrap is not None:
                return np.asarray(wrap(r.tolist()), dtype=float)
            return np.asarray(
                cell.displacement([0.0, 0.0, 0.0], r.tolist()), dtype=float
            )
        rows = np.atleast_2d(r)
        wrap_many = getattr(cell, "wrap_many", None)
        if wrap_many is not None:
            return np.asarray(wrap_many(rows), dtype=float)
        out = np.empty_like(rows, dtype=float)
        zero = [0.0, 0.0, 0.0]
        for i, row in enumerate(rows):
            out[i] = cell.displacement(zero, row.tolist())
        return out
    except ImportError:
        return _pbc_numpy(r, box, ibox)


def pbc_eon_legacy(
    r, box, ibox: np.ndarray | None = None, periodic: PeriodicSpec | None = None
) -> np.ndarray:
    """eOn numpy wrap. Kept so tests can compare it to minimage."""
    r = np.asarray(r, dtype=float)
    box = np.asarray(box, dtype=float)
    return _pbc_numpy(r, box, ibox, periodic)


def per_atom_norm(
    v,
    box,
    ibox: np.ndarray | None = None,
    periodic: PeriodicSpec | None = None,
) -> np.ndarray:
    """Per-row Euclidean norm after PBC (shape ``(N,)``)."""
    diff = pbc(v, box, ibox, periodic)
    return np.sqrt(np.sum(diff**2.0, axis=1))


def per_atom_norm_gen(
    v,
    box,
    ibox: np.ndarray | None = None,
    periodic: PeriodicSpec | None = None,
):
    """Yield per-row norms after PBC (legacy generator API)."""
    diff = pbc(v, box, ibox, periodic)
    for d in diff:
        yield np.linalg.norm(d)
