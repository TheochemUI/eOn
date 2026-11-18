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


_ORTHOGONAL_MEMO: dict = {}


def _wrap_orthogonal(rows: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """One image shift per axis, as ``minimage`` does in an orthogonal cell.

    A component at or above half the edge moves down one edge, one below
    minus half moves up one, and a component exactly at minus half stays.
    A row more than one edge out is moved once, not reduced, as there.
    """
    half = lengths / 2.0
    return np.where(
        rows >= half, rows - lengths, np.where(rows < -half, rows + lengths, rows)
    )


def _probe_rows(lengths: np.ndarray) -> np.ndarray:
    """Rows that separate one wrap rule from another: scattered rows up to
    one and a half edges out, and rows on and beside every half-edge tie."""
    rng = np.random.RandomState(0)
    rows = [rng.uniform(-1.5, 1.5, size=(256, 3)) * lengths]
    half = lengths / 2.0
    for axis in range(3):
        for sign in (1.0, -1.0):
            for edge in (half[axis], np.nextafter(half[axis], 0.0),
                         np.nextafter(half[axis], np.inf)):
                row = rng.uniform(-0.4, 0.4, size=3) * lengths
                row[axis] = sign * edge
                rows.append(row[None, :])
    rows.append(np.zeros((1, 3)))
    return np.vstack(rows)


def _orthogonal_lengths(cell, box: np.ndarray) -> np.ndarray | None:
    """Edge lengths that make :func:`_wrap_orthogonal` agree with *cell*, or None.

    minimage takes its edge lengths from the cell's volume and the face
    areas, which differs from the diagonal of the box in the last bit for
    some boxes. The candidate is therefore checked against the cell on the
    probe rows, bit for bit, once per box; a box that fails it, or one that
    is not orthogonal, wraps through the cell row by row.
    """
    key = box.tobytes()
    if key in _ORTHOGONAL_MEMO:
        return _ORTHOGONAL_MEMO[key]
    lengths = None
    diag = np.diagonal(box)
    is_ortho = getattr(cell, "is_ortho", None)
    if bool((diag > 0.0).all()) and (not callable(is_ortho) or is_ortho()):
        volume = diag[0] * diag[1] * diag[2]
        lengths = np.array(
            [
                volume / (diag[1] * diag[2]),
                volume / (diag[0] * diag[2]),
                volume / (diag[0] * diag[1]),
            ]
        )
        rows = _probe_rows(lengths)
        zero = [0.0, 0.0, 0.0]
        want = np.array([cell.displacement(zero, row.tolist()) for row in rows])
        got = _wrap_orthogonal(rows, lengths)
        if not (
            np.array_equal(got, want)
            and np.array_equal(np.signbit(got), np.signbit(want))
        ):
            lengths = None
    if len(_ORTHOGONAL_MEMO) >= 64:
        _ORTHOGONAL_MEMO.clear()
    _ORTHOGONAL_MEMO[key] = lengths
    return lengths


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
        lengths = _orthogonal_lengths(cell, box)
        if lengths is not None and bool(np.isfinite(rows).all()):
            return _wrap_orthogonal(rows, lengths)
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
