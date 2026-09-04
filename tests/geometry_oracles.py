"""Independent MIC oracles for eOn, LAMMPS, and GROMACS.

These are the published wrap rules, not calls into those codes.
minimage and linkcell are checked against this file. A disagreement
is a real kernel bug or a documented half-box convention split.
"""

from __future__ import annotations

import numpy as np


def eon_numpy_wrap(diff: np.ndarray, box: np.ndarray) -> np.ndarray:
    """eOn Python legacy: ``(frac % 1 + 1.5) % 1 - 0.5`` then H."""
    ibox = np.linalg.inv(box)
    frac = np.dot(diff, ibox)
    frac = (frac % 1.0 + 1.5) % 1.0 - 0.5
    return np.dot(frac, box)


def eon_cpp_wrap(diff: np.ndarray, box: np.ndarray) -> np.ndarray:
    """eOn C++ ``pbc::apply``: ``frac -= floor(frac + 0.5)`` then H.

    Matter.h: ``diff * cellInverse`` with lattice vectors as rows.
    """
    ibox = np.linalg.inv(box)
    frac = np.dot(np.asarray(diff, dtype=float), ibox)
    frac = frac - np.floor(frac + 0.5)
    return np.dot(frac, box)


def lammps_ortho_wrap(diff: np.ndarray, lx: float, ly: float, lz: float) -> np.ndarray:
    """LAMMPS ``Domain::minimum_image`` for an orthogonal box.

    If ``|d| > L/2``, subtract ``copysign(L, d)``.
    """
    out = np.asarray(diff, dtype=float).copy()
    half = np.array([0.5 * lx, 0.5 * ly, 0.5 * lz])
    span = np.array([lx, ly, lz])
    for i in range(3):
        if abs(out[i]) > half[i]:
            out[i] -= np.copysign(span[i], out[i])
    return out


def lammps_triclinic_wrap(
    diff: np.ndarray,
    xprd: float,
    yprd: float,
    zprd: float,
    xy: float,
    xz: float,
    yz: float,
) -> np.ndarray:
    """LAMMPS restricted-triclinic MIC via lamda coordinates.

    H columns are ``(xprd,0,0)``, ``(xy,yprd,0)``, ``(xz,yz,zprd)``.
    Wrap each lamda component to ``[-0.5, 0.5)``, then ``H @ lamda``.
    """
    h = np.array(
        [[xprd, xy, xz], [0.0, yprd, yz], [0.0, 0.0, zprd]], dtype=float
    )
    lamda = np.linalg.solve(h, np.asarray(diff, dtype=float))
    lamda = lamda - np.floor(lamda + 0.5)
    return h @ lamda


def gromacs_pbc_dx(x1: np.ndarray, x2: np.ndarray, box_vectors: np.ndarray) -> np.ndarray:
    """GROMACS ``pbc_dx`` triclinic walk (box vector ``d`` is row ``d``).

    From the last dimension to the first, subtract or add the whole
    box vector while that Cartesian component sits outside half the
    diagonal length. Orthorhombic boxes reduce to the usual ``L/2`` wrap.
    """
    dx = np.asarray(x1, dtype=float) - np.asarray(x2, dtype=float)
    box = np.asarray(box_vectors, dtype=float)
    for d in (2, 1, 0):
        hbox = 0.5 * box[d, d]
        while dx[d] > hbox:
            dx = dx - box[d]
        while dx[d] < -hbox:
            dx = dx + box[d]
    return dx
