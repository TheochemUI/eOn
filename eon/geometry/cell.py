"""Cell matrix helpers (row-vector lattice, ASE/vesin convention)."""

from __future__ import annotations

import numpy as np


def length_angle_to_box(boxlengths, angles) -> np.ndarray:
    """Convert cell lengths + CON angles (degrees) to a 3×3 row-vector box.

    ``angles`` is ``(alpha, beta, gamma)`` as CON header line 4 specifies:
    alpha = angle(b, c), beta = angle(a, c), gamma = angle(a, b). Vector
    ``a`` lies on x and ``b`` in the xy plane.
    """
    box = np.zeros((3, 3), dtype=float)
    ang = np.asarray(angles, dtype=float) * (np.pi / 180.0)
    lengths = np.asarray(boxlengths, dtype=float)
    alpha, beta, gamma = ang[0], ang[1], ang[2]
    box[0][0] = 1.0
    box[1][0] = np.cos(gamma)
    box[1][1] = np.sin(gamma)
    box[2][0] = np.cos(beta)
    box[2][1] = (np.cos(alpha) - box[1][0] * box[2][0]) / box[1][1]
    box[2][2] = np.sqrt(1.0 - box[2][0] ** 2 - box[2][1] ** 2)
    box[0, :] *= lengths[0]
    box[1, :] *= lengths[1]
    box[2, :] *= lengths[2]
    return box


def box_to_length_angle(box) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`length_angle_to_box`: lengths and (alpha, beta, gamma)."""
    box = np.asarray(box, dtype=float)
    lengths = np.zeros(3, dtype=float)
    lengths[0] = np.linalg.norm(box[0, :])
    lengths[1] = np.linalg.norm(box[1, :])
    lengths[2] = np.linalg.norm(box[2, :])
    angles = np.zeros(3, dtype=float)
    # alpha = angle(b, c), beta = angle(a, c), gamma = angle(a, b)
    angles[0] = np.arccos(
        np.dot(box[1, :] / lengths[1], box[2, :] / lengths[2])
    )
    angles[1] = np.arccos(
        np.dot(box[0, :] / lengths[0], box[2, :] / lengths[2])
    )
    angles[2] = np.arccos(
        np.dot(box[0, :] / lengths[0], box[1, :] / lengths[1])
    )
    angles *= 180.0 / np.pi
    return lengths, angles
