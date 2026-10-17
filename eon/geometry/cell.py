"""Cell matrix helpers (row-vector lattice, ASE/vesin convention)."""

from __future__ import annotations

import math

import numpy as np

# Same constant as eonc::helpers::pi. The CON writer divides by this.
_PI = 3.14159265358979323846


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
    """Inverse of :func:`length_angle_to_box`: lengths and (alpha, beta, gamma).

    The cosine is the raw dot product over the product of the lengths,
    clamped to [-1, 1], then ``libm`` ``acos``. That is the same map the
    C++ writer uses, so both sides spell the same angle bytes.
    """
    box = np.asarray(box, dtype=float)
    lengths = np.empty(3, dtype=float)

    def row_norm(i: int) -> float:
        return math.sqrt(
            float(box[i, 0]) ** 2 + float(box[i, 1]) ** 2 + float(box[i, 2]) ** 2
        )

    lengths[0] = row_norm(0)
    lengths[1] = row_norm(1)
    lengths[2] = row_norm(2)

    def angle_deg(i: int, j: int) -> float:
        denom = float(lengths[i]) * float(lengths[j])
        if abs(denom) < 1e-300:
            cosine = 0.0
        else:
            cosine = (
                float(box[i, 0]) * float(box[j, 0])
                + float(box[i, 1]) * float(box[j, 1])
                + float(box[i, 2]) * float(box[j, 2])
            ) / denom
        if cosine > 1.0:
            cosine = 1.0
        elif cosine < -1.0:
            cosine = -1.0
        return math.acos(cosine) * 180.0 / _PI

    # alpha = angle(b, c), beta = angle(a, c), gamma = angle(a, b)
    angles = np.array(
        [angle_deg(1, 2), angle_deg(0, 2), angle_deg(0, 1)], dtype=float
    )
    return lengths, angles
