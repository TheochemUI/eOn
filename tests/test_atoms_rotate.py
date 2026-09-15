"""Rotation helpers used to call bare cos/sin/acos (NameError)."""

import numpy as np
import pytest

pytest.importorskip("vesin")

from eon.atoms import get_mappings, get_rotation_matrix, rotm
from eon.structure import Structure


def test_rotm_identity_for_zero_angle():
    R = rotm(np.array([0.0, 0.0, 1.0]), 0.0)
    np.testing.assert_allclose(R, np.eye(3))


def test_get_rotation_matrix_z_90():
    R = get_rotation_matrix(np.array([0.0, 0.0, 1.0]), np.pi / 2)
    # rotate() left-multiplies row vectors: new_r = r @ R
    got = np.array([1.0, 0.0, 0.0]) @ R
    np.testing.assert_allclose(got, [0.0, 1.0, 0.0], atol=1e-12)


def test_get_mappings_identity():
    p = Structure(3)
    p.r = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]])
    p.box = np.eye(3) * 20.0
    p.names = ["H", "H", "H"]
    m = get_mappings(p, p.copy(), 1e-6, 2.0)
    assert m is not None
    assert len(m) == 3
