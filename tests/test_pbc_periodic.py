"""pbc() honors per-axis periodic flags."""

import numpy as np
import pytest

pytest.importorskip("vesin")

from eon.geometry.pbc import pbc, pbc_eon_legacy


def test_all_periodic_wraps_x():
    box = np.eye(3) * 10.0
    got = pbc_eon_legacy(np.array([9.0, 0.0, 0.0]), box)
    np.testing.assert_allclose(got, [-1.0, 0.0, 0.0], atol=1e-12)


def test_open_x_leaves_span():
    box = np.eye(3) * 10.0
    r = np.array([9.0, 0.0, 0.0])
    got = pbc_eon_legacy(r, box, periodic=(False, True, True))
    np.testing.assert_allclose(got, r, atol=1e-12)


def test_false_is_identity():
    box = np.eye(3) * 10.0
    r = np.array([[9.0, 8.0, 7.0], [-11.0, 0.0, 3.0]])
    got = pbc(r, box, periodic=False)
    np.testing.assert_allclose(got, r, atol=1e-12)


def test_slab_wraps_xy_not_z():
    box = np.eye(3) * 10.0
    r = np.array([9.0, 9.0, 9.0])
    got = pbc_eon_legacy(r, box, periodic=(True, True, False))
    np.testing.assert_allclose(got, [-1.0, -1.0, 9.0], atol=1e-12)
