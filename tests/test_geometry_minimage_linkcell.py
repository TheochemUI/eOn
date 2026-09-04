"""minimage wrap and linkcell neighbors must agree with the vesin path."""

from __future__ import annotations

import numpy as np
import pytest

from eon.geometry.neighbors import neighbor_list, neighbor_list_linkcell
from eon.geometry.pbc import _pbc_numpy, pbc
from eon.structure import Structure

minimage = pytest.importorskip("minimage")
linkcell = pytest.importorskip("linkcell")


def _pair() -> Structure:
    p = Structure(2)
    p.r = np.array([[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]], dtype=float)
    p.box = np.eye(3) * 10.0
    p.names = ["Cu", "Cu"]
    p.mass = np.full(2, 63.5)
    p.free = np.ones(2)
    return p


def test_pbc_matches_minimage_and_numpy():
    p = _pair()
    disp = p.r[1] - p.r[0]
    via_pbc = pbc(disp, p.box)
    via_numpy = _pbc_numpy(disp, p.box, None)
    cell = minimage.Cell.from_vesin(p.box.tolist())
    via_mi = np.asarray(cell.displacement([0.0, 0.0, 0.0], disp.tolist()))
    np.testing.assert_allclose(via_pbc, via_mi, atol=1e-12)
    np.testing.assert_allclose(via_numpy, via_mi, atol=1e-12)
    assert abs(float(np.dot(via_mi, via_mi)) - 0.64) < 1e-12


def test_linkcell_neighbors_match_vesin_on_wrap():
    p = _pair()
    vesin = neighbor_list(p, cutoff=1.0)
    lc = neighbor_list_linkcell(p, cutoff=1.0)
    assert vesin == lc
    assert vesin[0] == [1]
    assert vesin[1] == [0]
