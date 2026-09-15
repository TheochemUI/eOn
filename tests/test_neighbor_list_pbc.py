"""neighbor_list honors per-axis periodic flags."""

import numpy as np
import pytest

pytest.importorskip("vesin")

from eon.geometry.neighbors import neighbor_list
from eon.structure import Structure


def _pair_across_x(box=10.0):
    p = Structure(2)
    p.r = np.array([[0.1, 5.0, 5.0], [9.9, 5.0, 5.0]])
    p.box = np.eye(3) * box
    p.names = ["Al", "Al"]
    return p


def test_structure_defaults_all_periodic():
    p = Structure(1)
    assert p.periodic.tolist() == [True, True, True]


def test_wrap_pair_is_neighbor_when_periodic():
    p = _pair_across_x()
    nl = neighbor_list(p, cutoff=0.5)
    assert 1 in nl[0]
    assert 0 in nl[1]


def test_wrap_pair_dropped_when_x_open():
    p = _pair_across_x()
    p.periodic = np.array([False, True, True])
    nl = neighbor_list(p, cutoff=0.5)
    assert nl[0] == []
    assert nl[1] == []


def test_explicit_periodic_overrides_structure():
    p = _pair_across_x()
    p.periodic = np.array([False, True, True])
    nl = neighbor_list(p, cutoff=0.5, periodic=True)
    assert 1 in nl[0]
