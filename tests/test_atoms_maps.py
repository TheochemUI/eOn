"""Atom matching is bijective, and internal motion translates atom 0 the right way."""

from __future__ import annotations

import numpy as np

from eon.atoms import identical, internal_motion, point_energy_match
from eon.structure import Structure


def _pair(positions):
    p = Structure(len(positions))
    p.r = np.asarray(positions, dtype=float)
    p.box = np.eye(3) * 40.0
    p.names = ["Cu"] * len(positions)
    p.mass = np.ones(len(positions))
    p.free = np.ones((len(positions), 3))
    return p


def test_identical_rejects_a_map_that_reuses_one_atom():
    a = _pair([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    b = _pair([[0.01, 0.0, 0.0], [0.02, 0.0, 0.0]])
    assert identical(a, b, epsilon_r=0.1) is False


def test_identical_accepts_a_swap():
    a = _pair([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    b = _pair([[3.01, 0.0, 0.0], [0.01, 0.0, 0.0]])
    assert identical(a, b, epsilon_r=0.1) is True


def test_internal_motion_puts_atom_zero_on_the_reference():
    a = _pair([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    b = _pair([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    moved = internal_motion(a, b)
    assert np.allclose(moved.r[0], a.r[0])


def test_point_energy_match_forwards_use_identical(monkeypatch):
    seen = {}

    def fake_load(_path):
        return _pair([[0.0, 0.0, 0.0]])

    def fake_match(*args, **kwargs):
        seen["indistinguishable"] = args[4]
        seen["use_identical"] = kwargs["use_identical"]
        return True

    monkeypatch.setattr("eon.fileio.loadcon", fake_load)
    monkeypatch.setattr("eon.atoms.match", fake_match)
    assert point_energy_match("a.con", 1.0, "b.con", 1.0, 0.1, 0.1, 3.0, use_identical=True)
    assert seen == {"indistinguishable": True, "use_identical": True}
