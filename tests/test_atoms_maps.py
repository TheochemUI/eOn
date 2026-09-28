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


def test_identical_swaps_unlike_elements_on_the_same_sites():
    a = _pair([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]])
    b = _pair([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]])
    a.names = ["Cu", "Au"]
    b.names = ["Au", "Cu"]
    assert identical(a, b, epsilon_r=0.2) is True


def test_internal_motion_aligns_a_reversed_bond_and_a_twist():
    a = _pair([[1.0, 2.0, 3.0], [2.0, 2.0, 3.0], [1.0, 3.0, 3.0]])
    reversed_bond = _pair([[1.0, 2.0, 3.0], [0.0, 2.0, 3.0], [1.0, 3.0, 3.0]])
    twist = _pair([[1.0, 2.0, 3.0], [2.0, 2.0, 3.0], [1.0, 2.0, 4.0]])
    assert np.allclose(internal_motion(a, reversed_bond).r, a.r)
    assert np.allclose(internal_motion(a, twist).r, a.r)


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


def _structure(r, names):
    s = Structure(len(r))
    s.r = np.asarray(r, dtype=float)
    s.box = np.eye(3) * 50.0
    s.names = list(names)
    s.mass = np.ones(len(r))
    s.free = np.ones((len(r), 3))
    return s


def _brute_identical(a, b, eps):
    from itertools import permutations

    from eon.atoms import per_atom_norm

    ibox = np.linalg.inv(a.box)
    n = len(a)
    ok = np.array(
        [
            (per_atom_norm(a.r - b.r[i], a.box, ibox) <= eps)
            & (np.asarray(a.names) == b.names[i])
            for i in range(n)
        ]
    )
    return any(all(ok[i, j] for i, j in enumerate(p)) for p in permutations(range(n)))


def test_identical_agrees_with_every_permutation_on_random_clusters():
    rng = np.random.default_rng(20260928)
    eps = 0.3
    seen = {True: 0, False: 0}
    for _ in range(300):
        n = int(rng.integers(2, 7))
        r = rng.uniform(0.0, 2.0, size=(n, 3))
        names = rng.choice(["Pt", "Au"], size=n)
        perm = rng.permutation(n)
        r2 = r[perm] + rng.normal(0.0, 0.15, size=(n, 3))
        names2 = names[perm]
        if rng.random() < 0.3:
            k = int(rng.integers(n))
            names2 = names2.copy()
            names2[k] = "Au" if names2[k] == "Pt" else "Pt"
        a = _structure(r, names)
        b = _structure(r2, names2)
        want = _brute_identical(a, b, eps)
        assert identical(a, b, eps) == want
        seen[want] += 1
    # Both outcomes occur, so the comparison is not vacuous.
    assert seen[True] > 20 and seen[False] > 20
