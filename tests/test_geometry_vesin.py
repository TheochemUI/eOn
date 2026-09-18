"""Geometry kernels: PBC + vesin neighbor lists + process atoms."""

from __future__ import annotations

import numpy as np
import pytest

from eon.structure import Structure
from eon.geometry import (
    get_process_atoms,
    least_coordinated,
    neighbor_list,
    neighbor_list_pairs,
    neighbor_list_vectors,
    pbc,
    per_atom_norm,
    brute_neighbor_list,
)
from vesin import NeighborList as VesinNL


def _fcc(n=3, a=2.0) -> Structure:
    pts = [[i * a, j * a, k * a] for i in range(n) for j in range(n) for k in range(n)]
    p = Structure(len(pts))
    p.r = np.asarray(pts, dtype=float)
    p.box = np.eye(3) * (n * a)
    p.free = np.ones(len(pts))
    p.names = ["Cu"] * len(pts)
    p.mass = np.full(len(pts), 63.5)
    return p


def test_neighbor_list_uses_vesin_and_is_symmetric():
    p = _fcc(3, a=2.0)
    cutoff = 2.1
    nl = neighbor_list(p, cutoff)
    # Compare pair set to a direct vesin call
    calc = VesinNL(cutoff=cutoff, full_list=True)
    i, j = calc.compute(p.r, p.box, periodic=True, quantities="ij")
    pairs = {(int(a), int(b)) for a, b in zip(i, j) if a != b}
    built = set()
    for a, neigh in enumerate(nl):
        for b in neigh:
            built.add((a, b))
            assert a in nl[b]
            assert type(b) is int
    assert built == pairs


def test_brute_alias_matches_neighbor_list():
    p = _fcc(2, a=2.0)
    assert neighbor_list(p, 3.0) == brute_neighbor_list(p, 3.0)


@pytest.mark.parametrize(
    "free, expected",
    [
        ([1, 1, 1, 0], [0, 2]),
        ([[0, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 0]], [1]),
        ([0, 0, 0, 0], []),
    ],
)
def test_least_coordinated_selects_mobile_atoms(free, expected):
    p = Structure(4)
    p.r = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [8, 0, 0]], dtype=float)
    p.box = np.eye(3) * 20.0
    p.pbc = [False, False, False]
    p.free = free

    assert least_coordinated(p, 1.1) == expected


def test_neighbor_list_vectors_pbc():
    p = _fcc(2, a=2.0)
    vecs = neighbor_list_vectors(p, 2.1)
    assert len(vecs) == len(p)
    for center, vlist in enumerate(vecs):
        for v in vlist:
            assert v.shape == (3,)


def test_neighbor_list_pairs_matches_vesin_ijs():
    p = _fcc(3, a=2.0)
    cutoff = 2.1
    i, j, S = neighbor_list_pairs(p, cutoff)
    calc = VesinNL(cutoff=cutoff, full_list=True)
    vi, vj, vS = calc.compute(p.r, p.box, periodic=True, quantities="ijS")
    np.testing.assert_array_equal(i, vi)
    np.testing.assert_array_equal(j, vj)
    np.testing.assert_array_equal(S, vS)
    D = p.r[j] - p.r[i] + S.astype(float) @ p.box
    assert np.all(np.sum(D * D, axis=1) < cutoff * cutoff + 1e-12)


def test_pbc_packed_matches_rowwise_numpy():
    box = np.diag([10.0, 11.0, 12.0])
    diffs = np.array(
        [
            [0.2, 0.0, 0.0],
            [9.2, 0.0, 0.0],
            [-9.2, 0.1, 0.0],
            [4.9, 5.4, -5.9],
        ]
    )
    packed = pbc(diffs, box)
    rows = np.stack([pbc(row, box) for row in diffs])
    np.testing.assert_allclose(packed, rows, atol=1e-12)


def test_get_process_atoms_plain_ints_and_mobile():
    r = _fcc(3, a=3.0)
    p = r.copy()
    p.r[0] = p.r[0] + np.array([1.0, 0.0, 0.0])
    idxs = get_process_atoms(r, p, epsilon_r=0.2, nshells=1)
    assert 0 in idxs
    assert all(type(i) is int for i in idxs)
    # recycling metadata path
    line = "Indices of 'process' atoms = %s\n" % (repr(idxs),)
    assert "np.int64" not in line
    assert eval(line.split("=")[1].strip()) == idxs


def test_per_atom_norm_pbc():
    p = _fcc(2, a=3.0)
    v = p.r - p.r[0]
    nrm = per_atom_norm(v, p.box)
    assert nrm.shape == (len(p),)
    assert nrm[0] == pytest.approx(0.0, abs=1e-12)
