"""Tests that drive the real eon.atoms hot paths (vectorized implementations)."""

from __future__ import annotations

import numpy as np
import pytest

from eon import atoms as at


def _fcc_cell(n=4, a=4.0):
    """Simple cubic lattice of n^3 atoms for reproducible geometry tests."""
    pts = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                pts.append([i * a, j * a, k * a])
    p = at.Atoms(len(pts))
    p.r = np.asarray(pts, dtype=float)
    p.box = np.eye(3) * (n * a)
    p.free = np.ones(len(pts))
    p.names = ["Cu"] * len(pts)
    p.mass = np.full(len(pts), 63.5)
    return p


def test_free_r_matches_mask():
    p = _fcc_cell(3)
    p.free[::2] = 0
    fr = p.free_r()
    expected = p.r[p.atom_is_free()]
    assert fr.shape == expected.shape
    np.testing.assert_allclose(fr, expected)


def test_brute_neighbor_list_symmetric_and_cutoff():
    p = _fcc_cell(3, a=2.0)
    cutoff = 2.1
    nl = at.brute_neighbor_list(p, cutoff)
    assert len(nl) == len(p)
    # Self never listed
    for i, neigh in enumerate(nl):
        assert i not in neigh
    # Pair symmetry
    for i, neigh in enumerate(nl):
        for j in neigh:
            assert i in nl[j]
    # Distance check against real pbc path
    ibox = np.linalg.inv(p.box)
    for i, neigh in enumerate(nl):
        for j in neigh:
            d = np.linalg.norm(at.pbc(p.r[i] - p.r[j], p.box, ibox))
            assert d < cutoff


def test_get_process_atoms_includes_moved():
    r = _fcc_cell(3, a=3.0)
    p = r.copy()
    # Move atom 0 significantly
    p.r[0] = p.r[0] + np.array([1.0, 0.0, 0.0])
    idxs = at.get_process_atoms(r, p, epsilon_r=0.2, nshells=1)
    assert 0 in idxs
    assert len(idxs) >= 1


def test_get_process_atoms_plain_python_ints():
    """Indices must be built-in int so recycling metadata repr/eval works."""
    r = _fcc_cell(3, a=3.0)
    p = r.copy()
    p.r[0] = p.r[0] + np.array([1.0, 0.0, 0.0])
    idxs = at.get_process_atoms(r, p, epsilon_r=0.2, nshells=1)
    assert idxs
    assert all(type(i) is int for i in idxs)
    assert not any(isinstance(i, np.integer) for i in idxs)


def test_get_process_atoms_recycling_metadata_roundtrip():
    """Mirrors recycling.write/read_recycling_metadata process_atoms path."""
    r = _fcc_cell(3, a=3.0)
    p = r.copy()
    p.r[0] = p.r[0] + np.array([1.0, 0.0, 0.0])
    process_atoms = at.get_process_atoms(r, p, epsilon_r=0.2, nshells=1)
    line = "Indices of 'process' atoms = %s\n" % (repr(process_atoms),)
    assert "np.int64" not in line
    assert "numpy" not in line
    parsed = eval(line.split("=")[1].strip())
    assert parsed == process_atoms
    assert all(type(i) is int for i in parsed)


def test_get_process_atoms_fallback_argmax_plain_int():
    """When no atom exceeds epsilon_r, argmax path still returns plain ints."""
    r = _fcc_cell(2, a=3.0)
    p = r.copy()
    # Tiny displacement below epsilon
    p.r[1] = p.r[1] + np.array([0.01, 0.0, 0.0])
    idxs = at.get_process_atoms(r, p, epsilon_r=0.2, nshells=1)
    assert idxs
    assert all(type(i) is int for i in idxs)


def test_brute_neighbor_list_plain_ints():
    p = _fcc_cell(3, a=2.0)
    nl = at.brute_neighbor_list(p, 2.1)
    for neigh in nl:
        assert all(type(i) is int for i in neigh)


def test_per_atom_norm_shape():
    p = _fcc_cell(2)
    v = p.r - p.r[0]
    nrm = at.per_atom_norm(v, p.box)
    assert nrm.shape == (len(p),)
    assert nrm[0] == pytest.approx(0.0, abs=1e-12)
