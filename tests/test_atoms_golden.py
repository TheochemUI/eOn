"""Golden masters: shipped vectorized atoms helpers vs pre-#368 scalar code.

Reference implementations live in tests/reference/atoms_pre_vectorization.py
and were extracted from git commit ee3012498 (last main before #368).
"""

from __future__ import annotations

import numpy as np
import pytest

from eon import atoms as at
from tests.reference.atoms_pre_vectorization import (
    brute_neighbor_list_scalar,
    free_r_scalar,
    get_process_atoms_scalar,
)


def _fcc_cell(n=4, a=4.0):
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


@pytest.mark.parametrize("n,free_stride", [(2, 2), (3, 3), (4, 1)])
def test_free_r_matches_pre_vectorization(n, free_stride):
    p = _fcc_cell(n)
    p.free[::free_stride] = 0
    live = p.free_r()
    ref = free_r_scalar(p)
    assert live.shape == ref.shape
    np.testing.assert_allclose(live, ref, rtol=0, atol=0)


@pytest.mark.parametrize("n,cutoff", [(2, 4.1), (3, 2.1), (3, 4.1)])
def test_brute_neighbor_list_matches_pre_vectorization(n, cutoff):
    p = _fcc_cell(n, a=2.0)
    live = at.brute_neighbor_list(p, cutoff)
    ref = brute_neighbor_list_scalar(p, cutoff)
    assert len(live) == len(ref) == len(p)
    for i in range(len(p)):
        assert sorted(live[i]) == sorted(ref[i]), f"neighbors differ for atom {i}"
        assert all(type(j) is int for j in live[i])


@pytest.mark.parametrize(
    "disp,eps,nshells",
    [
        (np.array([1.0, 0.0, 0.0]), 0.2, 1),
        (np.array([0.5, 0.5, 0.0]), 0.2, 1),
        (np.array([0.01, 0.0, 0.0]), 0.2, 1),  # argmax fallback path
        (np.array([2.0, 0.0, 0.0]), 0.2, 2),
    ],
)
def test_get_process_atoms_matches_pre_vectorization(disp, eps, nshells):
    r = _fcc_cell(3, a=3.0)
    p = r.copy()
    p.r[0] = p.r[0] + disp
    live = at.get_process_atoms(r, p, epsilon_r=eps, nshells=nshells)
    ref = get_process_atoms_scalar(r, p, epsilon_r=eps, nshells=nshells)
    assert sorted(live) == sorted(ref)
    assert all(type(i) is int for i in live)


def test_per_atom_norm_matches_pbc_scalar_loop():
    """per_atom_norm was not rewritten in #368; pin against explicit pbc norm."""
    p = _fcc_cell(3, a=2.5)
    v = p.r - p.r[0]
    live = at.per_atom_norm(v, p.box)
    ibox = np.linalg.inv(p.box)
    ref = np.array(
        [np.linalg.norm(at.pbc(v[i], p.box, ibox)) for i in range(len(p))]
    )
    np.testing.assert_allclose(live, ref, rtol=1e-12, atol=1e-12)
