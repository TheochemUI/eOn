"""minimage / linkcell vs eOn, LAMMPS, and GROMACS wrap oracles."""

from __future__ import annotations

import numpy as np
import pytest

import sys
from pathlib import Path

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

from geometry_oracles import (  # noqa: E402
    eon_cpp_wrap,
    eon_numpy_wrap,
    gromacs_pbc_dx,
    lammps_ortho_wrap,
    lammps_triclinic_wrap,
)

minimage = pytest.importorskip("minimage")
linkcell = pytest.importorskip("linkcell")

from eon.geometry.neighbors import neighbor_list, neighbor_list_linkcell
from eon.geometry.pbc import pbc
from eon.structure import Structure


def _ortho_box(lx=10.0, ly=11.0, lz=12.0) -> np.ndarray:
    return np.diag([lx, ly, lz]).astype(float)


def _sheared_rows() -> np.ndarray:
    # a=(10,0,0), b=(5, 5*sqrt(3), 0), c=(0,0,10) — 60 deg in xy.
    return np.array(
        [[10.0, 0.0, 0.0], [5.0, 8.660254037844386, 0.0], [0.0, 0.0, 10.0]]
    )


@pytest.mark.parametrize(
    "diff",
    [
        np.array([0.2, 0.0, 0.0]),
        np.array([9.2, 0.0, 0.0]),
        np.array([-9.2, 0.1, 0.0]),
        np.array([4.9, 5.4, -5.9]),
        np.array([0.0, 0.0, 0.0]),
    ],
)
def test_ortho_minimage_matches_eon_lammps_gromacs(diff):
    box = _ortho_box()
    mi = minimage.Cell.ortho(10.0, 11.0, 12.0)
    got = np.asarray(mi.displacement([0.0, 0.0, 0.0], diff.tolist()))
    np.testing.assert_allclose(got, eon_numpy_wrap(diff, box), atol=1e-12)
    np.testing.assert_allclose(got, eon_cpp_wrap(diff, box), atol=1e-12)
    np.testing.assert_allclose(
        got, lammps_ortho_wrap(diff, 10.0, 11.0, 12.0), atol=1e-12
    )
    np.testing.assert_allclose(
        got, gromacs_pbc_dx(diff, np.zeros(3), box), atol=1e-12
    )
    np.testing.assert_allclose(pbc(diff, box), got, atol=1e-12)


def test_sheared_minimage_matches_lammps_triclinic_and_eon():
    rows = _sheared_rows()
    # Restricted triclinic: xprd=10, yprd=8.660..., zprd=10, xy=5.
    diff = np.array([9.7, 0.1, 1.0]) - np.array([0.2, 0.1, 1.0])
    mi = minimage.Cell.from_vesin(rows.tolist())
    got = np.asarray(mi.displacement([0.2, 0.1, 1.0], [9.7, 0.1, 1.0]))
    lam = lammps_triclinic_wrap(diff, 10.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0)
    np.testing.assert_allclose(got, lam, atol=1e-12)
    np.testing.assert_allclose(got, eon_cpp_wrap(diff, rows), atol=1e-12)
    np.testing.assert_allclose(got, eon_numpy_wrap(diff, rows), atol=1e-12)
    assert abs(float(mi.dist2([0.2, 0.1, 1.0], [9.7, 0.1, 1.0])) - 0.25) < 1e-12


def test_gromacs_walk_matches_fractional_on_mild_shear():
    rows = _sheared_rows()
    x1 = np.array([9.7, 0.1, 1.0])
    x2 = np.array([0.2, 0.1, 1.0])
    gmx = gromacs_pbc_dx(x1, x2, rows)
    eon = eon_cpp_wrap(x1 - x2, rows)
    np.testing.assert_allclose(gmx, eon, atol=1e-12)


def test_linkcell_knn_matches_vesin_sorted_by_minimage():
    rng = np.random.default_rng(0)
    n = 24
    box = _ortho_box(12.0, 12.0, 12.0)
    r = rng.random((n, 3)) * np.array([12.0, 12.0, 12.0])
    p = Structure(n)
    p.r = r
    p.box = box
    p.names = ["Cu"] * n
    p.mass = np.full(n, 63.5)
    p.free = np.ones(n)

    k = 4
    cutoff = 20.0
    vesin = neighbor_list(p, cutoff)
    cell = minimage.Cell.ortho(12.0, 12.0, 12.0)
    for i in range(n):
        d2 = []
        for j in vesin[i]:
            d2.append((float(cell.dist2(r[i].tolist(), r[j].tolist())), j))
        d2.sort()
        expect = [j for _d, j in d2[:k]]
        lc = neighbor_list_linkcell(p, cutoff=cutoff, k=k)
        assert set(expect) == set(lc[i])


def test_linkcell_wrap_pair_matches_vesin():
    p = Structure(2)
    p.r = np.array([[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]])
    p.box = np.eye(3) * 10.0
    p.names = ["Cu", "Cu"]
    p.mass = np.full(2, 63.5)
    vesin = neighbor_list(p, 1.0)
    lc = neighbor_list_linkcell(p, 1.0)
    assert vesin == lc == [[1], [0]]
    assert abs(float(minimage.Cell.ortho(10, 10, 10).dist2(
        [0.2, 0, 0], [9.4, 0, 0]
    )) - 0.64) < 1e-12
