"""eOn-native NEB path init (no ASE IDPP)."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def _lj_endpoints():
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.random_seed = 1
    pot = pyec.make_potential(pyec.PotType.LJ, params)

    def mk(x2: float):
        m = pyec.Matter(pot, params)
        m.resize(2)
        m.positions = np.array([[0.0, 0.0, 0.0], [x2, 0.0, 0.0]], dtype=np.float64)
        m.cell = np.eye(3) * 20.0
        m.atomic_numbers = np.array([1, 1], dtype=np.int64)
        m.masses = np.ones(2)
        m.fixed = np.zeros(2, dtype=np.int64)
        m.periodic = False
        return m

    return mk(1.2), mk(1.8), pot, params


def test_neb_linear_and_idpp_path():
    initial, final, pot, params = _lj_endpoints()
    n = 4
    lin = pyec.neb_linear_path(initial, final, n)
    assert len(lin) == n + 2  # endpoints + intermediates
    idpp = pyec.neb_idpp_path(initial, final, n, params)
    assert len(idpp) == n + 2
    # endpoints preserved (positions)
    assert np.allclose(np.asarray(lin[0].positions), np.asarray(initial.positions))
    assert np.allclose(np.asarray(idpp[0].positions), np.asarray(initial.positions))


def test_neb_initial_path_dispatches_idpp():
    initial, final, pot, params = _lj_endpoints()
    params.neb_init_method = pyec.NEBInit.IDPP
    path = pyec.neb_initial_path(initial, final, 3, params)
    assert len(path) == 5


def test_endpoint_neb_uses_init_method():
    initial, final, pot, params = _lj_endpoints()
    params.neb_images = 3
    params.neb_init_method = pyec.NEBInit.LINEAR
    params.neb_minimize_endpoints = False
    params.opt_max_iterations = 5
    params.opt_converged_force = 1.0
    neb = pyec.NudgedElasticBand(initial, final, params, pot)
    assert neb.n_path >= 2


def _slab_endpoints_mover_near_fixed():
    """Two frozen atoms plus one mover that starts next to a frozen site.

    IDPP pair forces between the mover and the nearby frozen atom pull that
    frozen atom if the optimizer is allowed to write its coordinates.
    """
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.random_seed = 1
    pot = pyec.make_potential(pyec.PotType.LJ, params)

    def mk(mover_xy):
        m = pyec.Matter(pot, params)
        m.resize(3)
        m.positions = np.array(
            [
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [mover_xy[0], mover_xy[1], 0.0],
            ],
            dtype=np.float64,
        )
        m.cell = np.eye(3) * 20.0
        m.atomic_numbers = np.array([1, 1, 1], dtype=np.int64)
        m.masses = np.ones(3)
        m.fixed = np.array([1, 1, 0], dtype=np.int64)
        m.periodic = False
        return m

    return mk((0.9, 0.0)), mk((0.9, 1.6)), pot, params


def _assert_fixed_unmoved(path, initial):
    frozen0 = np.asarray(initial.positions)[0]
    frozen1 = np.asarray(initial.positions)[1]
    for img in path:
        pos = np.asarray(img.positions)
        assert np.allclose(pos[0], frozen0, atol=1e-10), pos[0]
        assert np.allclose(pos[1], frozen1, atol=1e-10), pos[1]
        assert np.array_equal(np.asarray(img.fixed), np.asarray(initial.fixed))


def test_idpp_path_keeps_fixed_atoms():
    initial, final, pot, params = _slab_endpoints_mover_near_fixed()
    path = pyec.neb_idpp_path(initial, final, 3, params)
    assert len(path) == 5
    _assert_fixed_unmoved(path, initial)


def test_idpp_collective_path_keeps_fixed_atoms():
    initial, final, pot, params = _slab_endpoints_mover_near_fixed()
    path = pyec.neb_idpp_collective_path(initial, final, 3, params)
    assert len(path) == 5
    _assert_fixed_unmoved(path, initial)


def test_sidpp_path_keeps_fixed_atoms():
    initial, final, pot, params = _slab_endpoints_mover_near_fixed()
    path = pyec.neb_sidpp_path(initial, final, 3, params)
    assert len(path) == 5
    _assert_fixed_unmoved(path, initial)
