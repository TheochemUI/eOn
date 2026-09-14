"""First-class min-mode / saddle / hessian (class-first Dimer facade)."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def _lj_h2(params=None):
    params = params or pyec.Parameters()
    params.potential = pyec.PotType.LJ
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    m.cell = np.eye(3) * 20.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.masses = np.ones(2)
    m.fixed = np.zeros(2, dtype=np.int64)
    m.periodic = False
    return m, pot, params


def test_minmode_symbols_bound():
    for name in (
        "Dimer",
        "ClassicDimer",
        "ImprovedDimer",
        "Lanczos",
        "Davidson",
        "MinModeSaddleSearch",
        "SaddleStatus",
        "Hessian",
        "get_prefactors",
        "all_free_atoms",
        "moved_atoms",
        "built_with_gprd",
        "NEB",
    ):
        assert hasattr(pyec, name), name


def test_dimer_default_is_improved():
    m, pot, params = _lj_h2()
    params.dimer_max_iterations = 50
    params.dimer_rotations_max = 5
    d = pyec.Dimer(m, params, pot)  # method="improved" default
    direction = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], dtype=np.float64)
    d.compute(m, direction)
    assert np.isfinite(d.eigenvalue)
    assert np.asarray(d.eigenvector).shape == (2, 3)


def test_dimer_classic_via_method():
    m, pot, params = _lj_h2()
    d = pyec.Dimer(m, params, pot, method="classic")
    direction = np.random.default_rng(0).normal(size=(2, 3))
    d.compute(m, direction)
    assert np.isfinite(d.eigenvalue)


def test_classic_dimer_low_level():
    m, pot, params = _lj_h2()
    d = pyec.ClassicDimer(m, params, pot)
    d.compute(m, np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]))
    assert np.isfinite(d.eigenvalue)


def test_improved_dimer_low_level():
    m, pot, params = _lj_h2()
    params.dimer_rotations_max = 5
    dimer = pyec.ImprovedDimer(m, params, pot)
    dimer.compute(m, np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]))
    assert np.isfinite(dimer.eigenvalue)


def test_dimer_accelerant_gp_gated():
    m, pot, params = _lj_h2()
    if pyec.built_with_gprd():
        d = pyec.Dimer(m, params, pot, method="improved", accelerant="gp")
        d.compute(m, np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]))
        assert np.isfinite(d.eigenvalue)
    else:
        with pytest.raises(RuntimeError, match="gprd|WITH_GPRD|accelerant"):
            pyec.Dimer(m, params, pot, accelerant="gp")


def test_lanczos_and_davidson():
    m, pot, params = _lj_h2()
    for Cls in (pyec.Lanczos, pyec.Davidson):
        solver = Cls(m, params, pot)
        solver.compute(m, np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]]))
        assert np.isfinite(solver.eigenvalue)
    d = pyec.Dimer(m, params, pot, method="lanczos")
    d.compute(m, np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]]))
    assert np.isfinite(d.eigenvalue)


def test_minmode_saddle_search():
    m, pot, params = _lj_h2()
    params.saddle_minmode_method = "dimer"
    params.saddle_max_iterations = 20
    params.dimer_rotations_max = 3
    pos = np.asarray(m.positions).copy()
    pos[1, 0] += 0.15
    m.positions = pos
    mode = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    ss = pyec.MinModeSaddleSearch(m, mode, float(m.potential_energy), params, pot)
    status = ss.run()
    assert isinstance(status, int)
    assert pyec.saddle_status_message(status)


def test_hessian_and_free_atoms():
    m, pot, params = _lj_h2()
    atoms = pyec.all_free_atoms(m)
    assert atoms.shape[0] == 2
    hess = pyec.Hessian(params, m)
    H = np.asarray(hess.get_hessian(m, atoms))
    assert H.shape[0] == 3 * len(atoms)


def test_parameters_dimer_saddle_attrs():
    p = pyec.Parameters()
    p.dimer_improved = False
    assert p.dimer_improved is False
    p.saddle_minmode_method = "lanczos"
    assert p.saddle_minmode_method == "lanczos"


def test_resolve_mobile_atoms_all_is_free():
    m, pot, params = _lj_h2()
    free = pyec.free_atom_indices(m)
    all_m = pyec.resolve_mobile_atoms(m, "All")
    np.testing.assert_array_equal(free, all_m)
    np.testing.assert_array_equal(free, pyec.all_free_atoms(m))


def test_resolve_mobile_atoms_intersects_fixed():
    m, pot, params = _lj_h2()
    m.fixed = np.array([1, 0], dtype=np.int64)
    mobile = pyec.resolve_mobile_atoms(m, np.array([0, 1], dtype=np.int64))
    np.testing.assert_array_equal(mobile, np.array([1], dtype=np.int64))
    # free mask unchanged
    assert int(m.n_free) == 1


def test_lanczos_active_atoms_zeros_environment():
    """PHVA mobile set: Krylov only on atoms=; free mask stays full."""
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(4)
    m.positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [0.0, 1.2, 0.0],
            [1.2, 1.2, 0.0],
        ],
        dtype=np.float64,
    )
    m.cell = np.eye(3) * 20.0
    m.atomic_numbers = np.array([1, 1, 1, 1], dtype=np.int64)
    m.masses = np.ones(4)
    m.fixed = np.zeros(4, dtype=np.int64)
    m.periodic = False
    assert int(m.n_free) == 4

    direction = np.zeros((4, 3), dtype=np.float64)
    direction[0, 0] = 1.0
    direction[1, 0] = -0.5
    atoms = np.array([0, 1], dtype=np.int64)

    for Cls in (pyec.Lanczos, pyec.Davidson):
        solver = Cls(m, params, pot)
        solver.compute(m, direction, atoms=atoms)
        assert np.isfinite(solver.eigenvalue)
        ev = np.asarray(solver.eigenvector)
        assert ev.shape == (4, 3)
        # Environment (atoms 2,3) must be zero in the mode
        assert np.linalg.norm(ev[2]) < 1e-12
        assert np.linalg.norm(ev[3]) < 1e-12
        # Active block nonzero
        assert np.linalg.norm(ev[:2]) > 1e-12
        # free mask not rewritten
        assert int(m.n_free) == 4


def test_lanczos_phva_atoms_param_restricts_krylov():
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.lanczos_phva_atoms = "0,1"
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(3)
    m.positions = np.array(
        [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.6, 1.0, 0.0]], dtype=np.float64
    )
    m.cell = np.eye(3) * 20.0
    m.atomic_numbers = np.ones(3, dtype=np.int64)
    m.masses = np.ones(3)
    m.fixed = np.zeros(3, dtype=np.int64)
    direction = np.zeros((3, 3), dtype=np.float64)
    direction[0, 0] = 1.0
    solver = pyec.Lanczos(m, params, pot)
    solver.compute(m, direction)  # uses params.lanczos_phva_atoms
    ev = np.asarray(solver.eigenvector)
    assert np.linalg.norm(ev[2]) < 1e-12
    assert int(m.n_free) == 3


def test_lanczos_default_matches_explicit_free():
    m, pot, params = _lj_h2()
    direction = np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]], dtype=np.float64)
    free = pyec.free_atom_indices(m)
    a = pyec.Lanczos(m, params, pot)
    a.compute(m, direction)
    b = pyec.Lanczos(m, params, pot)
    b.compute(m, direction, atoms=free)
    assert a.eigenvalue == pytest.approx(b.eigenvalue, abs=1e-6)
    eva = np.asarray(a.eigenvector).ravel()
    evb = np.asarray(b.eigenvector).ravel()
    cos = abs(float(np.dot(eva, evb)) / (np.linalg.norm(eva) * np.linalg.norm(evb) + 1e-30))
    assert cos == pytest.approx(1.0, abs=1e-5)
