"""Seamless ASE Calculator → Matter Potential (zero-copy hot path)."""

from __future__ import annotations

import time

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")
pytest.importorskip("ase")
from ase import Atoms
from ase.calculators.lj import LennardJones


def test_potential_from_ase_lj_energy_forces():
    calc = LennardJones(epsilon=1.0, sigma=1.0, rc=10.0, ro=0.0, smooth=False)
    pot = pyec.potential_from_ase(calc)
    assert pot is not None
    # two atoms
    pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
    z = np.array([13, 13], dtype=np.int64)  # Al Z
    cell = np.eye(3) * 10.0
    E, F = pot.get_ef(pos, z, cell)
    assert np.isfinite(E)
    assert F.shape == (2, 3)
    # opposite forces along x
    assert F[0, 0] > 0 and F[1, 0] < 0


def test_from_ase_uses_atoms_calc():
    atoms = Atoms(
        "H2",
        positions=[[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
        cell=np.eye(3) * 12.0,
        pbc=True,
    )
    atoms.calc = LennardJones(epsilon=1.0, sigma=1.0, rc=5.0, smooth=False)
    # seamless: no explicit pot/params
    m = pyec.from_ase(atoms)
    assert m.n_atoms == 2
    e = float(m.potential_energy)
    assert np.isfinite(e)
    # force path via relax step or max_force
    f = float(m.max_force)
    assert np.isfinite(f)


def test_from_ase_without_calc_requires_pot():
    atoms = Atoms("H", positions=[[0, 0, 0]], cell=np.eye(3) * 5, pbc=False)
    with pytest.raises(ValueError, match="atoms.calc|potential"):
        pyec.from_ase(atoms)


def test_ase_hot_path_many_force_calls():
    """Repeated force() must stay bulk-fast (cached Atoms + view buffers).

    The old path built O(N) Python lists per call and was unusable for MD.
    400 evals on 48 atoms should finish well under a second on a workstation.
    """
    rng = np.random.default_rng(0)
    n = 48
    pos = rng.normal(scale=1.5, size=(n, 3)).astype(np.float64)
    # keep atoms from overlapping too hard
    pos[:, 0] += np.arange(n) * 1.2
    z = np.full(n, 13, dtype=np.int64)
    cell = np.eye(3) * (n * 1.5 + 10.0)
    calc = LennardJones(epsilon=0.1, sigma=1.0, rc=3.0, smooth=False)
    pot = pyec.potential_from_ase(calc)
    # warm-up (module import / Atoms cache)
    pot.get_ef(pos, z, cell)
    t0 = time.perf_counter()
    last_E = None
    for k in range(400):
        pos2 = pos + 1e-4 * k
        E, F = pot.get_ef(pos2, z, cell)
        last_E = E
        assert F.shape == (n, 3)
        assert np.isfinite(E)
        assert np.isfinite(F).all()
    elapsed = time.perf_counter() - t0
    assert last_E is not None
    # Loose ceiling: catches accidental O(N) Python list rebuild regressions.
    assert elapsed < 5.0, f"ASE hot path too slow: {elapsed:.3f}s for 400 evals"


def test_ase_force_matches_direct_calculator():
    """eOn AseCalcPotential must match ASE calculator results (bulk path)."""
    pos = np.array(
        [[0.0, 0.0, 0.0], [1.4, 0.1, 0.0], [0.2, 1.3, 0.0]], dtype=np.float64
    )
    z = np.array([1, 1, 1], dtype=np.int64)
    cell = np.eye(3) * 12.0
    calc = LennardJones(epsilon=1.0, sigma=1.0, rc=6.0, ro=0.0, smooth=False)
    atoms = Atoms(numbers=z, positions=pos, cell=cell, pbc=True)
    atoms.calc = calc
    E_ref = float(atoms.get_potential_energy())
    F_ref = np.asarray(atoms.get_forces(), dtype=np.float64)

    pot = pyec.potential_from_ase(
        LennardJones(epsilon=1.0, sigma=1.0, rc=6.0, ro=0.0, smooth=False)
    )
    E, F = pot.get_ef(pos, z, cell)
    assert E == pytest.approx(E_ref, rel=1e-10, abs=1e-10)
    assert np.allclose(F, F_ref, rtol=1e-10, atol=1e-10)


def test_matter_from_ase_core_roundtrip():
    """C++ matter_from_ase / matter_to_ase preserve geometry."""
    atoms = Atoms(
        "H2",
        positions=[[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]],
        cell=np.eye(3) * 10.0,
        pbc=True,
    )
    atoms.calc = LennardJones(epsilon=1.0, sigma=1.0, rc=5.0, smooth=False)
    m = pyec.from_ase(atoms)
    assert m.n_atoms == 2
    assert np.allclose(m.positions, atoms.get_positions())
    back = pyec.to_ase(m)
    assert len(back) == 2
    assert np.allclose(back.get_positions(), atoms.get_positions())
    assert np.allclose(back.get_cell(), atoms.get_cell())


def test_from_ase_uses_cpp_entry():
    assert hasattr(pyec, "matter_from_ase")
    assert callable(pyec.matter_from_ase)


def test_ase_potential_drives_matter_relax():
    """ASE calc is Matter's pot; relax/energy go through pyeonclient only."""
    handle = pyec.ase_potential(
        LennardJones(epsilon=1.0, sigma=1.0, rc=5.0, smooth=False)
    )
    params = pyec.Parameters()
    params.opt_max_iterations = 50
    params.opt_converged_force = 0.1
    m = pyec.Matter(handle.potential, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 12.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int32)
    m.masses = np.ones(2)
    m.periodic = False
    handle.bind_matter(m)
    assert handle.is_bound
    e0 = float(m.potential_energy)
    assert np.isfinite(e0)
    # second energy without geometry change
    e1 = float(m.potential_energy)
    assert e1 == pytest.approx(e0, rel=1e-12, abs=1e-12)
    out, ok = m.relax(inplace=False)
    assert np.isfinite(float(out.potential_energy))


def test_attach_ase_calculator_api():
    params = pyec.Parameters()
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 10.0
    m.atomic_numbers = np.array([13, 13], dtype=np.int32)
    m.masses = np.ones(2)
    pyec.attach_ase_calculator(
        m, LennardJones(epsilon=1.0, sigma=1.0, rc=6.0, smooth=False)
    )
    assert np.isfinite(float(m.potential_energy))
    assert np.isfinite(float(m.max_force))
