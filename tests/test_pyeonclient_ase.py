"""ASE ↔ Matter / Structure helpers (skip without pyeonclient and/or ase)."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")
ase = pytest.importorskip("ase")
from ase import Atoms
from ase.constraints import FixAtoms


@pytest.fixture
def pot_params():
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.quiet = True
    params.write_log = False
    pot = pyec.make_potential(params)
    return pot, params


def test_matter_to_ase_roundtrip(pot_params):
    pot, params = pot_params
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.cell = np.eye(3) * 15.0
    m.positions = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]])
    m.masses = np.array([1.0, 12.0])
    m.atomic_numbers = np.array([1, 6], dtype=np.int64)
    m.fixed = np.array([1, 0], dtype=np.int64)
    m.periodic = True

    atoms = pyec.to_ase(m)
    assert len(atoms) == 2
    assert list(atoms.get_atomic_numbers()) == [1, 6]
    np.testing.assert_allclose(atoms.get_positions(), m.positions)
    # FixAtoms present
    assert any(isinstance(c, FixAtoms) for c in atoms.constraints)

    m2 = pyec.from_ase(atoms, pot, params)
    assert m2.n_atoms == 2
    np.testing.assert_allclose(m2.positions, m.positions)
    np.testing.assert_array_equal(m2.fixed, m.fixed)
    np.testing.assert_array_equal(m2.atomic_numbers, m.atomic_numbers)


def test_structure_to_ase_roundtrip():
    from eon.structure import Structure

    s = Structure(2)
    s.r = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    s.box = np.eye(3) * 10.0
    s.mass = np.array([63.5, 63.5])
    s.free = np.array([1.0, 0.0])
    s.names = ["Cu", "Cu"]

    atoms = s.to_ase()
    assert len(atoms) == 2
    s2 = pyec.ase_to_structure(atoms)
    assert s2.names[0] == "Cu"
    np.testing.assert_allclose(s2.r, s.r)
    np.testing.assert_array_equal(s2.free[1], [0.0, 0.0, 0.0])


def test_matter_to_conframe_roundtrip(pot_params):
    pot, params = pot_params
    m = pyec.Matter(pot, params)
    m.resize(1)
    m.cell = np.eye(3) * 8.0
    m.positions = np.array([[1.0, 2.0, 3.0]])
    m.masses = np.array([1.0])
    m.atomic_numbers = np.array([1], dtype=np.int64)
    m.fixed = np.array([0], dtype=np.int64)

    frame = pyec.matter_to_conframe(m)
    m2 = pyec.conframe_to_matter(frame, pot, params)
    np.testing.assert_allclose(m2.positions, m.positions, atol=1e-6)
