"""Zero-copy Matter positions/forces views (binding-local, public Matter API)."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def _lj_h2():
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    m = pyec.Matter(pot, params)
    m.resize(2)
    m.positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], dtype=np.float64)
    m.cell = np.eye(3) * 20.0
    m.atomic_numbers = np.array([1, 1], dtype=np.int32)
    m.masses = np.ones(2, dtype=np.float64)
    m.fixed = np.zeros(2, dtype=np.int32)
    m.periodic = False
    return m, pot, params


def test_positions_view_shares_storage():
    m, pot, params = _lj_h2()
    a = m.positions
    b = m.positions
    assert a.shape == (2, 3)
    assert a.flags["C_CONTIGUOUS"]
    assert a.__array_interface__["data"][0] == b.__array_interface__["data"][0]


def test_positions_view_is_not_writeable():
    m, pot, params = _lj_h2()
    a = m.positions
    assert not a.flags["WRITEABLE"]
    with pytest.raises((ValueError, TypeError)):
        a += 0.1
    m.positions = np.array(a, copy=True) + 0.1
    assert m.needs_force_update


def test_bulk_set_positions_invalidates():
    m, pot, params = _lj_h2()
    _ = float(m.potential_energy)
    assert not m.needs_force_update
    m.positions = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]], dtype=np.float64)
    assert m.needs_force_update
    e = float(m.potential_energy)
    assert np.isfinite(e)


def test_forces_view_after_eval():
    m, pot, params = _lj_h2()
    f = m.forces
    assert f.shape == (2, 3)
    assert np.isfinite(f).all()
    f2 = m.forces
    assert f.__array_interface__["data"][0] == f2.__array_interface__["data"][0]


def test_cell_masses_z_fixed():
    m, pot, params = _lj_h2()
    assert m.cell.shape == (3, 3)
    assert m.masses.shape == (2,)
    assert m.atomic_numbers.shape == (2,)
    m.fixed = np.array([1, 0], dtype=np.int32)
    assert int(m.fixed[0]) == 1
    assert m.n_fixed == 1
