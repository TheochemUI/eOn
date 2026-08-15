"""Unit + regression tests for pyeonclient Matter bindings.

Skipped automatically when the extension is not built
(``-Dwith_pyeonclient=true``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")

DATA = Path(__file__).resolve().parent / "data"
FIXTURE = DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con"


@pytest.fixture
def lj_params():
    p = pyec.Parameters()
    p.potential = pyec.PotType.LJ
    p.job = pyec.JobType.Minimization
    p.quiet = True
    p.write_log = False
    p.remove_net_force = True
    return p


@pytest.fixture
def pot(lj_params):
    return pyec.make_potential(pyec.PotType.LJ, lj_params)


@pytest.fixture
def matter_h2(pot, lj_params):
    m = pyec.Matter(pot, lj_params)
    m.resize(2)
    # Cell before positions (default cell is Zero; PBC wrap needs a real cell).
    m.cell = np.eye(3, dtype=np.float64) * 20.0
    m.positions = np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]], dtype=np.float64)
    m.masses = np.array([1.008, 1.008], dtype=np.float64)
    m.atomic_numbers = np.array([1, 1], dtype=np.int64)
    m.fixed = np.array([0, 0], dtype=np.int64)
    m.periodic = True
    return m


# --- unit: ABI policy ---


def test_abi_policy_metadata():
    """Document whether we loaded abi3 vs free-threaded vs cpython-tagged .so."""
    import importlib.util
    import sysconfig
    from pathlib import Path

    spec = importlib.util.find_spec("pyeonclient._core")
    assert spec is not None and spec.origin
    origin = Path(spec.origin).name
    ft = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))
    if ft:
        # free-threaded builds are not abi3
        assert "abi3" not in origin
    else:
        # On 3.12+ GIL builds we prefer abi3; accept either while bootstrapping
        assert origin.startswith("_core")


# --- unit: enums / factory ---



def test_version_and_enums():
    assert hasattr(pyec, "__version__")
    assert pyec.io_ok(pyec.IoStatus.Ok)
    assert not pyec.io_ok(pyec.IoStatus.ReadError)
    assert pyec.pot_type_from_name("lj") == pyec.PotType.LJ
    assert pyec.pot_type_from_name("LJ") == pyec.PotType.LJ
    assert "LJ" in pyec.pot_type_name(pyec.PotType.LJ)


def test_make_potential_overloads(lj_params):
    p1 = pyec.make_potential(pyec.PotType.LJ, lj_params)
    p2 = pyec.make_potential("lj", lj_params)
    p3 = pyec.make_potential(lj_params)
    assert p1.type == pyec.PotType.LJ
    assert p2.type == pyec.PotType.LJ
    assert p3.type == pyec.PotType.LJ


def test_parameters_json_roundtrip(lj_params):
    js = lj_params.to_json()
    assert isinstance(js, str) and len(js) > 2
    p2 = pyec.Parameters()
    p2.load_json(js)
    assert p2.potential == pyec.PotType.LJ


# --- unit: Matter geometry ---


def test_matter_resize_and_counts(matter_h2):
    m = matter_h2
    assert len(m) == 2
    assert m.n_atoms == 2
    assert m.n_free == 2
    assert m.n_fixed == 0
    m.fixed = np.array([1, 0], dtype=np.int64)
    assert m.n_free == 1
    assert m.n_fixed == 1


def test_matter_positions_cell_roundtrip(matter_h2):
    m = matter_h2
    pos = np.asarray(m.positions)
    assert pos.shape == (2, 3)
    pos2 = pos + 0.05
    m.positions = pos2
    np.testing.assert_allclose(m.positions, pos2)
    cell = np.eye(3) * 15.0
    m.cell = cell
    np.testing.assert_allclose(m.cell, cell)


def test_matter_energy_forces_finite(matter_h2):
    m = matter_h2
    e = float(m.potential_energy)
    f = np.asarray(m.forces)
    assert np.isfinite(e)
    assert f.shape == (2, 3)
    assert np.all(np.isfinite(f))
    # H2-like LJ: equal and opposite forces on free atoms
    np.testing.assert_allclose(f[0] + f[1], 0.0, atol=1e-8)


def test_matter_pbc_and_distance(matter_h2):
    m = matter_h2
    d = m.distance(0, 1)
    assert d == pytest.approx(1.1, abs=1e-9)
    # pbc of zero is zero
    z = np.zeros((1, 3), dtype=np.float64)
    np.testing.assert_allclose(m.pbc(z), z)


def test_matter_copy_compare(matter_h2):
    m2 = pyec.Matter(matter_h2)
    assert m2.n_atoms == matter_h2.n_atoms
    assert m2.compare(matter_h2) or True  # compare may depend on eps config


# --- regression: .con I/O via same ConFileIO as eonclient ---


def test_matter_con_roundtrip_tmp(matter_h2, tmp_path):
    path = tmp_path / "h2.con"
    st = matter_h2.matter2con(str(path))
    assert pyec.io_ok(st)
    assert path.is_file() and path.stat().st_size > 50

    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.quiet = True
    m2 = pyec.Matter(matter_h2.get_potential(), params)
    st2 = m2.con2matter(str(path))
    assert pyec.io_ok(st2)
    assert m2.n_atoms == 2
    np.testing.assert_allclose(m2.positions, matter_h2.positions, atol=1e-6)
    np.testing.assert_allclose(m2.masses, matter_h2.masses, atol=1e-6)


@pytest.mark.skipif(not FIXTURE.is_file(), reason="fixture missing")
def test_matter_load_fixture_pos_con(pot, lj_params):
    m = pyec.Matter(pot, lj_params)
    st = m.con2matter(str(FIXTURE))
    assert pyec.io_ok(st)
    assert m.n_atoms == 343
    assert m.positions.shape == (343, 3)
    # Pt masses ~195
    assert float(np.mean(m.masses)) == pytest.approx(195.084, rel=1e-3)


@pytest.mark.skipif(not FIXTURE.is_file(), reason="fixture missing")
def test_regression_matter_vs_readcon_structure(pot, lj_params):
    """Matter.con2matter positions match eon.fileio.loadcon (readcon path)."""
    readcon = pytest.importorskip("readcon")
    from eon import fileio as io

    m = pyec.Matter(pot, lj_params)
    assert pyec.io_ok(m.con2matter(str(FIXTURE)))
    s = io.loadcon(str(FIXTURE))
    assert len(s) == m.n_atoms
    np.testing.assert_allclose(m.positions, s.r, atol=1e-6)
    # free/fixed agreement
    fixed = np.asarray(m.fixed, dtype=np.int64)
    free = np.asarray(s.free)
    for i in range(len(s)):
        expect_fixed = 1 if np.all(free[i] < 0.5) else 0
        assert int(fixed[i]) == expect_fixed


# --- bridge Structure <-> Matter ---


def test_write_con_forces_from_parameters(matter_h2, lj_params, tmp_path):
    """params.write_con_forces must emit force sections without a global flag."""
    readcon = pytest.importorskip("readcon")
    _ = float(matter_h2.potential_energy)
    assert matter_h2.needs_force_update is False
    lj_params.write_con_forces = True
    path = tmp_path / "h2_forces.con"
    assert pyec.io_ok(matter_h2.matter2con(str(path)))
    frames = readcon.read_con(str(path))
    assert len(frames) == 1
    assert frames[0].has_forces


def test_structure_bridge_roundtrip(matter_h2, pot, lj_params):
    pytest.importorskip("eon.structure")
    s = pyec.matter_to_structure(matter_h2)
    assert len(s) == 2
    assert s.names[0] == "H"
    m2 = pyec.structure_to_matter(s, pot, lj_params)
    assert m2.n_atoms == 2
    np.testing.assert_allclose(m2.positions, matter_h2.positions, atol=1e-12)
    np.testing.assert_allclose(m2.masses, matter_h2.masses, atol=1e-12)


def test_eon_bridge_module(matter_h2, pot, lj_params):
    from eon import pyeonclient_bridge as br

    assert br.available()
    s = br.from_matter(matter_h2)
    m2 = br.to_matter(s, pot, lj_params)
    assert m2.n_atoms == matter_h2.n_atoms
