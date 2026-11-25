"""[RgpotPot] table_path has to reach the in-process rgpot EAM potential.

A config that names a setfl file and backend = eam must evaluate that file.
The fixture is a synthetic one-element table whose functions are written in
its header, so the test recomputes energy and forces from those formulas
without going through rgpot.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")

TABLE = (
    Path(__file__).resolve().parents[1]
    / "client"
    / "unit_tests"
    / "data"
    / "systems"
    / "rgpot_eam_table"
    / "synthetic_one_element.setfl"
)
ATOMIC_NUMBER = 18
CUTOFF = 4.0
POSITIONS = np.array(
    [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.4, 1.3, 0.2]], dtype=np.float64
)


def _taper(r: float) -> tuple[float, float]:
    if r >= CUTOFF:
        return 0.0, 0.0
    s = 1.0 - r / CUTOFF
    return s**3, -3.0 * s**2 / CUTOFF


def _density(r: float) -> tuple[float, float]:
    t, dt = _taper(r)
    e = math.exp(-1.2 * r)
    return e * t, e * (dt - 1.2 * t)


def _pair(r: float) -> tuple[float, float]:
    t, dt = _taper(r)
    e = 2.0 * math.exp(-1.8 * r)
    return e * t, e * (dt - 1.8 * t)


def _embedding(rho: float) -> tuple[float, float]:
    return -1.5 * (1.0 - math.exp(-rho)), -1.5 * math.exp(-rho)


def _reference(positions: np.ndarray) -> tuple[float, np.ndarray]:
    n = len(positions)
    rho = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                rho[i] += _density(float(np.linalg.norm(positions[i] - positions[j])))[0]
    energy = sum(_embedding(r)[0] for r in rho)
    forces = np.zeros((n, 3))
    for i in range(n):
        for j in range(i + 1, n):
            d = positions[i] - positions[j]
            r = float(np.linalg.norm(d))
            dphi = _pair(r)[1]
            dd = _density(r)[1]
            slope = dphi + (_embedding(rho[i])[1] + _embedding(rho[j])[1]) * dd
            energy += _pair(r)[0]
            forces[i] -= slope * d / r
            forces[j] += slope * d / r
    return energy, forces


def _ini(table: Path | str, table_format: str = "setfl") -> str:
    return f"""[Potential]
potential = rgpot

[RgpotPot]
backend = eam
table_path = {table}
table_format = {table_format}
neighbor_skin = 0.5
"""


def _matter(params, positions: np.ndarray):
    pot = pyec.make_potential(params)
    matter = pyec.Matter(pot, params)
    matter.resize(len(positions))
    matter.cell = np.eye(3, dtype=np.float64) * 20.0
    matter.positions = positions
    matter.masses = np.full(len(positions), 39.948, dtype=np.float64)
    matter.atomic_numbers = np.full(len(positions), ATOMIC_NUMBER, dtype=np.int64)
    matter.fixed = np.zeros(len(positions), dtype=np.int64)
    matter.periodic = True
    return matter


@pytest.mark.parametrize("table_format", ["setfl", "eam/alloy", "alloy"])
def test_rgpot_eam_table_path_is_evaluated(table_format: str) -> None:
    params = pyec.Parameters()
    params.load_ini_text(_ini(TABLE, table_format))
    matter = _matter(params, POSITIONS)

    expected_energy, expected_forces = _reference(POSITIONS)
    energy = float(matter.potential_energy)
    forces = np.asarray(matter.forces)

    assert math.isfinite(energy)
    assert energy != 0.0
    assert energy == pytest.approx(expected_energy, rel=1e-4)
    np.testing.assert_allclose(forces, expected_forces, rtol=1e-3, atol=1e-5)
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1e-10)


def test_rgpot_eam_forces_are_the_energy_gradient() -> None:
    params = pyec.Parameters()
    params.load_ini_text(_ini(TABLE))
    h = 1e-5
    matter = _matter(params, POSITIONS)
    forces = np.asarray(matter.forces).copy()
    for atom, axis in ((0, 0), (1, 0), (2, 1)):
        plus = POSITIONS.copy()
        plus[atom, axis] += h
        minus = POSITIONS.copy()
        minus[atom, axis] -= h
        ep = float(_matter(params, plus).potential_energy)
        em = float(_matter(params, minus).potential_energy)
        assert forces[atom, axis] == pytest.approx(-(ep - em) / (2 * h), rel=1e-4, abs=1e-6)


def test_rgpot_eam_table_path_is_required(tmp_path: Path) -> None:
    missing = tmp_path / "missing.setfl"
    refused = pyec.Parameters()
    refused.load_ini_text(_ini(missing))
    with pytest.raises(ValueError, match="cannot open"):
        pyec.make_potential(refused)


def test_rgpot_eam_unknown_table_format_is_refused() -> None:
    refused = pyec.Parameters()
    refused.load_ini_text(_ini(TABLE, "not-a-format"))
    with pytest.raises(RuntimeError, match="unknown table format"):
        pyec.make_potential(refused)


@pytest.mark.parametrize(
    "ini",
    [
        pytest.param(_ini(""), id="empty"),
        pytest.param(
            "[Potential]\npotential = rgpot\n\n[RgpotPot]\nbackend = eam\n",
            id="omitted",
        ),
    ],
)
def test_rgpot_eam_empty_or_omitted_table_path_is_named(ini: str) -> None:
    refused = pyec.Parameters()
    refused.load_ini_text(ini)
    with pytest.raises(
        RuntimeError, match=r"\[RgpotPot\] table_path is empty or missing"
    ):
        pyec.make_potential(refused)


def test_rgpot_eam_relative_table_path_follows_the_config_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    job = tmp_path / "job"
    (job / "tables").mkdir(parents=True)
    (job / "tables" / "one.setfl").write_bytes(TABLE.read_bytes())
    (job / "config.ini").write_text(_ini("tables/one.setfl"))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    params = pyec.Parameters()
    params.load(str(job / "config.ini"))
    positions = POSITIONS
    energy = float(_matter(params, positions).potential_energy)
    assert energy == pytest.approx(_reference(positions)[0], rel=1e-4)

    # Text has no config directory, so the same path is relative to the cwd.
    from_text = pyec.Parameters()
    from_text.load_ini_text(_ini("tables/one.setfl"))
    with pytest.raises(ValueError, match="cannot open"):
        pyec.make_potential(from_text)
    monkeypatch.chdir(job)
    energy = float(_matter(from_text, positions).potential_energy)
    assert energy == pytest.approx(_reference(positions)[0], rel=1e-4)


def test_rgpot_eam_relative_table_path_is_fixed_when_the_config_loads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    job = tmp_path / "job"
    (job / "tables").mkdir(parents=True)
    (job / "tables" / "one.setfl").write_bytes(TABLE.read_bytes())
    (job / "config.ini").write_text(_ini("tables/one.setfl"))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()

    monkeypatch.chdir(tmp_path)
    params = pyec.Parameters()
    params.load("job/config.ini")
    monkeypatch.chdir(elsewhere)
    energy = float(_matter(params, POSITIONS).potential_energy)
    assert energy == pytest.approx(_reference(POSITIONS)[0], rel=1e-4)
