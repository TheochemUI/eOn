"""NudgedElasticBand first-class bindings (not Job.run black box)."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")

MODEL = os.environ.get("EON_PET_MAD_MODEL", "").strip()
_COOK_ROOT = os.environ.get("EON_PET_NEB_ROOT", "").strip()
COOK_ROOT = Path(_COOK_ROOT) if _COOK_ROOT else None


def test_neb_symbols_bound():
    assert hasattr(pyec, "NudgedElasticBand")
    assert hasattr(pyec, "NEBStatus")
    assert hasattr(pyec, "neb_read_file_paths")
    assert hasattr(pyec, "neb_write_results")
    assert hasattr(pyec, "neb_linear_path")
    assert hasattr(pyec, "pot_registry_total_force_calls")
    assert pyec.NEBStatus.GOOD is not None


def test_neb_linear_path_lj():
    params = pyec.Parameters()
    params.potential = pyec.PotType.LJ
    params.job = pyec.JobType.Nudged_Elastic_Band
    params.neb_images = 3
    pot = pyec.make_potential(pyec.PotType.LJ, params)
    a = pyec.Matter(pot, params)
    b = pyec.Matter(pot, params)
    a.resize(2)
    b.resize(2)
    a.positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    b.positions = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    a.cell = np.eye(3) * 20.0
    b.cell = np.eye(3) * 20.0
    a.atomic_numbers = np.array([1, 1], dtype=np.int64)
    b.atomic_numbers = np.array([1, 1], dtype=np.int64)
    a.masses = np.ones(2)
    b.masses = np.ones(2)
    path = pyec.neb_linear_path(a, b, 3)
    assert len(path) == 5  # 3 intermediates + 2 endpoints
    neb = pyec.NudgedElasticBand(a, b, params, pot)
    assert neb.num_images == 3
    assert neb.n_path == 5
    # one force update cycle
    neb.update_forces(False)
    assert neb.image(0).n_atoms == 2
    e0 = neb.image_energy(0)
    assert np.isfinite(e0)


@pytest.mark.skipif(
    not MODEL
    or not Path(MODEL).is_file()
    or COOK_ROOT is None
    or not (COOK_ROOT / "idppPath.dat").is_file(),
    reason="set EON_PET_MAD_MODEL and EON_PET_NEB_ROOT for cookbook NEB",
)
def test_neb_compute_cookbook_metatomic(tmp_path):
    """Full energy-weighted CI+OCI NEB via NudgedElasticBand.compute (not Job)."""
    import shutil

    work = tmp_path / "neb"
    work.mkdir()
    for name in ("reactant.con", "product.con", "idppPath.dat", "config.ini"):
        shutil.copy(COOK_ROOT / name, work / name)
    # ensure model path absolute in config
    cfg = (work / "config.ini").read_text()
    if "model_path" not in cfg or "pet-mad" not in cfg:
        (work / "config.ini").write_text(
            (COOK_ROOT / "config.ini").read_text()
        )
    files = pyec.neb_workdir(work)
    results = work / "results.dat"
    assert results.is_file()
    text = results.read_text()
    assert "GOOD" in text
    assert "time_seconds" in text
    assert (work / "neb.dat").is_file()
    assert (work / "neb.con").is_file(), "neb.con required for cookbook plots"
    assert (work / "sp.con").is_file()
    assert (work / "_potcalls.json").is_file()
    # peaks written when OCI extrema exist (cookbook golden has peak00/01)
    peaks = list(work.glob("peak*_pos.con"))
    assert len(peaks) >= 1, "expected MMF peak pos files for oxadiazole NEB"
    # parity vs golden
    golden = COOK_ROOT / "results.dat"
    if golden.is_file():
        import re

        def img(p):
            d = {}
            for line in Path(p).read_text().splitlines():
                m = re.match(r"([-\d.]+)\s+image(\d+)_energy", line)
                if m:
                    d[int(m.group(2))] = float(m.group(1))
            return d

        ne, ge = img(results), img(golden)
        md = max(abs(ne[k] - ge[k]) for k in set(ne) & set(ge))
        assert md <= 1e-3
