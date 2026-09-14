"""RGPOT backend=metatomic through pyeonclient (dlopen libmetatomic_engine)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

pyec = pytest.importorskip("pyeonclient")

MODEL = os.environ.get("EON_PET_MAD_MODEL", "").strip()
ENGINE = os.environ.get("RGPOT_METATOMIC_ENGINE", "").strip()
_COOK_ROOT = os.environ.get("EON_PET_NEB_ROOT", "").strip()
COOK = Path(_COOK_ROOT) if _COOK_ROOT else None

need = (
    not MODEL
    or not Path(MODEL).is_file()
    or not ENGINE
    or not Path(ENGINE).is_file()
    or COOK is None
    or not (COOK / "min_reactant" / "pos.con").is_file()
)


def test_rgpot_params_bound():
    p = pyec.Parameters()
    p.rgpot_backend = "metatomic"
    p.rgpot_engine_path = "/tmp/libmetatomic_engine.so"
    p.rgpot_model_path = "/tmp/model.pt"
    p.rgpot_device = "cpu"
    assert p.rgpot_backend == "metatomic"
    assert "metatomic" in p.rgpot_engine_path
    assert pyec.PotType.RGPOT is not None


@pytest.mark.skipif(need, reason="needs MODEL + RGPOT_METATOMIC_ENGINE + cookbook pos")
def test_rgpot_metatomic_minimize(tmp_path):
    import shutil

    work = tmp_path / "rgpot_min"
    work.mkdir()
    shutil.copy(COOK / "min_reactant" / "pos.con", work / "pos.con")
    (work / "config.ini").write_text(
        "[Main]\n"
        "job = minimization\n"
        "random_seed = 706253457\n"
        "[Potential]\n"
        "potential = RGPOT\n"
        "[RgpotPot]\n"
        "backend = metatomic\n"
        f"engine_path = {ENGINE}\n"
        f"model_path = {MODEL}\n"
        "device = cpu\n"
        "[Optimizer]\n"
        "max_iterations = 2000\n"
        "opt_method = lbfgs\n"
        "max_move = 0.1\n"
        "converged_force = 0.01\n"
    )
    files = pyec.run_job_in_directory(work)
    text = (work / "results.dat").read_text()
    assert "GOOD" in text or "potential_energy" in text
    assert "RGPOT" in text or "rgpot" in text.lower() or "METATOMIC" in text
    # energy near native fat path
    energy = None
    for line in text.splitlines():
        sp = line.split()
        if len(sp) == 2 and sp[1] == "potential_energy":
            energy = float(sp[0])
    assert energy is not None
    assert abs(energy - (-56.6036)) < 5e-2  # same science path via engine
    assert "time_seconds" in text
    assert (work / "_potcalls.json").is_file()


@pytest.mark.skipif(need, reason="needs MODEL + engine + cookbook NEB inputs")
def test_rgpot_metatomic_neb(tmp_path):
    import shutil

    work = tmp_path / "rgpot_neb"
    work.mkdir()
    for name in ("reactant.con", "product.con", "idppPath.dat"):
        shutil.copy(COOK / name, work / name)
    (work / "config.ini").write_text(
        "[Main]\n"
        "job = nudged_elastic_band\n"
        "random_seed = 706253457\n"
        "[Potential]\n"
        "potential = RGPOT\n"
        "[RgpotPot]\n"
        "backend = metatomic\n"
        f"engine_path = {ENGINE}\n"
        f"model_path = {MODEL}\n"
        "device = cpu\n"
        "[Nudged Elastic Band]\n"
        "images = 10\n"
        "initializer = file\n"
        "initial_path_in = idppPath.dat\n"
        "minimize_endpoints = false\n"
        "climbing_image_method = true\n"
        "climbing_image_converged_only = true\n"
        "ci_after = 0.5\n"
        "ci_after_rel = 0.8\n"
        "energy_weighted = true\n"
        "ew_ksp_min = 0.972\n"
        "ew_ksp_max = 9.72\n"
        "ci_mmf = true\n"
        "ci_mmf_after = 0.1\n"
        "ci_mmf_after_rel = 0.5\n"
        "ci_mmf_penalty_strength = 1.5\n"
        "ci_mmf_penalty_base = 0.4\n"
        "ci_mmf_angle = 0.9\n"
        "ci_mmf_nsteps = 1000\n"
        "[Optimizer]\n"
        "max_iterations = 1000\n"
        "opt_method = lbfgs\n"
        "max_move = 0.1\n"
        "converged_force = 0.01\n"
        "[Debug]\n"
        "write_movies = true\n"
    )
    files = pyec.rgpot_metatomic_neb_workdir(work, engine_path=ENGINE, model_path=MODEL)
    text = (work / "results.dat").read_text()
    assert "GOOD" in text
    assert (work / "neb.con").is_file() or "neb.dat" in text or (work / "neb.dat").is_file()
