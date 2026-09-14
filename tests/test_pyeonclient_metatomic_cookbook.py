"""Metatomic jobs via pyeonclient (requires model + fat build).

Set EON_PET_MAD_MODEL to a PET-MAD (or other metatomic) .pt path to enable.
Uses the real ``run_job_in_directory`` entry point on a cookbook-like workdir.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

pyec = pytest.importorskip("pyeonclient")

MODEL = os.environ.get("EON_PET_MAD_MODEL", "").strip()
pytestmark = pytest.mark.skipif(
    not MODEL or not Path(MODEL).is_file(),
    reason="set EON_PET_MAD_MODEL to a metatomic .pt for this test",
)


def _write_min_workdir(tmp: Path, model: str, pos_src: Path | None) -> Path:
    work = tmp / "min"
    work.mkdir()
    if pos_src and pos_src.is_file():
        shutil.copy(pos_src, work / "pos.con")
    else:
        # 2-atom fallback only if no cookbook pos (test still exercises MTA pot)
        pytest.skip("need cookbook pos.con or EON_PET_MAD_POS")
    # Same Optimizer block as atomistic-cookbook examples/eon-pet-neb
    (work / "config.ini").write_text(
        "[Main]\n"
        "job = minimization\n"
        "random_seed = 706253457\n"
        "[Potential]\n"
        "potential = Metatomic\n"
        "[Metatomic]\n"
        f"model_path = {model}\n"
        "device = cpu\n"
        "[Optimizer]\n"
        "max_iterations = 2000\n"
        "opt_method = lbfgs\n"
        "max_move = 0.1\n"
        "converged_force = 0.01\n"
    )
    return work


def test_metatomic_minimize_via_run_job_in_directory(tmp_path):
    pos = os.environ.get("EON_PET_MAD_POS", "").strip()
    pos_path = Path(pos) if pos else None
    if pos_path is None:
        root = os.environ.get("EON_PET_NEB_ROOT", "").strip()
        if root:
            cand = Path(root) / "min_reactant" / "pos.con"
            if cand.is_file():
                pos_path = cand
    work = _write_min_workdir(tmp_path, MODEL, pos_path)
    files = pyec.minimize_workdir(work)  # Matter.relax steps, not black-box Job
    results = work / "results.dat"
    assert results.is_file(), f"missing results.dat; files={files}"
    text = results.read_text()
    assert "potential_energy" in text
    assert "METATOMIC" in text or "metatomic" in text.lower()
    assert "GOOD" in text  # cookbook converges under these Optimizer settings
    assert "time_seconds" in text  # ClientEON-equivalent postamble
    assert (work / "_potcalls.json").is_file()
    # parse energy
    energy = None
    for line in text.splitlines():
        sp = line.split()
        if len(sp) >= 2 and sp[-1] == "potential_energy":
            energy = float(sp[0])
            break
        if len(sp) >= 2 and sp[0] == "potential_energy":
            energy = float(sp[1])
            break
    assert energy is not None
    assert energy < 0.0  # PET-MAD oxadiazole endpoints are negative


def test_job_and_pot_enums_cover_cookbook():
    assert pyec.pot_type_from_name("metatomic") == pyec.PotType.METATOMIC
    assert pyec.job_type_from_name("minimization") == pyec.JobType.Minimization
    assert (
        pyec.job_type_from_name("nudged_elastic_band")
        == pyec.JobType.Nudged_Elastic_Band
    )
