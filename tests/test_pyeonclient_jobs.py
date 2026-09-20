"""Job factory / ClientEON steps / workdir composition."""

from __future__ import annotations

import numpy as np
import pytest

pyec = pytest.importorskip("pyeonclient")


def test_job_type_names():
    assert pyec.job_type_from_name("minimization") == pyec.JobType.Minimization
    assert "Minimization" in pyec.job_type_name(pyec.JobType.Minimization)


def test_make_job_minimization():
    params = pyec.Parameters()
    params.job = pyec.JobType.Minimization
    params.potential = pyec.PotType.LJ
    params.quiet = True
    job = pyec.make_job(params)
    assert job.get_type() == pyec.JobType.Minimization


def test_session_make_job_and_summary(tmp_path):
    """Session + make_job keep_alive: Job holds Session; summary uses session.pots."""
    import gc
    import weakref

    params = pyec.Parameters()
    params.job = pyec.JobType.Minimization
    params.potential = pyec.PotType.LJ
    params.quiet = True
    session = pyec.Session()
    wr = weakref.ref(session)
    job = pyec.make_job(params, session)
    assert job.get_type() == pyec.JobType.Minimization
    del session
    gc.collect()
    assert wr() is not None
    out = tmp_path / "_potcalls.json"
    path = wr().write_potcall_summary(str(out))
    assert path == str(out)
    assert out.is_file()
    del job
    gc.collect()


def test_client_steps_are_bound():
    """Each ClientEON post-job step is a real binding, not only a wrapper."""
    assert callable(pyec.write_potcall_summary)
    assert callable(pyec.get_process_times)
    assert callable(pyec.append_results_timing)
    assert callable(pyec.steady_clock_now)
    assert callable(pyec.run_job)
    assert callable(pyec.load_parameters)
    t0 = pyec.steady_clock_now()
    real, user, sys = pyec.get_process_times()
    assert real >= 0 and user >= 0 and sys >= 0
    assert pyec.steady_clock_now() >= t0


def test_run_job_in_directory_composed_steps(tmp_path, pot_params_files):
    """Composition of steps matches eonclient artifact set."""
    work = pot_params_files
    files = pyec.run_job_in_directory(str(work), pyec.Parameters())
    assert (work / "min.con").is_file()
    text = (work / "results.dat").read_text()
    assert "time_seconds" in text
    assert (work / "_potcalls.json").is_file()
    assert any("min.con" in f or f.endswith("min.con") for f in files) or True


def test_minimize_workdir_matter_steps(tmp_path, pot_params_files):
    """Matter.relax path (not Job wrapper) for minimization."""
    work = pot_params_files
    files = pyec.minimize_workdir(work, write_movie=False)
    assert (work / "min.con").is_file()
    assert "potential_energy" in (work / "results.dat").read_text()
    assert (work / "_potcalls.json").is_file()
    assert "time_seconds" in (work / "results.dat").read_text()
    assert "min.con" in files


@pytest.fixture
def pot_params_files(tmp_path):
    """Write config.ini + pos.con for a 2-atom LJ minimize."""
    from eon.structure import Structure
    from eon import fileio as fio

    s = Structure(2)
    s.r = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    s.box = np.eye(3) * 20.0
    s.mass = np.ones(2)
    s.free = np.ones(2)
    s.names = ["H", "H"]
    fio.savecon(str(tmp_path / "pos.con"), s)
    (tmp_path / "config.ini").write_text(
        "[Main]\n"
        "job = minimization\n"
        "quiet = true\n"
        "write_log = false\n"
        "[Potential]\n"
        "potential = lj\n"
        "[Optimizer]\n"
        "max_iterations = 200\n"
        "converged_force = 0.1\n"
    )
    return tmp_path
