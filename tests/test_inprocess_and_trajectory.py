"""In-process minimization, a followed trajectory, and dynamics confidence."""

import runpy
import shutil
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from eon.akmc import kmc_step, main as akmc_main
from eon.akmcstate import AKMCState
from eon.amsel_superbasin_gate import amsel_discover_exit
from eon.communicator_inprocess import LocalInProcess
from tests.test_catalog_and_mains import _ini
from tests.test_library_bodies import _Comm, _atoms, _config, _states


def test_inprocess_minimizes_a_pair(tmp_path):
    cfg = _config(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    local = LocalInProcess(str(scratch), 1, config=cfg)
    assert local.cancel_state(0) == 0
    cluster = _atoms(2.5)
    local.submit_jobs([{"id": "pair", "structure": cluster}], {})
    finished = local.get_results()
    assert finished[0]["id"] == "pair"
    assert np.isfinite(finished[0]["_energy"])
    assert local.get_results() == []


def test_followed_trajectory_selects_the_matching_saddle(tmp_path):
    cfg, states, state, _product, proc, _reactant = _states(tmp_path)
    cfg.akmc_max_kmc_steps = 1
    for _ in range(100):
        state.inc_repeats()
    target = tmp_path / "target"
    procdata = target / "states" / "0" / "procdata"
    procdata.mkdir(parents=True)
    shutil.copy(state.proc_saddle_path(proc), procdata / ("saddle_%d.con" % proc))
    shutil.copy(state.proc_product_path(proc), procdata / ("product_%d.con" % proc))
    (target / "dynamics.txt").write_text(
        "\n".join(
            [
                "step reactant process product steptime total barrier rate energy",
                "----",
                "0 0 %d 1 1.0e-12 1.0e-12 0.30000 1.0e12 -1.00000" % proc,
                "",
            ]
        )
    )
    cfg.debug_target_trajectory = str(target)
    current, previous, time, steps = kmc_step(
        state, states, 0.0, states.kT, None, config=cfg
    )
    assert steps == 1
    assert previous.number == 0
    assert current.number == 1
    assert time > 0.0


def test_dynamics_confidence_counts_dynamics_saddles(tmp_path):
    cfg, _states_obj, state, _product, proc, _reactant = _states(tmp_path)
    cfg.akmc_confidence_scheme = "dynamics"
    cfg.recycling_on = True
    cfg.disp_moved_only = True
    from eon import fileio as io

    jobs = io.Table(str(Path(cfg.path_root) / "jobs.tbl"), ["state", "wuid", "type"])
    jobs.add_row({"state": 0, "wuid": 1, "type": "recycling"})
    assert state.get_confidence() == 0.0
    cfg.recycling_on = False
    state.increment_time(2.0, 300)
    path = Path(state.search_result_path)
    path.write_text(
        "header\n"
        "rule\n"
        "       1   dynamics      0.30000      0.20000          1          1          1    good-%d\n"
        % proc
    )
    confidence = state.get_confidence()
    assert 0.0 <= confidence <= 1.0
    assert isinstance(state, AKMCState)


def test_declined_restart_leaves_the_dynamics_file(tmp_path, monkeypatch, capsys):
    root, ini = _ini(tmp_path, "akmc")
    (root / "dynamics.txt").write_text("kept\n")
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-r"])
    monkeypatch.setattr("builtins.input", lambda prompt: "n")
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 1
    assert (root / "dynamics.txt").read_text() == "kept\n"
    assert "Not restarting" in capsys.readouterr().out


def test_displace_script_needs_a_reactant(tmp_path, monkeypatch):
    from eon import displace

    monkeypatch.setattr(sys, "argv", ["displace"])
    with pytest.raises(SystemExit) as caught:
        runpy.run_path(displace.__file__, run_name="__main__")
    assert caught.value.code == 1
    con = tmp_path / "reactant.con"
    from eon import fileio as io

    io.savecon(str(con), _atoms(2.5))
    monkeypatch.setattr(sys, "argv", ["displace", str(con), str(tmp_path / "out.con")])
    with pytest.raises(TypeError):
        runpy.run_path(displace.__file__, run_name="__main__")


def test_full_queue_makes_no_server_jobs(tmp_path, monkeypatch):
    cfg, states, state, _product, _proc, _reactant = _states(tmp_path)
    cfg.akmc_server_side_process_search = True
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_type = "local_inprocess"
    cfg.comm_job_buffer_size = 5
    cfg.comm_job_max_size = 1
    Path(cfg.path_incomplete).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)

    class _Full(_Comm):
        def queued_search_count(self):
            return 0

        def in_progress_search_count(self):
            return 1

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Full())
    from eon.explorer import ServerMinModeExplorer

    explorer = ServerMinModeExplorer(states, state, state, config=cfg)
    explorer.make_jobs()
    assert explorer.comm.submitted == []


def test_discover_exit_uses_the_stand_in_kernels(tmp_path, monkeypatch):
    cfg, states, state, _product, proc, _reactant = _states(tmp_path)
    cfg.amsel_e_min_init = 0.01
    cfg.debug_use_mean_time = False
    module = ModuleType("amsel")

    def discover_decide_status(*_args, **_kwargs):
        return ("accepted", [0], [1], [])

    def fpta(_transient, absorbing, _rates, _entry, _draw):
        return 1.25, [1.0] * len(absorbing)

    def mrm(_transient, absorbing, _rates, _entry):
        return 3.5, [1.0] * len(absorbing), None

    module.discover_decide_status = discover_decide_status
    module.fpta = fpta
    module.mrm = mrm
    monkeypatch.setitem(sys.modules, "amsel", module)
    sampled = amsel_discover_exit(state, states.get_state, cfg, uniform=lambda: 0.2)
    assert sampled is not None
    assert sampled.kernel == "fpta"
    assert sampled.mean_time == pytest.approx(1.25)
    assert sampled.proc_id == proc
    assert sampled.time_is_sample is True

    class _Problem:
        def __init__(self, transient, absorbing, rates):
            self.absorbing = absorbing

        def mrm(self, entry):
            return SimpleNamespace(tau_total=4.0, rate_to_absorbing=[1.0])

        def fpta(self, entry, draw):
            return SimpleNamespace(t_exit=0.5, weights=[1.0])

    del module.fpta
    del module.mrm
    module.AmcProblem = _Problem
    cfg.debug_use_mean_time = True
    mean = amsel_discover_exit(state, states.get_state, cfg, uniform=lambda: 0.2)
    assert mean.kernel == "mrm"
    assert mean.mean_time == pytest.approx(4.0)
    assert mean.time_is_sample is False
