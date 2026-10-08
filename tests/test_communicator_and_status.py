"""Bundled jobs, the local client, and the aKMC status and continuous paths."""

import os
import stat
import subprocess
import sys
from io import StringIO
from pathlib import Path

import pytest

from eon.akmc import main as akmc_main
from eon.communicator import Communicator, CommunicatorError, Local, Script
from eon.config import ConfigClass
from eon.explorer import ServerMinModeExplorer
from eon.movie import make_movie
from eon.recycling import SB_Recycling
from eon.superbasinscheme import TransitionCounting
from tests.test_catalog_and_mains import _ini
from tests.test_library_bodies import _Comm, _atoms, _buf, _config, _states


class _StatusComm(_Comm):
    def get_queue_size(self):
        return 2


def _bundle_config(tmp_path):
    cfg = _config(tmp_path)
    pot = tmp_path / "pot"
    pot.mkdir()
    cfg.path_pot = str(pot)
    cfg.debug_keep_all_results = False
    cfg.debug_results_path = "debug-results"
    return cfg


def test_bundles_split_results_and_a_local_client_runs(tmp_path):
    cfg = _bundle_config(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    comm = Communicator(str(scratch), 2, config=cfg)
    note = StringIO("hello\n")
    jobs = [
        {"id": "7", "note.txt": StringIO("a\n")},
        {"id": "8", "note.txt": StringIO("b\n")},
    ]
    paths = list(comm.make_bundles(jobs, {"shared.txt": (note, 0o644)}))
    assert len(paths) == 1
    bundle = Path(paths[0])
    assert (bundle / "shared.txt").is_file()
    assert (bundle / "note_0.txt").is_file()
    assert (bundle / "note_1.txt").is_file()
    (bundle / "results_0.dat").write_text("0 termination_reason\n")
    (bundle / "results_1.dat").write_text("1 termination_reason\n")
    (bundle / "mode_0.dat").write_text("0 0 0\n")
    size, bundled = comm.get_bundle_size(str(bundle))
    assert size == 2 and bundled is True
    assert comm.get_bundle_size(["results.dat"]) == (1, False)
    harvested = []
    for group in comm.unbundle(scratch, lambda name: name == "7"):
        harvested.extend(group)
    assert len(harvested) == 2
    assert harvested[0]["results.dat"].read().startswith("0")

    local_scratch = tmp_path / "local-scratch"
    local_scratch.mkdir()
    local = Local(str(local_scratch), "true", 1, 1, config=cfg)
    assert Path(local.client).name == "true"
    with pytest.raises(CommunicatorError, match="client"):
        Local(str(local_scratch), "not-a-real-eon-client", 1, 1, config=cfg)
    out = tmp_path / "local-out"
    out.mkdir()
    local.submit_jobs([{"id": "job", "pos.con": StringIO("Pt\n")}], {})
    assert list(local.get_results(str(out), lambda name: True)) == []
    failed = subprocess.Popen(["/usr/bin/false"])
    failed.wait()
    jobdir = local_scratch / "failed"
    jobdir.mkdir()
    (jobdir / "stderr.dat").write_text("boom\n")
    local.check_job((failed, str(jobdir), open(os.devnull), open(os.devnull)))


def test_script_communicator_submits_and_cancels(tmp_path):
    cfg = _bundle_config(tmp_path)
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    submit = scripts / "submit.sh"
    submit.write_text("#!/bin/sh\necho 42\n")
    queued = scripts / "queued.sh"
    queued.write_text("#!/bin/sh\nif [ -f %s ]; then exit 0; fi\necho 42\n" % (tmp_path / "done"))
    cancel = scripts / "cancel.sh"
    cancel.write_text("#!/bin/sh\nexit 0\n")
    for path in (submit, queued, cancel):
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    comm = Script(
        str(scratch),
        1,
        "eon",
        str(scripts),
        "queued.sh",
        "cancel.sh",
        "submit.sh",
        config=cfg,
    )
    comm.submit_jobs([{"id": "3", "pos.con": StringIO("Pt\n")}], {})
    assert comm.get_queue_size() == 1
    assert comm.cancel_state(3) == 1
    (tmp_path / "done").write_text("x\n")
    assert comm.get_queue_size() == 0
    with pytest.raises(CommunicatorError):
        comm.check_command(1, "no", "queued.sh")


def test_status_names_the_basin_and_continuous_stops(tmp_path, monkeypatch, capsys):
    root, ini = _ini(
        tmp_path,
        "akmc",
        "\n".join(
            [
                "[Coarse Graining]",
                "use_mcamc = true",
                "superbasin_scheme = transition_counting",
                "[AKMC]",
                "max_kmc_steps = 1",
                "[Communicator]",
                "type = local",
            ]
        ),
    )
    cfg = ConfigClass()
    cfg.init(str(ini))
    cfg.comp_use_identical = False
    cfg.sb_max_size = 0
    from eon import fileio as io
    from eon.akmcstatelist import AKMCStateList
    from tests.test_library_bodies import _result

    reactant = _atoms(2.5)
    io.savecon(str(root / "pos.con"), reactant)
    states = AKMCStateList(
        300.0 / 11604.5, 20.0, 40.0, initial_state=str(root / "pos.con"), config=cfg
    )
    state = states.get_state(0)
    proc = state.add_process(_result(reactant, _atoms(2.7), _atoms(3.1)))
    product = states.get_state(states.get_product_state(0, proc).number)
    Path(cfg.sb_path).mkdir(parents=True, exist_ok=True)
    (Path(cfg.sb_path) / "storage").mkdir(exist_ok=True)
    scheme = TransitionCounting(cfg.sb_path, states, states.kT, 1, config=cfg)
    scheme.register_transition(state, product)
    monkeypatch.chdir(root)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _StatusComm())
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-s", "-q"])
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 0
    assert "superbasin" in capsys.readouterr().out

    monkeypatch.setattr("eon.akmc.akmc", lambda config, steps=0: 1)
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-C", "-q"])
    akmc_main()


def test_server_registers_a_saddle_result(tmp_path, monkeypatch):
    cfg, states, state, _product, _proc, reactant = _states(tmp_path)
    cfg.akmc_server_side_process_search = True
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_type = "local_inprocess"
    cfg.comm_job_buffer_size = 1
    cfg.process_search_minimization_offset = 0.1
    Path(cfg.path_incomplete).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    comm = _Comm()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ServerMinModeExplorer(states, state, state, config=cfg)
    explorer.make_jobs()
    wuid = next(iter(explorer.wuid_to_search_id))
    comm.results = [
        {
            "name": "0_%s" % wuid,
            "results": {
                "job_type": "saddle_search",
                "termination_reason": 1,
                "potential_energy_saddle": -0.2,
                "potential_energy_reactant": -1.0,
                "barrier_reactant_to_product": 0.8,
                "total_force_calls": 1,
            },
            "saddle.con": _buf(_atoms(2.7)),
            "mode.dat": _buf(_atoms(2.7)),
        }
    ]
    registered = explorer.register_results()
    assert registered == 1


def test_movie_frames_and_a_resumed_recycling_list(tmp_path, monkeypatch):
    cfg, states, state, product, _proc, _reactant = _states(tmp_path)
    monkeypatch.chdir(tmp_path)
    frames = make_movie("processes,0,1", str(tmp_path), states, separate_files=True)
    assert frames is None
    assert list((tmp_path / "movies").glob("processes_0.poscar.*"))
    path = tmp_path / "sb"
    path.mkdir()
    (path / "current_sb_states").write_text("header\n[0, 1]\n")
    recycling = SB_Recycling(states, state, product, 0.2, True, path, "askmc", None)
    assert recycling.sb_state_nums[0][0] == 0
