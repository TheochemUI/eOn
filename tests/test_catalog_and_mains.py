"""Catalog storage, entry points, and the explorer branches still left open."""

import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from eon.akmc import kmc_step, main as akmc_main
from eon.basinhopping import main as hop_main
from eon.config import ConfigClass
from eon.escaperate import main as escape_main
from eon.explorer import ClientMinModeExplorer
from eon.parallelreplica import main as replica_main
from eon.process_catalog import insert, make_suggestion, query, was_queried
from eon.superbasin import Superbasin
from tests.test_library_bodies import _Comm, _atoms, _buf, _config, _states


class _Record:
    def __init__(self, barrier, prefactor, **kwargs):
        self.barrier_ev = barrier
        self.prefactor = prefactor
        for key, value in kwargs.items():
            setattr(self, key, value)


class _Store:
    _by_path = {}

    def __init__(self, path):
        self.path = str(path)
        self.rows = self._by_path.setdefault(self.path, {})

    def lookup(self, env):
        return list(self.rows.get(bytes(env), []))

    def insert(self, env, record):
        self.rows.setdefault(bytes(env), []).append(record)


def _install_amsel(monkeypatch):
    module = ModuleType("amsel")

    def _process(env_hash, name, barrier, prefactor, **kwargs):
        return _Record(barrier, prefactor, **kwargs)

    module.KdbProcess = _process
    module.KdbStore = _Store
    monkeypatch.setitem(sys.modules, "amsel", module)


def _ini(directory, job, extra=""):
    root = directory / "run"
    root.mkdir(parents=True)
    from eon import fileio as io

    io.savecon(str(root / "pos.con"), _atoms(2.5))
    path = root / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                f"job = {job}",
                "temperature = 300",
                "[Paths]",
                f"main_directory = {root}",
                f"results = {root}",
                f"states = {root / 'states'}",
                extra,
                "",
            ]
        )
    )
    return root, path


def test_catalog_stores_and_offers_one_saddle(tmp_path, monkeypatch):
    _install_amsel(monkeypatch)
    cfg, _states_obj, state, _product, proc, _reactant = _states(tmp_path)
    cfg.kdb_on = True
    cfg.kdb_path = str(tmp_path / "kdb")
    cfg.kdb_scratch_path = str(tmp_path / "scratch")
    Path(cfg.kdb_path).mkdir()
    assert insert(state, proc, cfg) is True
    assert insert(state, proc, cfg) is True
    assert query(state, cfg) is True
    assert was_queried(state, cfg) is True
    displacement, mode = make_suggestion(cfg, state)
    assert len(displacement) == len(state.get_reactant())
    assert mode.shape == state.get_reactant().r.shape
    assert make_suggestion(cfg, state) == (None, None)


def test_broken_target_trajectory_raises(tmp_path):
    cfg, states, state, _product, _proc, _reactant = _states(tmp_path)
    cfg.akmc_max_kmc_steps = 1
    cfg.debug_target_trajectory = str(tmp_path / "target")
    Path(cfg.debug_target_trajectory).mkdir()
    for _ in range(100):
        state.inc_repeats()
    with pytest.raises(RuntimeError, match="target trajectory"):
        kmc_step(state, states, 0.0, states.kT, None, config=cfg)


def test_recycling_explorer_queues_a_saddle(tmp_path, monkeypatch):
    cfg, states, previous, current, _proc, _reactant = _states(tmp_path)
    cfg.recycling_on = True
    cfg.recycling_save_sugg = True
    cfg.disp_moved_only = True
    cfg.kdb_on = True
    cfg.kdb_only = True
    cfg.kdb_path = str(tmp_path / "kdb")
    cfg.kdb_scratch_path = str(tmp_path / "scratch")
    Path(cfg.kdb_path).mkdir()
    cfg.comm_type = "local_inprocess"
    cfg.comm_job_buffer_size = 1
    cfg.debug_keep_all_results = True
    cfg.debug_results_path = "debug-results"
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    comm = _Comm(
        [
            {
                "name": "1_4",
                "results.dat": _buf(_atoms(2.5)),
                "min.con": _buf(_atoms(2.5)),
            }
        ]
    )
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ClientMinModeExplorer(states, previous, current, config=cfg)
    registered = explorer.register_results()
    assert registered == 0
    assert (Path(cfg.path_root) / "debug-results" / "1_4" / "min.con").is_file()
    explorer.make_jobs()
    assert comm.submitted
    assert "displacement.con" in comm.submitted[0]
    cfg.recycling_on = False
    cfg.disp_moved_only = False
    before = len(comm.submitted)
    empty = ClientMinModeExplorer(states, previous, current, config=cfg)
    empty.make_jobs()
    assert len(comm.submitted) == before


def test_confident_state_cancels_and_superbasin_gate_runs(tmp_path, monkeypatch):
    cfg, states, state, _product, proc, _reactant = _states(tmp_path)
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_job_buffer_size = 1
    cfg.amsel_discover_decide = True
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    for _ in range(100):
        state.inc_repeats()
    explorer = ClientMinModeExplorer(states, state, state, config=cfg)
    explorer.explore()
    assert comm.submitted == []
    basin_root = tmp_path / "basins"
    basin_root.mkdir()
    basin = Superbasin(basin_root, 2, state_list=[state], config=cfg)
    mean_time, exit_state, product_state, exit_proc, basin_id = basin.step(
        state, states.get_product_state
    )
    assert mean_time > 0.0
    assert exit_state.number == 0
    assert product_state.number == 1
    assert exit_proc == proc
    assert basin_id == 2


def test_job_mains_reset_report_and_refuse(tmp_path, monkeypatch, capsys):
    root, ini = _ini(tmp_path, "basin_hopping")
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["hop", str(ini), "-q", "-n"])
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Comm())
    hop_main()
    assert (root / "wuid.dat").is_file()

    wide = tmp_path / "wide"
    wide_root, wide_ini = _ini(wide, "basin_hopping", "[Communicator]\njobs_per_bundle = 2\n")
    monkeypatch.chdir(wide_root)
    monkeypatch.setattr(sys, "argv", ["hop", str(wide_ini), "-q"])
    with pytest.raises(SystemExit) as caught:
        hop_main()
    assert caught.value.code == 1

    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["hop", str(ini), "-R", "-q"])
    monkeypatch.setattr("builtins.input", lambda prompt: "n")
    with pytest.raises(SystemExit) as caught:
        hop_main()
    assert caught.value.code == 1

    escape_root, escape_ini = _ini(tmp_path / "escape", "escape_rate")
    monkeypatch.chdir(escape_root)
    monkeypatch.setattr(sys, "argv", ["escape", str(escape_ini), "-q", "-n"])
    escape_main()
    assert (escape_root / "info.txt").is_file() or list(escape_root.glob("info*.txt"))

    replica_root, replica_ini = _ini(tmp_path / "replica", "parallel_replica")
    monkeypatch.chdir(replica_root)
    monkeypatch.setattr(sys, "argv", ["replica", str(replica_ini), "-q", "-n"])
    replica_main()
    assert (replica_root / "pos.con").is_file()

    movie_root, movie_ini = _ini(tmp_path / "movie", "akmc")
    monkeypatch.chdir(movie_root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(movie_ini), "-m", "processes,0,1", "-q"])
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 0
    assert "Saved" in capsys.readouterr().out

    monkeypatch.setattr(sys, "argv", ["akmc", str(movie_ini), "-m", "graph", "-s", "-q"])
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 2

    both = tmp_path / "both"
    both_root, both_ini = _ini(
        both,
        "akmc",
        "[Coarse Graining]\nuse_mcamc = true\nuse_askmc = true\n",
    )
    monkeypatch.chdir(both_root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(both_ini), "-q"])
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 1

    restart = tmp_path / "restart"
    restart_root, restart_ini = _ini(
        restart,
        "akmc",
        "[Coarse Graining]\nuse_mcamc = true\n",
    )
    (restart_root / "dynamics.txt").write_text("x\n")
    Path(restart_root / "superbasins").mkdir()
    monkeypatch.chdir(restart_root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(restart_ini), "-r", "-f", "-q"])
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 0
    assert not (restart_root / "dynamics.txt").exists()

    locked = tmp_path / "locked"
    locked_root, locked_ini = _ini(locked, "akmc")
    (locked_root / "lockfile").write_text("%i\n" % os.getpid())
    monkeypatch.chdir(locked_root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(locked_ini), "-q"])
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 1


def test_config_in_the_working_directory_can_be_declined(tmp_path, monkeypatch):
    root, ini = _ini(tmp_path, "akmc")
    other = tmp_path / "elsewhere"
    other.mkdir()
    monkeypatch.chdir(other)
    (other / "config.ini").write_text(ini.read_text())
    monkeypatch.setattr("builtins.input", lambda prompt: "n")
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init("")
    assert caught.value.code == 3
