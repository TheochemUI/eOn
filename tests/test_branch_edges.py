"""Branches the coverage table still left inside otherwise covered functions."""

import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.basinhopping import main as hop_main
from eon.escaperate import get_pr_metadata, main as escape_main
from eon.explorer import ProcessSearch
from eon.geometry.neighbors import neighbor_list_pairs
from eon.parallelreplica import main as replica_main
from tests.test_catalog_and_mains import _ini
from tests.test_library_bodies import _atoms, _buf, _config, _mode, _states


def _search(tmp_path, reactant):
    tmp_path.mkdir(parents=True, exist_ok=True)
    cfg = _config(tmp_path)
    cfg.process_search_minimization_offset = 0.1
    Path(cfg.path_incomplete).mkdir(parents=True, exist_ok=True)
    ini = Path(cfg.config_path)
    ini.write_text(ini.read_text() + "\n[Saddle Search]\nmethod = min_mode\n")
    displacement = reactant.copy()
    displacement.r[1, 0] += 0.2
    mode = np.zeros_like(reactant.r)
    mode[1, 0] = 1.0
    return cfg, ProcessSearch(reactant, displacement, mode, "random", 1, 0, config=cfg)


def _minimized(name, atoms_, reason):
    return {
        "name": name,
        "results": {
            "job_type": "minimization",
            "termination_reason": reason,
            "potential_energy": -1.0,
            "total_force_calls": 2,
        },
        "min.con": _buf(atoms_),
    }


def test_unconnected_minima_and_a_minimization_that_does_not_finish(tmp_path):
    reactant = _atoms(2.5)
    cfg, search = _search(tmp_path, reactant)
    search.process_result(
        {
            "name": "0_1",
            "results": {
                "job_type": "saddle_search",
                "termination_reason": 0,
                "potential_energy_saddle": -0.7,
                "potential_energy_reactant": -1.0,
                "barrier_reactant_to_product": 0.3,
                "total_force_calls": 1,
            },
            "saddle.con": _buf(_atoms(2.7)),
            "mode.dat": _mode(_atoms(2.7)),
        }
    )
    search.get_job(0)
    other = _atoms(4.0)
    search.process_result(_minimized("0_2", other, 1))
    search.get_job(0)
    stalled = search.process_result(_minimized("0_3", other, 1))
    assert stalled["results"]["termination_reason"] == 9

    cfg2, search2 = _search(tmp_path / "open", reactant)
    search2.process_result(
        {
            "name": "0_1",
            "results": {
                "job_type": "saddle_search",
                "termination_reason": 0,
                "potential_energy_saddle": -0.7,
                "potential_energy_reactant": -1.0,
                "barrier_reactant_to_product": 0.3,
                "total_force_calls": 1,
            },
            "saddle.con": _buf(_atoms(2.7)),
            "mode.dat": _mode(_atoms(2.7)),
        }
    )
    search2.get_job(0)
    search2.process_result(_minimized("0_2", other, 0))
    search2.get_job(0)
    disconnected = search2.process_result(_minimized("0_3", other, 0))
    assert disconnected["results"]["termination_reason"] == 6


def test_hop_reset_removes_the_state_table(tmp_path, monkeypatch, capsys):
    root, ini = _ini(tmp_path, "basin_hopping")
    (root / "states").mkdir()
    (root / "wuid.dat").write_text("4\n")
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["hop", str(ini), "-R"])
    monkeypatch.setattr("builtins.input", lambda prompt: "y")
    with pytest.raises(SystemExit) as caught:
        hop_main()
    assert caught.value.code == 0
    assert not (root / "wuid.dat").exists()
    assert "Reset." in capsys.readouterr().out

    quiet, quiet_ini = _ini(tmp_path / "talk", "basin_hopping")
    monkeypatch.chdir(quiet)
    monkeypatch.setattr(sys, "argv", ["hop", str(quiet_ini)])
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: type("C", (), {
        "get_results": lambda self, path, keep: iter(()),
        "queued_search_count": lambda self: 0,
        "submit_jobs": lambda self, searches, invariants: None,
    })())
    hop_main()
    assert (quiet / "wuid.dat").is_file()


def test_partial_replica_metadata_and_a_spoken_escape(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    info = Path(cfg.path_results) / "info.txt"
    info.write_text("[Simulation Information]\ncurrent_state = later\n")
    assert get_pr_metadata(cfg) == (0, 0.0, 0)

    root, ini = _ini(tmp_path / "escape", "escape_rate")
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["escape", str(ini), "-n"])
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: type("C", (), {
        "get_results": lambda self, path, keep: iter(()),
        "queued_search_count": lambda self: 0,
        "submit_jobs": lambda self, searches, invariants: None,
        "cancel_state": lambda self, number: 0,
    })())
    escape_main()
    assert (root / "info.txt").is_file()


def test_retune_rates_when_the_temperature_changes(tmp_path):
    cfg, _states_obj, state, _product, _proc, _reactant = _states(tmp_path)
    cfg.akmc_eq_rate = 1.0
    state.info.set("MetaData", "kT", 1.0)
    state.load_process_table(force=True)
    rate = next(iter(state.get_process_table().values()))["rate"]
    assert rate > 0.0


def test_nonperiodic_neighbor_pairs_and_a_failed_atomic_write(tmp_path):
    cluster = _atoms(2.5)
    pairs = neighbor_list_pairs(cluster, 3.5, periodic=False)
    assert pairs[0].size >= 1
    path = tmp_path / "half.txt"
    with pytest.raises(RuntimeError, match="stop"):
        with io.atomic_write(str(path)) as handle:
            handle.write("partial\n")
            raise RuntimeError("stop")
    assert not path.exists()
