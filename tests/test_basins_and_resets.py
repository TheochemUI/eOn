"""Energy-level basins, reset entry points, and internal motion."""

import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from eon import atoms
from eon.akmc import get_superbasin_scheme
from eon.communicator_inprocess import _params_from_invariants, _require_pyeonclient
from eon.config import ConfigClass
from eon.escaperate import main as escape_main
from eon.explorer import ClientMinModeExplorer
from eon.parallelreplica import main as replica_main
from eon.recycling import Recycling
from eon.structure import Structure
from eon.superbasinscheme import EnergyLevel, RateThreshold
from tests.test_catalog_and_mains import _ini
from tests.test_library_bodies import _Comm, _atoms, _config, _states


def _trio():
    atoms_ = Structure(3)
    atoms_.names = ["Cu", "Cu", "Cu"]
    atoms_.mass[:] = 63.5
    atoms_.box = np.diag([12.0, 12.0, 12.0])
    atoms_.r[1] = [2.2, 0.0, 0.0]
    atoms_.r[2] = [1.1, 1.9, 0.0]
    return atoms_


def test_energy_level_merges_when_the_level_passes_the_saddle(tmp_path):
    cfg, states, start, end, _proc, _reactant = _states(tmp_path)
    cfg.comp_use_identical = False
    cfg.sb_max_size = 0
    root = tmp_path / "levels"
    root.mkdir()
    (root / "storage").mkdir()
    scheme = EnergyLevel(str(root), states, states.kT, 10.0, config=cfg)
    scheme.register_transition(start, start)
    scheme.register_transition(start, end)
    basin = scheme.get_containing_superbasin(start)
    assert basin is not None
    assert set(basin.state_numbers) == {0, 1}
    assert scheme.get_statelist(start)[0].number in (0, 1)
    assert scheme.levels[end] > end.get_energy()

    rates = tmp_path / "rates"
    rates.mkdir()
    (rates / "storage").mkdir()
    threshold = RateThreshold(str(rates), states, states.kT, 0.0, config=cfg)
    threshold.register_transition(start, end)
    assert threshold.get_containing_superbasin(start) is not None
    with pytest.raises(ValueError, match="Unknown superbasin"):
        get_superbasin_scheme(states, _bad_scheme(cfg))


def _bad_scheme(cfg):
    cfg.sb_scheme = "not-a-scheme"
    return cfg


def test_resets_discard_the_run_directory(tmp_path, monkeypatch):
    root, ini = _ini(tmp_path / "escape", "escape_rate")
    (root / "dynamics.txt").write_text("x\n")
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["escape", str(ini), "-R"])
    monkeypatch.setattr("builtins.input", lambda prompt: "y")
    with pytest.raises(SystemExit) as caught:
        escape_main()
    assert caught.value.code == 0
    assert not (root / "dynamics.txt").exists()

    replica, replica_ini = _ini(tmp_path / "replica", "parallel_replica")
    (replica / "dynamics.txt").write_text("x\n")
    monkeypatch.chdir(replica)
    monkeypatch.setattr(sys, "argv", ["replica", str(replica_ini), "-R"])
    monkeypatch.setattr("builtins.input", lambda prompt: "n")
    with pytest.raises(SystemExit) as caught:
        replica_main()
    assert caught.value.code == 0
    assert (replica / "dynamics.txt").is_file()

    wide, wide_ini = _ini(
        tmp_path / "wide",
        "parallel_replica",
        "[Communicator]\njobs_per_bundle = 2\n",
    )
    monkeypatch.chdir(wide)
    monkeypatch.setattr(sys, "argv", ["replica", str(wide_ini), "-q"])
    with pytest.raises(SystemExit) as caught:
        replica_main()
    assert caught.value.code == 1


def test_dynamics_recycling_weights_the_moved_atoms(tmp_path, monkeypatch):
    cfg, states, previous, current, _proc, _reactant = _states(tmp_path)
    cfg.saddle_method = "dynamics"
    cfg.recycling_on = True
    cfg.disp_moved_only = True
    cfg.kdb_on = False
    cfg.comm_type = "local"
    cfg.comm_job_buffer_size = 1
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "pos.con").write_text("placeholder\n")
    ini = Path(cfg.config_path)
    ini.write_text(ini.read_text() + "\n[Saddle Search]\nmethod = dynamics\n")
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ClientMinModeExplorer(states, previous, current, config=cfg)
    explorer.make_jobs()
    assert comm.submitted
    meta = Path(current.path) / "recycling_info"
    assert meta.is_file()
    again = Recycling(
        states,
        previous,
        current,
        cfg.recycling_move_distance,
        False,
        config=cfg,
    )
    assert again.process_number >= 0


def test_internal_motion_and_inprocess_parameters(tmp_path):
    left = _trio()
    right = left.copy()
    right.r[1] = [0.0, 2.2, 0.0]
    right.r[2] = [-1.9, 1.1, 0.0]
    moved = atoms.internal_motion(left, right)
    assert len(moved) == 3
    assert np.all(np.isfinite(moved.r))
    pc = _require_pyeonclient()
    params = _params_from_invariants(
        pc, {"config.ini": (StringIO("[Main]\njob = akmc\n"), 0o644)}
    )
    assert params is not None
    bare = _params_from_invariants(pc, {})
    assert bare.quiet is True
    cfg = _config(tmp_path)
    from eon import fileio as io

    nested = tmp_path / "a" / "b"
    nested.mkdir(parents=True)
    (nested / "c").mkdir()
    io.remove_tree_and_empty_parents(nested / "c")
    assert not (tmp_path / "a").exists() or (tmp_path / "config.ini").is_file()
