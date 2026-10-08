"""Reset answers, a process that returns to its reactant, and foreign formats."""

import sys
from pathlib import Path

import pytest

from eon import fileio as io
from eon.akmc import main as akmc_main
from eon.config import ConfigClass
from eon.parallelreplica import main as replica_main
from tests.test_catalog_and_mains import _ini
from tests.test_library_bodies import _atoms, _states


def test_replica_reset_removes_dynamics_and_speaks(tmp_path, monkeypatch, capsys):
    root, ini = _ini(tmp_path, "parallel_replica")
    (root / "dynamics.txt").write_text("old\n")
    (root / "states").mkdir()
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["replica", str(ini), "-R"])
    monkeypatch.setattr("builtins.input", lambda prompt: "y")
    with pytest.raises(SystemExit) as caught:
        replica_main()
    assert caught.value.code == 0
    assert not (root / "dynamics.txt").exists()
    assert "Reset" in capsys.readouterr().out

    spoken, spoken_ini = _ini(tmp_path / "spoken", "parallel_replica")
    monkeypatch.chdir(spoken)
    monkeypatch.setattr(sys, "argv", ["replica", str(spoken_ini), "-n"])
    monkeypatch.setattr(
        "eon.communicator.get_communicator",
        lambda config: type(
            "C",
            (),
            {
                "get_results": lambda self, path, keep: iter(()),
                "queued_search_count": lambda self: 0,
                "submit_jobs": lambda self, searches, invariants: None,
                "cancel_state": lambda self, number: 0,
            },
        )(),
    )
    replica_main()
    assert (spoken / "info.txt").is_file()


def test_restart_also_removes_superbasin_files(tmp_path, monkeypatch, capsys):
    root, ini = _ini(
        tmp_path,
        "akmc",
        "[Coarse Graining]\nuse_mcamc = true\n",
    )
    cfg = ConfigClass()
    cfg.init(str(ini))
    state_dir = Path(cfg.path_states) / "0"
    state_dir.mkdir(parents=True)
    marker = state_dir / cfg.sb_state_file
    marker.write_text("1 1\n")
    Path(cfg.sb_path).mkdir(parents=True, exist_ok=True)
    (root / "dynamics.txt").write_text("old\n")
    answers = iter(["y", "y"])
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", ["akmc", str(ini), "-r"])
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))
    with pytest.raises(SystemExit) as caught:
        akmc_main()
    assert caught.value.code == 0
    assert not marker.exists()
    assert "removed" in capsys.readouterr().out


def test_a_process_that_returns_to_its_reactant(tmp_path):
    _cfg, states, state, _product, _proc, reactant = _states(tmp_path)
    extra = state.allocate_process_id(b"return", b"self")
    forward = state.get_process(_proc)
    state.append_process_table(
        id=extra,
        saddle_energy=forward["saddle_energy"],
        prefactor=forward["prefactor"],
        product=-1,
        product_energy=state.get_energy(),
        product_prefactor=forward["product_prefactor"],
        barrier=forward["barrier"],
        rate=forward["rate"],
        repeats=0,
    )
    io.savecon(state.proc_product_path(extra), reactant)
    found = states.get_product_state(0, extra)
    assert found.number == 0


def test_foreign_xyz_without_chemfiles_is_refused(tmp_path):
    path = tmp_path / "pair.xyz"
    path.write_text("2\n\nPt 0 0 0\nPt 2.5 0 0\n")
    try:
        loaded = io.loadxyz(str(path))
    except RuntimeError as exc:
        assert "chemfiles" in str(exc)
    else:
        assert len(loaded) == 2
