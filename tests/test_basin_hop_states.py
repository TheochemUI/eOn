"""Basin hopping records a new minimum and can draw it back."""

import shutil
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.basinhopping import BHStates, register_results
from eon.config import ConfigClass
from eon.structure import Structure


def _config(directory):
    path = directory / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = basin_hopping",
                "temperature = 300",
                "[Paths]",
                f"main_directory = {directory}",
                f"results = {directory}",
                f"states = {directory / 'states'}",
                "",
            ]
        )
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    cfg.bh_initial_state_pool_size = 5
    return cfg


def _min():
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    atoms.r[0] = [1.0, 2.0, 3.0]
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    return atoms, buf


def test_new_minimum_can_be_drawn(tmp_path):
    atoms, buf = _min()
    states = BHStates(_config(tmp_path))
    assert states.get_random_minimum() is None
    states.add_state({"min.con": buf}, {"minimum_energy": -1.25})
    drawn = states.get_random_minimum()
    assert drawn is not None
    loaded = io.loadcon(drawn)
    assert loaded.names == ["Pt"]
    assert np.allclose(loaded.r, atoms.r)
    again = StringIO()
    io.savecon(again, atoms)
    again.seek(0)
    states.add_state({"min.con": again}, {"minimum_energy": -1.25})
    assert states.energy_table.rows[0]["repeats"] == 1


def test_add_state_accepts_an_existing_directory(tmp_path):
    atoms, buf = _min()
    states = BHStates(_config(tmp_path))
    state_dir = Path(tmp_path) / "states" / "0"
    state_dir.mkdir(parents=True)
    states.add_state({"min.con": buf}, {"minimum_energy": -1.25})
    assert (state_dir / "minimum.con").is_file()
    loaded = io.loadcon(str(state_dir / "minimum.con"))
    assert np.allclose(loaded.r, atoms.r)


class _Killed(BaseException):
    pass


def _kill_after_table_write(monkeypatch):
    write = io.Table.write

    def write_then_die(table):
        write(table)
        raise _Killed

    monkeypatch.setattr(io.Table, "write", write_then_die)


def _killed_pass(cfg, monkeypatch):
    _, buf = _min()
    with monkeypatch.context() as m:
        _kill_after_table_write(m)
        with pytest.raises(_Killed):
            BHStates(cfg).add_state({"min.con": buf}, {"minimum_energy": -1.25})


def test_draw_after_a_kill_at_the_row_write(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    _killed_pass(cfg, monkeypatch)
    states = BHStates(cfg)
    assert len(states.energy_table) == 1
    drawn = states.get_random_minimum()
    assert io.loadcon(drawn).names == ["Pt"]


def test_energy_match_after_a_kill_at_the_row_write(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    _killed_pass(cfg, monkeypatch)
    _, buf = _min()
    states = BHStates(cfg)
    assert not states.add_state({"min.con": buf}, {"minimum_energy": -1.25})
    assert states.energy_table.rows[0]["repeats"] == 1


def test_add_state_replaces_a_directory_without_a_row(tmp_path):
    atoms, buf = _min()
    state_dir = Path(tmp_path) / "states" / "0"
    state_dir.mkdir(parents=True)
    (state_dir / "minimum.con").write_text("torn\n")
    states = BHStates(_config(tmp_path))
    assert states.add_state({"min.con": buf}, {"minimum_energy": -1.25})
    loaded = io.loadcon(str(state_dir / "minimum.con"))
    assert np.allclose(loaded.r, atoms.r)
    assert [r["state"] for r in states.energy_table.rows] == [0]


class _EmptyComm:
    def get_results(self, path, keep):
        return iter(())


def test_register_results_accepts_an_existing_job_directory(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    jobs = Path(cfg.path_jobs_in)
    jobs.mkdir(parents=True)
    monkeypatch.setattr(shutil, "rmtree", lambda path: None)
    register_results(_EmptyComm(), BHStates(cfg), cfg)
    assert jobs.is_dir()
