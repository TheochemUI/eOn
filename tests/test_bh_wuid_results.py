"""Basin hopping stores wuid.dat with the results, and a killed rewrite keeps it."""

import os
from pathlib import Path

import numpy as np

from eon import fileio as io
from eon.basinhopping import basinhopping
from eon.config import ConfigClass
from eon.structure import Structure


class _Comm:
    def get_results(self, path, keep):
        return iter(())

    def queued_search_count(self):
        return 0

    def submit_jobs(self, searches, invariants):
        return None


def _config(root, results):
    path = root / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = basin_hopping",
                "temperature = 300",
                "[Paths]",
                f"main_directory = {root}",
                f"results = {results}",
                f"states = {root / 'states'}",
                f"jobs_in = {root / 'jobs_in'}",
                f"scratch = {root / 'scratch'}",
                "",
            ]
        )
    )
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    io.savecon(str(root / "pos.con"), atoms)
    cfg = ConfigClass()
    cfg.init(str(path))
    cfg.bh_initial_state_pool_size = 0
    cfg.comm_job_buffer_size = 1
    return cfg


def test_wuid_is_stored_under_results(tmp_path, monkeypatch):
    root = tmp_path / "root"
    results = tmp_path / "results"
    root.mkdir()
    results.mkdir()
    cfg = _config(root, results)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Comm())
    basinhopping(cfg)
    assert (results / "wuid.dat").read_text().strip() == "1"
    assert not (tmp_path / "wuid.dat").exists()


def test_killed_wuid_rewrite_keeps_the_previous_id(tmp_path, monkeypatch):
    root = tmp_path / "root"
    results = tmp_path / "results"
    root.mkdir()
    results.mkdir()
    (results / "wuid.dat").write_text("7\n")
    cfg = _config(root, results)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Comm())

    def killed(src, dst):
        raise OSError("killed")

    monkeypatch.setattr(os, "replace", killed)
    try:
        basinhopping(cfg)
    except OSError as exc:
        assert "killed" in str(exc)
    else:
        raise AssertionError("wuid rewrite completed without os.replace")
    assert (results / "wuid.dat").read_text().strip() == "7"
    assert not (tmp_path / "wuid.dat").exists()
