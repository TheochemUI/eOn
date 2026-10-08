"""Basin hopping, the catalog, escape-rate searches, and a reloaded table."""

from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from eon import atoms
from eon import fileio as io
from eon.basinhopping import basinhopping, make_searches
from eon.config import ConfigClass
from eon.escaperate import make_searches as pr_make_searches
from eon.process_catalog import make_suggestion, query
from eon.server import _warn_pos_con_in_potfiles, select_job_runner
from eon.structure import Structure


class _Comm:
    def __init__(self):
        self.submitted = []

    def get_results(self, path, keep):
        return iter(())

    def queued_search_count(self):
        return 0

    def submit_jobs(self, searches, invariants):
        self.submitted.extend(searches)


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
    return cfg


def test_missing_hop_root_raises(tmp_path):
    cfg = SimpleNamespace(path_root=str(tmp_path / "absent"))
    try:
        basinhopping(cfg)
    except FileNotFoundError as exc:
        assert "root directory does not exist" in str(exc)
        return
    raise AssertionError("basinhopping returned without raising")


def test_negative_pool_raises(tmp_path):
    cfg = _config(tmp_path)
    (tmp_path / "pos.con").write_text("placeholder\n")
    cfg.bh_initial_state_pool_size = -1
    cfg.comm_job_buffer_size = 1
    try:
        make_searches(_Comm(), 0, SimpleNamespace(), cfg)
    except ValueError as exc:
        assert "negative" in str(exc)
        return
    raise AssertionError("make_searches returned without raising")


def test_hop_pass_queues_the_reactant(tmp_path, monkeypatch):
    root = tmp_path / "run"
    root.mkdir()
    atoms_ = Structure(1)
    atoms_.names = ["Pt"]
    atoms_.mass[:] = 195.084
    atoms_.box = np.diag([12.0, 12.0, 12.0])
    io.savecon(str(root / "pos.con"), atoms_)
    cfg = _config(root)
    cfg.bh_initial_state_pool_size = 0
    cfg.comm_job_buffer_size = 1
    cfg.recycling_on = False
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    comm = _Comm()
    monkeypatch.chdir(root)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    basinhopping(cfg)
    assert comm.submitted
    assert comm.submitted[0]["pos.con"].getvalue()
    assert Path("wuid.dat").read_text().strip() == "1"


def test_catalog_without_a_store_returns_nothing(tmp_path):
    cfg = _config(tmp_path)
    cfg.kdb_path = str(tmp_path / "kdb")
    Path(cfg.kdb_path).mkdir()
    state = SimpleNamespace(number=0, reactant_path=str(tmp_path / "reactant.con"))
    assert query(state, cfg) is False
    assert make_suggestion(cfg, state) == (None, None)


def test_escape_searches_use_the_reactant(tmp_path):
    atoms_ = Structure(1)
    atoms_.names = ["Cu"]
    atoms_.mass[:] = 63.5
    atoms_.box = np.diag([10.0, 10.0, 10.0])
    state = SimpleNamespace(number=4, get_reactant=lambda: atoms_)
    cfg = _config(tmp_path)
    cfg.comm_job_buffer_size = 2
    comm = _Comm()
    wuid = pr_make_searches(comm, state, 8, cfg)
    assert wuid == 10
    assert [item["id"] for item in comm.submitted] == ["4_8", "4_9"]


def test_runner_names_and_potfile_warning(tmp_path, capsys):
    assert select_job_runner("akmc").__name__ == "main"
    assert select_job_runner("basin_hopping").__name__ == "main"
    assert select_job_runner("not-a-job") is None
    pot = tmp_path / "pot"
    pot.mkdir()
    (pot / "pos.con").write_text("x\n")
    _warn_pos_con_in_potfiles(SimpleNamespace(path_pot=str(pot)))
    assert "pos.con" in capsys.readouterr().out


def test_identical_clusters_and_a_reloaded_table(tmp_path):
    left = Structure(2)
    left.names = ["Pt", "Pt"]
    left.box = np.diag([15.0, 15.0, 15.0])
    left.r[1] = [2.5, 0.0, 0.0]
    right = left.copy()
    assert atoms.identical(left, right, 0.1)
    right.r[1, 0] = 4.0
    assert atoms.identical(left, right, 0.1) is False
    assert atoms.symbol_for_z(atoms.atomic_number("Cu")) == "Cu"
    path = tmp_path / "rates.tbl"
    table = io.Table(str(path), ["state", "rate"])
    table.add_row({"state": 1, "rate": 2.5})
    again = io.Table(str(path), ["state", "rate"])
    assert len(again) == 1
    assert again.rows[0]["rate"] == 2.5
