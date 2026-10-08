"""One aKMC pass with no finished search returns zero steps."""

from pathlib import Path

import numpy as np

from eon import fileio as io
from eon.akmc import akmc
from eon.config import ConfigClass
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

    def cancel_state(self, number):
        return 0


def test_one_pass_returns_zero_steps(tmp_path, monkeypatch):
    root = tmp_path / "run"
    root.mkdir()
    atoms = Structure(2)
    atoms.names = ["Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[1] = [2.5, 0.0, 0.0]
    io.savecon(str(root / "pos.con"), atoms)
    ini = root / "config.ini"
    ini.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "[Paths]",
                f"main_directory = {root}",
                f"results = {root}",
                f"states = {root / 'states'}",
                "",
            ]
        )
    )
    cfg = ConfigClass()
    cfg.init(str(ini))
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.sb_on = False
    cfg.comm_job_buffer_size = 1
    cfg.comm_type = "local_inprocess"
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    steps = akmc(cfg)
    assert steps == 0
    assert comm.submitted
    assert "displacement.con" in comm.submitted[0]
