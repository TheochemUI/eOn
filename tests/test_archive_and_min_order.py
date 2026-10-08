"""Debug archives, the swapped minimization order, and a Lambert value."""

import math
from pathlib import Path

import numpy as np
import pytest

from eon import atoms
from eon.akmcstate import lambertw
from eon.communicator import get_communicator
from eon.explorer import ProcessSearch, ServerMinModeExplorer, _archive_debug_result
from tests.test_library_bodies import _Comm, _atoms, _buf, _config, _mode, _states


class _Blob:
    def __init__(self, text):
        self.text = text
        self.pos = 0

    def tell(self):
        return self.pos

    def seek(self, pos):
        self.pos = pos

    def read(self):
        return self.text


def test_archive_reads_a_stream_and_rejects_a_missing_communicator(tmp_path):
    cfg = _config(tmp_path)
    cfg.debug_keep_all_results = True
    cfg.debug_results_path = "kept"
    _archive_debug_result(
        cfg,
        {
            "name": "stream",
            "results.dat": _Blob("0 termination_reason\n"),
            "min.con": _Blob("Pt\n"),
        },
    )
    kept = Path(cfg.path_root) / "kept" / "stream"
    assert (kept / "results.dat").read_text().startswith("0")
    assert (kept / "min.con").read_text() == "Pt\n"
    assert "not-an-element" not in atoms.elements
    assert "Cu" in atoms.elements
    with pytest.raises(TypeError):
        get_communicator(None)
    cfg.comm_type = "mpi"
    with pytest.raises(ModuleNotFoundError):
        get_communicator(cfg)
    near = -0.3678794411714423 + 1.0e-5
    assert math.isfinite(lambertw(near))


def test_second_minimum_can_be_the_reactant(tmp_path, monkeypatch):
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
    ini = Path(cfg.config_path)
    ini.write_text(ini.read_text() + "\n[Saddle Search]\nmethod = min_mode\n")
    comm = _Comm()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ServerMinModeExplorer(states, state, state, config=cfg)
    explorer.make_jobs()
    wuid = next(iter(explorer.wuid_to_search_id))
    search = next(iter(explorer.process_searches.values()))
    assert isinstance(search, ProcessSearch)
    search.process_result(
        {
            "name": "0_%s" % wuid,
            "results": {
                "job_type": "saddle_search",
                "termination_reason": 0,
                "potential_energy_saddle": -0.7,
                "potential_energy_reactant": -1.0,
                "barrier_reactant_to_product": 0.3,
                "total_force_calls": 1,
            },
            "saddle.con": _buf(_atoms(2.7)),
            "mode.dat": _mode(reactant),
        }
    )
    before = len(comm.submitted)
    explorer.make_jobs()
    assert len(comm.submitted) > before
    product = _atoms(3.4)
    search.process_result(
        {
            "name": "0_b",
            "results": {
                "job_type": "minimization",
                "termination_reason": 0,
                "potential_energy": -0.8,
                "total_force_calls": 2,
            },
            "min.con": _buf(product),
        }
    )
    _job, kind = search.get_job(0)
    assert kind == "min2"
    finished = search.process_result(
        {
            "name": "0_c",
            "results": {
                "job_type": "minimization",
                "termination_reason": 0,
                "potential_energy": -1.0,
                "total_force_calls": 2,
            },
            "min.con": _buf(reactant),
        }
    )
    assert finished["results"]["potential_energy_reactant"] == pytest.approx(-1.0)
    assert np.isfinite(finished["results"]["barrier_reactant_to_product"])
