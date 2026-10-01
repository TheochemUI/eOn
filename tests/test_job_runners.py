"""Unit coverage for the basin-hopping and parallel-replica runners."""
from __future__ import annotations

import configparser
from io import StringIO
from types import SimpleNamespace

import numpy as np
import pytest

from eon import atoms
from eon import fileio as io
from eon.basinhopping import BHStates, make_searches, register_results
from eon.parallelreplica import (
    get_pr_metadata,
    make_searches as pr_make_searches,
    step,
    write_pr_metadata,
)
from eon.structure import Structure


def _structure(pos):
    s = Structure(1)
    s.r[0] = np.asarray(pos, dtype=float)
    s.free[0] = 1.0
    s.names[0] = "H"
    s.mass[0] = 1.0
    s.box = np.eye(3) * 20.0
    return s


def _con_buf(structure):
    buf = StringIO()
    io.savecon(buf, structure)
    buf.seek(0)
    return buf


def _bh_config(tmp_path):
    states = tmp_path / "states"
    root = tmp_path / "root"
    root.mkdir()
    (root / "pos.con").write_text("reactant\n")
    ini = tmp_path / "config.ini"
    ini.write_text("[Main]\nrandom_seed = 1\n")
    return SimpleNamespace(
        path_states=str(states),
        path_root=str(root),
        path_jobs_in=str(tmp_path / "jobs_in"),
        config_path=str(ini),
        comp_eps_e=1.0e-6,
        comp_eps_r=0.1,
        comp_neighbor_cutoff=3.0,
        comp_check_rotation=False,
        comp_use_identical=False,
        comp_remove_translation=False,
        bh_initial_state_pool_size=1,
        comm_job_buffer_size=2,
    )


def _result(structure, energy, reason=0):
    return {
        "min.con": _con_buf(structure),
        "results.dat": StringIO(
            f"{reason} termination_reason\n{energy} minimum_energy\n"
        ),
    }


class _Comm:
    def __init__(self, results=(), queued=0):
        self._results = list(results)
        self.queued = queued
        self.submitted = []

    def get_results(self, path, keep_result):
        assert keep_result("anything") is True
        return self._results

    def queued_search_count(self):
        return self.queued

    def submit_jobs(self, searches, invariants):
        self.submitted.extend(searches)


def test_bh_states_add_repeat_and_pool(tmp_path):
    cfg = _bh_config(tmp_path)
    states = BHStates(cfg)
    first = _structure([0.0, 0.0, 0.0])
    assert states.get_random_minimum() is None
    info = {"minimum_energy": -1.5, "termination_reason": 0}
    assert states.add_state(_result(first, -1.5), info) is True
    assert states.add_state(_result(first, -1.5), info) is False
    drawn = states.get_random_minimum()
    assert drawn is not None
    loaded = io.loadcon(drawn)
    assert atoms.match(
        first, loaded, cfg.comp_eps_r, cfg.comp_neighbor_cutoff, False,
        check_rotation=False, use_identical=False,
    )
    other = _structure([4.0, 0.0, 0.0])
    info2 = {"minimum_energy": -0.2, "termination_reason": 0}
    assert states.add_state(_result(other, -0.2), info2) is True
    assert len(states.energy_table) == 2


def test_bh_register_and_make_searches(tmp_path):
    cfg = _bh_config(tmp_path)
    states = BHStates(cfg)
    struct = _structure([1.0, 2.0, 3.0])
    good = _result(struct, -2.0, reason=0)
    skipped = {
        "min.con": _con_buf(struct),
        "results.dat": StringIO("1 termination_reason\n"),
    }
    comm = _Comm([good, skipped])
    register_results(comm, states, cfg)
    assert len(states.energy_table) == 1

    comm.queued = 5
    cfg.comm_job_buffer_size = 5
    assert make_searches(comm, 4, states, cfg) == 4

    comm.queued = 0
    cfg.comm_job_buffer_size = 2
    cfg.bh_initial_state_pool_size = 0
    assert make_searches(comm, 4, states, cfg) == 6
    assert len(comm.submitted) == 2
    assert comm.submitted[0]["id"] == "4"

    cfg.bh_initial_state_pool_size = -1
    with pytest.raises(SystemExit):
        make_searches(comm, 9, states, cfg)


def test_pr_metadata_roundtrip_and_step(tmp_path):
    results = tmp_path / "results"
    cfg = SimpleNamespace(
        path_results=str(results),
        path_root=str(tmp_path),
        main_job="parallel_replica",
        config_path=str(tmp_path / "config.ini"),
        comm_job_buffer_size=1,
    )
    (tmp_path / "config.ini").write_text("[Main]\njob = parallel_replica\n")
    assert get_pr_metadata(cfg) == (0, 0, 0)

    parser = configparser.RawConfigParser()
    write_pr_metadata(parser, 3, 1.5e-9, 8)
    io.write_info_txt(cfg, parser)
    assert get_pr_metadata(cfg) == (3, pytest.approx(1.5e-9), 8)

    # Missing keys fall back to the zero state.
    (results / "info.txt").write_text("[Simulation Information]\n")
    assert get_pr_metadata(cfg) == (0, 0.0, 0)

    class _State:
        def __init__(self, number):
            self.number = number
            self.energy = -3.0

        def zero_time(self):
            self.time = 0.0

        def get_time(self):
            return 0.0

        def get_energy(self):
            return self.energy

        def get_process(self, process_id):
            return process_id

        def get_reactant(self):
            return _structure([0.0, 0.0, 0.0])

    class _States:
        def get_product_state(self, number, process_id):
            assert number == 1
            assert process_id == 4
            return _State(2)

    current, previous = step(
        10.0, _State(1), _States(), {"process_id": 4, "time": 2.5}, cfg,
    )
    assert current.number == 2
    assert previous.number == 1
    assert (results / "dynamics.txt").is_file()

    comm = _Comm(queued=3)
    cfg.comm_job_buffer_size = 1
    assert pr_make_searches(comm, _State(1), 5, cfg) == 5
    comm.queued = 0
    cfg.comm_job_buffer_size = 2
    assert pr_make_searches(comm, _State(1), 5, cfg) == 7
    assert comm.submitted[-1]["id"].startswith("1_")


def test_runner_helpers_require_config():
    with pytest.raises(TypeError):
        get_pr_metadata()
    with pytest.raises(TypeError):
        step(0.0, None, None, {})
    with pytest.raises(TypeError):
        pr_make_searches(None, None, 0)
