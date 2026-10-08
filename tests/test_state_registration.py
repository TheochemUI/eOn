"""Registering a saddle builds a product state and a finite rate."""

from io import StringIO

import numpy as np
import pytest

from eon import fileio as io
from eon.akmcstatelist import AKMCStateList
from eon.config import ConfigClass
from eon.mcamc.mcamc import estimate_condition, guess_precision, np_mcamc
from eon.structure import Structure


def _atoms(x1):
    atoms = Structure(2)
    atoms.names = ["Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[0] = [0.0, 0.0, 0.0]
    atoms.r[1] = [x1, 0.0, 0.0]
    return atoms


def _buf(atoms):
    buf = StringIO()
    io.savecon(buf, atoms)
    buf.seek(0)
    return buf


def _mode(atoms):
    mode = np.zeros_like(atoms.r)
    mode[1, 0] = 0.1
    buf = StringIO()
    io.save_mode(buf, mode)
    buf.seek(0)
    return buf


def _config(directory):
    path = directory / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
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


def _result(reactant, saddle, product, energy_saddle, wuid=1):
    results = StringIO()
    io.save_results_dat(
        results,
        {
            "potential_energy_saddle": energy_saddle,
            "potential_energy_reactant": -1.0,
            "potential_energy_product": -0.8,
            "prefactor_reactant_to_product": 1.0e12,
            "prefactor_product_to_reactant": 1.0e12,
            "barrier_reactant_to_product": energy_saddle - (-1.0),
            "displacement_saddle_distance": 0.2,
            "force_calls_saddle": 4,
            "force_calls_minimization": 3,
            "force_calls_prefactors": 1,
        },
    )
    results.seek(0)
    return {
        "wuid": wuid,
        "type": "random",
        "results": {
            "potential_energy_saddle": energy_saddle,
            "potential_energy_reactant": -1.0,
            "potential_energy_product": -0.8,
            "prefactor_reactant_to_product": 1.0e12,
            "prefactor_product_to_reactant": 1.0e12,
            "barrier_reactant_to_product": energy_saddle - (-1.0),
            "displacement_saddle_distance": 0.2,
            "force_calls_saddle": 4,
            "force_calls_minimization": 3,
            "force_calls_prefactors": 1,
        },
        "reactant.con": _buf(reactant),
        "saddle.con": _buf(saddle),
        "product.con": _buf(product),
        "mode.dat": _mode(saddle),
        "results.dat": results,
    }


def test_new_saddle_creates_a_product_state(tmp_path):
    reactant = _atoms(2.5)
    saddle = _atoms(2.7)
    product = _atoms(3.1)
    con = tmp_path / "reactant.con"
    io.savecon(str(con), reactant)
    cfg = _config(tmp_path)
    kT = 300.0 / 11604.5
    states = AKMCStateList(kT, 20.0, 40.0, initial_state=str(con), config=cfg)
    state = states.get_state(0)
    proc_id = state.add_process(_result(reactant, saddle, product, -0.7))
    assert proc_id is not None
    assert state.get_good_saddle_count() == 1
    assert state.get_unique_saddle_count() == 1
    table = state.get_ratetable()
    assert table[0][0] == proc_id
    assert table[0][1] > 0.0
    product_state = states.get_product_state(0, proc_id)
    assert product_state.number == 1
    assert states.get_num_states() == 2
    assert state.get_process(proc_id)["product"] == 1
    reverse = product_state.get_process_table()
    assert reverse
    assert next(iter(reverse.values()))["product"] == 0
    assert 0.0 <= state.get_confidence() <= 1.0


def test_repeated_saddle_does_not_add_a_process(tmp_path):
    reactant = _atoms(2.5)
    saddle = _atoms(2.7)
    product = _atoms(3.1)
    con = tmp_path / "reactant.con"
    io.savecon(str(con), reactant)
    cfg = _config(tmp_path)
    states = AKMCStateList(300.0 / 11604.5, 20.0, 40.0, initial_state=str(con), config=cfg)
    state = states.get_state(0)
    first = state.add_process(_result(reactant, saddle, product, -0.7, wuid=1))
    second = state.add_process(_result(reactant, saddle, product, -0.7, wuid=2))
    assert first is not None
    assert second is None
    assert state.get_process(first)["repeats"] == 1


class _Comm:
    def __init__(self, results):
        self.results = results
        self.submitted = []

    def get_results(self, path, keep):
        for result in self.results:
            if keep(result["name"]):
                yield result

    def queued_search_count(self):
        return 0

    def submit_jobs(self, searches, invariants):
        self.submitted.extend(searches)

    def cancel_state(self, number):
        return 0


def test_explorer_records_a_bad_saddle_and_queues_a_search(tmp_path, monkeypatch):
    from eon.explorer import ClientMinModeExplorer

    reactant = _atoms(2.5)
    con = tmp_path / "reactant.con"
    io.savecon(str(con), reactant)
    cfg = _config(tmp_path)
    cfg.comm_job_buffer_size = 1
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_type = "local_inprocess"
    from pathlib import Path
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    states = AKMCStateList(
        300.0 / 11604.5, 20.0, 40.0, initial_state=str(con), config=cfg
    )
    state = states.get_state(0)
    comm = _Comm(
        [
            {
                "name": "0_7",
                "results": {
                    "termination_reason": 1,
                    "barrier_reactant_to_product": 0.0,
                    "displacement_saddle_distance": 0.0,
                    "force_calls_saddle": 1,
                    "force_calls_minimization": 0,
                    "force_calls_prefactors": 0,
                },
            }
        ]
    )
    monkeypatch.setattr(
        "eon.communicator.get_communicator", lambda config: comm
    )
    explorer = ClientMinModeExplorer(states, None, state, config=cfg)
    explorer.job_table.add_row({"state": 0, "wuid": 7, "type": "random"})
    registered = explorer.register_results()
    assert registered == 1
    assert state.get_bad_saddle_count() == 1
    explorer.make_jobs()
    assert comm.submitted
    assert "displacement.con" in comm.submitted[0]


def test_catalog_insert_without_the_store_returns_false(tmp_path):
    from pathlib import Path

    from eon.process_catalog import insert

    reactant = _atoms(2.5)
    saddle = _atoms(2.7)
    product = _atoms(3.1)
    con = tmp_path / "reactant.con"
    io.savecon(str(con), reactant)
    cfg = _config(tmp_path)
    cfg.kdb_path = str(tmp_path / "kdb")
    Path(cfg.kdb_path).mkdir()
    states = AKMCStateList(
        300.0 / 11604.5, 20.0, 40.0, initial_state=str(con), config=cfg
    )
    state = states.get_state(0)
    proc_id = state.add_process(_result(reactant, saddle, product, -0.7))
    assert insert(state, proc_id, cfg) is False


def test_numpy_markov_times_are_finite():
    q_matrix = np.array([[0.0, 0.1], [0.2, 0.0]])
    exit_matrix = np.array([[0.9], [0.8]])
    costs = np.array([1.0, 1.0])
    times, branching, residual = np_mcamc(q_matrix, exit_matrix, costs)
    assert times.shape == (2,)
    assert branching.shape == (2, 1)
    assert np.isfinite(times).all()
    assert np.isfinite(residual)
    assert guess_precision(q_matrix, exit_matrix) in {"f", "d", "dd", "qd"}
    assert estimate_condition(q_matrix, exit_matrix) > 0.0
