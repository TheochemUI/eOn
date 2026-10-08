"""Library bodies the coverage table still left unexecuted."""

import ast
import math
import sys
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from eon import atoms
from eon import approach
from eon import communicator
from eon import concorpus
from eon import fileio as io
from eon.akmc import akmc, kmc_step
from eon.akmcstate import lambertw
from eon.akmcstatelist import AKMCStateList
from eon.askmc import ASKMC
from eon.basinhopping import basinhopping
from eon.config import ConfigClass
from eon.contracts.load_outcome import load_outcome
from eon.escaperate import parallelreplica as escape_pass
from eon.escaperate import register_results as escape_register
from eon.explorer import ProcessSearch, ServerMinModeExplorer
from eon.geometry.process import get_process_atoms
from eon.mcamc.mcamc import c_mcamc, guess_precision
from eon.mpiwait import signal_handler
from eon.parallelreplica import parallelreplica, register_results
from eon.process_catalog import (
    _accepts,
    _mark_consumed,
    _mark_queried,
    _max_distance,
    _mode_cosine,
    env_hash,
    pack_frame_key,
    unpack_frame_key,
    was_queried,
)
from eon.prstatelist import PRStateList
from eon.recycling import SB_Recycling
from eon.server import _fallback_single_job
from eon.structure import Structure
from eon.superbasin import Superbasin
from eon.superbasinscheme import TransitionCounting


class _Comm:
    def __init__(self, results=()):
        self.results = list(results)
        self.submitted = []

    def get_results(self, path, keep):
        for result in self.results:
            name = result.get("name", "kept")
            if keep(name):
                yield result

    def queued_search_count(self):
        return 0

    def in_progress_search_count(self):
        return 0

    def submit_jobs(self, searches, invariants):
        self.submitted.extend(searches)

    def cancel_state(self, number):
        return 0


def _atoms(x1, names=("Pt", "Pt")):
    atoms_ = Structure(len(names))
    atoms_.names = list(names)
    atoms_.mass[:] = 195.084
    atoms_.box = np.diag([20.0, 20.0, 20.0])
    if len(names) > 1:
        atoms_.r[1] = [x1, 0.0, 0.0]
    return atoms_


def _buf(atoms_):
    buf = StringIO()
    io.savecon(buf, atoms_)
    buf.seek(0)
    return buf


def _mode(atoms_):
    mode = np.zeros_like(atoms_.r)
    mode[min(1, len(atoms_) - 1), 0] = 0.1
    buf = StringIO()
    io.save_mode(buf, mode)
    buf.seek(0)
    return buf


def _config(directory, job="akmc"):
    path = directory / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                f"job = {job}",
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


def _result(reactant, saddle, product, energy_saddle=-0.7):
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
        "wuid": 1,
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


def _states(tmp_path):
    reactant = _atoms(2.5)
    con = tmp_path / "reactant.con"
    io.savecon(str(con), reactant)
    cfg = _config(tmp_path)
    kT = 300.0 / 11604.5
    states = AKMCStateList(kT, 20.0, 40.0, initial_state=str(con), config=cfg)
    state = states.get_state(0)
    proc = state.add_process(_result(reactant, _atoms(2.7), _atoms(3.1)))
    created = states.get_product_state(0, proc)
    # The list caches the state it builds while linking the reverse process.
    product = states.get_state(created.number)
    return cfg, states, state, product, proc, reactant


def _dat(**fields):
    buf = StringIO()
    io.save_results_dat(buf, fields)
    buf.seek(0)
    return buf


def test_catalog_keys_queries_and_a_rejected_frame(tmp_path):
    key = pack_frame_key(7, 3)
    assert unpack_frame_key(key) == (7, 3)
    with pytest.raises(ValueError, match="12"):
        unpack_frame_key(b"short")
    cluster = _atoms(2.5)
    digest = env_hash(cluster)
    assert len(digest) == 16
    assert env_hash(cluster) == digest
    cfg = _config(tmp_path)
    cfg.kdb_scratch_path = str(tmp_path / "scratch")
    state = SimpleNamespace(number=4)
    assert was_queried(state, cfg) is False
    _mark_queried(state, cfg)
    _mark_queried(state, cfg)
    assert was_queried(state, cfg) is True
    _mark_consumed(cfg, state, 2)
    assert _max_distance(cluster, cluster) == 0.0
    assert _mode_cosine(np.zeros((2, 3)), cluster, cluster) == 0.0
    moved = cluster.copy()
    moved.r[1, 0] = 4.0
    cosine = _mode_cosine(moved.r - cluster.r, cluster, moved)
    assert cosine == pytest.approx(1.0)
    empty = Structure(0)
    assert _max_distance(empty, empty) == 0.0
    process = SimpleNamespace(
        saddle_frame_key=pack_frame_key(1, 0),
        reactant_frame_key=pack_frame_key(1, 1),
    )
    assert _accepts(process, cluster, tmp_path / "missing-corpus", 0.0, 0.2) is None


def test_hop_registers_a_new_minimum_and_a_repeat(tmp_path, monkeypatch):
    root = tmp_path / "run"
    root.mkdir()
    cluster = _atoms(2.5)
    io.savecon(str(root / "pos.con"), cluster)
    cfg = _config(root, job="basin_hopping")
    cfg.bh_initial_state_pool_size = 1
    cfg.comm_job_buffer_size = 1
    cfg.recycling_on = False
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    minimum = _buf(cluster)
    comm = _Comm(
        [
            {"results.dat": _dat(minimum_energy=-1.5, termination_reason=0), "min.con": minimum},
            {
                "results.dat": _dat(minimum_energy=-1.5, termination_reason=0),
                "min.con": _buf(cluster),
            },
        ]
    )
    monkeypatch.chdir(root)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    basinhopping(cfg)
    table = (Path(cfg.path_states) / "state_table").read_text()
    assert "0" in table
    assert (Path(cfg.path_states) / "0" / "minimum.con").is_file()
    assert comm.submitted
    assert Path("wuid.dat").read_text().strip() == "1"


def test_superbasin_step_leaves_through_the_only_exit(tmp_path):
    cfg, states, state, product, proc, _reactant = _states(tmp_path)
    for _ in range(100):
        state.inc_repeats()
    basin_root = tmp_path / "basins"
    basin_root.mkdir()
    basin = Superbasin(basin_root, 3, state_list=[state], config=cfg)

    def get_product(number, process_id):
        assert number == 0
        assert process_id == proc
        return product

    mean_time, exit_state, product_state, exit_proc, basin_id = basin.step(
        state, get_product
    )
    assert mean_time > 0.0
    assert exit_state.number == 0
    assert product_state.number == 1
    assert exit_proc == proc
    assert basin_id == 3
    assert basin.contains_state(state) is True
    assert basin.contains_state(SimpleNamespace()) is False
    assert 0.0 <= basin.get_confidence() <= 1.0
    assert basin.get_lowest_confidence_state().number == 0
    again = Superbasin(basin_root, 3, get_state=states.get_state, config=cfg)
    assert again.state_numbers == [0]
    storage = tmp_path / "kept"
    storage.mkdir()
    again.delete(storage)
    assert (storage / "3").is_file()


def test_kmc_step_crosses_a_confident_superbasin(tmp_path):
    cfg, states, state, _product, proc, _reactant = _states(tmp_path)
    cfg.akmc_max_kmc_steps = 1
    cfg.sb_on = True
    cfg.askmc_on = False
    for _ in range(100):
        state.inc_repeats()
    basin_root = tmp_path / "basins"
    basin_root.mkdir()
    basin = Superbasin(basin_root, 1, state_list=[state], config=cfg)

    class _Scheme:
        def get_containing_superbasin(self, current):
            if current.number == state.number:
                return basin
            return None

        def register_transition(self, start, end):
            return None

        def write_data(self):
            return None

    current, previous, time, steps = kmc_step(
        state, states, 0.0, states.kT, _Scheme(), config=cfg
    )
    assert steps == 1
    assert time > 0.0
    assert previous.number == 0
    assert current.number == 1
    text = (Path(cfg.path_results) / "dynamics.txt").read_text()
    assert text


def test_askmc_raises_a_barrier_both_ways(tmp_path):
    cfg, states, state, product, proc, _reactant = _states(tmp_path)
    kT = states.kT
    ask = ASKMC(
        kT,
        states,
        0.5,
        1.5,
        2.0,
        True,
        True,
        True,
        str(tmp_path),
        20.0,
        str(tmp_path / "recycle"),
    )
    forward = dict(state.get_process(proc))
    forward["view_count"] = 10
    reverse_id = next(iter(product.get_process_table()))
    backward = dict(product.get_process(reverse_id))
    backward["view_count"] = 10
    ask.save_modified_process_table(state, {proc: forward})
    ask.save_modified_process_table(product, {reverse_id: backward})
    old_rate = forward["rate"]
    ask.raiseup(state, product, 0, 0)
    raised = ask.get_modified_process_table(state)[proc]["rate"]
    assert raised == pytest.approx(old_rate / 1.5)
    assert (Path(ask.recycle_path) / "current_sb_states").is_file()
    checks, changes = ask.get_askmc_metadata()
    assert checks == 0
    assert changes == 1
    with pytest.raises(KeyError, match="no process"):
        ask.get_process_id({}, 9, "find")
    with pytest.raises(ValueError, match="unknown process lookup"):
        ask.get_process_id({}, 9, "other")


def test_confidence_schemes_and_a_stored_bad_saddle(tmp_path):
    cfg, states, state, _product, proc, reactant = _states(tmp_path)
    state.inc_proc_random_count(proc)
    for _ in range(12):
        state.inc_proc_random_count(proc)
    cfg.akmc_confidence_correction = True
    cfg.akmc_confidence_scheme = "new"
    new_conf = state.get_confidence()
    assert 0.0 <= new_conf <= 1.0
    cfg.akmc_confidence_scheme = "sampling"
    sampling = state.get_confidence()
    assert 0.0 <= sampling <= 1.0
    cfg.akmc_confidence_scheme = "dynamics"
    cfg.disp_moved_only = False
    state.increment_time(5.0, 300)
    state.increment_time(1.0, 300)
    by_temp = state.get_time_by_temp()
    assert by_temp[300] > 0.0
    dynamics = state.get_confidence()
    assert 0.0 <= dynamics <= 1.0
    assert lambertw(0.0) == 0.0
    assert math.isfinite(lambertw(2.0))
    with pytest.raises(ValueError):
        lambertw(-1.0)
    bad = _result(reactant, _atoms(2.7), _atoms(3.1), energy_saddle=-0.2)
    bad["results"]["termination_reason"] = 15
    bad["results"]["simulation_time"] = 2.0
    bad["results"]["md_temperature"] = 400
    bad["wuid"] = 9
    state.register_bad_saddle(bad, store=True)
    assert state.get_bad_saddle_count() >= 1
    assert state.get_total_saddle_count() >= 1
    assert list(Path(state.bad_procdata_path).glob("saddle_9.con"))
    assert state.get_process_reactant(proc)
    assert state.get_process_saddle(proc) is not None
    assert state.get_process_mode(proc).shape[1] == 3
    fresh = states.get_state(0)
    assert fresh.get_number_of_searches() >= 1


def test_reverse_candidate_links_an_unassigned_process(tmp_path):
    cfg, states, state, product, proc, reactant = _states(tmp_path)
    extra = product.allocate_process_id(b"candidate", b"extra")
    forward = state.get_process(proc)
    product.append_process_table(
        id=extra,
        saddle_energy=forward["saddle_energy"],
        prefactor=forward["product_prefactor"],
        product=-1,
        product_energy=forward["product_energy"],
        product_prefactor=forward["prefactor"],
        barrier=forward["barrier"],
        rate=forward["rate"],
        repeats=0,
    )
    io.savecon(product.proc_product_path(extra), reactant)
    cfg.akmc_eq_rate = 1.0
    states.register_process(0, 1, proc)
    linked = product.get_process_table()[extra]
    assert linked["product"] == 0
    assert linked["rate"] > 0.0
    extra2 = product.allocate_process_id(b"connect", b"extra")
    product.append_process_table(
        id=extra2,
        saddle_energy=forward["saddle_energy"],
        prefactor=forward["prefactor"],
        product=-1,
        product_energy=state.get_energy(),
        product_prefactor=forward["product_prefactor"],
        barrier=forward["barrier"],
        rate=forward["rate"],
        repeats=0,
    )
    io.savecon(product.proc_saddle_path(extra2), _atoms(2.8))
    io.savecon(product.proc_reactant_path(extra2), product.get_reactant())
    io.savecon(product.proc_product_path(extra2), reactant)
    io.save_mode(product.proc_mode_path(extra2), np.zeros((2, 3)))
    Path(product.proc_results_path(extra2)).write_text("0 termination_reason\n")
    states.connect_states([state, product])
    assert product.get_process(extra2)["product"] == 0


def test_server_search_finishes_both_minima(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg, states, state, _product, _proc, reactant = _states(tmp_path)
    cfg.akmc_server_side_process_search = True
    cfg.comm_job_buffer_size = 1
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.process_search_minimization_offset = 0.1
    Path(cfg.path_incomplete).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ServerMinModeExplorer(states, state, state, config=cfg)
    assert explorer.register_results() == 0
    explorer.make_jobs()
    assert comm.submitted
    assert Path("searches.log").is_file()
    displacement = reactant.copy()
    displacement.r[1, 0] += 0.3
    mode = np.zeros_like(reactant.r)
    mode[1, 0] = 1.0
    search = ProcessSearch(reactant, displacement, mode, "random", 4, state.number, config=cfg)
    job, kind = search.get_job(state.number)
    assert kind == "saddle_search"
    assert "displacement.con" in job
    saddle = _atoms(2.7)
    search.process_result(
        {
            "name": "0_4",
            "type": "random",
            "results": {
                "job_type": "saddle_search",
                "termination_reason": 0,
                "potential_energy_saddle": -0.7,
                "potential_energy_reactant": -1.0,
                "potential_energy_product": -0.8,
                "barrier_reactant_to_product": 0.3,
                "total_force_calls": 2,
            },
            "saddle.con": _buf(saddle),
            "mode.dat": _mode(saddle),
        }
    )
    assert search.get_saddle() is not None
    _job, kind = search.get_job(state.number)
    assert kind == "min1"
    search.process_result(
        {
            "name": "0_5",
            "results": {
                "job_type": "minimization",
                "termination_reason": 0,
                "potential_energy": -1.0,
                "total_force_calls": 3,
            },
            "min.con": _buf(reactant),
        }
    )
    _job, kind = search.get_job(state.number)
    assert kind == "min2"
    finished = search.process_result(
        {
            "name": "0_6",
            "results": {
                "job_type": "minimization",
                "termination_reason": 0,
                "potential_energy": -0.8,
                "total_force_calls": 3,
            },
            "min.con": _buf(_atoms(3.1)),
        }
    )
    assert finished["results"]["barrier_reactant_to_product"] == pytest.approx(0.3)
    assert search.get_job(state.number) == (None, None)
    assert search.get_job(99) == (None, None)


def test_replica_and_escape_register_a_transition(tmp_path, monkeypatch):
    root = tmp_path / "run"
    root.mkdir()
    cluster = _atoms(2.5)
    product = _atoms(3.1)
    io.savecon(str(root / "pos.con"), cluster)
    cfg = _config(root, job="parallel_replica")
    cfg.comm_job_buffer_size = 1
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    fields = {
        "speedup": 2.5,
        "transition_found": 1,
        "transition_time_s": 3.0,
        "simulation_time_s": 1.0,
        "potential_energy_reactant": -1.0,
        "potential_energy_product": -0.8,
    }
    result = {
        "name": "0_3",
        "reactant.con": _buf(cluster),
        "product.con": _buf(product),
        "results.dat": _dat(**fields),
    }
    comm = _Comm([result])
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    monkeypatch.chdir(root)
    parallelreplica(cfg)
    assert (root / "states" / "1").is_dir()
    text = (root / "dynamics.txt").read_text()
    assert "0" in text

    idle = {
        "name": "0_4",
        "results.dat": _dat(speedup=1.0, transition_found=0, simulation_time_s=4.0),
    }
    states = PRStateList(str(root / "pos.con"), config=cfg)
    count, transition, speedup = register_results(
        _Comm([idle]), states.get_state(0), states, cfg
    )
    assert count == 1
    assert transition is None
    assert speedup == pytest.approx(1.0)

    escape_root = tmp_path / "escape"
    escape_root.mkdir()
    io.savecon(str(escape_root / "pos.con"), cluster)
    escape_cfg = _config(escape_root, job="escape_rate")
    escape_cfg.comm_job_buffer_size = 0
    Path(escape_cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    registered, found, gained = escape_register(
        _Comm([result]),
        PRStateList(str(escape_root / "pos.con"), config=escape_cfg).get_state(0),
        PRStateList(str(escape_root / "pos.con"), config=escape_cfg),
        escape_cfg,
    )
    assert registered == 1
    assert found["time"] > 0.0
    assert gained == pytest.approx(2.5)
    with pytest.raises(FileNotFoundError, match="root directory does not exist"):
        escape_pass(SimpleNamespace(path_root=str(tmp_path / "absent")))


def test_transition_counting_opens_a_basin(tmp_path):
    cfg, states, state, product, _proc, _reactant = _states(tmp_path)
    cfg.comp_use_identical = False
    cfg.sb_max_size = 0
    root = tmp_path / "scheme"
    root.mkdir()
    (root / "storage").mkdir()
    (root / "notes.txt").write_text("debris\n")
    scheme = TransitionCounting(root, states, states.kT, 1, config=cfg)
    scheme.register_transition(state, state)
    scheme.register_transition(state, product)
    assert scheme.get_containing_superbasin(state) is not None
    assert scheme.get_count(state)[product] >= 1


def test_recycling_resumes_and_starts_from_askmc(tmp_path):
    _cfg, states, state, product, _proc, _reactant = _states(tmp_path)
    path = tmp_path / "sb-recycle"
    path.mkdir()
    (path / "recycling_data.txt").write_text(
        "Superbasin Recycling Metadata\n"
        "'Prev_Superbasin_state_list' = [[0, 1]]\n"
        "'in_progress' = True\n"
    )
    resumed = SB_Recycling(states, state, product, 0.2, False, path, None, None)
    assert resumed.in_progress is True
    assert resumed.sb_states[0][0].number == 0
    suggestion = resumed.make_suggestion()
    assert suggestion == (None, None) or suggestion[0] is not None

    fresh = tmp_path / "fresh"
    (fresh).mkdir()
    started = SB_Recycling(states, state, product, 0.2, False, fresh, "askmc", None)
    assert started.in_progress is False
    listed = tmp_path / "listed"
    listed.mkdir()
    (listed / "current_sb_states").write_text("header\n[4, 5]\n")
    missed = SB_Recycling(states, state, product, 0.2, False, listed, "askmc", None)
    assert missed.in_progress is False


def test_table_poscar_atoms_and_small_helpers(tmp_path):
    path = tmp_path / "rates.tbl"
    table = io.Table(str(path), ["state", "rate"])
    table.add_row({"state": 1, "rate": 2.5})
    table.add_row({"state": 2, "rate": 0.5})
    assert table.min_value("rate") == 0.5
    assert table.max_row("rate")["state"] == 1
    assert table.get_rows("state", 1)[0]["rate"] == 2.5
    assert table.get_column("state") == [1, 2]
    assert table.delete_row("state", 2) == 1
    assert table.delete_row_func("state", lambda value: value == 1) == 1
    with pytest.raises(io.TableException, match="Mismatched"):
        table.add_row({"state": 3})
    bare = io.Table(str(path))
    assert len(bare) == 0
    assert bare.columns == ["state", "rate"]
    assert "state" in repr(bare)
    with pytest.raises(io.TableException, match="not optional"):
        empty = io.Table(str(tmp_path / "new.tbl"))
        len(empty)

    cluster = _atoms(2.5)
    poscar = tmp_path / "POSCAR"
    io.saveposcar(str(poscar), cluster)
    loaded = io.loadposcar(str(poscar))
    assert loaded.names[0] == "Pt"
    assert loaded.r.shape == (2, 3)
    selective = tmp_path / "selective.poscar"
    selective.write_text(
        "\n".join(
            [
                "Pt",
                "1.0",
                "10 0 0",
                "0 10 0",
                "0 0 10",
                "1",
                "Selective dynamics",
                "Cartesian",
                "0 0 0 T T F",
                "",
            ]
        )
    )
    one = io.loadposcar(str(selective))
    assert one.free[0, 2] == 0.0
    direct = tmp_path / "direct.poscar"
    io.saveposcar(str(direct), cluster, direct=True)
    assert io.loadposcar(str(direct)).r.shape == (2, 3)

    other = cluster.copy()
    assert atoms.match(
        cluster, other, 0.1, 3.5, True, check_rotation=True, use_identical=True
    ) == {0: 0, 1: 1}
    shifted = cluster.copy()
    shifted.r += 0.2
    assert atoms.remove_net_translation(cluster, shifted) is not None
    io.savecon(str(tmp_path / "a.con"), cluster)
    io.savecon(str(tmp_path / "b.con"), other)
    assert atoms.point_energy_match(
        str(tmp_path / "a.con"), -1.0, str(tmp_path / "b.con"), -1.0, 0.1, 0.2, 3.5
    )
    assert (
        atoms.points_energies_match(
            str(tmp_path / "a.con"),
            -1.0,
            [str(tmp_path / "b.con")],
            [-1.0],
            0.1,
            0.2,
            3.5,
        )
        == 0
    )
    assert atoms.points_energies_match(
        str(tmp_path / "a.con"), -1.0, [], [], 0.1, 0.2, 3.5
    ) is None
    assert len(atoms.cna(cluster, 3.5)) == 2
    assert len(atoms.cnat(cluster, 3.5)) == 2
    assert len(atoms.cnar(cluster, 3.5)) == 2
    assert 0 in atoms.not_TCP(cluster, 3.5)
    assert 0 in atoms.not_TCP_or_BCC(cluster, 3.5)
    assert 0 in atoms.not_HCP_or_FCC(cluster, 3.5)
    matrix = atoms.rotm(np.array([0.0, 0.0, 1.0]), 0.3)
    assert matrix.shape == (3, 3)
    turned = atoms.rotate(cluster.r, np.array([0.0, 0.0, 1.0]), np.zeros(3), 0.2)
    assert turned.shape == cluster.r.shape
    assert atoms.crystal_spacegroup(cluster)
    mobile = get_process_atoms(cluster, _atoms(3.1), epsilon_r=0.2, nshells=0)
    assert mobile == [1]

    rate = approach.arrhenius_rate(0.1, 300.0)
    assert rate > 0.0
    with pytest.raises(ValueError):
        approach.arrhenius_rate(-0.1, 300.0)
    tau = approach.prepared_imbalance(rate, rate)
    assert approach.survival(0.0, tau) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        approach.prepared_imbalance(0.0, 0.0)
    with pytest.raises(ValueError):
        approach.survival(-1.0, tau)
    assert approach.recurrence_time(2.0) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        approach.recurrence_time(0.0)

    outcome = load_outcome('{"jobType": "akmc", "statusCode": 0, "statusText": "ok"}')
    assert outcome["statusCode"] == 0
    with pytest.raises(ValueError, match="trajectories"):
        load_outcome({"jobType": "akmc", "statusCode": 0, "statusText": "ok", "trajectory": []})
    with pytest.raises(KeyError):
        load_outcome({})
    with pytest.raises(TypeError):
        load_outcome([1, 2])

    con = tmp_path / "frame.con"
    io.savecon(str(con), cluster)
    assert concorpus.store_frame_text(con, "  ") is None
    blob = con.read_text()
    stored = concorpus.store_frame_text(con, blob)
    if stored is None:
        assert concorpus.stored_frame_text(con) is None
    else:
        assert stored[1] == 0
        text = concorpus.load_frame_text(
            concorpus.corpus_dir(con), stored[0], stored[1]
        )
        assert "Pt" in text
        assert concorpus.stored_frame_text(con)
    concorpus.mirror_con_path(con)
    assert tmp_path in concorpus.directories_with_con(tmp_path) or any(
        path.name == tmp_path.name for path in concorpus.directories_with_con(tmp_path)
    )

    lock_path = tmp_path / "lockfile"
    lock_path.write_text("99999999\n")
    from eon.locking import LockFile

    lock = LockFile(lock_path)
    assert lock.islocked() is False
    assert lock.aquirelock() is True
    assert lock.islocked() is True
    lock.removelock()

    signal_handler(15, None)
    from eon import mpiwait

    assert mpiwait.QUIT is True
    mpiwait.QUIT = False

    Q = np.zeros((2, 2))
    R = np.ones((2, 1))
    assert guess_precision(Q, R) == "d"
    with pytest.raises(ValueError, match="Unknown prec"):
        c_mcamc(Q, R, np.ones(2), prec="nope")

    import eon

    eon.__dict__.pop("server", None)
    assert callable(eon.server)
    with pytest.raises(AttributeError):
        getattr(eon, "missing_name")
    assert "server" in dir(eon)
    import eon.__main__ as entry

    assert callable(entry.main)


def test_each_displacement_kind(tmp_path):
    from eon.displace import DisplacementManager

    cfg = _config(tmp_path)
    reactant = _atoms(2.5)
    water = Structure(3)
    water.names = ["H", "H", "O"]
    water.mass[:] = [1.0, 1.0, 16.0]
    water.box = np.diag([12.0, 12.0, 12.0])
    water.r[0] = [0.0, 0.8, 0.0]
    water.r[1] = [0.8, -0.2, 0.0]
    water.r[2] = [0.0, 0.0, 0.0]
    weights = [
        "displace_random_weight",
        "displace_listed_atom_weight",
        "displace_listed_type_weight",
        "displace_under_coordinated_weight",
        "displace_least_coordinated_weight",
        "displace_not_FCC_HCP_weight",
        "displace_not_TCP_BCC_weight",
        "displace_not_TCP_weight",
        "displace_water_weight",
    ]
    cfg.disp_listed_atoms = [0]
    cfg.disp_listed_types = ["Pt"]
    cfg.molecule_list = []
    cfg.disp_at_random = 0
    for name in weights:
        for other in weights:
            setattr(cfg, other, 0.0)
        setattr(cfg, name, 1.0)
        subject = water if name == "displace_water_weight" else reactant
        made, mode = DisplacementManager(subject, None, config=cfg).make_displacement()
        assert len(made) == len(subject)
        assert np.all(np.isfinite(mode))


def test_config_parser_rejections_and_communicator_branches(tmp_path, monkeypatch):
    missing = tmp_path / "absent.ini"
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(missing))
    assert caught.value.code == 2

    unknown = tmp_path / "unknown.ini"
    unknown.write_text("[Not A Section]\nvalue = 1\n")
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(unknown))
    assert caught.value.code == 1

    bad_int = tmp_path / "bad.ini"
    bad_int.write_text("[Main]\nrandom_seed = no\n")
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(bad_int))
    assert caught.value.code == 1

    server_side = tmp_path / "server.ini"
    server_side.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "[AKMC]",
                "server_side_process_search = true",
                "[Prefactor]",
                "default_value = 0",
                "[Communicator]",
                "jobs_per_bundle = 1",
                "[Paths]",
                f"main_directory = {tmp_path}",
                f"results = {tmp_path}",
                f"states = {tmp_path / 'states'}",
                "",
            ]
        )
    )
    with pytest.raises(SystemExit):
        ConfigClass().init(str(server_side))

    bundle = tmp_path / "bundle.ini"
    text = server_side.read_text().replace("default_value = 0", "default_value = 1e12")
    text = text.replace("jobs_per_bundle = 1", "jobs_per_bundle = 2")
    bundle.write_text(text)
    with pytest.raises(SystemExit):
        ConfigClass().init(str(bundle))

    mpi_ini = tmp_path / "mpi.ini"
    mpi_ini.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "random_seed = 1",
                "[Communicator]",
                "type = mpi",
                "[Paths]",
                f"main_directory = {tmp_path}",
                f"results = {tmp_path}",
                f"states = {tmp_path / 'states'}",
                "",
            ]
        )
    )
    previous = sys.excepthook
    cfg = ConfigClass()
    cfg.init(str(mpi_ini))
    assert sys.excepthook is not previous
    sys.excepthook = previous
    assert cfg.comm_type == "mpi"

    local_ini = tmp_path / "local.ini"
    local_ini.write_text(
        mpi_ini.read_text().replace("type = mpi", "type = local")
    )
    local = ConfigClass()
    local.init(str(local_ini))
    assert local.comm_local_ncpus >= 1

    cluster_ini = tmp_path / "cluster.ini"
    cluster_ini.write_text(
        mpi_ini.read_text().replace("type = mpi", "type = cluster")
    )
    cluster = ConfigClass()
    cluster.init(str(cluster_ini))
    assert cluster.comm_script_name_prefix

    job = tmp_path / "job"
    job.mkdir()
    (job / "results.dat").write_text("0 termination_reason\n")
    (job / "pos.con").write_text("Pt\n")
    (job / "return_files.dat").write_text("pos.con\nmissing.con\npos.con\n")
    harvested = communicator.harvest_job_files(job)
    assert "pos.con" in harvested
    assert "missing.con" not in harvested
    assert "results.dat" in harvested
    plain = tmp_path / "plain"
    plain.mkdir()
    (plain / "mode.dat").write_text("0 0 0\n")
    assert "mode.dat" in communicator.harvest_job_files(plain)
    assert communicator.client_environment({"PATH": "/usr/bin"})["PATH"] == "/usr/bin"

    monkeypatch.chdir(tmp_path)
    (tmp_path / "note.txt").write_text("hello\n")
    (tmp_path / "output").mkdir()
    cfg = _config(tmp_path)
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    _fallback_single_job(cfg)
    assert comm.submitted[0]["id"] == "output"
    assert "note.txt" in comm.submitted[0]
    assert (tmp_path / "output_old").is_dir()


def test_akmc_drops_stale_superbasins_when_temperature_changes(tmp_path, monkeypatch):
    root = tmp_path / "run"
    root.mkdir()
    io.savecon(str(root / "pos.con"), _atoms(2.5))
    cfg = _config(root)
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.sb_on = False
    cfg.comm_job_buffer_size = 1
    cfg.comm_type = "local_inprocess"
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    comm = _Comm()
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    assert akmc(cfg) == 0
    cfg.main_temperature = 500.0
    Path(cfg.sb_path).mkdir(parents=True, exist_ok=True)
    (Path(cfg.sb_path) / "marker").write_text("x\n")
    state_dir = Path(cfg.path_states) / "0"
    marker = state_dir / cfg.sb_state_file
    marker.write_text("1 1\n")
    assert akmc(cfg) == 0
    assert not Path(cfg.sb_path).exists()
    assert not marker.exists()


def test_movie_fastest_path_and_pr_state_list(tmp_path, capsys):
    from eon.movie import fastest_path, get_fastest_process_rate

    _cfg, states, state, product, proc, _reactant = _states(tmp_path)
    frames = fastest_path(tmp_path, states, full=True)
    assert len(frames) >= 2
    assert get_fastest_process_rate(state, product) > 0.0
    out = capsys.readouterr().out
    assert "0" in out
    with pytest.raises(ValueError, match="no process"):
        get_fastest_process_rate(product, product)

    root = tmp_path / "pr"
    root.mkdir()
    io.savecon(str(root / "pos.con"), _atoms(2.5))
    cfg = _config(root)
    listed = PRStateList(str(root / "pos.con"), config=cfg)
    first = listed.get_state(0)
    first.add_process(
        {
            "reactant.con": _buf(_atoms(2.5)),
            "product.con": _buf(_atoms(3.1)),
            "results.dat": _dat(
                potential_energy_reactant=-1.0,
                potential_energy_product=-0.8,
                transition_time_s=1.5,
            ),
            "results": {
                "potential_energy_reactant": -1.0,
                "potential_energy_product": -0.8,
                "transition_time_s": 1.5,
            },
        }
    )
    process_id = next(iter(first.get_process_table()))
    made = listed.get_product_state(0, process_id)
    assert made.number == 1
    assert first.get_process(process_id)["product"] == 1
