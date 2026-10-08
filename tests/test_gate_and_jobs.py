"""Gate decisions, in-process job types, and the remaining state branches."""

import math
from io import StringIO

import numpy as np
from pathlib import Path
from types import SimpleNamespace

import pytest

from eon import fileio as io
from eon.amsel_superbasin_gate import (
    AmselSuperbasinReject,
    apply_gate_to_superbasin,
    basin_view_from_state,
    discover_decide_for_superbasin,
    load_persisted_split,
    persist_split,
    status_pair,
)
from eon.communicator import CommunicatorError
from eon.communicator_inprocess import (
    LocalInProcess,
    _require_pyeonclient,
    _run_inprocess_job,
)
from eon.config import ConfigClass
from eon.displace import DisplacementManager
from eon.explorer import ProcessSearch, ServerMinModeExplorer
from eon.movie import priorityDictionary
from eon.superbasinscheme import EnergyLevel
from tests.test_library_bodies import _Comm, _atoms, _buf, _config, _mode, _result, _states


class _Basin:
    def __init__(self, path, states):
        self.path = str(path)
        self.states = list(states)
        self.state_numbers = [state.number for state in states]
        self.state_dict = {state.number: state for state in states}
        self.wrote = 0

    def write_data(self):
        self.wrote += 1


def test_gate_splits_rejects_and_remembers(tmp_path):
    cfg, states, first, second, _proc, _reactant = _states(tmp_path)
    basin = _Basin(tmp_path / "basin-1", [first, second])
    assert load_persisted_split(SimpleNamespace(path=None)) is None
    assert load_persisted_split(basin) is None
    persist_split(basin, [0])
    assert load_persisted_split(basin) == [0]
    (Path(str(basin.path) + ".amsel_split.json")).write_text("{")
    assert load_persisted_split(basin) is None
    persist_split(basin, [])
    assert load_persisted_split(basin) is None

    assert apply_gate_to_superbasin(
        basin, first, {"available": True, "status": "accepted"}
    ) == "accepted"
    assert apply_gate_to_superbasin(
        basin, first, {"available": True, "status": "split_required"}
    ) == "split_required"
    assert apply_gate_to_superbasin(
        _Basin(tmp_path / "basin-2", [first, second]),
        first,
        {"available": True, "status": "split_required", "primary_transient": [99]},
    ) == "split_required"

    splitting = _Basin(tmp_path / "basin-3", [first, second])
    assert apply_gate_to_superbasin(
        splitting,
        first,
        {"available": True, "status": "split_required", "primary_transient": [0]},
    ) == "split_required"
    assert splitting.state_numbers == [0]
    assert splitting.wrote == 1
    remembered = _Basin(tmp_path / "basin-3", [first, second])
    assert apply_gate_to_superbasin(
        remembered,
        first,
        {"available": True, "status": "split_required", "primary_transient": [0, 1]},
    ) == "split_cached"
    assert remembered.state_numbers == [0]

    with pytest.raises(AmselSuperbasinReject) as caught:
        apply_gate_to_superbasin(
            basin,
            first,
            {"available": True, "status": "rejected", "reason": "no basin"},
        )
    assert caught.value.status == "rejected"
    assert status_pair("unavailable", False) is True
    assert status_pair("unavailable", True) is False
    assert status_pair("not-a-status", True) is False

    view = basin_view_from_state(first, states.get_state, 0.01)
    assert 0 in view.state_numbers
    assert 1 in view.state_numbers

    class _Broken:
        def discover_decide_status(self, *args, **kwargs):
            raise RuntimeError("amsel down")

    # discover_decide imports the name from the module.
    import sys
    from types import ModuleType

    module = ModuleType("amsel")

    def discover_decide_status(*_args, **_kwargs):
        raise RuntimeError("amsel down")

    module.discover_decide_status = discover_decide_status
    sys.modules["amsel"] = module
    try:
        with pytest.raises(RuntimeError, match="amsel down"):
            discover_decide_for_superbasin(view, first, on_error="raise")
        fallback = discover_decide_for_superbasin(view, first, on_error="fallback_single")
        assert fallback["status"] == "fallback_single"
        assert fallback["available"] is True
        missing = discover_decide_for_superbasin(view, first, on_error="unavailable_mcamc")
        assert missing["status"] == "unavailable"
        assert missing["available"] is False
    finally:
        sys.modules.pop("amsel", None)


def test_inprocess_job_types_and_con_text(tmp_path):
    cfg = _config(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    local = LocalInProcess(str(scratch), 1, config=cfg)
    cluster = _atoms(2.5)
    text = StringIO()
    io.savecon(text, cluster)
    text.seek(0)
    local.submit_jobs([{"id": "text", "pos.con": text}], {})
    done = local.get_results()
    assert done[0]["id"] == "text"
    assert done[0]["min.con"].read().startswith("Generated") or "Pt" in done[0]["min.con"].getvalue()
    raw = StringIO()
    io.savecon(raw, cluster)
    local.submit_jobs([{"id": "bytes", "pos.con": raw.getvalue().encode()}], {})
    assert local.get_results()[0]["id"] == "bytes"
    with pytest.raises(CommunicatorError, match="could not read"):
        local.submit_jobs([{"id": "bad", "pos.con": StringIO("not a con\n")}], {})
    with pytest.raises(CommunicatorError, match="needs a Structure"):
        local.submit_jobs([{"id": "none"}], {})

    pc = _require_pyeonclient()
    params = pc.Parameters()
    params.potential = pc.PotType.LJ
    params.quiet = True
    params.write_log = False
    pot = pc.make_potential(params)
    from pyeonclient.bridge import structure_to_matter

    matter = structure_to_matter(cluster, pot, params)
    expected = [
        (pc.JobType.Hessian, "hessian"),
        (pc.JobType.Prefactor, "prefactor"),
        (pc.JobType.Finite_Difference, "finite_difference"),
        (pc.JobType.Point, "point"),
    ]
    for kind, name in expected:
        payload = _run_inprocess_job(pc, kind, matter, pot, params, {})
        assert payload["job_type"] == name
        assert math.isfinite(payload["energy"])
    with pytest.raises(CommunicatorError, match="no dispatch"):
        _run_inprocess_job(pc, "not-a-job", matter, pot, params, {})


def test_repeat_saddle_high_barrier_and_a_listed_script(tmp_path, monkeypatch):
    cfg, states, state, _product, proc, reactant = _states(tmp_path)
    again = _result(reactant, _atoms(2.7), _atoms(3.1))
    again["wuid"] = 2
    assert state.add_process(again) is None
    assert state.get_process(proc)["repeats"] >= 1

    high = _result(reactant, _atoms(2.9), _atoms(3.2), energy_saddle=8.0)
    high["wuid"] = 3
    assert state.add_process(high) is None

    broken = _result(reactant, _atoms(2.2), _atoms(3.3))
    broken["mode.dat"] = StringIO("not a mode\n")
    broken["wuid"] = 4
    assert state.add_process(broken) is None

    cfg.displace_atom_kmc_state_script = "ids.py"
    cfg.saddle_method = "dynamics"
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_type = "local_inprocess"
    cfg.comm_job_buffer_size = 0
    state.info.set("Saddle Search", "displace_atom_list", "0, 1")
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Comm())
    from eon.explorer import ClientMinModeExplorer

    explorer = ClientMinModeExplorer(states, state, state, config=cfg)
    assert cfg.disp_listed_from_script is True
    assert explorer.generate_displacement() == (None, None, "dynamics")

    class _Boom(_Comm):
        def submit_jobs(self, searches, invariants):
            raise RuntimeError("queue down")

    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Boom())
    cfg.akmc_server_side_process_search = True
    cfg.comm_job_buffer_size = 1
    cfg.saddle_method = "min_mode"
    Path(cfg.path_incomplete).mkdir(parents=True, exist_ok=True)
    server = ServerMinModeExplorer(states, state, state, config=cfg)
    server.make_jobs()
    bare = ProcessSearch(
        reactant, reactant, np.zeros_like(reactant.r), "random", 3, 0, config=cfg
    )
    assert bare.get_saddle() is None
    assert bare.get_saddle_file() is None


def test_water_picks_one_molecule_and_a_capped_basin(tmp_path):
    cfg = _config(tmp_path)
    from eon.structure import Structure

    water = Structure(3)
    water.names = ["H", "H", "O"]
    water.mass[:] = [1.0, 1.0, 16.0]
    water.box = np.diag([12.0, 12.0, 12.0])
    water.r[0] = [0.0, 0.8, 0.0]
    water.r[1] = [0.8, -0.2, 0.0]
    water.r[2] = [0.0, 0.0, 0.0]
    for name in (
        "displace_random_weight",
        "displace_listed_atom_weight",
        "displace_listed_type_weight",
        "displace_under_coordinated_weight",
        "displace_least_coordinated_weight",
        "displace_not_FCC_HCP_weight",
        "displace_not_TCP_BCC_weight",
        "displace_not_TCP_weight",
        "displace_water_weight",
    ):
        setattr(cfg, name, 0.0)
    cfg.displace_water_weight = 1.0
    cfg.disp_at_random = 1
    cfg.molecule_list = []
    made, mode = DisplacementManager(water, None, config=cfg).make_displacement()
    assert len(made) == 3
    assert np.all(np.isfinite(mode))

    level_root = tmp_path / "levels"
    level_root.mkdir()
    cfg2, states, first, second, _proc, _reactant = _states(level_root)
    cfg2.sb_max_size = 1
    root = level_root / "scheme"
    root.mkdir()
    (root / "storage").mkdir()
    scheme = EnergyLevel(str(root), states, states.kT, 0.05, config=cfg2)
    scheme.make_basin([first, second])
    assert scheme.get_containing_superbasin(first) is None
    cfg2.sb_max_size = 0
    scheme.make_basin([first, second])
    assert scheme.get_containing_superbasin(first) is not None
    scheme.make_basin([first, second])
    assert len(scheme.superbasins) == 1


def test_priority_queue_drops_a_stale_key():
    queue = priorityDictionary()
    queue["late"] = 5.0
    queue["late"] = 0.1
    queue["mid"] = 1.0
    queue.setdefault("last", 3.0)
    assert list(queue) == ["late", "mid", "last"]
    with pytest.raises(IndexError):
        queue.smallest()


def test_config_rejects_a_bad_float_and_an_unknown_option(tmp_path):
    path = tmp_path / "bad.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "temperature = no",
                "not_a_key = 1",
                "[Coarse Graining]",
                "use_mcamc = maybe",
                "",
            ]
        )
    )
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(path))
    assert caught.value.code == 1
