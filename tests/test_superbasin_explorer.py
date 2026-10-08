"""Superbasin exploration, coarse-grained recycling, and the precision script."""

import sys
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from eon import fileio as io
from eon.akmc import akmc
from eon.explorer import ClientMinModeExplorer
from eon.process_catalog import insert
from eon.recycling import SB_Recycling
from eon.superbasin import Superbasin
from eon.superbasinscheme import TransitionCounting
from tests.test_catalog_and_mains import _Record, _Store
from tests.test_library_bodies import _Comm, _atoms, _buf, _dat, _states


def test_akmc_explores_the_least_confident_basin_state(tmp_path, monkeypatch):
    cfg, states, state, product, proc, _reactant = _states(tmp_path)
    cfg.sb_on = True
    cfg.sb_scheme = "transition_counting"
    cfg.sb_tc_ntrans = 1
    cfg.comp_use_identical = False
    cfg.sb_max_size = 0
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_job_buffer_size = 1
    cfg.comm_type = "local_inprocess"
    Path(cfg.sb_path).mkdir(parents=True, exist_ok=True)
    (Path(cfg.sb_path) / "storage").mkdir(exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    scheme = TransitionCounting(cfg.sb_path, states, states.kT, 1, config=cfg)
    scheme.register_transition(state, product)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Comm())
    assert akmc(cfg) == 0
    assert scheme.get_containing_superbasin(state) is not None
    assert (Path(cfg.path_results) / "info.txt").is_file()


def test_coarse_grained_recycling_builds_the_next_reactants(tmp_path):
    cfg, states, previous, current, _proc, _reactant = _states(tmp_path)
    root = tmp_path / "basins"
    root.mkdir()
    basin = Superbasin(root, 4, state_list=[previous, current], config=cfg)
    superbasining = SimpleNamespace(
        get_containing_superbasin=lambda state: basin if state.number == previous.number else None
    )
    recycling = SB_Recycling(
        states,
        previous,
        current,
        0.2,
        False,
        tmp_path / "recycle",
        "mcacm",
        superbasining,
    )
    assert recycling.sb_state_nums
    assert recycling.sb_state_nums[0][0] == previous.number


def test_client_registers_a_basin_result_and_a_catalog_mode(tmp_path, monkeypatch):
    from types import ModuleType

    cfg, states, state, product, proc, reactant = _states(tmp_path)
    cfg.recycling_on = False
    cfg.kdb_on = True
    cfg.kdb_only = False
    cfg.kdb_path = str(tmp_path / "kdb")
    cfg.kdb_scratch_path = str(tmp_path / "scratch")
    Path(cfg.kdb_path).mkdir()
    cfg.comm_type = "local"
    cfg.comm_job_buffer_size = 1
    cfg.debug_keep_bad_saddles = True
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    module = ModuleType("amsel")

    def _process(*args, **kwargs):
        return _Record(args[2], args[3], **kwargs)

    module.KdbProcess = _process
    module.KdbStore = _Store
    monkeypatch.setitem(sys.modules, "amsel", module)
    assert insert(state, proc, cfg) is True
    root = tmp_path / "basins"
    root.mkdir()
    basin = Superbasin(root, 2, state_list=[state], config=cfg)
    results = StringIO()
    io.save_results_dat(
        results,
        {
            "termination_reason": 0,
            "potential_energy_saddle": -0.7,
            "potential_energy_reactant": -1.0,
            "potential_energy_product": -0.8,
            "prefactor_reactant_to_product": 1.0e12,
            "prefactor_product_to_reactant": 1.0e12,
            "barrier_reactant_to_product": 0.3,
            "displacement_saddle_distance": 0.2,
            "force_calls_saddle": 1,
            "force_calls_minimization": 1,
            "force_calls_prefactors": 1,
        },
    )
    results.seek(0)
    comm = _Comm(
        [
            {
                "name": "0_3",
                "results.dat": results,
                "reactant.con": _buf(reactant),
                "saddle.con": _buf(_atoms(2.7)),
                "product.con": _buf(_atoms(3.1)),
                "mode.dat": _buf(_atoms(2.7)),
            }
        ]
    )
    monkeypatch.chdir(tmp_path)
    (tmp_path / "pos.con").write_text("x\n")
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ClientMinModeExplorer(states, state, state, superbasin=basin, config=cfg)
    explorer.job_table.add_row({"state": 0, "wuid": 3, "type": "random"})
    registered = explorer.register_results()
    assert registered == 1
    displacement, mode, kind = explorer.generate_displacement()
    assert kind == "kdb"
    assert displacement is not None
    assert mode.shape[1] == 3


def test_precision_script_runs_both_tables(capsys):
    from eon.mcamc.test import main, random_chain

    rates, exits, _absorb = random_chain(3, 1, 1e-3)
    assert rates.shape == (3, 3)
    assert exits.shape == (3, 1)
    assert main() is None
    printed = capsys.readouterr().out
    assert "PRECISION TESTING" in printed
    assert "PERFORMANCE TESTING" in printed
