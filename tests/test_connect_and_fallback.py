"""Unlinked processes, a rejected basin, and a similar recycled saddle."""

import sys
from pathlib import Path
from types import ModuleType

import numpy as np

from eon import fileio as io
from eon.akmc import kmc_step
from eon.recycling import SB_Recycling
from eon.superbasin import Superbasin
from tests.test_library_bodies import _atoms, _states


def test_unlinked_processes_join_the_other_state(tmp_path):
    cfg, states, state, product, proc, reactant = _states(tmp_path)
    forward = state.get_process(proc)
    extra = product.allocate_process_id(b"join", b"back")
    product.append_process_table(
        id=extra,
        saddle_energy=forward["saddle_energy"],
        prefactor=forward["prefactor"],
        product=-1,
        product_energy=state.get_energy(),
        product_prefactor=forward["product_prefactor"],
        barrier=forward["barrier"],
        rate=forward["rate"],
        repeats=0,
    )
    io.savecon(product.proc_product_path(extra), reactant)
    io.savecon(product.proc_saddle_path(extra), _atoms(2.7))
    io.savecon(product.proc_reactant_path(extra), product.get_reactant())
    io.save_mode(product.proc_mode_path(extra), np.zeros((2, 3)))
    Path(product.proc_results_path(extra)).write_text("0 termination_reason\n")
    states.connect_state_sets([product], [state])
    assert product.get_process(extra)["product"] == state.number


def test_a_rejected_basin_falls_back_to_one_state(tmp_path, monkeypatch):
    cfg, states, state, _product, _proc, _reactant = _states(tmp_path)
    for _ in range(100):
        state.inc_repeats()
    cfg.sb_on = True
    cfg.akmc_max_kmc_steps = 1
    cfg.amsel_discover_decide = True
    cfg.askmc_on = False
    root = tmp_path / "basins"
    root.mkdir()
    basin = Superbasin(root, 1, state_list=[state], config=cfg)
    module = ModuleType("amsel")

    def discover_decide_status(*_args, **_kwargs):
        return ("rejected", [], [], [])

    module.discover_decide_status = discover_decide_status
    monkeypatch.setitem(sys.modules, "amsel", module)

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


def test_similar_rates_compare_the_product_geometry(tmp_path):
    cfg, states, previous, current, proc, _reactant = _states(tmp_path)
    previous.load_process_table()
    current.load_process_table()
    previous.procs[proc]["rate"] = 1.0
    previous.procs[proc]["barrier"] = 0.2
    for row in current.procs.values():
        row["rate"] = 1.0
        row["barrier"] = 0.2
    root = tmp_path / "basins"
    root.mkdir()
    basin = Superbasin(root, 2, state_list=[previous, current], config=cfg)
    recycling = SB_Recycling(
        states,
        previous,
        current,
        0.2,
        False,
        tmp_path / "recycle",
        "mcacm",
        type("S", (), {"get_containing_superbasin": lambda self, state: basin})(),
    )
    assert recycling.sb_state_nums[0][0] == previous.number
