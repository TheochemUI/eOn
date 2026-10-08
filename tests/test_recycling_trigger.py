"""Superbasin recycling stops when the trigger process is missing or has no rate."""

import numpy as np

from eon.recycling import SB_Recycling
from eon.structure import Structure


def _atoms():
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([12.0, 12.0, 12.0])
    return atoms


class _State:
    def __init__(self, number, procs):
        self.number = number
        self.procs = procs
        self.reactant = _atoms()
        self.path = "."

    def get_reactant(self):
        return self.reactant

    def load_process_table(self):
        return None


def _recycler(procs):
    previous = _State(0, procs)
    current = _State(1, {})
    other = _State(2, {7: {"rate": 0.0, "barrier": 0.2, "product": 3}})
    recycler = SB_Recycling.__new__(SB_Recycling)
    recycler.in_progress = True
    recycler.previous_state = previous
    recycler.current_state = current
    recycler.move_distance = 5.0
    recycler.sb_states = [[other, None]]
    recycler.sb_state_nums = [[2, None]]
    recycler.states = None
    return recycler


def test_missing_trigger_process_stops_recycling():
    recycler = _recycler({})
    recycler.generate_corresponding_states()
    assert recycler.in_progress is False


def test_zero_rate_does_not_divide(tmp_path):
    recycler = _recycler({4: {"rate": 1.0e6, "barrier": 0.2, "product": 1}})
    other = recycler.sb_states[0][0]
    product = tmp_path / "product_7.con"
    # The zero-rate process is skipped before any product file is opened.
    other.path = tmp_path
    (tmp_path / "procdata").mkdir()
    recycler.generate_corresponding_states()
    assert recycler.in_progress is False
    assert not product.exists()
