"""A zero-probability superbasin exit is not drawn."""

from types import SimpleNamespace

import numpy as np

from eon.superbasin import Superbasin


class _State:
    def __init__(self, number, procs):
        self.number = number
        self.procs = procs

    def get_process_table(self):
        return self.procs


def test_zero_exit_probability_is_not_drawn(monkeypatch):
    basin = Superbasin.__new__(Superbasin)
    basin.config = SimpleNamespace(
        amsel_discover_decide=False,
        sb_amsel_discover_decide=False,
    )
    basin.id = 1
    basin.state_numbers = [0]
    state = _State(
        0,
        {
            1: {"product": 2, "rate": 0.0},
            2: {"product": 3, "rate": 1.0},
        },
    )
    basin.state_dict = {0: state}
    basin.states = [state]
    monkeypatch.setattr(
        "eon.superbasin.mcamc",
        lambda q_matrix, exit_matrix, costs: ([1.0], np.array([[0.0, 1.0]]), 0.0),
    )
    monkeypatch.setattr(np.random, "sample", lambda: 0.0)

    def product_of(_number, proc_id):
        return SimpleNamespace(number=9, proc_id=proc_id)

    _mean, _exit_state, product, proc_id, basin_id = basin.step(state, product_of)
    assert proc_id == 2
    assert product.proc_id == 2
    assert basin_id == 1
