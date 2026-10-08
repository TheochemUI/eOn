"""AS-KMC substitutes a modified rate and counts a first visit."""

from pathlib import Path

import numpy as np

from eon.askmc import ASKMC


class _State:
    def __init__(self, number, path, procs):
        self.number = number
        self.path = path
        self.procs = procs

    def load_process_table(self):
        return None


def _ask(path):
    return ASKMC(
        300.0 / 11604.5,
        None,
        0.9,
        1.5,
        2.0,
        False,
        False,
        False,
        path,
        20.0,
        None,
    )


def _proc():
    return {
        3: {
            "saddle_energy": -0.5,
            "prefactor": 1.0e12,
            "product": 1,
            "product_energy": -1.0,
            "product_prefactor": 1.0e12,
            "barrier": 0.2,
            "rate": 1.0e6,
        }
    }


def test_modified_rate_replaces_the_table_rate(tmp_path):
    state_dir = tmp_path / "state0"
    state_dir.mkdir()
    state = _State(0, state_dir, _proc())
    ask = _ask(tmp_path)
    compiled = ask.compile_process_table(state)
    assert compiled[3]["rate"] == 1.0e6
    ask.append_modified_process_table(
        state, 3, -0.5, 1.0e12, 1, -1.0, 1.0e12, 0.2, 4.0, 2
    )
    compiled = ask.compile_process_table(state)
    assert compiled[3]["rate"] == 4.0
    assert compiled[3]["view_count"] == 2
    table = ask.get_ratetable(state)
    assert table == [(3, 4.0)]


def test_first_visit_records_a_view_count(tmp_path):
    left_dir = tmp_path / "left"
    right_dir = tmp_path / "right"
    left_dir.mkdir()
    right_dir.mkdir()
    left = _State(0, left_dir, _proc())
    right = _State(1, right_dir, {})
    ask = _ask(tmp_path)
    ask.register_transition(left, right)
    modified = ask.get_modified_process_table(left)
    assert modified[3]["view_count"] == 1
    assert modified[3]["product"] == 1
    ask.register_transition(left, left)
    assert ask.get_modified_process_table(left)[3]["view_count"] == 1


def test_edge_helpers_use_numpy_equality(tmp_path):
    ask = _ask(tmp_path)
    ask.edgelist = [(0, 1), (1, 2), (0, 1)]
    assert ask.edgelist_to_statelist() == [0, 1, 2]
    assert ask.is_equal(1.0, 1.0 + 1.0e-8) is True
    assert ask.is_equal(1.0, 1.1) is False
    assert ask.in_array([0, 1], np.array([[2, 3], [0, 1]])) == 1
    assert ask.in_array([4, 5], np.array([[2, 3]])) == 0
