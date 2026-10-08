"""State graphs name the fastest process and a shortest path."""

from types import SimpleNamespace

import numpy as np
import pytest

from eon.movie import (
    Graph,
    fastest_path,
    get_fastest_process_id,
    get_fastest_process_rate,
    make_graph,
)
from eon.structure import Structure


class _State:
    def __init__(self, number, processes):
        self.number = number
        self.processes = processes
        self.reactant = Structure(1)
        self.reactant.r[0] = [float(number), 0.0, 0.0]
        self.saddles = {}

    def get_process_table(self):
        return self.processes

    def get_reactant(self):
        return self.reactant

    def get_process_saddle(self, pid):
        return self.saddles[pid]


def test_fastest_process_is_the_larger_rate():
    left = _State(
        0,
        {
            1: {"product": 1, "rate": 1.0},
            2: {"product": 1, "rate": 4.0},
            3: {"product": 2, "rate": 9.0},
        },
    )
    right = _State(1, {})
    assert get_fastest_process_id(left, right) == 2
    assert get_fastest_process_rate(left, right) == 4.0


def test_missing_process_names_the_states():
    left = _State(0, {1: {"product": 2, "rate": 1.0}})
    right = _State(1, {})
    with pytest.raises(ValueError, match="no process from state 0 to 1"):
        get_fastest_process_id(left, right)


def test_shortest_path_visits_the_bridge():
    a = _State(0, {1: {"product": 1, "rate": 2.0}})
    b = _State(1, {2: {"product": 2, "rate": 3.0}})
    c = _State(2, {})
    a.saddles[1] = Structure(1)
    b.saddles[2] = Structure(1)
    states = SimpleNamespace(
        get_num_states=lambda: 3,
        get_state=lambda n: (a, b, c)[n],
    )
    graph = make_graph(states)
    assert "0 -- 1" in graph.dot()
    assert "1 -- 2" in graph.dot()
    path = graph.shortest_path(a, c)
    assert [node.number for node in path] == [0, 1, 2]
    frames = fastest_path("unused", states)
    assert len(frames) == 5
    assert np.isfinite(frames[0].r).all()
