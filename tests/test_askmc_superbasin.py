"""AS-KMC names a missing product, and transition counting opens a basin."""

from __future__ import annotations

import textwrap

import pytest

from eon.askmc import ASKMC
from eon.config import ConfigClass
from eon.superbasinscheme import TransitionCounting


class _State:
    def __init__(self, number, path, procs):
        self.number = number
        self.path = path
        self.procs = procs

    def load_process_table(self):
        return None


class _States:
    def __init__(self, mapping):
        self.mapping = mapping
        self.connected = None

    def get_state(self, number):
        return self.mapping[int(number)]

    def connect_states(self, states):
        self.connected = list(states)

    def connect_state_sets(self, start_states, end_states):
        self.sets = (set(start_states), set(end_states))


def _config(directory):
    path = directory / "config.ini"
    path.write_text(
        textwrap.dedent(
            f"""
            [Main]
            job = akmc
            temperature = 300

            [Paths]
            main_directory = {directory}
            results = {directory}
            states = {directory / "states"}
            """
        ).lstrip()
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def _ask(directory):
    cfg = _config(directory)
    state_dir = directory / "state0"
    state_dir.mkdir()
    state = _State(
        0,
        state_dir,
        {3: {"barrier": 0.2, "rate": 1.0e6, "product": 1, "prefactor": 1.0e12,
             "saddle_energy": -0.5, "product_energy": -1.0,
             "product_prefactor": 1.0e12}},
    )
    ask = ASKMC(
        300.0 / 11604.5,
        _States({0: state}),
        0.9,
        1.5,
        2.0,
        True,
        False,
        False,
        directory,
        20.0,
        None,
    )
    return ask, state


def test_find_without_a_product_names_the_state(tmp_path):
    ask, state = _ask(tmp_path)
    with pytest.raises(KeyError, match="no process to state 9"):
        ask.get_process_id(state.procs, 9, "find")
    assert ask.get_process_id(state.procs, 1, "try") == 3
    assert ask.get_process_id(state.procs, 9, "try") is None


def test_rate_table_and_metadata_round_trip(tmp_path):
    ask, state = _ask(tmp_path)
    table = ask.get_ratetable(state)
    assert table == [(3, 1.0e6)]
    assert ask.get_askmc_metadata() == (0, 0)
    ask.save_askmc_metadata(2, 4)
    assert ask.get_askmc_metadata() == (2, 4)


def test_two_transitions_open_a_basin(tmp_path):
    cfg = _config(tmp_path)
    left_dir = tmp_path / "states" / "0"
    right_dir = tmp_path / "states" / "1"
    left_dir.mkdir(parents=True)
    right_dir.mkdir()
    left = _State(0, left_dir, {})
    right = _State(1, right_dir, {})
    states = _States({0: left, 1: right})
    scheme = TransitionCounting(
        tmp_path / "superbasin",
        states,
        300.0 / 11604.5,
        2,
        config=cfg,
    )
    scheme.register_transition(left, right)
    assert scheme.get_containing_superbasin(left) is None
    scheme.register_transition(left, right)
    basin = scheme.get_containing_superbasin(left)
    assert basin is not None
    assert basin.contains_state(right)
    assert {state.number for state in states.connected} == {0, 1}
