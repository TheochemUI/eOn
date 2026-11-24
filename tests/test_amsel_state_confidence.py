"""[AKMC] confidence_scheme = amsel and the catalogue it hands to amsel."""

from __future__ import annotations

import sys
import textwrap
import types

import pytest

from eon.amsel_confidence import STATE_CONFIDENCE_SCHEMES, amsel_state_confidence
from eon.config import ConfigClass


class _State:
    def __init__(self, number, procs, counts, states):
        self.number, self.procs, self._counts, self.statelist = (
            number,
            procs,
            counts,
            states,
        )

    def load_process_table(self):
        pass

    def get_proc_random_count(self):
        return dict(self._counts)


class _List(list):
    def get_state(self, number):
        return self[number]

    def get_num_states(self):
        return len(self)


def _network():
    states = _List()
    # 0 -> 1 was found by searching 0; 1 holds its reverse but never found it.
    states.append(
        _State(
            0, {7: {"product": 1, "rate": 3.2e11, "barrier": 0.0895}}, {7: 5}, states
        )
    )
    states.append(
        _State(
            1,
            {
                1: {"product": -1, "rate": 3.4e9, "barrier": 0.383},
                2: {"product": 0, "rate": 2.4e11, "barrier": 0.105},
            },
            {1: 400},
            states,
        )
    )
    return states


def _config(directory, extra=""):
    path = directory / "config.ini"
    path.write_text(
        textwrap.dedent(
            f"""
            [Main]
            job = akmc
            temperature = 600

            [AKMC]
            confidence_scheme = amsel

            [Paths]
            main_directory = {directory}
            results = {directory}
            states = {directory / "states"}
            """
        ).lstrip()
        + textwrap.dedent(extra)
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def test_the_configuration_reaches_the_controller(tmp_path):
    cfg = _config(tmp_path)
    assert cfg.akmc_confidence_scheme == "amsel"
    assert cfg.amsel_state_confidence == "heuristic"
    cfg = _config(
        tmp_path,
        """
        [amsel]
        state_confidence = two-sided
        reverse_barrier_tolerance = 0.01
        """,
    )
    assert cfg.amsel_state_confidence == "two-sided"
    assert cfg.amsel_reverse_barrier_tolerance == 0.01


def test_far_side_saddles_and_own_finds_are_separated(monkeypatch):
    seen = []

    def state_confidence(scheme, **kw):
        seen.append((scheme, kw))
        return 0.25, (0, len(kw["implied"]))

    monkeypatch.setitem(
        sys.modules, "amsel", types.SimpleNamespace(state_confidence=state_confidence)
    )
    states = _network()
    value = amsel_state_confidence(
        states[1], scheme="two-sided", barrier_tolerance=0.02
    )
    assert value == 0.25
    scheme, kw = seen[-1]
    assert scheme == "two-sided"
    assert kw["hits"] == [400] and kw["barriers"] == [0.383]
    assert kw["implied"] == [(2.4e11, 0.105)]
    assert kw["network"] == (0, 2)


def test_amsel_scores_a_state_that_missed_its_entry_saddle():
    pytest.importorskip("amsel", reason="amsel is optional")
    import amsel

    if not hasattr(amsel, "state_confidence"):
        pytest.skip("amsel without state_confidence")
    states = _network()
    score = {
        s: amsel_state_confidence(states[1], scheme=s, barrier_tolerance=0.02)
        for s in STATE_CONFIDENCE_SCHEMES
    }
    assert score["good-turing"] == 1.0
    assert score["reverse-capped"] == pytest.approx(3.4e9 / (3.4e9 + 2.4e11))
    assert score["two-sided"] <= score["reverse-capped"]
