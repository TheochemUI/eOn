"""A draw a few ulps past the rate sum still steps; a real mismatch does not."""

import numpy as np

from eon.akmc import kmc_step
from tests.test_kmc_rate_table import _State, _States, _config


class _Short(float):
    """Rate that adds as ``value`` and divides to ``fraction``.

    The running sum of the divided terms can then end well short of a
    draw that is still below 1.
    """

    def __new__(cls, value, fraction):
        obj = float.__new__(cls, value)
        obj.fraction = fraction
        return obj

    def __truediv__(self, other):
        return self.fraction


def _step(tmp_path, state, states, monkeypatch, draw):
    monkeypatch.setattr(np.random, "random_sample", lambda: draw)
    return kmc_step(state, states, 0.0, 0.025, None, config=_config(tmp_path))


def test_short_rate_sum_takes_the_last_process(tmp_path, monkeypatch):
    rates = [1.0] * 7
    ratesum = sum(rates, 0.0)
    running = 0.0
    for rate in rates:
        running += rate / ratesum
    draw = float(np.nextafter(running, 1.0))
    procs = {
        pid: {"rate": 1.0, "product": 10 + pid, "barrier": 0.2, "prefactor": 1.0}
        for pid in range(1, 8)
    }
    states = {0: _State(0, procs)}
    for pid in range(1, 8):
        states[10 + pid] = _State(10 + pid, {})
    landed, _previous, time, steps = _step(
        tmp_path, states[0], _States(states), monkeypatch, draw
    )
    assert steps == 1
    assert time > 0.0
    assert landed.number == 17


def test_draw_past_the_sum_does_not_count_a_step(tmp_path, monkeypatch):
    procs = {
        1: {"rate": _Short(1.0, 0.1), "product": 1, "barrier": 0.2, "prefactor": 1.0},
        2: {"rate": _Short(1.0, 0.1), "product": 1, "barrier": 0.2, "prefactor": 1.0},
    }
    state0 = _State(0, procs)
    state1 = _State(1, {})
    landed, _previous, time, steps = _step(
        tmp_path,
        state0,
        _States({0: state0, 1: state1}),
        monkeypatch,
        0.5,
    )
    assert steps == 0
    assert time == 0.0
    assert landed.number == 0


def test_zero_rate_table_does_not_count_a_step(tmp_path, monkeypatch):
    def _no_draw():
        raise AssertionError("a zero-rate table must not draw")

    monkeypatch.setattr(np.random, "random_sample", _no_draw)
    state0 = _State(
        0,
        {3: {"rate": 0.0, "product": 1, "barrier": 0.2, "prefactor": 1.0}},
    )
    state1 = _State(1, {})
    landed, _previous, time, steps = kmc_step(
        state0,
        _States({0: state0, 1: state1}),
        0.0,
        0.025,
        None,
        config=_config(tmp_path),
    )
    assert steps == 0
    assert time == 0.0
    assert landed.number == 0
