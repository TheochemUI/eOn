"""A draw that matches no rate does not advance the KMC step count."""

import numpy as np

from tests.test_kmc_rate_table import _State, _States, _config
from eon.akmc import kmc_step


def test_missed_rate_draw_does_not_count_a_step(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    state0 = _State(
        0,
        {3: {"rate": 1.0, "product": 0, "barrier": 0.2, "prefactor": 1.0e12}},
    )
    # random_sample is half-open at 1, so a value of 1 matches no row
    # under the strict cumulative test.
    monkeypatch.setattr(np.random, "random_sample", lambda: 1.0)
    state, _previous, time, steps = kmc_step(
        state0,
        _States({0: state0}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 0
    assert time == 0.0
    assert state.number == 0
