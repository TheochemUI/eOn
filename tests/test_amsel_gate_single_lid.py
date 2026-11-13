"""One lid runs from discover through the exit, and None does not draw."""
from __future__ import annotations

import textwrap

import numpy as np
import pytest

from eon.akmc import kmc_step
from eon.amsel_superbasin_gate import amsel_discover_exit
from eon.config import ConfigClass


class _State:
    def __init__(self, number, procs):
        self.number = number
        self.procs = procs

    def get_process_table(self):
        return dict(self.procs)

    def get_process(self, pid):
        return self.procs[int(pid)]

    def get_ratetable(self, superbasin=None):
        return [
            (pid, proc["rate"], proc.get("prefactor", 0.0))
            for pid, proc in self.procs.items()
        ]

    def get_confidence(self, superbasin=None):
        return 1.0

    def get_energy(self):
        return 0.0


class _States:
    def __init__(self, mapping):
        self.mapping = mapping

    def get_state(self, number):
        return self.mapping[int(number)]

    def get_product_state(self, reactant, proc_id):
        product = self.mapping[int(reactant)].procs[int(proc_id)]["product"]
        return self.mapping[int(product)]


class _Gate:
    amsel_e_min_init = 0.5
    amsel_e_min_step = 0.05
    amsel_e_min_floor = 0.05
    amsel_cv_threshold = 10.0
    amsel_on_error = "fallback_single"
    debug_use_mean_time = False


def _one_reactant():
    entry = _State(
        0,
        {
            0: {"product": 1, "rate": 5.0e9, "barrier": 0.05165, "prefactor": 1.0e13},
            1: {"product": 2, "rate": 3.0e9, "barrier": 0.08, "prefactor": 1.0e13},
            2: {"product": 3, "rate": 2.0e9, "barrier": 0.12, "prefactor": 1.0e13},
        },
    )
    states = _States({0: entry, 1: _State(1, {}), 2: _State(2, {}), 3: _State(3, {})})
    return entry, states


def _two_minima(escape):
    left = _State(
        0,
        {
            0: {"product": 1, "rate": 1.0e11, "barrier": 0.032, "prefactor": 1.0e13},
            1: {"product": 2, "rate": 1.0, "barrier": escape, "prefactor": 1.0e13},
        },
    )
    right = _State(1, {0: {"product": 0, "rate": 1.0e11, "barrier": 0.032, "prefactor": 1.0e13}})
    states = {0: left, 1: right, 2: _State(2, {})}
    return left, states


def _config(directory, discover):
    flag = "true" if discover else "false"
    text = textwrap.dedent(
        f"""
        [Main]
        job = akmc
        temperature = 300
        random_seed = -1

        [AKMC]
        confidence = 0.0
        confidence_scheme = old
        max_kmc_steps = 1

        [Paths]
        main_directory = {directory}
        results = {directory}

        [Debug]
        use_mean_time = false
        stop_criterion = 1e8

        [amsel]
        discover_decide = {flag}
        e_min_init = 0.5
        e_min_step = 0.05
        e_min_floor = 0.05
        """
    ).lstrip()
    path = directory / "config.ini"
    path.write_text(text)
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def test_one_reactant_below_the_lid_does_not_draw():
    entry, states = _one_reactant()
    rng = np.random.default_rng(1)
    before = rng.bit_generator.state
    got = amsel_discover_exit(
        entry, states.get_state, _Gate(), uniform=rng.random
    )
    assert got is None
    assert rng.bit_generator.state == before


def test_one_reactant_step_matches_the_plain_clock(tmp_path):
    entry, states = _one_reactant()
    plain_dir = tmp_path / "plain"
    amsel_dir = tmp_path / "amsel"
    plain_dir.mkdir()
    amsel_dir.mkdir()
    plain = _config(plain_dir, discover=False)
    amsel = _config(amsel_dir, discover=True)
    np.random.seed(9000)
    _current, _previous, plain_time, plain_steps = kmc_step(
        entry, states, 0.0, 0.025, None, config=plain
    )
    np.random.seed(9000)
    _current, _previous, amsel_time, amsel_steps = kmc_step(
        entry, states, 0.0, 0.025, None, config=amsel
    )
    assert plain_steps == 1 and amsel_steps == 1
    assert amsel_time == plain_time


def test_two_minima_escape_above_the_lid_is_fpta():
    entry, states = _two_minima(0.69)
    got = amsel_discover_exit(
        entry, lambda number: states[int(number)], _Gate(), uniform=np.random.default_rng(7).random
    )
    assert got is not None
    assert got.kernel == "fpta"
    assert got.time_is_sample
    assert got.barrier == pytest.approx(0.69)


def test_escape_below_the_initial_lid_uses_the_settled_lid():
    entry, states = _two_minima(0.20)
    got = amsel_discover_exit(
        entry, lambda number: states[int(number)], _Gate(), uniform=np.random.default_rng(3).random
    )
    assert got is not None
    assert got.kernel == "fpta"
    assert got.barrier == pytest.approx(0.20)
