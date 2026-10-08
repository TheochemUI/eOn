"""A confident state with no rate does not exit, and the default trajectory flag stays off."""

from __future__ import annotations

import textwrap

from eon.akmc import kmc_step
from eon.config import ConfigClass


class _State:
    def __init__(self, number, procs, energy=0.0, confidence=1.0):
        self.number = number
        self.procs = procs
        self.energy = energy
        self.confidence = confidence

    def get_process(self, pid):
        return self.procs[pid]

    def get_ratetable(self):
        return [
            (pid, proc["rate"], proc.get("prefactor", 1.0))
            for pid, proc in self.procs.items()
        ]

    def get_confidence(self, superbasin=None):
        return self.confidence

    def get_energy(self):
        return self.energy


class _States:
    def __init__(self, mapping):
        self.mapping = mapping

    def get_state(self, number):
        return self.mapping[int(number)]

    def get_product_state(self, reactant, proc_id):
        product = self.mapping[int(reactant)].procs[proc_id]["product"]
        return self.mapping[int(product)]


def _config(directory):
    path = directory / "config.ini"
    path.write_text(
        textwrap.dedent(
            f"""
            [Main]
            job = akmc
            temperature = 300

            [AKMC]
            confidence = 0.9
            confidence_scheme = old
            max_kmc_steps = 1

            [Paths]
            main_directory = {directory}
            results = {directory}

            [Debug]
            use_mean_time = true
            stop_criterion = 1e8
            """
        ).lstrip()
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def test_default_trajectory_flag_does_not_abort(tmp_path):
    cfg = _config(tmp_path)
    assert cfg.debug_target_trajectory in (False, "False")
    state0 = _State(
        0,
        {3: {"rate": 1.0e8, "product": 1, "barrier": 0.2, "prefactor": 1.0e12}},
    )
    state1 = _State(1, {})
    current, _previous, _time, steps = kmc_step(
        state0,
        _States({0: state0, 1: state1}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 1
    assert current.number == 1


def test_empty_rate_table_does_not_count_a_step(tmp_path):
    cfg = _config(tmp_path)
    state0 = _State(0, {})
    current, _previous, time, steps = kmc_step(
        state0,
        _States({0: state0}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 0
    assert time == 0.0
    assert current.number == 0
