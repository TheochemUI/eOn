"""A rate sum that stops short of one still takes a KMC step."""

from __future__ import annotations

import textwrap

import numpy as np

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


def test_short_rate_sum_still_steps(tmp_path, monkeypatch):
    rates = [1.0] * 7
    running = 0.0
    for rate in rates:
        running += rate / sum(rates)
    monkeypatch.setattr(
        np.random,
        "random_sample",
        lambda: float(np.nextafter(running, 1.0)),
    )
    procs = {
        pid: {"rate": 1.0, "product": 1, "barrier": 0.2, "prefactor": 1.0}
        for pid in range(1, 8)
    }
    state0 = _State(0, procs)
    state1 = _State(1, {})
    cfg = _config(tmp_path)
    current, _previous, time, steps = kmc_step(
        state0,
        _States({0: state0, 1: state1}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 1
    assert time > 0.0
    assert current.number == 1
