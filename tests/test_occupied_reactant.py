"""A superbasin exit is registered from the state the trajectory occupies."""
from __future__ import annotations

import textwrap

from eon.akmc import kmc_step
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
        return [(pid, proc["rate"], proc.get("prefactor", 0.0)) for pid, proc in self.procs.items()]

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


class _Exit:
    mean_time = 1.0
    exit_number = 0
    proc_id = 32
    time_is_sample = True


class _Scheme:
    def __init__(self):
        self.seen = []

    def get_containing_superbasin(self, state):
        return None

    def register_transition(self, start, end):
        self.seen.append((start.number, end.number))

    def write_data(self):
        return None


def _config(directory):
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
        discover_decide = true
        e_min_init = 0.5
        e_min_step = 0.05
        e_min_floor = 0.05
        """
    ).lstrip()
    path = directory / "config.ini"
    path.write_text(text)
    cfg = ConfigClass()
    cfg.init(str(path))
    cfg.sb_on = True
    return cfg


def test_exit_on_another_minimum_registers_the_occupied_state(tmp_path, monkeypatch):
    """The process may belong to the exit minimum. The transition starts at the occupied state."""
    occupied = _State(11, {0: {"product": 0, "rate": 1.0e12, "barrier": 0.01, "prefactor": 1.0e13}})
    owner = _State(0, {32: {"product": 33, "rate": 1.0, "barrier": 0.5, "prefactor": 1.0e13}})
    product = _State(33, {})
    states = _States({11: occupied, 0: owner, 33: product})
    monkeypatch.setattr("eon.amsel_superbasin_gate.amsel_discover_exit", lambda *_args, **_kwargs: _Exit())
    scheme = _Scheme()
    kmc_step(occupied, states, 0.0, 0.025, scheme, config=_config(tmp_path))
    assert scheme.seen == [(11, 33)]
