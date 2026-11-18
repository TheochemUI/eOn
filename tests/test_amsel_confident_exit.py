"""A superbasin exit waits for every transient state to reach the configured confidence.

amsel reads the edges it is given as each transient state's complete
outgoing catalog. A state with only some exits found has a smaller escape
rate than the real one, so the exit time is too long and the missing exits
carry no weight. The network below is the smallest case. Two states
exchange population quickly. State 0 has its exit, state 1 has not yet had
its exit found.

The confidence is eOn's own: ``AKMCState.get_confidence`` on states that
carry a real process table, compared with ``config.akmc_confidence``.
"""
from __future__ import annotations

import math
import textwrap
from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest

from eon import akmc as akmc_module
from eon import amsel_superbasin_gate as gate
from eon.akmcstate import AKMCState
from eon.config import ConfigClass

CONFIDENCE = 0.99
IN_BASIN_BARRIER = 0.10
ENTRY_EXIT_BARRIER = 0.60
PARTNER_EXIT_BARRIER = 0.65
SWAP_RATE = 1.0e6
ENTRY_EXIT_RATE = 0.6
PARTNER_EXIT_RATE = 1.257
SEARCHED = 1000
UNSEARCHED = 20


def _config(directory, confidence=CONFIDENCE, discover=True):
    text = textwrap.dedent(
        f"""
        [Main]
        job = akmc
        temperature = 300
        random_seed = -1

        [AKMC]
        confidence = {confidence}
        confidence_scheme = old
        max_kmc_steps = 1

        [Paths]
        main_directory = {directory}
        results = {directory}

        [Debug]
        use_mean_time = true
        stop_criterion = 1e8

        [amsel]
        discover_decide = {"true" if discover else "false"}
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


class _Info:
    def __init__(self, values):
        self.values = dict(values)

    def get(self, section, key, default=None):
        return self.values.get(key, default)

    def set(self, section, key, value):
        self.values[key] = value


class _EonState(AKMCState):
    """A state with an in-memory process table and eOn's own bookkeeping."""

    def __init__(self, number, procs, config, repeats):
        self.config = config
        self.number = number
        self.procs = procs
        self.statelist = SimpleNamespace(
            kT=config.main_temperature / 11604.5,
            thermal_window=config.akmc_thermal_window,
        )
        self.info = _Info({"repeats": repeats})

    def load_process_table(self, force=False):
        return None

    def get_energy(self):
        return 0.0


class _States:
    def __init__(self, mapping):
        self.mapping = mapping

    def get_state(self, number):
        return self.mapping[int(number)]

    def get_product_state(self, reactant, proc_id):
        product = self.mapping[int(reactant)].procs[proc_id]["product"]
        return self.mapping[int(product)]


def _proc(cfg, product, rate, barrier):
    kT = cfg.main_temperature / 11604.5
    return {
        "product": product,
        "rate": rate,
        "barrier": barrier,
        "prefactor": rate * math.exp(barrier / kT),
    }


def _entry_procs(cfg):
    return {
        0: _proc(cfg, 1, SWAP_RATE, IN_BASIN_BARRIER),
        1: _proc(cfg, 2, ENTRY_EXIT_RATE, ENTRY_EXIT_BARRIER),
    }


def _partner_partial(cfg):
    return {0: _proc(cfg, 0, SWAP_RATE, IN_BASIN_BARRIER)}


def _partner_complete(cfg):
    procs = _partner_partial(cfg)
    procs[1] = _proc(cfg, 3, PARTNER_EXIT_RATE, PARTNER_EXIT_BARRIER)
    return procs


def _network(cfg, partner_procs, partner_repeats, entry_repeats=SEARCHED):
    entry = _EonState(0, _entry_procs(cfg), cfg, entry_repeats)
    partner = _EonState(1, partner_procs, cfg, partner_repeats)
    mapping = {
        0: entry,
        1: partner,
        2: _EonState(2, {}, cfg, 0),
        3: _EonState(3, {}, cfg, 0),
    }
    return entry, partner, _States(mapping)


def _mean_exit_time(swap_01, swap_10, exit_0, exit_1):
    """First-passage time from state 0 for two states, by the linear solve."""
    generator = np.array(
        [
            [swap_01 + exit_0, -swap_01],
            [-swap_10, swap_10 + exit_1],
        ]
    )
    return float(np.linalg.solve(generator, np.ones(2))[0])


def _outcome(entry, states, cfg, seed=5):
    """('exit', AmselExit) or ('needs_search', error) from the real gate."""
    rng = np.random.default_rng(seed)
    try:
        return "exit", gate.amsel_discover_exit(
            entry, states.get_state, cfg, uniform=rng.random
        )
    except gate.AmselNeedsSearch as need:
        return "needs_search", need


def test_partial_catalog_tau_is_far_from_the_true_tau():
    partial = _mean_exit_time(SWAP_RATE, SWAP_RATE, ENTRY_EXIT_RATE, 0.0)
    true = _mean_exit_time(SWAP_RATE, SWAP_RATE, ENTRY_EXIT_RATE, PARTNER_EXIT_RATE)
    assert partial == pytest.approx(3.33, rel=2e-3)
    assert true == pytest.approx(1.077, rel=1e-3)
    assert partial / true > 3.0


def test_unconfident_partner_blocks_the_exit(tmp_path):
    cfg = _config(tmp_path)
    entry, partner, states = _network(cfg, _partner_partial(cfg), UNSEARCHED)
    assert entry.get_confidence() >= cfg.akmc_confidence
    assert partner.get_confidence() < cfg.akmc_confidence
    kind, got = _outcome(entry, states, cfg)
    assert kind == "needs_search", (
        "exit taken with tau=%r from a catalog missing state 1's exit" % got.mean_time
    )
    assert [number for number, _confidence in got.states] == [1]
    assert got.states[0][1] == pytest.approx(partner.get_confidence())


def test_no_draw_when_the_gate_sends_states_back_to_search(tmp_path):
    cfg = _config(tmp_path)
    entry, _partner, states = _network(cfg, _partner_partial(cfg), UNSEARCHED)
    rng = np.random.default_rng(31)
    before = rng.bit_generator.state
    try:
        got = gate.amsel_discover_exit(
            entry, states.get_state, cfg, uniform=rng.random
        )
    except gate.AmselNeedsSearch:
        got = None
    assert got is None
    assert rng.bit_generator.state == before


def test_exit_time_matches_the_analytic_tau_once_every_state_is_confident(tmp_path):
    cfg = _config(tmp_path)
    entry, partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    assert partner.get_confidence() >= cfg.akmc_confidence
    kind, got = _outcome(entry, states, cfg)
    assert kind == "exit"
    assert got.kernel == "mrm"
    true = _mean_exit_time(SWAP_RATE, SWAP_RATE, ENTRY_EXIT_RATE, PARTNER_EXIT_RATE)
    assert got.mean_time == pytest.approx(true, rel=1e-6)
    assert got.mean_time == pytest.approx(1.077, rel=1e-3)


def test_confidence_equal_to_the_threshold_is_enough(tmp_path):
    cfg = _config(tmp_path)
    entry, partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    cfg.akmc_confidence = min(entry.get_confidence(), partner.get_confidence())
    kind, _got = _outcome(entry, states, cfg)
    assert kind == "exit"
    cfg.akmc_confidence = math.nextafter(cfg.akmc_confidence, math.inf)
    kind, got = _outcome(entry, states, cfg)
    assert kind == "needs_search"
    below = {
        number
        for number, state in ((0, entry), (1, partner))
        if state.get_confidence() < cfg.akmc_confidence
    }
    assert below
    assert {number for number, _confidence in got.states} == below


class _NaNPartner(_EonState):
    def get_confidence(self, superbasin=None):
        return float("nan")


class _NonePartner(_EonState):
    def get_confidence(self, superbasin=None):
        return None


class _TextPartner(_EonState):
    def get_confidence(self, superbasin=None):
        return "high"


class _InfPartner(_EonState):
    def get_confidence(self, superbasin=None):
        return float("inf")


class _NoConfidencePartner:
    def __init__(self, number, procs):
        self.number = number
        self.procs = procs

    def get_process_table(self):
        return dict(self.procs)


@pytest.mark.parametrize("kind", ["nan", "none", "text", "inf", "absent"])
def test_missing_or_nan_confidence_needs_search(tmp_path, kind):
    cfg = _config(tmp_path)
    entry, _partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    procs = _partner_complete(cfg)
    if kind == "absent":
        replacement = _NoConfidencePartner(1, procs)
    else:
        cls = {
            "nan": _NaNPartner,
            "none": _NonePartner,
            "text": _TextPartner,
            "inf": _InfPartner,
        }[kind]
        replacement = cls(1, procs, cfg, SEARCHED)
    states.mapping[1] = replacement
    outcome, got = _outcome(entry, states, cfg)
    assert outcome == "needs_search", "exit taken with %s confidence" % kind
    (number, confidence), = got.states
    assert number == 1
    assert confidence is None


def test_a_missing_configured_confidence_is_an_error(tmp_path):
    cfg = _config(tmp_path)
    entry, _partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    bare = SimpleNamespace(
        amsel_e_min_init=0.5,
        amsel_e_min_step=0.05,
        amsel_e_min_floor=0.05,
        amsel_cv_threshold=10.0,
        amsel_on_error="fallback_single",
        debug_use_mean_time=True,
    )
    with pytest.raises((AttributeError, TypeError, ValueError)):
        gate.amsel_discover_exit(
            entry, states.get_state, bare, uniform=np.random.default_rng(1).random
        )


@pytest.mark.parametrize("threshold", [float("nan"), None, "0.99"])
def test_a_malformed_configured_confidence_is_an_error(tmp_path, threshold):
    cfg = _config(tmp_path)
    entry, _partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    cfg.akmc_confidence = threshold
    with pytest.raises((TypeError, ValueError)):
        gate.amsel_discover_exit(
            entry, states.get_state, cfg, uniform=np.random.default_rng(1).random
        )


def test_kmc_step_does_not_hop_until_every_state_is_confident(tmp_path):
    cfg = _config(tmp_path)
    entry, partner, states = _network(cfg, _partner_partial(cfg), UNSEARCHED)
    np.random.seed(77)
    current, _previous, time, steps = akmc_module.kmc_step(
        entry, states, 0.0, 0.025, None, config=cfg
    )
    assert steps == 0
    assert time == 0.0
    assert current is entry
    assert not (tmp_path / "dynamics.txt").exists()

    partner.procs = _partner_complete(cfg)
    partner.info.values["repeats"] = SEARCHED
    np.random.seed(77)
    current, _previous, time, steps = akmc_module.kmc_step(
        entry, states, 0.0, 0.025, None, config=cfg
    )
    assert steps == 1
    true = _mean_exit_time(SWAP_RATE, SWAP_RATE, ENTRY_EXIT_RATE, PARTNER_EXIT_RATE)
    assert time == pytest.approx(true, rel=1e-6)


class _RecordingExplorer:
    calls: list = []

    def __init__(self, states, previous_state, state, superbasin=None, config=None):
        self.record = (state, superbasin)

    def explore(self):
        _RecordingExplorer.calls.append(self.record)


def _drive_step(monkeypatch, tmp_path, cfg, states, entry):
    _RecordingExplorer.calls = []
    monkeypatch.setattr(
        akmc_module,
        "get_akmc_metadata",
        lambda config: (entry.number, 0.0, -1, True, config.main_temperature),
    )
    monkeypatch.setattr(akmc_module, "get_statelist", lambda kT, config=None: states)
    monkeypatch.setattr(
        akmc_module.explorer, "get_minmodexplorer", lambda config: _RecordingExplorer
    )
    monkeypatch.setattr(akmc_module.io, "write_info_txt", lambda config, parser: None)
    monkeypatch.setattr(akmc_module.io, "save_prng_state", lambda path: None)
    np.random.seed(5)
    return akmc_module.akmc(cfg)


def test_unconfident_partner_is_sent_back_to_the_search_loop(monkeypatch, tmp_path):
    cfg = _config(tmp_path)
    entry, partner, states = _network(cfg, _partner_partial(cfg), UNSEARCHED)
    steps = _drive_step(monkeypatch, tmp_path, cfg, states, entry)
    assert steps == 0
    (searched, superbasin), = _RecordingExplorer.calls
    assert searched is partner
    assert superbasin is None


def test_confident_basin_is_searched_from_the_current_state(monkeypatch, tmp_path):
    cfg = _config(tmp_path)
    entry, _partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    _drive_step(monkeypatch, tmp_path, cfg, states, entry)
    (searched, _superbasin), = _RecordingExplorer.calls
    assert searched is entry


class _Table:
    def __init__(self, number, procs):
        self.number = number
        self.procs = procs

    def get_process_table(self):
        return self.procs


def _graph_basin(state_numbers, cfg):
    states = {
        0: _Table(0, _entry_procs(cfg)),
        1: _Table(1, _partner_complete(cfg)),
    }
    return SimpleNamespace(
        state_numbers=list(state_numbers),
        state_dict=states,
        states=[states[n] for n in state_numbers],
    )


def test_a_revisited_state_adds_its_processes_once(tmp_path):
    cfg = _config(tmp_path)
    once = gate.build_graph_from_superbasin(_graph_basin([0, 1], cfg), 0)
    again = gate.build_graph_from_superbasin(_graph_basin([0, 1, 0, 1, 0], cfg), 0)
    assert again == once
    candidates, rates, barriers = again
    assert candidates == [0, 1]
    assert len(barriers) == len(rates)
    assert max(Counter((a, b, r) for a, b, r in rates).values()) == 1


def test_a_revisited_state_gives_amsel_one_catalog(tmp_path):
    cfg = _config(tmp_path)
    entry = SimpleNamespace(number=0)
    once = gate.discover_decide_for_superbasin(_graph_basin([0, 1], cfg), entry)
    again = gate.discover_decide_for_superbasin(_graph_basin([0, 1, 0], cfg), entry)
    assert once["status"] in ("accepted", "retightened", "split_required")
    assert again["status"] == once["status"]
    assert sorted(again["raw"][3]) == sorted(once["raw"][3])


def test_a_resolved_again_state_is_not_counted_twice(tmp_path):
    cfg = _config(tmp_path)
    entry, _partner, states = _network(cfg, _partner_complete(cfg), SEARCHED)
    fresh = {
        0: lambda: _EonState(0, _entry_procs(cfg), cfg, SEARCHED),
        1: lambda: _EonState(1, _partner_complete(cfg), cfg, SEARCHED),
        2: lambda: _EonState(2, {}, cfg, 0),
        3: lambda: _EonState(3, {}, cfg, 0),
    }
    rng = np.random.default_rng(5)
    got = gate.amsel_discover_exit(
        entry, lambda number: fresh[int(number)](), cfg, uniform=rng.random
    )
    true = _mean_exit_time(SWAP_RATE, SWAP_RATE, ENTRY_EXIT_RATE, PARTNER_EXIT_RATE)
    assert got.mean_time == pytest.approx(true, rel=1e-6)
