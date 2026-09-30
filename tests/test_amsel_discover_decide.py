"""discover_decide runs from [amsel] alone and exits through MRM or FPTA."""
from __future__ import annotations

import logging
import textwrap

import pytest

from eon.akmc import kmc_step
from eon.config import ConfigClass


class _State:
    def __init__(self, number, procs, energy=0.0, confidence=1.0):
        self.number = number
        self.procs = procs
        self.energy = energy
        self.confidence = confidence

    def get_process_table(self):
        return self.procs

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
        if int(number) not in self.mapping:
            raise KeyError(number)
        return self.mapping[int(number)]

    def get_product_state(self, reactant, proc_id):
        product = self.mapping[int(reactant)].procs[proc_id]["product"]
        return self.mapping[int(product)]


def _config(directory, extra, use_mean_time="true", confidence=0.0, max_kmc_steps=1):
    text = textwrap.dedent(
        f"""
        [Main]
        job = akmc
        temperature = 300

        [AKMC]
        confidence = {confidence}
        confidence_scheme = old
        max_kmc_steps = {max_kmc_steps}

        [Paths]
        main_directory = {directory}
        results = {directory}

        [Debug]
        use_mean_time = {use_mean_time}
        stop_criterion = 1e8
        """
    ).lstrip()
    text += textwrap.dedent(extra).strip() + "\n"
    path = directory / "config.ini"
    path.write_text(text)
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def _isomer1():
    # Isomer 1: 0.20 eV back to state 0, 0.30 eV flip out to state 2.
    state1 = _State(
        1,
        {
            7: {"rate": 1.0e8, "product": 0, "barrier": 0.20, "prefactor": 1.0e12},
            8: {"rate": 1.0e2, "product": 2, "barrier": 0.30, "prefactor": 1.0e12},
        },
    )
    state0 = _State(
        0,
        {
            3: {"rate": 1.0e8, "product": 1, "barrier": 0.20, "prefactor": 1.0e12},
        },
    )
    state2 = _State(2, {})
    return state1, _States({0: state0, 1: state1, 2: state2})


def test_status_vocabulary_rejects_the_disagreed_pair():
    from eon.amsel_superbasin_gate import status_pair

    assert status_pair("unavailable", False)
    assert status_pair("fallback_single", True)
    assert not status_pair("unavailable", True)
    assert not status_pair("fallback_single", False)
    assert status_pair("accepted", True)
    assert status_pair("rejected_no_metastable_basin", True)


def test_ini_without_use_mcamc_logs_discover_decide_status(tmp_path, monkeypatch, caplog):
    """An [amsel] ini with discover_decide and no use_mcamc logs the status."""
    monkeypatch.setitem(__import__("sys").modules, "amsel", None)
    cfg = _config(
        tmp_path,
        """
        [amsel]
        discover_decide = true
        """,
    )
    assert cfg.amsel_discover_decide is True
    assert cfg.sb_on is False
    state0 = _State(
        0,
        {1: {"rate": 1.0, "product": 1, "barrier": 0.40, "prefactor": 1.0}},
    )
    state1 = _State(1, {})
    caplog.set_level(logging.INFO, logger="superbasin.amsel_gate")
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
    assert any("amsel discover_decide status" in message for message in caplog.messages)
    assert any("status=unavailable" in message for message in caplog.messages)


def test_isomer1_exit_comes_from_mrm(tmp_path, monkeypatch, caplog):
    """0.20 eV stays in the basin. The 0.30 eV flip is the MRM exit."""
    import sys
    import types

    seen = {}

    def decide(entry, candidates, rates, barriers, e_init, e_step, e_floor, cv):
        seen["e_init"] = e_init
        seen["barriers"] = list(barriers)
        seen["entry"] = entry
        return (
            "accepted",
            [1, 0],
            [2],
            [(1, 0, 1.0e8), (0, 1, 1.0e8), (1, 2, 1.0e2)],
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            (0.25, 0, False),
            (1, 1, []),
        )

    def mrm(transient, absorbing, rates, entry):
        seen["kernel"] = "mrm"
        seen["absorbing"] = list(absorbing)
        seen["transient"] = list(transient)
        assert 0 not in absorbing
        assert 2 in absorbing
        return (4.2e-3, [1.0], [0.5, 0.5])

    def fpta(*_args, **_kwargs):
        raise AssertionError("mean time selects MRM, not FPTA")

    mod = types.ModuleType("amsel")
    mod.discover_decide_status = decide
    mod.mrm = mrm
    mod.fpta = fpta
    monkeypatch.setitem(sys.modules, "amsel", mod)

    cfg = _config(
        tmp_path,
        """
        [amsel]
        discover_decide = true
        e_min_init = 0.25
        """,
    )
    assert cfg.sb_on is False
    entry, states = _isomer1()
    caplog.set_level(logging.INFO, logger="superbasin.amsel_gate")
    current, previous, time, steps = kmc_step(
        entry, states, 0.0, 0.025, None, config=cfg
    )
    assert steps == 1
    assert previous.number == 1
    assert current.number == 2
    assert time == pytest.approx(4.2e-3)
    assert seen["e_init"] == 0.25
    assert seen["entry"] == 1
    assert min(seen["barriers"]) < 0.25
    assert max(seen["barriers"]) > 0.25
    assert seen["kernel"] == "mrm"
    assert any(
        "amsel discover_decide status=accepted" in message for message in caplog.messages
    )
    assert any(
        "exit kernel=mrm" in message and "product=2" in message
        for message in caplog.messages
    )


def test_isomer1_sampled_exit_comes_from_fpta(tmp_path, monkeypatch, caplog):
    """With the mean clock off, the 0.30 eV exit time is one FPTA sample."""
    import sys
    import types

    def decide(entry, candidates, rates, barriers, e_init, e_step, e_floor, cv):
        return (
            "accepted",
            [1, 0],
            [2],
            [(1, 0, 1.0e8), (0, 1, 1.0e8), (1, 2, 1.0e2)],
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            (0.25, 0, False),
            (1, 1, []),
        )

    def mrm(*_args, **_kwargs):
        raise AssertionError("a sampled step selects FPTA, not MRM")

    def fpta(transient, absorbing, rates, entry, draw):
        assert 2 in absorbing
        assert 0.0 < float(draw) < 1.0
        return (7.5e-4, [1.0])

    mod = types.ModuleType("amsel")
    mod.discover_decide_status = decide
    mod.mrm = mrm
    mod.fpta = fpta
    monkeypatch.setitem(sys.modules, "amsel", mod)

    cfg = _config(
        tmp_path,
        """
        [amsel]
        discover_decide = true
        e_min_init = 0.25
        """,
        use_mean_time="false",
    )
    entry, states = _isomer1()
    caplog.set_level(logging.INFO, logger="superbasin.amsel_gate")
    current, previous, time, steps = kmc_step(
        entry, states, 0.0, 0.025, None, config=cfg
    )
    assert steps == 1
    assert previous.number == 1
    assert current.number == 2
    assert time == pytest.approx(7.5e-4)
    assert any("exit kernel=fpta" in message for message in caplog.messages)


def _install_mrm(monkeypatch, decide, mrm):
    import sys
    import types

    def fpta(*_args, **_kwargs):
        raise AssertionError("mean time selects MRM, not FPTA")

    mod = types.ModuleType("amsel")
    mod.discover_decide_status = decide
    mod.mrm = mrm
    mod.fpta = fpta
    monkeypatch.setitem(sys.modules, "amsel", mod)


def _accepted_basin(entry, candidates, rates, barriers, e_init, e_step, e_floor, cv):
    return (
        "accepted",
        [1, 0],
        [2],
        [(1, 0, 1.0e8), (0, 1, 1.0e8), (1, 2, 1.0e2)],
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        (0.25, 0, False),
        (1, 1, []),
    )


def test_discover_decide_steps_when_old_confidence_is_zero(tmp_path, monkeypatch, caplog):
    """Scheme old at confidence 0 still leaves through the 0.30 eV MRM exit."""
    seen = {}

    def mrm(transient, absorbing, rates, entry):
        seen["kernel"] = "mrm"
        seen["transient"] = list(transient)
        assert 2 in absorbing
        return (4.2e-3, [1.0], [0.5, 0.5])

    _install_mrm(monkeypatch, _accepted_basin, mrm)
    cfg = _config(
        tmp_path,
        """
        [amsel]
        discover_decide = true
        e_min_init = 0.25
        """,
        confidence=0.6,
    )
    assert cfg.akmc_confidence == pytest.approx(0.6)
    assert cfg.akmc_confidence_scheme == "old"
    entry, states = _isomer1()
    for state in states.mapping.values():
        state.confidence = 0.0
    caplog.set_level(logging.INFO, logger="superbasin.amsel_gate")
    current, previous, time, steps = kmc_step(
        entry, states, 0.0, 0.025, None, config=cfg
    )
    assert steps == 1
    assert previous.number == 1
    assert current.number == 2
    assert time == pytest.approx(4.2e-3)
    assert seen["kernel"] == "mrm"
    assert any("status=accepted" in message for message in caplog.messages)


def test_old_confidence_holds_the_step_when_discover_decide_is_off(tmp_path):
    """Without discover_decide, confidence 0 does not take a kinetic Monte Carlo step."""
    cfg = _config(tmp_path, "", confidence=0.6)
    assert cfg.amsel_discover_decide is False
    state0 = _State(
        0,
        {1: {"rate": 1.0, "product": 1, "barrier": 0.40, "prefactor": 1.0}},
        confidence=0.0,
    )
    state1 = _State(1, {}, confidence=0.0)
    current, previous, time, steps = kmc_step(
        state0,
        _States({0: state0, 1: state1}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 0
    assert current.number == 0
    assert previous.number == 0
    assert time == 0.0


def test_lone_barrier_is_a_direct_mrm_exit(tmp_path, monkeypatch, caplog):
    """A 0.30 eV saddle with no faster edge leaves. The transient set is the entry."""
    seen = {}

    def decide(entry, candidates, rates, barriers, e_init, e_step, e_floor, cv):
        seen["decide_entry"] = entry
        return (
            "rejected_no_metastable_basin",
            [entry],
            [],
            [],
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            (0.25, 0, False),
            (1, 0, []),
        )

    def mrm(transient, absorbing, rates, entry):
        seen.setdefault("calls", []).append(
            (list(transient), list(absorbing), list(rates), int(entry))
        )
        return (1.0e-2, [1.0], [1.0])

    _install_mrm(monkeypatch, decide, mrm)
    cfg = _config(
        tmp_path,
        """
        [amsel]
        discover_decide = true
        e_min_init = 0.25
        """,
        confidence=0.6,
        max_kmc_steps=0,
    )
    state0 = _State(
        0,
        {4: {"rate": 1.0e2, "product": 1, "barrier": 0.30, "prefactor": 1.0e12}},
        confidence=0.0,
    )
    state1 = _State(
        1,
        {5: {"rate": 1.0e2, "product": 2, "barrier": 0.30, "prefactor": 1.0e12}},
        confidence=0.0,
    )
    state2 = _State(2, {}, confidence=0.0)
    caplog.set_level(logging.INFO, logger="superbasin.amsel_gate")
    current, previous, time, steps = kmc_step(
        state0,
        _States({0: state0, 1: state1, 2: state2}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 1
    assert previous.number == 0
    assert current.number == 1
    assert time == pytest.approx(1.0e-2)
    assert seen["calls"] == [([0], [1], [(0, 1, 1.0e2)], 0)]
    assert any(
        "status=accepted" in message and "available=True" in message
        for message in caplog.messages
    )
    assert any("exit kernel=mrm" in message and "product=1" in message for message in caplog.messages)


def test_missing_amsel_does_not_step_below_confidence(tmp_path, monkeypatch, caplog):
    """A missing amsel package logs unavailable and does not hop at confidence 0."""
    monkeypatch.setitem(__import__("sys").modules, "amsel", None)
    cfg = _config(
        tmp_path,
        """
        [amsel]
        discover_decide = true
        e_min_init = 0.25
        """,
        confidence=0.6,
    )
    state0 = _State(
        0,
        {4: {"rate": 1.0e2, "product": 1, "barrier": 0.30, "prefactor": 1.0e12}},
        confidence=0.0,
    )
    state1 = _State(1, {}, confidence=0.0)
    caplog.set_level(logging.INFO, logger="superbasin.amsel_gate")
    current, _previous, time, steps = kmc_step(
        state0,
        _States({0: state0, 1: state1}),
        0.0,
        0.025,
        None,
        config=cfg,
    )
    assert steps == 0
    assert current.number == 0
    assert time == 0.0
    assert any("status=unavailable" in message for message in caplog.messages)
