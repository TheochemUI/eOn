"""Import + smoke tests for the eon.mcamc amsel shim.

Validates that bead amsel-yo8 Step 1 (re-export full kernel surface
through eon.mcamc) is live: every entry point importable + callable
on a closed-form two-state chain that every kernel has a closed
answer for.
"""
from __future__ import annotations

import math

import pytest

try:
    import amsel  # type: ignore
except ImportError:
    amsel = None

import eon.mcamc as mcamc_mod


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_reexports_match_bead_spec():
    for name in (
        "mcamc",
        "adaptive_clock",
        "discover_adaptive",
        "discover_adaptive_status",
        "discover_decide_status",
        "discover_decide_diagnostics_status",
        "fpta",
        "fpt_spectrum",
        "mrm_direct",
        "ngt",
        "reduced_kinetics",
        "as_kmc_rate",
        "as_kmc_nf_from_delta",
        "as_kmc_nf_kaiser",
        "stiffness_step",
        "discover_fichthorn",
    ):
        assert hasattr(mcamc_mod, name), f"eon.mcamc missing {name}"


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_mrm_direct_single_state_closed_form():
    tau_total, rates_to_abs, x = mcamc_mod.mrm_direct(
        transient=[0], absorbing=[10],
        rates=[(0, 10, 42.0)], entry=0,
    )
    assert tau_total == pytest.approx(1.0 / 42.0)
    assert rates_to_abs == pytest.approx([42.0])
    assert x == pytest.approx([1.0])


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_fpta_single_state_analytic():
    # One transient state, one absorbing state, rate = lambda. Exit
    # time is -ln(r) / lambda; weights concentrate on the single
    # absorbing slot.
    lam = 5.0
    r = 0.5
    t, w = mcamc_mod.fpta(
        transient=[0], absorbing=[10],
        rates=[(0, 10, lam)], entry=0, r=r,
    )
    assert t == pytest.approx(-math.log(r) / lam)
    assert w == pytest.approx([1.0])


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_reduced_kinetics_two_mode_flags_rank1_failure():
    reduced = mcamc_mod.reduced_kinetics(
        transient=[0, 1],
        absorbing=[10, 20],
        rates=[(0, 1, 1.0e-4), (1, 0, 1.0e-4), (0, 10, 1.0), (1, 20, 1.0e-4)],
        entry=0,
    )
    assert reduced.slow_subspace_rank == 2
    assert reduced.rank1_invalidity > 0.0
    assert not reduced.one_rate_clock_is_plausible(1.0e-6)


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_adaptive_clock_single_state_uses_mean():
    t_exit, weights, mode, reduced = mcamc_mod.adaptive_clock(
        transient=[0], absorbing=[10],
        rates=[(0, 10, 42.0)], entry=0,
    )
    assert mode == "mean"
    assert t_exit == pytest.approx(1.0 / 42.0)
    assert weights == pytest.approx([1.0])
    assert reduced.one_rate_clock_is_plausible(1.0e-6)


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_adaptive_clock_two_mode_samples():
    rng = type("DummyRng", (), {"random": staticmethod(lambda: 0.5)})()
    t_exit, weights, mode, reduced = mcamc_mod.adaptive_clock(
        transient=[0, 1],
        absorbing=[10, 20],
        rates=[(0, 1, 1.0e-4), (1, 0, 1.0e-4), (0, 10, 1.0), (1, 20, 1.0e-4)],
        entry=0,
        rng=rng,
    )
    assert mode == "sampled"
    assert t_exit > 0.0
    assert sum(weights) == pytest.approx(1.0)
    assert reduced.slow_subspace_rank == 2


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_ngt_two_state_chain_agrees_with_closed_form():
    # 0 -> 1 at rate 1e3, 1 -> absorbing at rate 1.0. Slow leak.
    # Expected MFPT = 1 + 1/1e3. rate = 1 / MFPT. committor = 1.
    rate, mfpt, committor = mcamc_mod.ngt(
        transient=[0, 1], absorbing=[10],
        rates=[(0, 1, 1.0e3), (1, 10, 1.0)],
        source=[0], target=[10],
    )
    expected_mfpt = 1.0 + 1.0 / 1.0e3
    assert mfpt == pytest.approx(expected_mfpt, rel=1e-10)
    assert rate == pytest.approx(1.0 / expected_mfpt, rel=1e-10)
    assert committor == pytest.approx(1.0, abs=1e-12)


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_as_kmc_bounds_round_trip():
    alpha, delta = 2.0, 1e-3
    assert mcamc_mod.as_kmc_rate(1.0e6, repeats=5, alpha=alpha, n_f=10) == pytest.approx(1.0e6)
    # Kaiser >> Ferasat on tight delta (documented divergence).
    assert (
        mcamc_mod.as_kmc_nf_kaiser(alpha, delta)
        > mcamc_mod.as_kmc_nf_from_delta(alpha, delta) * 10
    )


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_discover_fichthorn_partitions_chain():
    # 0 -> 1 (TS 0.1), 1 -> 2 (TS 0.5). e_min=0.3 puts 0,1 in
    # transient, 2 on the absorbing border.
    transient, absorbing, rates = mcamc_mod.discover_fichthorn(
        entry=0,
        candidate_states=[0, 1, 2],
        rates=[(0, 1, 1.0), (1, 2, 1.0)],
        ts_energies=[0.1, 0.5],
        e_min=0.3,
    )
    assert set(transient) == {0, 1}
    assert set(absorbing) == {2}
    assert len(rates) == 2


@pytest.mark.skipif(amsel is None, reason="real amsel wheel not importable in test env")
def test_stiffness_identity_when_no_firings():
    s_out, scaled = mcamc_mod.stiffness_step(
        s=[1.0, 1.0], n_fwd=[0, 0], n_rev=[0, 0],
        rates_unscaled=[1.0, 2.0],
        pe_tol=0.02, min_sep=500, down_limit=5.0,
    )
    assert s_out == pytest.approx([1.0, 1.0])
    assert scaled == pytest.approx([1.0, 2.0])


def test_adaptive_clock_uses_fake_amsel_object(monkeypatch):
    class FakeReduced:
        slow_subspace_rank = 2
        rank1_invalidity = 1.0

        @staticmethod
        def one_rate_clock_is_plausible(_tol):
            return False

    class FakeFpta:
        t_exit = 3.0
        weights = [0.25, 0.75]

    class FakeProblem:
        def __init__(self, transient, absorbing, rates):
            self.transient = transient
            self.absorbing = absorbing
            self.rates = rates

        @staticmethod
        def reduced_kinetics(entry):
            assert entry == 0
            return FakeReduced()

        @staticmethod
        def fpta(entry, r):
            assert entry == 0
            assert r == pytest.approx(0.5)
            return FakeFpta()

    fake_amsel = type("FakeAmsel", (), {"AmcProblem": FakeProblem})()
    monkeypatch.setattr(mcamc_mod, "_amsel", fake_amsel)

    rng = type("DummyRng", (), {"random": staticmethod(lambda: 0.5)})()
    t_exit, weights, mode, reduced = mcamc_mod.adaptive_clock(
        transient=[0, 1],
        absorbing=[10, 20],
        rates=[(0, 1, 1.0)],
        entry=0,
        rng=rng,
    )
    assert t_exit == pytest.approx(3.0)
    assert weights == pytest.approx([0.25, 0.75])
    assert mode == "sampled"
    assert reduced.slow_subspace_rank == 2
