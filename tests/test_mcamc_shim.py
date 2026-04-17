"""Import + smoke tests for the eon.mcamc amsel shim.

Validates that bead amsel-yo8 Step 1 (re-export full kernel surface
through eon.mcamc) is live: every entry point importable + callable
on a closed-form two-state chain that every kernel has a closed
answer for.
"""
from __future__ import annotations

import math

import pytest

amsel = pytest.importorskip("amsel")

import eon.mcamc as mcamc_mod


def test_reexports_match_bead_spec():
    for name in (
        "mcamc",
        "fpta",
        "mrm_direct",
        "ngt",
        "as_kmc_rate",
        "as_kmc_nf_from_delta",
        "as_kmc_nf_kaiser",
        "stiffness_step",
    ):
        assert hasattr(mcamc_mod, name), f"eon.mcamc missing {name}"


def test_mrm_direct_single_state_closed_form():
    tau_total, rates_to_abs, x = mcamc_mod.mrm_direct(
        transient=[0], absorbing=[10],
        rates=[(0, 10, 42.0)], entry=0,
    )
    assert tau_total == pytest.approx(1.0 / 42.0)
    assert rates_to_abs == pytest.approx([42.0])
    assert x == pytest.approx([1.0])


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


def test_as_kmc_bounds_round_trip():
    alpha, delta = 2.0, 1e-3
    assert mcamc_mod.as_kmc_rate(1.0e6, repeats=5, alpha=alpha, n_f=10) == pytest.approx(1.0e6)
    # Kaiser >> Ferasat on tight delta (documented divergence).
    assert (
        mcamc_mod.as_kmc_nf_kaiser(alpha, delta)
        > mcamc_mod.as_kmc_nf_from_delta(alpha, delta) * 10
    )


def test_stiffness_identity_when_no_firings():
    s_out, scaled = mcamc_mod.stiffness_step(
        s=[1.0, 1.0], n_fwd=[0, 0], n_rev=[0, 0],
        rates_unscaled=[1.0, 2.0],
        pe_tol=0.02, min_sep=500, down_limit=5.0,
    )
    assert s_out == pytest.approx([1.0, 1.0])
    assert scaled == pytest.approx([1.0, 2.0])
