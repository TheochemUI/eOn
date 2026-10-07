"""A prepared aKMC imbalance decays as 1/(k_f + k_b)."""

import math
from pathlib import Path

from eon.akmcstate import AKMCState
from eon.approach import arrhenius_rate, prepared_imbalance, recurrence_time, survival

ROOT = Path(__file__).resolve().parents[1]
K_B = 8.617333262145177e-5


def test_the_nitride_channel_gap_is_the_two_state_rate():
    source = (ROOT / "eon" / "akmcstate.py").read_text(encoding="utf-8")
    assert "prepared_imbalance(" in source
    guide = (ROOT / "docs" / "source" / "user_guide" / "akmc.md").read_text(encoding="utf-8")
    assert "A cluster has no box length" in guide
    assert "4 * pi" not in (ROOT / "eon" / "approach.py").read_text(encoding="utf-8")

    forward = arrhenius_rate(0.2447, 1050.0)
    reverse = arrhenius_rate(1.189, 1050.0)
    tau = prepared_imbalance(forward, reverse)
    assert math.isclose(tau, 14.945e-12, rel_tol=1e-3)
    assert math.isclose(survival(10e-12, tau), 0.5122, rel_tol=1e-3)
    assert math.isclose(survival(30e-12, tau), 0.1343, rel_tol=1e-3)
    assert math.isclose(tau * math.log(2), 10.359e-12, rel_tol=1e-3)
    cold = prepared_imbalance(arrhenius_rate(0.2447, 600.0), arrhenius_rate(1.189, 600.0))
    assert math.isclose(cold, 113.603e-12, rel_tol=1e-3)
    assert math.isclose(recurrence_time(arrhenius_rate(1.864, 1050.0)), 0.883e-3, rel_tol=1e-2)

    class _List:
        kT = K_B * 1050.0

    class _Config:
        akmc_eq_rate = 0.0

    state = AKMCState.__new__(AKMCState)
    state.number = 1
    state.statelist = _List()
    state.config = _Config()
    proc = {
        "barrier": 0.2447,
        "prefactor": 1.0e12,
        "product_prefactor": 1.0e12,
        "product_energy": 0.2447 - 1.189,
        "product": 2,
    }
    assert math.isclose(state.population_gap(proc, 0.0), tau, rel_tol=1e-12)
    proc["product"] = 1
    assert state.population_gap(proc, 0.0) is None
