"""Energy-level superbasin (energy lid) must not absorb an empty product.

The Mishin EAM Cu di-interstitial at 119 K, the ``eon_lid`` arm of the
amsel bench, dies in ``Superbasin.step`` once the lid has passed the
lowest barrier (5.04 kT). The product of that hop has an empty process
table: no escape has been found yet. The scheme still merges it, the
transient rate-sum for that state is zero, and the exit draw raises
``ValueError: Failed to select exit process``.

This graph is that situation with the atomistics removed. States 0 and
1 are ordinary minima and do merge. State 2 is the absorbing product,
energy below the lid, process table empty. State 3 is a second exit so
the absorption probabilities are not trivial.
"""
from __future__ import annotations

import math
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

import eon.superbasin as superbasin_mod
from eon.superbasinscheme import EnergyLevel

# eOn stores temperature in kelvin and converts with this factor.
_KT_PER_KELVIN = 1.0 / 11604.5
_TEMPERATURE_K = 119.0
_KT = _TEMPERATURE_K * _KT_PER_KELVIN
_LOWEST_BARRIER = 5.04 * _KT

# Rates on the merged pair. The exit branching and the mean escape time
# have a closed form (derived in ``_exact_absorbing_chain``).
_R01 = 50.0
_R10 = 40.0
_R0P = 4.0
_R0Q = 1.0


class _State:
    def __init__(self, number, energy, procs, path):
        self.number = number
        self.energy = energy
        self.procs = procs
        self.path = path

    def get_energy(self):
        return self.energy

    def get_process_table(self):
        return self.procs

    def __eq__(self, other):
        return isinstance(other, _State) and self.number == other.number

    def __hash__(self):
        return hash(self.number)


class _States:
    def __init__(self, states):
        self._states = list(states)

    def get_num_states(self):
        return len(self._states)

    def get_state(self, number):
        return self._states[number]

    def connect_states(self, states):
        return None


def _exact_absorbing_chain():
    """Mean absorption time and exit probabilities from either transient state.

    States 0 and 1 are the transient set. The only exits leave state 0,
    with rates ``_R0P`` and ``_R0Q``. State 1 returns to state 0 at
    ``_R10``. For an absorbing chain the exit probabilities from either
    start are the normalized exit rates, and the mean times are

        t0 = (r10 + r01) / (r10 * (r0P + r0Q))
        t1 = (1 / det) * (1 / c0 + 1 / c1)

    with ``c0 = r01 + r0P + r0Q``, ``c1 = r10`` and
    ``det = (r0P + r0Q) / c0``.
    """
    c0 = _R01 + _R0P + _R0Q
    c1 = _R10
    exit_sum = _R0P + _R0Q
    t0 = (_R10 + _R01) / (_R10 * exit_sum)
    det = exit_sum / c0
    t1 = (1.0 / det) * (1.0 / c0 + 1.0 / c1)
    p_product = _R0P / exit_sum
    p_other = _R0Q / exit_sum
    return t0, t1, p_product, p_other


def _build(root: Path):
    for number in range(4):
        (root / "states" / str(number)).mkdir(parents=True)
    paths = [str(root / "states" / str(number)) for number in range(4)]
    # Product energy is below the saddle (the lid). Its table is empty.
    metastable = _State(
        0,
        0.0,
        {
            10: {"product": 1, "barrier": _LOWEST_BARRIER, "rate": _R01},
            20: {"product": 2, "barrier": _LOWEST_BARRIER, "rate": _R0P},
            30: {"product": 3, "barrier": 3.0 * _LOWEST_BARRIER, "rate": _R0Q},
        },
        paths[0],
    )
    neighbor = _State(
        1,
        0.003,
        {11: {"product": 0, "barrier": _LOWEST_BARRIER - 0.003, "rate": _R10}},
        paths[1],
    )
    product = _State(2, -0.25, {}, paths[2])
    other = _State(3, 0.10, {}, paths[3])
    states = _States([metastable, neighbor, product, other])
    config = SimpleNamespace(
        sb_max_size=0,
        comp_use_identical=False,
        sb_state_file="superbasin",
        sb_superbasin_confidence=False,
        amsel_discover_decide=False,
    )
    scheme = EnergyLevel(
        str(root / "superbasins"),
        states,
        _KT,
        0.01,
        config=config,
    )
    return scheme, metastable, neighbor, product, other


def _fill_until_transient_basin(scheme, metastable, neighbor):
    for visit in range(500):
        start, end = (metastable, neighbor) if visit % 2 == 0 else (neighbor, metastable)
        scheme.register_transition(start, end)
        basin = scheme.get_containing_superbasin(metastable)
        if basin is not None and neighbor in basin.states:
            return basin
    raise AssertionError("transient minima were not merged under the energy lid")


def test_empty_product_below_the_lid_stays_an_exit(monkeypatch):
    t0, t1, p_product, p_other = _exact_absorbing_chain()
    assert math.isclose(p_product + p_other, 1.0)

    root = Path(tempfile.mkdtemp(prefix="energy-lid-"))
    scheme = None
    try:
        scheme, metastable, neighbor, product, other = _build(root)
        _fill_until_transient_basin(scheme, metastable, neighbor)

        # The hop that used to pull the empty product into the basin.
        scheme.register_transition(metastable, product)

        basin = scheme.get_containing_superbasin(metastable)
        assert basin is not None
        by_number = {0: metastable, 1: neighbor, 2: product, 3: other}

        def product_of(state_number, proc_id):
            table = by_number[state_number].get_process_table()
            return by_number[table[proc_id]["product"]]

        draw = {"u": 0.0}

        def sample():
            return draw["u"]

        monkeypatch.setattr(superbasin_mod.np.random, "sample", sample)

        def step_from(entry, u):
            draw["u"] = u
            mean_time, exit_state, landed, proc_id, _basin_id = basin.step(
                entry, product_of
            )
            return mean_time, exit_state.number, landed.number, proc_id

        # Unfixed, the empty product is inside the basin and this draw
        # raises ValueError: Failed to select exit process.
        mean0, exit_number, landed, proc = step_from(metastable, 0.0)
        assert {state.number for state in basin.states} == {0, 1}
        assert scheme.get_containing_superbasin(product) is None
        # The visit is recorded. Only the lid comparison is deferred,
        # so the product level has risen above its minimum energy.
        assert scheme.levels[product] > product.get_energy()
        assert exit_number == 0
        assert proc == 20
        assert landed == 2
        assert math.isclose(mean0, t0, rel_tol=1e-12, abs_tol=1e-12)

        # The selector's switch from the product exit to the other exit
        # is the absorption probability, to the resolution of the bisection.
        lo = 0.0
        hi = 1.0 - 1e-15
        _mean, _exit, landed_hi, proc_hi = step_from(metastable, hi)
        assert proc_hi == 30
        assert landed_hi == 3
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            _mean, _exit, _landed, proc_mid = step_from(metastable, mid)
            if proc_mid == 20:
                lo = mid
            else:
                hi = mid
        assert math.isclose(hi, p_product, rel_tol=0.0, abs_tol=1e-12)

        mean1, exit1, landed1, proc1 = step_from(neighbor, 0.0)
        assert exit1 == 0
        assert proc1 == 20
        assert landed1 == 2
        assert math.isclose(mean1, t1, rel_tol=1e-12, abs_tol=1e-12)
    finally:
        # write_data runs from EnergyLevel.__del__. Drop the scheme
        # while the state directories still exist.
        del scheme
        shutil.rmtree(root, ignore_errors=True)
