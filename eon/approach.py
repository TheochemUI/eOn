"""Prepared population imbalance on a two-state rate matrix.

Approach to equilibrium and adaptive kinetic Monte Carlo are the same
experiment on different generators. This module propagates the rate matrix.
A cluster has no box length, so the decay time is not turned into a
conductivity.
"""

import math

# eV/K. The same constant the tunnelling rates use.
K_B = 8.617333262145177e-5


def arrhenius_rate(barrier_eV, temperature_K, prefactor_hz=1.0e12):
    """Harmonic transition-state rate in Hz."""
    if temperature_K <= 0.0 or prefactor_hz <= 0.0 or barrier_eV < 0.0:
        raise ValueError("a rate needs a positive temperature, a positive prefactor, and a non-negative barrier")
    return prefactor_hz * math.exp(-barrier_eV / (K_B * temperature_K))


def prepared_imbalance(forward_rate, reverse_rate):
    """Decay time, in seconds, of a prepared population imbalance.

    The two-state generator relaxes as ``exp(-t/tau)`` with
    ``tau = 1/(k_f + k_b)``.
    """
    total = forward_rate + reverse_rate
    if not total > 0.0 or not math.isfinite(total):
        raise ValueError("a two-state gap needs a positive total rate")
    return 1.0 / total


def survival(time_s, decay_time_s):
    """Remaining imbalance ``exp(-t/tau)``."""
    if decay_time_s <= 0.0 or time_s < 0.0:
        raise ValueError("survival needs a positive decay time and a non-negative time")
    return math.exp(-time_s / decay_time_s)


def recurrence_time(rate):
    """Mean return time of a saddle that comes back to the same basin.

    That process has no population mode.
    """
    if not rate > 0.0 or not math.isfinite(rate):
        raise ValueError("a recurrence needs a positive rate")
    return 1.0 / rate
