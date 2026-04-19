"""MCAMC kernel entry points.

The canonical superbasin solver is `mcamc(Q, R, c)` -- the Ferasat-
corrected MRM that dispatches through amsel's fpta_bundle and falls
back to libqd/numpy when amsel rejects the problem. `Superbasin.step`
uses it directly.

The surrounding amsel kernel surface (fpta, mrm, ngt, askmc,
stiffness) is re-exported here so benchmarks and research scripts can
reach every literature method through a single import path without
having to know about amsel-python vs the legacy libmcamc.so layer.
Import shape intentionally matches the bead amsel-yo8 spec.
"""
import numpy

from .mcamc import mcamc

try:
    import amsel as _amsel
except ImportError:
    _amsel = None


def _require_amsel():
    if _amsel is None:
        raise ImportError(
            "amsel is not installed; run `pip install amsel` (or "
            "`maturin develop` from the amsel source tree) to enable "
            "the full kernel surface."
        )
    return _amsel


def fpta(transient, absorbing, rates, entry, r):
    """Ferasat-corrected FPTA. See amsel.fpta."""
    return _require_amsel().fpta(transient, absorbing, rates, entry, r)


def mrm_direct(transient, absorbing, rates, entry):
    """Ferasat-corrected MRM on a typed AmcProblem. See amsel.mrm.

    Distinct from `mcamc(Q, R, c)` above: mcamc takes dense (Q, R, c)
    matrices in eon's legacy convention, while `mrm_direct` takes the
    typed transient/absorbing/edge list amsel uses natively. For
    head-to-head benchmarks prefer the direct path -- it skips the
    dense matrix assembly.
    """
    return _require_amsel().mrm(transient, absorbing, rates, entry)


def ngt(transient, absorbing, rates, source, target, order=None):
    """Wales 2009 graph-transformation kernel. See amsel.ngt."""
    return _require_amsel().ngt(transient, absorbing, rates, source, target, order)


def as_kmc_rate(rate, repeats, alpha, n_f):
    """Chatterjee-Voter 2010 / Ferasat 2020 A.18 rate modification."""
    return _require_amsel().as_kmc_rate(rate, repeats, alpha, n_f)


def as_kmc_nf_from_delta(alpha, delta):
    """Ferasat 2020 A.17 N_f bound."""
    return _require_amsel().as_kmc_nf_from_delta(alpha, delta)


def as_kmc_nf_kaiser(alpha, delta):
    """Kaiser-Gosswein-Gagliardi 2020 Eq 15 N_f bound."""
    return _require_amsel().as_kmc_nf_kaiser(alpha, delta)


def stiffness_step(s, n_fwd, n_rev, rates_unscaled, pe_tol, min_sep, down_limit):
    """Prats-Li-Stamatakis 2025 adaptive stiffness scaling, one update."""
    return _require_amsel().stiffness_step(
        s, n_fwd, n_rev, rates_unscaled, pe_tol, min_sep, down_limit,
    )


def discover_fichthorn(entry, candidate_states, rates, ts_energies, e_min, max_size=0):
    """Fichthorn-Lin 2013 local-rule superbasin discovery.

    Build an AmcProblem-shaped partition `(transient, absorbing,
    rates)` from a KMC candidate pool by BFS under a TS-energy
    cutoff `e_min`. `max_size=0` disables the transient-set cap;
    any positive value bounds the partition size (overflow raises).
    """
    return _require_amsel().discover_fichthorn(
        entry, candidate_states, rates, ts_energies, e_min, max_size,
    )


def discover_adaptive(entry, candidate_states, rates, ts_energies, e_min_init,
                      e_min_step, e_min_floor, cv_threshold=10.0,
                      min_absorbing=1, max_transient_fraction=1.0, max_size=0):
    """Adaptive discover wrapper. See amsel.discover_adaptive."""
    return _require_amsel().discover_adaptive(
        entry, candidate_states, rates, ts_energies,
        e_min_init, e_min_step, e_min_floor, cv_threshold,
        min_absorbing, max_transient_fraction, max_size,
    )


def discover_adaptive_status(entry, candidate_states, rates, ts_energies, e_min_init,
                             e_min_step, e_min_floor, cv_threshold=10.0,
                             min_absorbing=1, max_transient_fraction=1.0, max_size=0):
    """Adaptive discover status wrapper. See amsel.discover_adaptive_status."""
    return _require_amsel().discover_adaptive_status(
        entry, candidate_states, rates, ts_energies,
        e_min_init, e_min_step, e_min_floor, cv_threshold,
        min_absorbing, max_transient_fraction, max_size,
    )


def discover_decide_status(entry, candidate_states, rates, ts_energies, e_min_init,
                           e_min_step, e_min_floor, cv_threshold=10.0,
                           min_absorbing=1, max_transient_fraction=1.0, max_size=0):
    """Split-aware discover decision. See amsel.discover_decide_status."""
    return _require_amsel().discover_decide_status(
        entry, candidate_states, rates, ts_energies,
        e_min_init, e_min_step, e_min_floor, cv_threshold,
        min_absorbing, max_transient_fraction, max_size,
    )


def discover_decide_diagnostics_status(entry, candidate_states, rates, ts_energies, e_min_init,
                                       e_min_step, e_min_floor, cv_threshold=10.0,
                                       min_absorbing=1, max_transient_fraction=1.0,
                                       max_size=0):
    """Split-aware discover decision with diagnostics."""
    return _require_amsel().discover_decide_diagnostics_status(
        entry, candidate_states, rates, ts_energies,
        e_min_init, e_min_step, e_min_floor, cv_threshold,
        min_absorbing, max_transient_fraction, max_size,
    )


def fpt_spectrum(transient, absorbing, rates, entry):
    """Typed FPT spectrum. See amsel.AmcProblem.fpt_spectrum."""
    problem = _require_amsel().AmcProblem(transient, absorbing, rates)
    return problem.fpt_spectrum(entry)


def reduced_kinetics(transient, absorbing, rates, entry):
    """Typed reduced kinetics. See amsel.AmcProblem.reduced_kinetics."""
    problem = _require_amsel().AmcProblem(transient, absorbing, rates)
    return problem.reduced_kinetics(entry)


def adaptive_clock(transient, absorbing, rates, entry, rank_tol=1.0e-8, rng=None):
    """Choose a mean or sampled clock from ``reduced_kinetics``.

    Returns ``(t_exit, weights, mode, reduced)`` where ``mode`` is either
    ``"mean"`` or ``"sampled"``.
    """
    amsel = _require_amsel()
    problem = amsel.AmcProblem(transient, absorbing, rates)
    reduced = problem.reduced_kinetics(entry)
    if reduced.one_rate_clock_is_plausible(rank_tol):
        mrm = problem.mrm(entry)
        weights = numpy.asarray(mrm.rate_to_absorbing, dtype=float) * float(mrm.tau_total)
        weights = weights / weights.sum()
        return float(mrm.tau_total), weights.tolist(), "mean", reduced

    if rng is None:
        r = float(numpy.random.random())
    else:
        r = float(rng.random())
    fpta_res = problem.fpta(entry, r)
    weights = numpy.asarray(fpta_res.weights, dtype=float)
    weights = weights / weights.sum()
    return float(fpta_res.t_exit), weights.tolist(), "sampled", reduced


__all__ = [
    "mcamc",
    "fpta",
    "fpt_spectrum",
    "mrm_direct",
    "ngt",
    "reduced_kinetics",
    "adaptive_clock",
    "as_kmc_rate",
    "as_kmc_nf_from_delta",
    "as_kmc_nf_kaiser",
    "stiffness_step",
    "discover_fichthorn",
    "discover_adaptive",
    "discover_adaptive_status",
    "discover_decide_status",
    "discover_decide_diagnostics_status",
]
