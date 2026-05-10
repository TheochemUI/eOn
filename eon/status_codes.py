"""Process-search termination codes shared with the C++ client.

Mirror of ``client/StatusTypes.h``. The integer values are wire format
written to ``results.dat`` by eonclient, so they MUST match the C++
``eonc::SaddleStatus`` enumerators byte-for-byte. The Python side
exposes them as ``IntEnum`` members so server code can compare against
``SaddleStatus.Good`` instead of the magic literal ``0``, while still
accepting raw ints transparently from ``parse_results``.

Whenever ``client/StatusTypes.h`` adds, removes, or renumbers a code,
update this module in the same commit and add a row to the label /
slug tables.
"""

from __future__ import annotations

from enum import IntEnum


class SaddleStatus(IntEnum):
    """Saddle-search termination reason. Wire-equivalent to
    ``eonc::SaddleStatus`` in ``client/StatusTypes.h``."""

    Good = 0
    Init = 1
    BadNoConvex = 2
    BadHighEnergy = 3
    BadMaxConcaveIterations = 4
    BadMaxIterations = 5
    BadNotConnected = 6
    BadPrefactor = 7
    BadHighBarrier = 8
    BadMinima = 9
    FailedPrefactor = 10
    PotentialFailed = 11
    NonnegativeAbort = 12
    NonlocalAbort = 13
    NegativeBarrier = 14
    BadMdTrajectoryTooShort = 15
    BadNoNegativeModeAtSaddle = 16
    BadNoBarrier = 17
    ZeromodeAbort = 18
    OptimizerError = 19
    DimerLostMode = 20
    DimerRestoredBest = 21
    BadArtnError = 22


class MinimizationStatus(IntEnum):
    """Minimizer termination as reported in ``results.dat`` for
    ``job_type == 'minimization'``. Distinct from the broader
    ``SaddleStatus`` -- the minimizer only emits these three."""

    Good = 0
    MaxIterations = 1
    PotentialFailed = 2


_UNKNOWN = "unknown_exit_code"


# Display labels (title case) used by akmc bad-saddle bookkeeping.
SADDLE_STATUS_LABELS: dict[SaddleStatus, str] = {
    SaddleStatus.Good: "Good",
    SaddleStatus.Init: "Init",
    SaddleStatus.BadNoConvex: "Saddle Search No Convex Region",
    SaddleStatus.BadHighEnergy: "Saddle Search Terminated High Energy",
    SaddleStatus.BadMaxConcaveIterations:
        "Saddle Search Terminated Concave Iterations",
    SaddleStatus.BadMaxIterations:
        "Saddle Search Terminated Total Iterations",
    SaddleStatus.BadNotConnected: "Not Connected",
    SaddleStatus.BadPrefactor: "Bad Prefactor",
    SaddleStatus.BadHighBarrier: "Bad Barrier",
    SaddleStatus.BadMinima: "Minimum Not Converged",
    SaddleStatus.FailedPrefactor: "Failed Prefactor Calculation",
    SaddleStatus.PotentialFailed: "Potential Failed",
    SaddleStatus.NonnegativeAbort: "Nonnegative Displacement Abort",
    SaddleStatus.NonlocalAbort: "Nonlocal Abort",
    SaddleStatus.NegativeBarrier: "Negative Barrier",
    SaddleStatus.BadMdTrajectoryTooShort: "MD Trajectory Too Short",
    SaddleStatus.BadNoNegativeModeAtSaddle: "No Negative Mode at Saddle",
    SaddleStatus.BadNoBarrier: "No Forward Barrier in Minimized Band",
    SaddleStatus.ZeromodeAbort: "MinMode Zero Mode Abort",
    SaddleStatus.OptimizerError: "Optimizer Error",
    SaddleStatus.DimerLostMode: "Dimer Lost Mode",
    SaddleStatus.DimerRestoredBest: "Dimer Restored Best Mode",
    SaddleStatus.BadArtnError: "ARTn Error",
}


# Slug-case labels (snake_case) used by explorer.py for routing
# decisions. Codes that the legacy table left as ``unknown_exit_code``
# stay that way to preserve downstream string-equality semantics.
SADDLE_STATUS_SLUGS: dict[SaddleStatus, str] = {
    SaddleStatus.Good: "good",
    SaddleStatus.Init: _UNKNOWN,
    SaddleStatus.BadNoConvex: "no_convex",
    SaddleStatus.BadHighEnergy: "high_energy",
    SaddleStatus.BadMaxConcaveIterations: "max_concave_iterations",
    SaddleStatus.BadMaxIterations: "max_iterations",
    SaddleStatus.BadNotConnected: _UNKNOWN,
    SaddleStatus.BadPrefactor: _UNKNOWN,
    SaddleStatus.BadHighBarrier: _UNKNOWN,
    SaddleStatus.BadMinima: _UNKNOWN,
    SaddleStatus.FailedPrefactor: _UNKNOWN,
    SaddleStatus.PotentialFailed: "potential_failed",
    SaddleStatus.NonnegativeAbort: "nonnegative_abort",
    SaddleStatus.NonlocalAbort: "nonlocal_abort",
    SaddleStatus.NegativeBarrier: _UNKNOWN,
    SaddleStatus.BadMdTrajectoryTooShort: "md_trajectory_too_short",
    SaddleStatus.BadNoNegativeModeAtSaddle: "no_negative_mode_at_saddle",
    SaddleStatus.BadNoBarrier: "no_barrier",
    SaddleStatus.ZeromodeAbort: "zero_mode_abort",
    SaddleStatus.OptimizerError: "optimizer_error",
    SaddleStatus.DimerLostMode: "dimer_lost_mode",
    SaddleStatus.DimerRestoredBest: "dimer_restored_best",
    SaddleStatus.BadArtnError: "artn_error",
}


MINIMIZATION_STATUS_SLUGS: dict[MinimizationStatus, str] = {
    MinimizationStatus.Good: "good",
    MinimizationStatus.MaxIterations: "max_iterations",
    MinimizationStatus.PotentialFailed: "potential_failed",
}


def saddle_status_label(code: int) -> str:
    """Human-readable label for a saddle-search termination code.

    Tolerant of unknown codes -- returns a generic placeholder rather
    than raising, so a future C++ rev cannot cause an IndexError or
    KeyError to crash long-running akmc jobs.
    """
    try:
        return SADDLE_STATUS_LABELS[SaddleStatus(code)]
    except (ValueError, KeyError):
        return f"Unknown termination code ({code})"


def saddle_status_slug(code: int) -> str:
    """Slug-case label for a saddle-search termination code, with
    the same unknown-code tolerance as :func:`saddle_status_label`."""
    try:
        return SADDLE_STATUS_SLUGS[SaddleStatus(code)]
    except (ValueError, KeyError):
        return _UNKNOWN


def minimization_status_slug(code: int) -> str:
    """Slug-case label for a minimization termination code."""
    try:
        return MINIMIZATION_STATUS_SLUGS[MinimizationStatus(code)]
    except (ValueError, KeyError):
        return _UNKNOWN
