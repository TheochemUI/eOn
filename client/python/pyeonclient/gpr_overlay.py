"""Optional gpr_optim overlay for pyeonclient Matter workflows.

Requires the ``gpr_optim`` package (nanobind step surface + GPNEBSession).
Does not change the default ``NEB(accelerant="gp")`` path (eOn GPSurrogateJob);
call this explicitly when you want gpr_optim's production GP-OIE / OCINEB stack.
"""

from __future__ import annotations

from typing import Any, Optional


def gpr_optim_available() -> bool:
    try:
        import gpr_optim  # noqa: F401

        return True
    except ImportError:
        return False


def gpneb_session_from_matters(
    initial: Any,
    final: Any,
    n_intermediate: int = 8,
    *,
    calculator: Any,
    gp: Any = None,
    path_method: str = "idpp",
    **session_kwargs: Any,
) -> Any:
    """Build a :class:`gpr_optim.session.GPNEBSession` from Matter endpoints.

    Parameters
    ----------
    initial, final
        pyeonclient.Matter (or duck-types with positions / atomic_numbers).
    n_intermediate
        Interior image count.
    calculator
        ASE Calculator used as the oracle (e.g. MetatomicCalculator).
    gp
        Optional pre-built ``gpr_optim.GaussianProcess``. Default: production-ish
        Laplace-slice chain without forcing board pins.
    path_method
        C++ path init name (idpp / linear / sidpp / ...).
    **session_kwargs
        Forwarded to GPNEBSession (acquisition, use_mmf, climb, ...).
    """
    if not gpr_optim_available():
        raise ImportError(
            "gpneb_session_from_matters requires gpr_optim "
            "(pip/pixi install the gpr_optim package)"
        )
    from gpr_optim import GaussianProcess
    from gpr_optim.session import GPNEBSession

    if gp is None:
        gp = GaussianProcess(
            uncertainty="laplace",
            psis=True,
            psis_deterministic=True,
            use_delta=True,
        )
    return GPNEBSession.from_matters(
        initial,
        final,
        n_intermediate,
        oracle=calculator,
        gp=gp,
        path_method=path_method,
        calculator=calculator,
        **session_kwargs,
    )


def compute_gpneb_from_matters(
    initial: Any,
    final: Any,
    n_intermediate: int = 8,
    *,
    calculator: Any,
    fmax: float = 0.05,
    steps: int = 500,
    gp: Any = None,
    path_method: str = "idpp",
    **session_kwargs: Any,
) -> Any:
    """Full prepare+compute from Matter endpoints; returns GPNEBSession."""
    sess = gpneb_session_from_matters(
        initial,
        final,
        n_intermediate,
        calculator=calculator,
        gp=gp,
        path_method=path_method,
        **session_kwargs,
    )
    sess.compute(fmax=fmax, steps=steps)
    return sess


__all__ = [
    "gpr_optim_available",
    "gpneb_session_from_matters",
    "compute_gpneb_from_matters",
]
