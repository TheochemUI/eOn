#!/usr/bin/env python
"""MCAMC solver: amsel primary path, legacy libqd fallback.

The original MCAMC solver in this directory was a ctypes wrapper around
libmcamc.so (C++, Eigen, libqd for adaptive double-double/quad-double
arithmetic). amsel (lode-org) is now the preferred path -- pure Rust,
f64 with proper scaling, zero-copy DLPack I/O, free-threaded wheels.
The legacy path is retained as a fallback for ill-conditioned cases
that exceed f64's reach; switching between them is automatic.
"""
from __future__ import annotations

import ctypes
import logging
import os
from os.path import abspath, dirname, join

import numpy as np

logger = logging.getLogger("mcamc")


# ---------------------------------------------------------------------------
# amsel (primary path)
# ---------------------------------------------------------------------------

try:
    import amsel as _amsel

    _AMSEL_AVAILABLE = True
except ImportError as _amsel_import_error:
    _amsel = None
    _AMSEL_AVAILABLE = False
    logger.debug("amsel not importable: %s; legacy MCAMC path will be used", _amsel_import_error)


def _mcamc_amsel(Q, R, c):
    """Run amsel's fpta_bundle on numpy inputs, return (t, B, residual)."""
    q_arr = np.ascontiguousarray(Q, dtype=np.float64)
    r_arr = np.ascontiguousarray(R, dtype=np.float64)
    c_arr = np.ascontiguousarray(c, dtype=np.float64)
    t_dl, b_dl, residual = _amsel.fpta_bundle(q_arr, r_arr, c_arr)
    t = np.array(np.from_dlpack(t_dl), dtype=np.float64, copy=True)
    B = np.array(np.from_dlpack(b_dl), dtype=np.float64, copy=True)
    return t, B, residual


# ---------------------------------------------------------------------------
# Legacy ctypes-backed libmcamc.so (fallback for ill-conditioned cases)
# ---------------------------------------------------------------------------


def estimate_condition(Q, R):
    Qflat = list(Q.ravel())
    Qflat = (ctypes.c_double * len(Qflat))(*Qflat)
    Rflat = list(R.ravel())
    Rflat = (ctypes.c_double * len(Rflat))(*Rflat)
    _libmcamc.estimate_condition.restype = ctypes.c_double
    return _libmcamc.estimate_condition(Q.shape[0], Qflat, R.shape[1], Rflat)


def guess_precision(Q, R):
    cond = estimate_condition(Q, R)
    k = np.log10(cond)
    if k < 11:
        return "d"
    if k < 26:
        return "dd"
    if k < 60:
        return "qd"
    return "-"


def _mcamc_legacy_c(Q, R, c, prec="dd"):
    Qflat = list(Q.ravel())
    Qflat = (ctypes.c_double * len(Qflat))(*Qflat)
    Rflat = list(R.ravel())
    Rflat = (ctypes.c_double * len(Rflat))(*Rflat)
    cflat = list(c)
    cflat = (ctypes.c_double * len(cflat))(*cflat)

    Bflat = (ctypes.c_double * len(Rflat))()
    tflat = (ctypes.c_double * len(cflat))()
    residual = (ctypes.c_double * 1)()

    if prec == "f":
        solve = _libmcamc.solve_float
    elif prec == "d":
        solve = _libmcamc.solve_double
    elif prec == "dd":
        solve = _libmcamc.solve_double_double
    elif prec == "qd":
        solve = _libmcamc.solve_quad_double
    else:
        raise ValueError(f'Unknown prec value "{prec}"')

    solve(Q.shape[0], Qflat, R.shape[1], Rflat, cflat, Bflat, tflat, residual)
    B = np.array(list(Bflat)).reshape(R.shape)
    t = list(tflat)
    residual = list(residual)[0]
    return t, B, residual


def _mcamc_numpy(Q, R, c):
    """Pure-numpy fallback (no libmcamc.so, no amsel)."""
    Q = Q.copy()
    R = R.copy()
    c = c.copy()
    for i in range(len(c)):
        c[i] = 1.0 / c[i]
        Q[i] *= c[i]
        R[i] *= c[i]
    A = np.identity(Q.shape[0]) - Q
    t = np.linalg.solve(A, c)
    B = np.linalg.solve(A, R)
    residual = (A.dot(t) - c).max()
    return t, B, residual


_libpath = join(dirname(abspath(__file__)), "libmcamc.so")
try:
    _libmcamc = ctypes.CDLL(_libpath)
    _LEGACY_C_AVAILABLE = True
except OSError:
    _libmcamc = None
    _LEGACY_C_AVAILABLE = False


def _mcamc_legacy(Q, R, c, prec=None):
    if _LEGACY_C_AVAILABLE:
        return _mcamc_legacy_c(Q, R, c, prec=prec or guess_precision(Q, R))
    return _mcamc_numpy(Q, R, c)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

_FORCE_LEGACY = os.environ.get("EON_MCAMC_LEGACY", "").lower() in ("1", "true", "yes")


def mcamc(Q, R, c, prec=None):
    """Solve the MCAMC system and return (t, B, residual).

    Dispatches to amsel when available (pure Rust, f64 + scaling),
    falls back to the legacy libqd-backed C path for ill-conditioned
    cases amsel rejects, and finally to a pure-numpy solve if neither
    is available. Set ``EON_MCAMC_LEGACY=1`` to bypass amsel entirely.
    """
    if _AMSEL_AVAILABLE and not _FORCE_LEGACY:
        try:
            return _mcamc_amsel(Q, R, c)
        except _amsel.AmselError as exc:
            logger.warning(
                "amsel rejected MCAMC problem (%s); falling back to legacy path", exc
            )
    return _mcamc_legacy(Q, R, c, prec=prec)
