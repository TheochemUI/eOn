#!/usr/bin/env python
import ctypes
import numpy as np
from os.path import join, abspath, dirname
import logging
logger = logging.getLogger('mcamc')


def _numpy_condition(Q, R):
    """Largest absolute row scale divided by that row's exit sum.

    The scale of a row is the sum of that row of Q and of R. The exit
    sum is that row of R. A zero scale or a zero exit makes the
    condition non-finite.
    """
    rates = np.array(Q, dtype=float, copy=True)
    exits = np.array(R, dtype=float, copy=True)
    scale = rates.sum(axis=1) + exits.sum(axis=1)
    exit_sum = exits.sum(axis=1)
    zero_scale = bool(np.any(scale == 0.0))
    zero_exit = bool(np.any(exit_sum == 0.0))
    if zero_scale or zero_exit:
        return float("inf")
    return float(np.max(np.abs(scale / exit_sum)))


def estimate_condition(Q, R):
    if libmcamc is None:
        return _numpy_condition(Q, R)
    Qflat = list(np.asarray(Q, dtype=float).ravel())
    Qflat = (ctypes.c_double * len(Qflat))(*Qflat)
    Rflat = list(np.asarray(R, dtype=float).ravel())
    Rflat = (ctypes.c_double * len(Rflat))(*Rflat)
    libmcamc.estimate_condition.restype = ctypes.c_double
    return libmcamc.estimate_condition(Q.shape[0], Qflat, R.shape[1], Rflat)


def guess_precision(Q, R):
    cond = estimate_condition(Q, R)
    k = np.log10(cond)
    if k < 11:
        prec = 'd'
    elif k < 26:
        prec = 'dd'
    elif k < 60:
        prec = 'qd'
    else:
        #problem is too ill-conditioned to solve accurately
        prec = '-'
    return prec


def c_mcamc(Q, R, c, prec='dd'):
    Qflat = list(Q.ravel())
    Qflat = (ctypes.c_double * len(Qflat))(*Qflat)
    Rflat = list(R.ravel())
    Rflat = (ctypes.c_double * len(Rflat))(*Rflat)
    cflat = list(c)
    cflat = (ctypes.c_double * len(cflat))(*cflat)

    Bflat = (ctypes.c_double * len(Rflat))()
    tflat = (ctypes.c_double * len(cflat))()

    residual = (ctypes.c_double * 1)()

    if prec == 'f':
        solve = libmcamc.solve_float
    elif prec == 'd':
        solve = libmcamc.solve_double
    elif prec == 'dd':
        solve = libmcamc.solve_double_double
    elif prec == 'qd':
        solve = libmcamc.solve_quad_double
    else:
        raise ValueError('Unknown prec value "%s"' % prec)
    # void solve(int Qsize, double *Qflat, int Rcols, double *Rflat,
    #           double *c_in, double *B, double *t)
    solve(Q.shape[0], Qflat, R.shape[1], Rflat, cflat, Bflat, tflat, residual)

    B = np.array(list(Bflat)).reshape(R.shape)
    t = list(tflat)
    residual = list(residual)[0]

    return t, B, residual


def np_mcamc(Q, R, c, prec="NA"):
    for i in range(len(c)):
        c[i] = 1.0/c[i]
        Q[i] *= c[i]
        R[i] *= c[i]
    A = np.identity(Q.shape[0]) - Q
    t = np.linalg.solve(A, c)
    B = np.linalg.solve(A, R)
    residual = (A.dot(t) - c).max()
    return t, B, residual


libpath = join(dirname(abspath(__file__)), 'libmcamc.so')
libmcamc = None
try:
    libmcamc = ctypes.CDLL(libpath)
    mcamc = c_mcamc
except OSError:
    logger.debug("Was unable to use libmcamc, using numpy instead.")
    mcamc = np_mcamc
