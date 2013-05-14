#!/usr/bin/env python
import ctypes
import numpy as np
from os.path import abspath

def mcamc(Q, R, c, prec='dd'):
    Qflat = list(Q.ravel())
    Qflat = (ctypes.c_double * len(Qflat))(*Qflat)
    Rflat = list(R.ravel())
    Rflat = (ctypes.c_double * len(Rflat))(*Rflat)
    cflat = list(c)
    cflat = (ctypes.c_double * len(cflat))(*cflat)

    Bflat = (ctypes.c_double * len(Rflat))()
    tflat = (ctypes.c_double * len(cflat))()

    residual = (ctypes.c_double * 1)()

    libmcamc = ctypes.CDLL(abspath('libmcamc.so'))

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
