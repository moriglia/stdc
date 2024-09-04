# SPDX-License-Identifier: GPL-3.0-or-later

from cython.view cimport array as cvarray
import numpy as np

cdef extern:
    void filteriir(double* b, long* Nb,
                   double* a, long* Na,
                   double* x, long* Nx,
                   double* y,
                   double* y_initial)
    # double[:] ext_filteriir(double* b, long* Nb,
    #                         double* a, long* Na,
    #                         double* x, long* Nx,
    #                         double* y_initial)

cpdef double[:] ext_filteriir_initial(double[:] b, double[:] a, double[:] x, double [:] y_initial):
    cdef long Nb, Na, Nx
    cdef double [:] y
    Nb = len(b)
    Na = len(a)
    Nx = len(x)
    y  = cvarray(shape=(Nx,), itemsize=sizeof(double), format="d")
    filteriir(&b[0], &Nb, &a[0], &Na, &x[0], &Nx, &y[0], &y_initial[0])
    return y


cpdef double[:] ext_filteriir(double[:] b, double[:] a, double[:] x):
    cdef long Nb, Na, Nx
    cdef double [:] y
    Nb = len(b)
    Na = len(a)
    Nx = len(x)
    y  = cvarray(shape=(Nx,), itemsize=sizeof(double), format="d")
    filteriir(&b[0], &Nb, &a[0], &Na, &x[0], &Nx, &y[0], NULL)
    return y

