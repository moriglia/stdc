# SPDX-License-Identifier: GPL-3.0-or-later

from cython.view cimport array as cvarray
import numpy as np

cdef extern:
    void filteriir_(double* b, long* Nb,
                   double* a, long* Na,
                   double* x, long* Nx,
                   double* y,
                   double* y_initial)
    void filter_butter2_impulseinvariance_iir_(double* w0T,
                                              double* x,
                                              long* Nx,
                                              double* y,
                                              double* y_initial)
    

cpdef double[:] filteriir_initial(double[:] b, double[:] a, double[:] x, double [:] y_initial):
    cdef long Nb, Na, Nx
    cdef double [:] y
    Nb = len(b)
    Na = len(a)
    Nx = len(x)
    y  = np.zeros(Nx, dtype=np.double) # cvarray(shape=(Nx,), itemsize=sizeof(double), format="d")
    filteriir_(&b[0], &Nb, &a[0], &Na, &x[0], &Nx, &y[0], &y_initial[0])
    return y


cpdef double[:] filteriir(double[:] b, double[:] a, double[:] x):
    cdef long Nb, Na, Nx
    cdef double [:] y
    Nb = len(b)
    Na = len(a)
    Nx = len(x)
    y  = np.zeros(Nx, dtype=np.double) # = cvarray(shape=(Nx,), itemsize=sizeof(double), format="d")
    filteriir_(&b[0], &Nb, &a[0], &Na, &x[0], &Nx, &y[0], NULL)
    return y


cpdef double[:] filter_butter2_impulseinvariance_iir_initial(double w0T, double [:] x, double [:] y_initial):
    # Mind that the size of y_initial must be 2, more than 2 is ok, but useless
    if len(y_initial) < 2:
        raise ValueError("y_initial should be of size 2")
    cdef long Nx
    cdef double [:] y
    Nx = len(x)
    y = np.zeros(Nx, dtype=np.double) # = cvarray(shape=(Nx,), itemsize=sizeof(double), format="d")
    filter_butter2_impulseinvariance_iir_(&w0T, &x[0], &Nx, &y[0], &y_initial[0])
    return y


cpdef double[:] filter_butter2_impulseinvariance_iir(double w0T, double [:] x):
    cdef long Nx
    cdef double [:] y
    Nx = len(x)
    y = np.zeros(Nx, dtype=np.double) #= cvarray(shape=(Nx,), itemsize=sizeof(double), format="d")
    filter_butter2_impulseinvariance_iir_(&w0T, &x[0], &Nx, &y[0], NULL)
    return y
