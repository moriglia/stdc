# SPDX-License-Identifier: GPL-3.0-or-later

cdef extern:
    void dummy()
    int mult_ints(int*, int*)
    double inner_prod(double*, double*)
    void print_array_(double *, long *) 

def dummy_from_fortran():
    dummy()
    return


cpdef int mult_ints_from_fortran(int a, int b):
    print("Hello")
    return mult_ints(&a,&b)


cpdef double inner_prod_f(double[:] a, double[:] b):
    return inner_prod(&a[0], &b[0])



