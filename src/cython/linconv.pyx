from cython.view cimport array as cvarray


cdef extern:
    void overlap_save_(double * x,
                       long   * Nx,
                       double * h,
                       long   * M,
                       long   * N,
                       double * y)
    
    void overlap_save_iteration_(double * x,
                                 double * H,
                                 long   * N,
                                 long   * M,
                                 double * y,
                                 long   * N_new)    

    void mult_hermitian_(double * X,
                         double * H,
                         long   * N,
                         double * Y)

    void rfft_(double * x, int * N, int * jsign)
    
    void ifft_hermitian_(double * Y,
                         long * N_half)


cpdef double [:] overlap_save(double [:] x, double [:] h, long N):
    cdef long Nx = len(x)
    cdef long M  = len(h)
    cdef double [:] y = cvarray(shape=(Nx+M-1,),
                                itemsize=sizeof(double),
                                format="d")
    overlap_save_(&x[0], &Nx, &h[0], &M, &N, &y[0])

    return y

cpdef double [:] overlap_save_iteration(double [:] x, double [:] H, long M):
    cdef long N = len(x),
    cdef long N_new
    cdef double [:] res

    if (len(H) != N):
        raise ValueError(f"Size of DFT of H does not match size of x: {N}")

    N_new = N - M + 1
    res = cvarray(shape=(N_new,),
                  itemsize=sizeof(double),
                  format="d")

    overlap_save_iteration_(&x[0],
                            &H[0],
                            &N,
                            &M,
                            &res[0],
                            &N_new)
    return res

    

cpdef double [:] mult_hermitian(double [:] X, double [:] H):
    cdef long N = len(X)
    cdef double [:] Y
    
    if (len(H) != N ):
        raise ValueError("Sizes do not correspond")

    Y = cvarray(shape=(N,),
                itemsize=sizeof(double),
                format="d")

    mult_hermitian_(&X[0], &H[0], &N, &Y[0])

    return Y


cpdef double [:] rfft(double [:] x):
    cdef int N = len(x)
    cdef int jsign = -1
    cdef double [:] res

    res = cvarray(shape=(N,),
                  itemsize=sizeof(double),
                  format="d")
    res[:] = x[:]

    rfft_(&res[0], &N, &jsign)
    return res


cpdef double [:] ifft_hermitian(double [:] Y):
    cdef long N = len(Y)
    cdef double [:] res

    res = cvarray(shape=(N,),
                  itemsize=sizeof(double),
                  format="d")

    res[:] = Y # Copy all elements

    ifft_hermitian_(&res[0], &N)

    return res
    
    
