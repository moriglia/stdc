# Simulation Techniques for Digital Communications

This repository contains the code for the simulations performed as part of the
exam for the course mentioned in the title

## Project structure

### Processing 
Processing functions are mainly coded in `fortran`, their code reside in `src/fortran/`.
These functions are then re-defined as extern functions in `cython` code, residing in `src/cython`


### Simulations
Simulation scripts are in the `sims/` folder. They can `import stdc`.
If they are launched as modules `./pywithlib.sh -m sims.sim_iir`
(see [this section](#running-python-for-this-project) for the expanation for `pywithlib.sh`)


### Compilation

Just `make` it. It should:
1. compile the object files into `build/**/*.o` from the `fortran` sources with `gfortran`
2. link the object files into libraries and move them into the `lib/` folder
3. compile the `cython` files, dynamically link them to the libraries and move them into `stdc/`


## Running python for this project

This project includes some shared objects imported from `./lib`.
In order to let `python3` know that shared libraries may reside there, you have some options:

1. Export the library path for the session (or add it to your `.bashrc`)
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/lib/
#then
python3 [etc...]
```
2. Modify the library path just for the current command, like this:
`LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/lib/ python3 [etc...]`
3. use the ready made script `./pywithlib.sh [etc...]`
For instance you can just type `./pywithlib.sh -m sims.sim_iir` 
instead of `python3 -m sims.sim_iir`

## Modules

### `stdc.linconv`

Performs linear convolution using the Overlap-Add algorithm.
The processing code is in `src/fortran/libstdc/linconv.f90`:
```fortran
subroutine overlap_save(x, Nx, h, M, N, y)
! ...
end subroutine overlap_save

subroutine overlap_save_iteration(...)
! ...
end subroutine overlap_save_iteration

! other subroutines ...
```
This code also calls an FFT subroutine from `src/fortran/libcern/fft.f`, so both of them need be compiled:

```shell
gfortran -O2 -c src/fortran/libstdc/linconv.f90 -fpic -o build/stdc/linconv.o
gfortran -O2 -c src/fortran/libcern/fft.f -fpic -o build/cern/fft.o
```

These subroutines are wrapped by a cython code in `src/cython/linconv.pyx`.
```cython
cdef extern:
    # note the final "_", added to the name by gfortran
    void overlap_save_(double* x, long* Nx, double* h, long* M, long*N, double* y)
    void overlap_save_iteration_(...)

cpdef double [:] overlap_save(double [:] x, double [:] h, long N):
    #  [ wrapper preparation ]
    # then call fortran function
    overlap_save_(&x[0], &Nx, &h[0], &M, &N, &y[0])
    return y

# [Other wrappers ...]
```
In order to compile an extension with these code we write a `setup.py` file with this code appended to the list of extensions:
```python
Extension(
    name = "linconv",
    sources = ["src/cython/linconv.pyx"],
    libraries = ["gfortran"],
    extra_objects = ["build/cern/fft.o",
                     "build/stdc/linconv.o"]
)
```
By calling: `python3 setup.py build_ext -b stdc/ --only linconv` we get this extension compiled.

This workflow is summarised by the command `make stdc/linconv.cpython-310-x86_64-linux-gnu.so`

From the project root directory, it is possible to import the module extension by adding `from stdc import linconv` to the preable of your python file.

### `stdc.tdomain`

This submodule offers the `filteriir` function, to perform iir filtering given the parameters of the $\mathcal{Z}$ transform of the unit sample response of the filter.

As for the `stdc.linconv` module, we the code in `src/fortran/libstdc/filter(butterworth2)iir.f90`, 
```fortran
subroutine filteriir(b, Nb, a, Na, x, Nx, y, y_initial)
  ! [body of the routine omitted]
end subroutine filteriir
```
which is compiled to `build/stdc/filter(butterworth2)iir.o` with 
```shell
gfortran -O2 -c src/fortran/libstdc/filter(butterworth2)iir.f90 -fpic -o build/stdc/filter(butterworth2)iir.o
```

The function symbols in the produced object files are referenced in the `src/cython/tdomain.pyx`, similarly as above:
```cython
cdef extern:
    void filteriir_(double* b, long* Nb,
                    double* a, long* Na,
                    double* x, long* Nx,
                    double* y,
                    double* y_initial)
    
cpdef double[:] filteriir(double[:] b, double[:] a, double[:] x):
    # [wrapper setup omitted]
    filteriir_(&b[0], &Nb, &a[0], &Na, &x[0], &Nx, &y[0], NULL)
    return y
```
and linked by the cython compiler, called from `setup.py`:
```python
Extension(
    name = "tdomain",
    sources = ["src/cython/tdomain.pyx"],
    extra_objects = ["build/stdc/filteriir.o",
                     "build/stdc/filterbutterworth2iir.o"]
)
```
then compiled with: `python3 setup.py build_ext -b stdc/ --only linconv`.
Again it is faster to use `make stdc/tdomain.cpython-310-x86_64-linux-gnu.so`

From the project root directory, it is possible to import the module extension by adding `from stdc import linconv` to the preable of your python file.

### `stdc.dummy`: a slightly different story

The sole purpose of this module is to demonstrate a different approach: instead of linking the objects files together into `stdc/tdomain.cpython-310-x86_64-linux-gnu.so`, this shared object will only contain the wrappers, while the referenced fortran functions will be retrieved from an external runtime library.

#### Preparation of the external library

First we compile the external library with:
1. `gfortran -c src/fortran/dummy.f90 -o build/dummy.o`
2. `gfortran -shared build/dummy.o -o lib/libdummy.so`

#### Linking the external module to the runtime library

In the `setup.py`:
```python
Extension(
    name = "dummy",
    sources = ["src/cython/dummy.pyx"],
    libraries = ["dummy"], # will assume a runtime library file named 'libdummy.so'
    library_dirs=["lib"]
)
```

#### Correctly import the module

In order to correctly import the module, python must be executed with the `lib` directory included in the `LD_LIBRARY_PATH` environment variable.
For instance compare:
1. `python3 -c "from stdc import dummy" # Yields an error` 
2. `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd) python3 -c "from stdc import dummy" # succeeds`

As already mentioned option 2 is wrapped by the shell script `pywithlib.sh`. For instance try:
```shell
$ ./pywithlib.sh -c "from stdc import dummy; import numpy as np; dummy.dummy_from_fortran(); print(dummy.mult_ints_from_fortran(5,7)); print(dummy.inner_prod_f(np.array([1.0,2,3]), np.array([3,2,1.])))"
 Hello from fortran code
Hello
35
10.0
```
