# SPDX-License-Identifier: GPL-3.0-or-later

from distutils.core import setup, Extension
from Cython.Build import cythonize
import sys



all_extensions = {
    "dummy" : Extension(
        name = "dummy",
        sources = ["src/cython/dummy.pyx"],
        libraries = ["dummy"],
        library_dirs=["lib"]
    ),
    "tdomain" :  Extension(
        name = "tdomain",
        sources = ["src/cython/tdomain.pyx"],
        extra_objects = ["build/stdc/filteriir.o",
                         "build/stdc/filterbutterworth2iir.o"]
    ),
    "linconv" : Extension(
        name = "linconv",
        sources = ["src/cython/linconv.pyx"],
        libraries = ["gfortran"],
        extra_objects = ["build/cern/fft.o",
                         "build/stdc/linconv.o"]
    )
}


# Parse the command-line argument
only_module = None
if '--only' in sys.argv:
    module_index = sys.argv.index('--only') + 1
    if module_index < len(sys.argv):
        only_module = sys.argv.pop(module_index)
        sys.argv.remove('--only')
        
# Select the specific module to build
if only_module:
    ext_modules = [all_extensions[only_module]]
else:
    ext_modules = list(all_extensions.values())

    
setup(
    name = "stdc",
    ext_modules = cythonize(ext_modules, language_level=3)
)
