FORTDIR := src/fortran
CPYSUFFIX ?= cpython-310-x86_64-linux-gnu.so

.PHONY: stdc libs

ALLMODULES = $(addprefix stdc/,$(addsuffix .$(CPYSUFFIX), dummy tdomain linconv))
stdc: $(ALLMODULES)

stdc/%.$(CPYSUFFIX): src/cython/%.pyx
	python3 setup.py build_ext -b stdc/ --only $(*F)

#####################
# Necessary folders #
#####################
%/:
	mkdir -p $@

#######################
# Module dependancies #
#######################
stdc/dummy.$(CPYSUFFIX)   : lib/libdummy.so
stdc/tdomain.$(CPYSUFFIX) : build/stdc/filteriir.o build/stdc/filterbutterworth2iir.o
stdc/linconv.$(CPYSUFFIX) : lib/libcern.so lib/libstdc.so # build/stdc/linconv.o build/cern/fft.o

#######################
# Recipes for objects #
#######################


build/stdc/%.o : $(FORTDIR)/libstdc/%.f* build/stdc/
	gfortran -O2 -c $< -fpic -o $@
build/cern/%.o : $(FORTDIR)/libcern/%.f* build/cern/
	gfortran -O2 -c $< -fpic -o $@

build/%.o : $(FORTDIR)/%.f* build/
	gfortran -O2 -c $< -fpic -o $@


###############################
# Rules for runtime libraries #
###############################

LIBNAMES = dummy stdc cern
LIBS = $(foreach libname,$(LIBNAMES),lib/lib$(libname).so)
libs : $(LIBS) lib/ build/

lib/libdummy.so: build/dummy.o lib/
	gfortran -shared $< -o $@
lib/libstdc.so : $(patsubst $(FORTDIR)/libstdc/%.f90, build/stdc/%.o, $(wildcard $(FORTDIR)/libstdc/*.f90))
# lib/libcern.so : $(patsubst $(FORTDIR)/libcern/%.f, build/cern/%.o, $(wildcard $(FORTDIR)/libcern/*.f))
lib/libcern.so : build/cern/fft.o


lib/lib%.so : lib/
	gfortran -shared build/$*/*.o -o $@


###########
# Cleanup #
###########

clean:
	rm -rf build/*
	rm -f stdc/*.$(CPYSUFFIX)
	rm -rf stdc/__pycache__ stdc/build
	rm -f lib/*.so
