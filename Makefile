FORTDIR := src/fortran
CPYSUFFIX ?= cpython-310-x86_64-linux-gnu.so

.PHONY: stdc libs

ALLMODULES = $(addprefix stdc/,$(addsuffix .$(CPYSUFFIX), dummy tdomain))
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
stdc/tdomain.$(CPYSUFFIX) : lib/libstdc.so


##############################
# Generic recipe for objects #
##############################

OBJRECIPE = gfortran -c $(wordlist 2, $(words $^),$^) -fpic -o $@
build/%.o : build/ $(FORTDIR)/%.f*
	$(OBJRECIPE)
build/%.o : build/ $(FORTDIR)/lib%/*.f*
	$(OBJRECIPE)


###############################
# Rules for runtime libraries #
###############################

LIBNAMES = dummy stdc cern
LIBS = $(foreach libname,$(LIBNAMES),lib/lib$(libname).so)

libs : $(LIBS) lib/ build/

lib/lib%.so :  build/%.o lib/
	gfortran -shared $< -o $@


###########
# Cleanup #
###########

clean:
	rm -rf build/*
	rm -f stdc/*.$(CPYSUFFIX)
	rm -rf stdc/__pycache__ stdc/build
	rm -f lib/*.so
