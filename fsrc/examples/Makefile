# This is the Makefile for the tests.
# Please change the variables FC and FFLAGS to match your fortran compiler.
#
# N.B.:
# The .F90 and .f90 files will be compiled in the enumeration order of
# the .o files. The order matters, because the compilation of each .o file 
# depends on the .o files (and the corresponding .mod files) preceding it.

FC = gfortran -Wall -Wextra -fcheck=all -pedantic
#FC = ifort -traceback -warn all -check all
#FC = flang -Wall -Wextra -Mstandard 
#FC = nagfor -C 
#FC = sunf95 -w3 -u -U
#FC = af95 -no-pie -m1
#FC = pgfortran -Minform=warn
FFLAGS = -g 
DEPS = ../common/*.h 
COMMON = ../common
NEWUOA = ../newuoa

.PHONY: test_newuoa clean

test_newuoa: clean newuoa_exmp
	./newuoa_exmp 
	rm -f newuoa_exmp
 
newuoa_exmp: clean newuoa_exmp.f90 \
	consts.o info.o debug.o memory.o infnan.o lina.o output.o pintrf.o \
	preproc.o initialize.o trustregion.o vlagbeta.o geometry.o shiftbase.o update.o newuob.o newuoa.o 
	$(FC) $(FFLAGS) -o newuoa_exmp newuoa_exmp.f90 *.o
	rm -f *.o *.mod *.dbg

# Compile the Fortran code providing generic modules
%.o: $(COMMON)/%.F90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $< 

# Compile the Fortran code providing solver-specific modules
%.o: $(NEWUOA)/%.f90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $< 

clean:
	rm -f *.o *.mod *.dbg
	rm -f newuoa_exmp
