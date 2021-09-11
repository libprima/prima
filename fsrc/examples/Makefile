# This is the Makefile for the tests.
# Please change the variables FC and FFLAGS to match your fortran compiler.
#
# N.B.:
# The .F90 and .f90 files will be compiled in the enumeration order of
# the .o files. The order matters, because the compilation of each .o file
# depends on the .o files (and the corresponding .mod files) preceding it.

FFLAGS = -g
DEPS = ../common/*.h
COMMON = ../common
NEWUOA = ../newuoa
TEST_NEWUOA = atest_newuoa ftest_newuoa gtest_newuoa itest_newuoa ntest_newuoa ptest_newuoa stest_newuoa xtest_newuoa

.PHONY: test_newuoa clean

test_newuoa: $(TEST_NEWUOA)

atest_newuoa: FC = af95 -no-pie -m1
ftest_newuoa: FC = flang -Wall -Wextra -Mstandard
gtest_newuoa: FC = gfortran -Wall -Wextra -fcheck=all -pedantic
itest_newuoa: FC = ifort -warn all -check all
ntest_newuoa: FC = nagfor -C -C=alias -C=dangling -C=intovf -C=undefined
ptest_newuoa: FC = pgfortran -Minform=warn
stest_newuoa: FC = sunf95 -w3 -u -U
xtest_newuoa: FC = ifx -warn all

%_newuoa: %_newuoa_exmp
	./$<

%_newuoa_exmp: newuoa_exmp.f90 \
	consts.o info.o debug.o memory.o infnan.o linalg.o \
	ratio.o resolution.o history.o selectx.o checkexit.o output.o preproc.o pintrf.o evaluate.o \
	initialize.o trustregion.o vlagbeta.o geometry.o shiftbase.o update.o newuob.o newuoa.o
	$(FC) $(FFLAGS) -o newuoa_exmp newuoa_exmp.f90 *.o

# Compile the Fortran code providing generic modules
%.o: $(COMMON)/%.F90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: $(COMMON)/%.f90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $<

# Compile the Fortran code providing solver-specific modules
%.o: $(NEWUOA)/%.f90 $(DEPS)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod *.dbg
	rm -f newuoa_exmp
