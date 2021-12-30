# This Makefile illustrates how to use the modern-Fortran version of Powell's solvers. We intend to
# try as many compilers as possible.
#
# Coded by Zaikun Zhang (www.zhangzk.net).
#
# Started: July 2020
#
# Last Modified: September 13, 2021
#
# N.B.:
# The .F90 and .f90 files will be compiled in the enumeration order of the .o files. The order
# matters, because the compilation of each .o file depends on the .o files (and the corresponding
# .mod files) preceding it.

.PHONY: test_newuoa clean

####################################################################################################
# Variables
# Fortran standard to follow. We aim to make the code compatible with F2003, F2008, and F2018.
FS = 03
FSTD = 20$(FS)
# Default options for all the compilers.
FFLAGS = -O3
# Common directories.
COMMON = ../common
# Headers.
HEADERS = $(COMMON)/*.h
# Solver directories.
UOBYQA = ../uobyqa
NEWUOA = ../newuoa
BOBYQA = ../bobyqa
LINCOA = ../lincoa
COBYLA = ../cobyla

TESTS = test_uobyqa test_newuoa test_bobyqa test_lincoa test_cobyla

$(foreach TST, $(TESTS), $(eval $(TST): SOLVER := $(shell echo $(TST) | sed "s/test_/\.\.\//")))


####################################################################################################
# All the tests
test_newuoa:
	make atest_newuoa
	make dtest_newuoa
	make ftest_newuoa
	make gtest_newuoa
	make itest_newuoa
	make ntest_newuoa
	make stest_newuoa
	make vtest_newuoa
	make xtest_newuoa

####################################################################################################
# Here are the compilers to test.

# Absoft af95
ates%: FC = af95 -no-pie

# AMD AOCC Flang
dtes%: FC = /opt/AMD/aocc-compiler-3.1.0/bin/flang

# LLVM Flang
ftes%: FC = flang

# GNU gfortran
gtes%: FC = gfortran

# Intel ifort
ites%: FC = ifort

# NAG nagfor
ntes%: FC = nagfor

# NVIDIA nvfortran (aka, pgfortran)
vtes%: FC = nvfortran

# Oracle sunf95
stes%: FC = sunf95

# Intel ifx
xtes%: FC = ifx

####################################################################################################
# Making a compiler-specific test
%_newuoa: %_newuoa_exmp
	./$<
	make clean

%_cobyla: %_cobyla_exmp
	./$<
	make clean

####################################################################################################
# Compile the binary needed for a compiler-specific test
%_newuoa_exmp: newuoa_exmp.f90 \
	consts.o info.o debug.o memory.o infnan.o linalg.o \
	ratio.o resolution.o history.o selectx.o checkexit.o output.o string.o preproc.o pintrf.o evaluate.o \
	initialize.o trustregion.o vlagbeta.o geometry.o shiftbase.o update.o newuob.o newuoa.o
	$(FC) $(FFLAGS) -o $@ newuoa_exmp.f90 *.o

%_cobyla_exmp: cobyla_exmp.f90 \
	consts.o info.o debug.o memory.o infnan.o linalg.o \
	ratio.o resolution.o history.o selectx.o checkexit.o output.o string.o preproc.o pintrf.o evaluate.o \
	initialize.o trustregion.o vlagbeta.o geometry.o shiftbase.o update.o newuob.o newuoa.o
	$(FC) $(FFLAGS) -o $@ cobyla_exmp.f90 *.o

# Compile the Fortran code providing generic modules
%.o: $(COMMON)/%.*90 $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

# Compile the Fortran code providing solver-specific modules
%.o: $(SOLVER)/%.f90 $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

####################################################################################################
# Cleaning up.
clean:
	rm -f *.o *.mod *.dbg
	rm -f *_*_exmp
