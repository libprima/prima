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

.PHONY: test clean

####################################################################################################
# Variables
SOLVER := $(shell basename $(PWD))
# Fortran standard to follow. We aim to make the code compatible with F2003, F2008, and F2018.
FS := 03
FSTD := 20$(FS)
# Default options for all the compilers.
FFLAGS := -O3
# Common directories.
COMMON := ../../common
# Headers.
HEADERS := $(COMMON)/*.h
# Solver source files.
SOLVER_SRC := ../../$(SOLVER)

####################################################################################################
# All the tests
test:
	$(MAKE) atest
	$(MAKE) dtest
	$(MAKE) ftest
	$(MAKE) gtest
	$(MAKE) itest
	$(MAKE) ntest
	$(MAKE) stest
	$(MAKE) vtest
	$(MAKE) xtest

####################################################################################################
# Here are the compilers to test.

# Absoft af95
AFORT :=
ifneq ("$(wildcard /etc/debian_version)","")
	AFORT := af95 -no-pie
else
	AFORT := af95
endif
atest: FC := $(AFORT) -m3

# AMD AOCC Flang
AFLANG := $(shell find /opt/AMD \( -type l -o -type f \) -executable -name flang -print \
	| sort | tail -n 1 2> /dev/null || echo AFLANG_NOT_FOUND)
dtest: FC := $(AFLANG) -Wall -Wextra -std=f$(FSTD) -Mstandard

# Flang
ftest: FC := flang -Wall -Wextra -std=f$(FSTD) -Mstandard

# GNU gfortran
gtest: FC := gfortran -Wall -Wextra -Wno-function-elimination -std=f$(FSTD) -fall-intrinsics

# Intel ifort
itest: FC := ifort -warn all -stand f$(FS)

# NAG nagfor
ntest: FC := nagfor -C -f$(FSTD)

# NVIDIA nvfortran (aka, pgfortran)
VFORT := $(shell find /opt/nvidia \( -type l -o -type f \) -executable -name nvfortran -print \
	| sort | tail -n 1 2> /dev/null || echo NVFORTRAN_NOT_FOUND)
vtest: FC := $(VFORT) -C -Wall -Wextra -Mstandard

# Oracle sunf95
stest: FC := sunf95 -w3 -u -U

# Intel ifx
xtest: FC := ifx -warn all -stand f$(FS)

####################################################################################################
# Compile the binary needed for a compiler-specific test
%test: $(SOLVER)_exmp.f90 \
	consts.o info.o debug.o memory.o infnan.o linalg.o \
	ratio.o resolution.o history.o selectx.o checkexit.o output.o string.o preproc.o pintrf.o evaluate.o \
	initialize.o trustregion.o geometry.o update.o cobylb.o cobyla.o
	$(FC) $(FFLAGS) -o $@ $(SOLVER)_exmp.f90 *.o
	./$@
	@$(MAKE) clean

# Compile the Fortran code providing generic modules
%.o: $(COMMON)/%.*90 $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

# Compile the Fortran code providing solver-specific modules
%.o: $(SOLVER_SRC)/%.f90 $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

####################################################################################################
# Cleaning up.
clean:
	@rm -f *.o *.mod *.dbg
	@rm -f *test
