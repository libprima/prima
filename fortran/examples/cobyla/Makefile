# This Makefile illustrates how to use the modern-Fortran version of Powell's solvers. We intend to
# try as many compilers as possible.
#
# Coded by Zaikun Zhang (www.zhangzk.net).
#
# Started: July 2020
#
# Last Modified: February 13, 2026
#
# N.B.:
# The .F90 and .f90 files will be compiled in the enumeration order of the .o files. The order
# matters, because the compilation of each .o file depends on the .o files (and the corresponding
# .mod files) preceding it.


.PHONY: test clean

# All the tests
test:
	$(MAKE) dtest
	$(MAKE) ftest
	$(MAKE) gtest
	$(MAKE) itest
	$(MAKE) mtest
	$(MAKE) ntest
	$(MAKE) vtest
	$(MAKE) xtest
	$(MAKE) rtest

# Variables
SOLVER := $(shell basename $(CURDIR))
# Define the Fortran standard to follow. We aim to make the code compatible with F2008 and above.
# Make sure that your compiler supports the selected standard.
FS ?= 18  # Set FS if it does not have a value.
FSTD := 20$(FS)
# Default options for all the compilers.
FFLAGS := -O -g
# Common directories.
COMMON := ../../common
# Headers.
HEADERS := $(COMMON)/*.h
# Solver source files.
SOLVER_SRC := ../../$(SOLVER)
# Overridable selection for the program to compile.
EXAMPLE_NUM ?= 1
EXAMPLE_SRC := $(SOLVER)_example_$(EXAMPLE_NUM).f90

# Include the common Makefile defining the compilers and their options.
include ../Makefile.common


####################################################################################################
# Compile the binary needed for a compiler-specific test
%test: $(EXAMPLE_SRC) \
	consts.o infos.o debug.o huge.o inf.o infnan.o memory.o string.o linalg.o powalg.o \
	ratio.o redrho.o history.o selectx.o checkexit.o fprint.o message.o preproc.o pintrf.o evaluate.o \
	update.o initialize.o trustregion.o geometry.o cobylb.o cobyla.o
	$(FC) $(FFLAGS) -o $@ $(EXAMPLE_SRC) *.o
	./$@
	@$(MAKE) clean

# Compile the Fortran code providing common modules
%.o: $(COMMON)/%.*90 $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

# Compile the Fortran code providing solver-specific modules
%.o: $(SOLVER_SRC)/%.f90 $(HEADERS)
	$(FC) $(FFLAGS) -c -o $@ $<

####################################################################################################
# Cleaning up.
clean:
	@rm -f *.o *.mod *.dbg *.cmdx *.cmod *.ilm *.stb
	@rm -f *test
