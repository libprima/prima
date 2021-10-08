# This Makefile intends to test the modern-Fortran version of Powell's solvers on as many compilers
# as possible.
#
# The following compilers are tested:
# a: Absoft af95
# d: AOCC flang
# f: LLVM flang
# g: GNU gfortran
# n: NAG nagfor
# i: Intel ifort
# s: Oracle sunf95
# v: NVIDIA nvfortran
# x: Intel ifx
#
# The following tests are available, where C = a, d, f, g, n, i, s, v, x, and N = 2, 4, 8:
# Ctest.S: test solver S with compiler X
# Ctest.S_c: test solver S with compiler X; compilation only (do not run the binary)
# Ctest_iN: test solver S with compiler X and integer kind INT(N*8)
# Ctest_iN_c: test solver S with compiler X and integer kind INT(N*8), compilation only
#
# Coded by Zaikun ZHANG (www.zhangzk.net).
#
# Started: September 2021
#
# Last Modified: October 4, 2021
####################################################################################################


.PHONY: all test clean *.cobyla *.uobyqa *.newuoa *.bobyqa *.lincoa

MAKEFILE_DIR = ./makefiles

%.cobyla: SOLVER = cobyla
%.uobyqa: SOLVER = uobyqa
%.newuoa: SOLVER = newuoa
%.bobyqa: SOLVER = bobyqa
%.lincoa: SOLVER = lincoa

%test:
	@make $@.cobyla $@.uobyqa $@.newuoa $@.bobyqa $@.lincoa

%test_c:
	@make $@.cobyla $@.uobyqa $@.newuoa $@.bobyqa $@.lincoa

%.cobyla %.uobyqa %.newuoa %.bobyqa %.lincoa:
	mkdir -p test.$(SOLVER)
	cp $(MAKEFILE_DIR)/Makefile.common test.$(SOLVER)/
	cp $(MAKEFILE_DIR)/Makefile.$(SOLVER) test.$(SOLVER)/Makefile
	cd test.$(SOLVER)/; make $(subst .$(SOLVER),,$@)

clean:
	@rm -rf log test.cobyla test.uobyqa test.newuoa test.bobyqa test.lincoa
