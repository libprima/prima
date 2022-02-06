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

%.cobyla: SOLVER := cobyla
%.uobyqa: SOLVER := uobyqa
%.newuoa: SOLVER := newuoa
%.bobyqa: SOLVER := bobyqa
%.lincoa: SOLVER := lincoa

# Root directory of this project.
ROOT_DIR := $(CURDIR)/../../..
# Temporary directory. It should not be a subdirectory of $(ROOT_DIR). The ?= assignment assigns the
# value to TEST_DIR if TEST_DIR is not an exported environment variable of the shell.
TEST_DIR ?= /tmp
# Root directory of this test. It is a subdirectory of $(TEST_DIR) and a partial copy of $(ROOT_DIR).
TEST_ROOT := $(TEST_DIR)/neupdfo
# The test directory
MTEST_DIR := $(TEST_ROOT)/matlab/mex_gateways/tests
# The makefile directory
MAKEFILE_DIR := $(MTEST_DIR)/makefiles

# GNU make manual: Recursive make commands should always use the variable MAKE, not the explicit
# command name 'make'.
%test:
	@$(MAKE) $@.cobyla $@.uobyqa $@.newuoa $@.bobyqa $@.lincoa

%test_c:
	@$(MAKE) $@.cobyla $@.uobyqa $@.newuoa $@.bobyqa $@.lincoa

%.cobyla %.uobyqa %.newuoa %.bobyqa %.lincoa:
	@$(eval SOLVER_MTEST_DIR := $(MTEST_DIR)/test.$(SOLVER))  # The test directory for the tested solver
	@mkdir -p $(TEST_ROOT)
	@cp -r $(ROOT_DIR)/fsrc $(TEST_ROOT)/
	@mkdir -p $(TEST_ROOT)/matlab/
	@cp -r $(ROOT_DIR)/matlab/mex_gateways $(TEST_ROOT)/matlab/
	@mkdir -p $(SOLVER_MTEST_DIR)
	@cp $(MAKEFILE_DIR)/Makefile.common $(SOLVER_MTEST_DIR)/
	@cp $(MAKEFILE_DIR)/Makefile.$(SOLVER) $(SOLVER_MTEST_DIR)/Makefile
	@cd $(SOLVER_MTEST_DIR) && $(MAKE) $(subst .$(SOLVER),,$@) --no-print-directory

clean:
	@rm -rf $(TEST_ROOT) || :
