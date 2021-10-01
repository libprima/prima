# This Makefile intends to test the modern-Fortran version of Powell's solvers on as many compilers
# as possible.
#
# Coded by Zaikun ZHANG (www.zhangzk.net).
#
# Started: September 2021
#
# Last Modified: September 13, 2021


#$(VERBOSE).SILENT:

SHELL = /bin/bash

####################################################################################################
# Variables

# Fortran standard to follow. We aim to make the code compatible with F2003, F2008, and F2018.
FS = 03
FSTD = 20$(FS)

# Default options for all the compilers.
FFLAGS = -g -O0

# Log file stamp.
RANDNUM := $(shell awk 'BEGIN{srand();printf("%d", 65536*rand())}')
TIME := $(shell date +%Y.%m.%d_%H.%M.%S)
STAMP := $(TIME)_$(RANDNUM)

# Log directory.
LOG_DIR := ./log/
# Trash directory.
TRASH_DIR := ./trash/

# Common directories.
COMMON_DIR = ./common/
# Solver directories.
NEWUOA_DIR = ./newuoa/
# Headers.
HEADERS = $(COMMON_DIR)/*.h

# Sources.
SRC_DIRS = ../common/ ../newuoa/
TESTSRC = $(shell cat ./ffiles.txt)
COMMONSRC = $(shell sed "s/^/\.\/common\//" ../common/ffiles.txt)
NEWUOASRC = $(shell sed "s/^/\.\/newuoa\//" ../newuoa/ffiles.txt)

# Tests.
TESTS = ntest atest dtest ftest gtest itest stest vtest xtest
TESTS_QP = atest gtest ntest itest stest
TESTS_NO_QP = dtest ftest vtest xtest
$(foreach TST, $(TESTS_QP), $(eval $(TST): SUBTESTS := \
	$(TST)_i2_r4_d1_tst $(TST)_i4_r4_d1_tst $(TST)_i8_r4_d1_tst \
	$(TST)_i2_r4_d0_tst $(TST)_i4_r4_d0_tst $(TST)_i8_r4_d0_tst \
	$(TST)_i2_r8_d1_tst $(TST)_i4_r8_d1_tst $(TST)_i8_r8_d1_tst \
	$(TST)_i2_r8_d0_tst $(TST)_i4_r8_d0_tst $(TST)_i8_r8_d0_tst \
	$(TST)_i2_r16_d0_tst $(TST)_i4_r16_d0_tst $(TST)_i8_r16_d0_tst \
	$(TST)_i2_r16_d1_tst $(TST)_i4_r16_d1_tst $(TST)_i8_r16_d1_tst \
	))
$(foreach TST, $(TESTS_NO_QP), $(eval $(TST): SUBTESTS := \
	$(TST)_i2_r4_d1_tst $(TST)_i4_r4_d1_tst $(TST)_i8_r4_d1_tst \
	$(TST)_i2_r4_d0_tst $(TST)_i4_r4_d0_tst $(TST)_i8_r4_d0_tst \
	$(TST)_i2_r8_d1_tst $(TST)_i4_r8_d1_tst $(TST)_i8_r8_d1_tst \
	$(TST)_i2_r8_d0_tst $(TST)_i4_r8_d0_tst $(TST)_i8_r8_d0_tst \
	))

#.PHONY: all test $(TESTS) $(SUBTESTS) source header_* clean*
.PHONY: all test source_* header_* clean*

# Operating system.
OSTYPE :=
ifeq ($(OS),Windows_NT)
	OSTYPE = WINDOWS
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		OSTYPE = MAC
	else  # We assume the OS is LINUX by default.
		OSTYPE = LINUX
	endif
endif

# Define SED.
# When calling "sed -i" on macOS, it is obligatory to specify a string (e.g., .bak) after -i as the
# extension for saving a backup. If the string is "", then no backup will be saved. If no string is
# specified, then an error will be raised, saying "invalid command code".
SED :=
ifeq ($(OSTYPE), MAC)
	SED = @sed -i ""
else
	SED = @sed -i
endif


####################################################################################################
# Here are the testing options for each compiler. We try to turn on as many debugging options as
# possible, in particular, the options regarding the following are included whenever available.
# 1. Standard conformity.
# 2. Runtime bound checking of arrays.
# 3. IEEE 754 Floating-Point Exceptions (FPE).

# Absoft af95
# For af95, -Rb imposes bound checking for arrays. But it is unhappy with zero-sized arrays as of
# Absoft Pro 21.0. Thus -Rb is not included in the following options.
AFORT = af95 -no-pie -m1 -en -et -Rc -Rs -Rp
atest_i2_r4_d1_tst atest_i4_r4_d1_tst atest_i8_r4_d1_tst atest_i2_r4_d0_tst atest_i4_r4_d0_tst atest_i8_r4_d0_tst: \
	FC = $(AFORT) -TENV:simd_zmask=off
atest_i2_r8_d1_tst atest_i4_r8_d1_tst atest_i8_r8_d1_tst atest_i2_r8_d0_tst atest_i4_r8_d0_tst atest_i8_r8_d0_tst: \
	FC = $(AFORT) \
	-TENV:simd_zmask=off -TENV:simd_omask=off -TENV:simd_imask=off #-TENV:simd_dmask=off -TENV:simd_umask=off
atest_i2_r16_d1_tst atest_i4_r16_d1_tst atest_i8_r16_d1_tst atest_i2_r16_d0_tst atest_i4_r16_d0_tst atest_i8_r16_d0_tst: \
	FC = $(AFORT) -TENV:simd_zmask=off

# AMD AOCC Flang
DFORT = /opt/AMD/aocc-compiler-3.1.0/bin/flang -std=f$(FSTD) -Mstandard -Wall -Wextra
dtest_i2_r4_d1_tst dtest_i4_r4_d1_tst dtest_i8_r4_d1_tst dtest_i2_r4_d0_tst dtest_i4_r4_d0_tst dtest_i8_r4_d0_tst: \
	FC = $(DFORT)
dtest_i2_r8_d1_tst dtest_i4_r8_d1_tst dtest_i8_r8_d1_tst dtest_i2_r8_d0_tst dtest_i4_r8_d0_tst dtest_i8_r8_d0_tst: \
	FC = $(DFORT) -ffp-exception-behavior=strict
dtest_i2_r16_d1_tst dtest_i4_r16_d1_tst dtest_i8_r16_d1_tst dtest_i2_r16_d0_tst dtest_i4_r16_d0_tst dtest_i8_r16_d0_tst: \
	FC = $(DFORT) -ffp-exception-behavior=strict

# LLVM Flang
FFORT = flang -std=f$(FSTD) -Mstandard -Wall -Wextra
ftest_i2_r4_d1_tst ftest_i4_r4_d1_tst ftest_i8_r4_d1_tst ftest_i2_r4_d0_tst ftest_i4_r4_d0_tst ftest_i8_r4_d0_tst: \
	FC = $(FFORT)
ftest_i2_r8_d1_tst ftest_i4_r8_d1_tst ftest_i8_r8_d1_tst ftest_i2_r8_d0_tst ftest_i4_r8_d0_tst ftest_i8_r8_d0_tst: \
	FC = $(FFORT)
ftest_i2_r16_d1_tst ftest_i4_r16_d1_tst ftest_i8_r16_d1_tst ftest_i2_r16_d0_tst ftest_i4_r16_d0_tst ftest_i8_r16_d0_tst: \
	FC = $(FFORT)

# GNU gfortran
# In the debug mode, our code includes the BACKTRACE function if the compiler is gfortran. This
# makes the code not standard-conforming. So -std=f$(FSTD) is excluded for the following options.
GFORT = gfortran -Wall -Wextra -pedantic -Wampersand -Wconversion  -Wuninitialized \
	-Wmaybe-uninitialized -Wsurprising -Waliasing  -Wimplicit-interface -Wimplicit-procedure \
	-Wintrinsics-std -Wunderflow -Wuse-without-only -Wrealloc-lhs -Wrealloc-lhs-all -Wdo-subscript \
	-Wunused-parameter -fPIC -fimplicit-none -fbacktrace -fcheck=all \
	-finit-real=nan -finit-integer=-9999999
gtest_i2_r4_d1_tst gtest_i4_r4_d1_tst gtest_i8_r4_d1_tst gtest_i2_r4_d0_tst gtest_i4_r4_d0_tst gtest_i8_r4_d0_tst: \
	FC = $(GFORT) -ffpe-trap=zero
gtest_i2_r8_d1_tst gtest_i4_r8_d1_tst gtest_i8_r8_d1_tst gtest_i2_r8_d0_tst gtest_i4_r8_d0_tst gtest_i8_r8_d0_tst: \
	FC = $(GFORT) -ffpe-trap=zero,invalid,overflow,#,underflow,denorm
gtest_i2_r16_d1_tst gtest_i4_r16_d1_tst gtest_i8_r16_d1_tst gtest_i2_r16_d0_tst gtest_i4_r16_d0_tst gtest_i8_r16_d0_tst: \
	FC = $(GFORT) -ffpe-trap=zero,invalid,overflow,#,underflow,denorm

# Intel ifort
IFORT = ifort -stand f$(FS) -warn all -check all -debug extended -fimplicit-none \
	-traceback -debug-parameters all -fp-stack-check
itest_i2_r4_d1_tst itest_i4_r4_d1_tst itest_i8_r4_d1_tst itest_i2_r4_d0_tst itest_i4_r4_d0_tst itest_i8_r4_d0_tst: \
	FC = $(IFORT) -fp-trap=divzero
itest_i2_r8_d1_tst itest_i4_r8_d1_tst itest_i8_r8_d1_tst itest_i2_r8_d0_tst itest_i4_r8_d0_tst itest_i8_r8_d0_tst: \
	FC = $(IFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict
itest_i2_r16_d1_tst itest_i4_r16_d1_tst itest_i8_r16_d1_tst itest_i2_r16_d0_tst itest_i4_r16_d0_tst itest_i8_r16_d0_tst: \
	FC = $(IFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict

# NAG nagfor
NFORT = nagfor -colour=error:red,warn:magenta,info:yellow \
	-f$(FSTD) -info -gline -u -C -C=alias -C=dangling -C=intovf -C=undefined -kind=unique \
	-Warn=allocation -Warn=constant_coindexing -Warn=subnormal
ntest_i2_r4_d1_tst ntest_i4_r4_d1_tst ntest_i8_r4_d1_tst ntest_i2_r4_d0_tst ntest_i4_r4_d0_tst ntest_i8_r4_d0_tst: \
	FC = $(NFORT) -nan -ieee=full
ntest_i2_r8_d1_tst ntest_i4_r8_d1_tst ntest_i8_r8_d1_tst ntest_i2_r8_d0_tst ntest_i4_r8_d0_tst ntest_i8_r8_d0_tst: \
	FC = $(NFORT) -nan -ieee=stop
ntest_i2_r16_d1_tst ntest_i4_r16_d1_tst ntest_i8_r16_d1_tst ntest_i2_r16_d0_tst ntest_i4_r16_d0_tst ntest_i8_r16_d0_tst: \
	FC = $(NFORT) -nan -ieee=stop

# NVIDIA nvfortran (aka, pgfortran)
VFORT = nvfortran -C -Mstandard -Minform=warn -Mbounds -Mchkstk -Mchkptr
vtest_i2_r4_d1_tst vtest_i4_r4_d1_tst vtest_i8_r4_d1_tst vtest_i2_r4_d0_tst vtest_i4_r4_d0_tst vtest_i8_r4_d0_tst: \
	FC = $(VFORT) -Kieee -Ktrap=divz
vtest_i2_r8_d1_tst vtest_i4_r8_d1_tst vtest_i8_r8_d1_tst vtest_i2_r8_d0_tst vtest_i4_r8_d0_tst vtest_i8_r8_d0_tst: \
	FC = $(VFORT) -Kieee -Ktrap=divz,ovf,inv#,unf,denorm
vtest_i2_r16_d1_tst vtest_i4_r16_d1_tst vtest_i8_r16_d1_tst vtest_i2_r16_d0_tst vtest_i4_r16_d0_tst vtest_i8_r16_d0_tst: \
	FC = $(VFORT) -Kieee -Ktrap=divz,ovf,inv#,unf,denorm

# Oracle sunf95
SFORT = sunf95 -w3 -u -U -ansi -xcheck=%all -C
stest_i2_r4_d1_tst stest_i4_r4_d1_tst stest_i8_r4_d1_tst stest_i2_r4_d0_tst stest_i4_r4_d0_tst stest_i8_r4_d0_tst: \
	FC = $(SFORT) -fnonstd -ftrap=division
stest_i2_r8_d1_tst stest_i4_r8_d1_tst stest_i8_r8_d1_tst stest_i2_r8_d0_tst stest_i4_r8_d0_tst stest_i8_r8_d0_tst: \
	FC = $(SFORT) \
	-fnonstd -ftrap=overflow,division,invalid#,underflow
stest_i2_r16_d1_tst stest_i4_r16_d1_tst stest_i8_r16_d1_tst stest_i2_r16_d0_tst stest_i4_r16_d0_tst stest_i8_r16_d0_tst: \
	FC = $(SFORT) \
	-fnonstd -ftrap=overflow,division,invalid#,underflow

# Intel ifx
XFORT = ifx -stand f$(FS) -warn all -check all -debug extended -fimplicit-none \
	-traceback -debug-parameters all -fp-stack-check
xtest_i2_r4_d1_tst xtest_i4_r4_d1_tst xtest_i8_r4_d1_tst xtest_i2_r4_d0_tst xtest_i4_r4_d0_tst xtest_i8_r4_d0_tst: \
	FC = $(XFORT) -fp-trap=divzero
xtest_i2_r8_d1_tst xtest_i4_r8_d1_tst xtest_i8_r8_d1_tst xtest_i2_r8_d0_tst xtest_i4_r8_d0_tst xtest_i8_r8_d0_tst: \
	FC = $(XFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict
xtest_i2_r16_d1_tst xtest_i4_r16_d1_tst xtest_i8_r16_d1_tst xtest_i2_r16_d0_tst xtest_i4_r16_d0_tst xtest_i8_r16_d0_tst: \
	FC = $(XFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict


####################################################################################################
# All the tests.
test:
	 make $(TESTS)

####################################################################################################
# Make a compiler-specific test.
# N.B.: Execute checktest by  "bash ./checktest" instead of just "./checktest", because
# "./checktest" may not be executable on all operating systems.
%test:
	@mkdir -p $(LOG_DIR)
	@printf "\n$@ starts.\n\n" | tee -a "$(LOG_DIR)/$@_$(STAMP).log"
	@make $(SUBTESTS)
	bash ./checktest --warning "$(LOG_DIR)/$@_$(STAMP).log"
	@printf "\n$@ ends!.\n\n" | tee -a "$(LOG_DIR)/$@_$(STAMP).log"
	@make clean

# Make a compiler-specific subtest.
# N.B.: The prerequisite should NOT include files that may not exist, e.g., $(COMMONSRC), because
# Make will ignore the rule if such files do not exist, resulting in "No rule to make ..." mistake.
%_tst: source_% header_%
	@mkdir -p $(LOG_DIR)
	@printf "\n$@ starts.\n\n" | tee -a "$(LOG_DIR)/$@_$(STAMP).log"
	$(FC) $(FFLAGS) -o $@ $(COMMONSRC) $(NEWUOASRC) $(TESTSRC) 2>&1 | tee -a "$(LOG_DIR)/$@_$(STAMP).log"
	#./$@ 2>&1 | tee -a "$(LOG_DIR)/$@_$(STAMP).log"
	bash ./checktest --error "$(LOG_DIR)/$@_$(STAMP).log"
	bash ./checktest --warning "$(LOG_DIR)/$@_$(STAMP).log"
	@printf "\n$@ ends at $(shell date +%Y.%m.%d_%H.%M.%S).\n\n" | tee -a "$(LOG_DIR)/$@_$(STAMP).log"
	@cat "$(LOG_DIR)/$@_$(STAMP).log" >> "$(LOG_DIR)/$(shell echo $@ | sed "s/_.*//")_$(STAMP).log"
	@rm "$(LOG_DIR)/$@_$(STAMP).log"
	@make cleanbin

# Make a local copy of the source code for the test.
source_%:
	@printf "\nMaking a copy of the source code for the test.\n"
	@for DIR in $(SRC_DIRS); do cp -r "$$DIR" ./; done
	@printf "\nDone.\n\n"

# Adapt the header file for the test.
# For REAL128:
# Flang and nvfortran do not support REAL128;
# AOCC Flang complains about a symbol lookup error: undefined symbol: "fort_rnumq_i8";
# ifx encounters an internal error when compiling with REAL128.
header_%: IK = $(shell expr 8 \* $$(echo $@ | sed "s/.*_i//" | sed "s/_.*//"))
header_%: RP = $(shell expr 8 \* $$(echo $@ | sed "s/.*_r//" | sed "s/_.*//"))
header_%: DBG = $(shell echo $@ | sed "s/.*_d//" | sed "s/_.*//")
header_%: QPAVLB = $(shell expr $$(echo $@ | sed "s/.*_r//" | sed "s/_.*//") / 16)
header_%: source_%
	@printf "\nAdapting the header file for the test.\n"
	$(SED) "0,/^#define __QP_AVAILABLE__ [0-9]*/s//#define __QP_AVAILABLE__ $(QPAVLB)/" $(HEADERS)
	$(SED) "0,/^#define __INTEGER_KIND__ [0-9]*/s//#define __INTEGER_KIND__ $(IK)/" $(HEADERS)
	$(SED) "0,/^#define __REAL_PRECISION__ [0-9]*/s//#define __REAL_PRECISION__ $(RP)/" $(HEADERS)
	$(SED) "0,/^#define __DEBUGGING__ [0-9]*/s//#define __DEBUGGING__ $(DBG)/" $(HEADERS)
	@printf "\nDone.\n\n"


####################################################################################################
# Cleaning up.
cleanbin:
	@printf "\nCleaning up the binaries ...\n"
	@rm -f *.o *.mod *.dbg *.cmdx *.cmod *.ilm *.stb
	@rm -f *_tst *_tst_cobyla *_tst_uobyqa *_tst_newuoa *_tst_bobyqa *_tst_lincoa *_output.txt
	@printf "\nDone.\n\n"

cleansrc:
	@printf "\nCleaning up the source code ...\n"
	@if [ `ls -1 $(COMMON_DIR) 2>/dev/null | wc -l ` -gt 0 ] ; \
		then mkdir -p $(TRASH_DIR); \
		rm -rf $(TRASH_DIR)/$(COMMON_DIR); \
		mv $(COMMON_DIR) $(TRASH_DIR) > /dev/null 2>&1; fi
	@if [ `ls -1 $(NEWUOA_DIR) 2>/dev/null | wc -l ` -gt 0 ] ; \
		then mkdir -p $(TRASH_DIR); rm -rf $(TRASH_DIR)/$(NEWUOA_DIR); \
		mv $(NEWUOA_DIR) $(TRASH_DIR) > /dev/null 2>&1; fi
	@printf "\nDone.\n\n"

clean: cleanbin cleansrc

cleanall:
	@make clean
	rm -rf $(TRASH_DIR) $(LOG_DIR)
