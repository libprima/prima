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
# Common directories.
COMMON = ./common/
# Headers.
HEADERS = $(COMMON)/*.h
# Solver directories.
NEWUOA = ./newuoa/
# Sources.
TESTSRC = $(shell cat ./ffiles.txt)
COMMONSRC = $(shell sed "s/^/\.\/common\//" ../common/ffiles.txt)
NEWUOASRC = $(shell sed "s/^/\.\/newuoa\//" ../newuoa/ffiles.txt)
# Tests.
TESTS = ntest atest dtest ftest gtest itest stest vtest xtest
# Log file.
RANDNUM := $(shell bash -c 'echo $$RANDOM')
$(foreach TST, $(TESTS), $(eval $(TST): \
	TESTLOG = $(TST)_$(shell date +%Y.%m.%d_%H.%M.%S)_$(RANDNUM).log))

#.PHONY: all test ntest atest dtest ftest gtest itest stest vtest xtest src header_* clean
.PHONY: all test src header_* clean


####################################################################################################
# All the tests
test:
	 make $(TESTS)


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
a_i2_r4_d1_tst a_i4_r4_d1_tst a_i8_r4_d1_tst a_i2_r4_d0_tst a_i4_r4_d0_tst a_i8_r4_d0_tst: \
	FC = $(AFORT) -TENV:simd_zmask=off
a_i2_r8_d1_tst a_i4_r8_d1_tst a_i8_r8_d1_tst a_i2_r8_d0_tst a_i4_r8_d0_tst a_i8_r8_d0_tst: \
	FC = $(AFORT) \
	-TENV:simd_zmask=off -TENV:simd_omask=off -TENV:simd_imask=off #-TENV:simd_dmask=off -TENV:simd_umask=off
a_i2_r16_d1_tst a_i4_r16_d1_tst a_i8_r16_d1_tst a_i2_r16_d0_tst a_i4_r16_d0_tst a_i8_r16_d0_tst: \
	FC = $(AFORT) -TENV:simd_zmask=off

# AMD AOCC Flang
DFORT = /opt/AMD/aocc-compiler-3.1.0/bin/flang -std=f$(FSTD) -Mstandard -Wall -Wextra
d_i2_r4_d1_tst d_i4_r4_d1_tst d_i8_r4_d1_tst d_i2_r4_d0_tst d_i4_r4_d0_tst d_i8_r4_d0_tst: \
	FC = $(DFORT)
d_i2_r8_d1_tst d_i4_r8_d1_tst d_i8_r8_d1_tst d_i2_r8_d0_tst d_i4_r8_d0_tst d_i8_r8_d0_tst: \
	FC = $(DFORT) -ffp-exception-behavior=strict
d_i2_r16_d1_tst d_i4_r16_d1_tst d_i8_r16_d1_tst d_i2_r16_d0_tst d_i4_r16_d0_tst d_i8_r16_d0_tst: \
	FC = $(DFORT) -ffp-exception-behavior=strict

# LLVM Flang
FFORT = flang -std=f$(FSTD) -Mstandard -Wall -Wextra
f_i2_r4_d1_tst f_i4_r4_d1_tst f_i8_r4_d1_tst f_i2_r4_d0_tst f_i4_r4_d0_tst f_i8_r4_d0_tst: \
	FC = $(FFORT)
f_i2_r8_d1_tst f_i4_r8_d1_tst f_i8_r8_d1_tst f_i2_r8_d0_tst f_i4_r8_d0_tst f_i8_r8_d0_tst: \
	FC = $(FFORT)
f_i2_r16_d1_tst f_i4_r16_d1_tst f_i8_r16_d1_tst f_i2_r16_d0_tst f_i4_r16_d0_tst f_i8_r16_d0_tst: \
	FC = $(FFORT)

# GNU gfortran
# In the debug mode, our code includes the BACKTRACE function if the compiler is gfortran. This
# makes the code not standard-conforming. So -std=f$(FSTD) is excluded for the following options.
GFORT = gfortran -Wall -Wextra -pedantic -Wampersand -Wconversion  -Wuninitialized \
	-Wmaybe-uninitialized -Wsurprising -Waliasing  -Wimplicit-interface -Wimplicit-procedure \
	-Wintrinsics-std -Wunderflow -Wuse-without-only -Wrealloc-lhs -Wrealloc-lhs-all -Wdo-subscript \
	-Wunused-parameter -fPIC -fimplicit-none -fbacktrace -fcheck=all \
	-finit-real=nan -finit-integer=-9999999
g_i2_r4_d1_tst g_i4_r4_d1_tst g_i8_r4_d1_tst g_i2_r4_d0_tst g_i4_r4_d0_tst g_i8_r4_d0_tst: \
	FC = $(GFORT) -ffpe-trap=zero
g_i2_r8_d1_tst g_i4_r8_d1_tst g_i8_r8_d1_tst g_i2_r8_d0_tst g_i4_r8_d0_tst g_i8_r8_d0_tst: \
	FC = $(GFORT) -ffpe-trap=zero,invalid,overflow,#,underflow,denorm
g_i2_r16_d1_tst g_i4_r16_d1_tst g_i8_r16_d1_tst g_i2_r16_d0_tst g_i4_r16_d0_tst g_i8_r16_d0_tst: \
	FC = $(GFORT) -ffpe-trap=zero,invalid,overflow,#,underflow,denorm

# Intel ifort
IFORT = ifort -stand f$(FS) -warn all -check all -debug extended -fimplicit-none \
	-traceback -debug-parameters all -fp-stack-check
i_i2_r4_d1_tst i_i4_r4_d1_tst i_i8_r4_d1_tst i_i2_r4_d0_tst i_i4_r4_d0_tst i_i8_r4_d0_tst: \
	FC = $(IFORT) -fp-trap=divzero
i_i2_r8_d1_tst i_i4_r8_d1_tst i_i8_r8_d1_tst i_i2_r8_d0_tst i_i4_r8_d0_tst i_i8_r8_d0_tst: \
	FC = $(IFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict
i_i2_r16_d1_tst i_i4_r16_d1_tst i_i8_r16_d1_tst i_i2_r16_d0_tst i_i4_r16_d0_tst i_i8_r16_d0_tst: \
	FC = $(IFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict

# NAG nagfor
NFORT = nagfor -colour=error:red,warn:magenta,info:yellow \
	-f$(FSTD) -info -gline -u -C -C=alias -C=dangling -C=intovf -C=undefined -kind=unique \
	-Warn=allocation -Warn=constant_coindexing -Warn=subnormal
n_i2_r4_d1_tst n_i4_r4_d1_tst n_i8_r4_d1_tst n_i2_r4_d0_tst n_i4_r4_d0_tst n_i8_r4_d0_tst: \
	FC = $(NFORT) -nan -ieee=full
n_i2_r8_d1_tst n_i4_r8_d1_tst n_i8_r8_d1_tst n_i2_r8_d0_tst n_i4_r8_d0_tst n_i8_r8_d0_tst: \
	FC = $(NFORT) -nan -ieee=stop
n_i2_r16_d1_tst n_i4_r16_d1_tst n_i8_r16_d1_tst n_i2_r16_d0_tst n_i4_r16_d0_tst n_i8_r16_d0_tst: \
	FC = $(NFORT) -nan -ieee=stop

# NVIDIA nvfortran (aka, pgfortran)
VFORT = nvfortran -C -Mstandard -Minform=warn -Mbounds -Mchkstk -Mchkptr
v_i2_r4_d1_tst v_i4_r4_d1_tst v_i8_r4_d1_tst v_i2_r4_d0_tst v_i4_r4_d0_tst v_i8_r4_d0_tst: \
	FC = $(VFORT) -Kieee -Ktrap=divz
v_i2_r8_d1_tst v_i4_r8_d1_tst v_i8_r8_d1_tst v_i2_r8_d0_tst v_i4_r8_d0_tst v_i8_r8_d0_tst: \
	FC = $(VFORT) -Kieee -Ktrap=divz,ovf,inv#,unf,denorm
v_i2_r16_d1_tst v_i4_r16_d1_tst v_i8_r16_d1_tst v_i2_r16_d0_tst v_i4_r16_d0_tst v_i8_r16_d0_tst: \
	FC = $(VFORT) -Kieee -Ktrap=divz,ovf,inv#,unf,denorm

# Oracle sunf95
SFORT = sunf95 -w3 -u -U -ansi -xcheck=%all -C
s_i2_r4_d1_tst s_i4_r4_d1_tst s_i8_r4_d1_tst s_i2_r4_d0_tst s_i4_r4_d0_tst s_i8_r4_d0_tst: \
	FC = $(SFORT) -fnonstd -ftrap=division
s_i2_r8_d1_tst s_i4_r8_d1_tst s_i8_r8_d1_tst s_i2_r8_d0_tst s_i4_r8_d0_tst s_i8_r8_d0_tst: \
	FC = $(SFORT) \
	-fnonstd -ftrap=overflow,division,invalid#,underflow
s_i2_r16_d1_tst s_i4_r16_d1_tst s_i8_r16_d1_tst s_i2_r16_d0_tst s_i4_r16_d0_tst s_i8_r16_d0_tst: \
	FC = $(SFORT) \
	-fnonstd -ftrap=overflow,division,invalid#,underflow

# Intel ifx
XFORT = ifx -stand f$(FS) -warn all -check all -debug extended -fimplicit-none \
	-traceback -debug-parameters all -fp-stack-check
x_i2_r4_d1_tst x_i4_r4_d1_tst x_i8_r4_d1_tst x_i2_r4_d0_tst x_i4_r4_d0_tst x_i8_r4_d0_tst: \
	FC = $(XFORT) -fp-trap=divzero
x_i2_r8_d1_tst x_i4_r8_d1_tst x_i8_r8_d1_tst x_i2_r8_d0_tst x_i4_r8_d0_tst x_i8_r8_d0_tst: \
	FC = $(XFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict
x_i2_r16_d1_tst x_i4_r16_d1_tst x_i8_r16_d1_tst x_i2_r16_d0_tst x_i4_r16_d0_tst x_i8_r16_d0_tst: \
	FC = $(XFORT) -ftrapuv -init=snan,array -fpe0 -fpe-all=0 -assume ieee_fpe_flags \
	-fp-trap=divzero,invalid,overflow#,underflow,denormal
	#-no-ftz -fp-model strict


####################################################################################################
# Making a compiler-specific test
%test: %_i2_r4_d1_tst %_i4_r4_d1_tst %_i8_r4_d1_tst \
	%_i2_r4_d0_tst %_i4_r4_d0_tst %_i8_r4_d0_tst \
	%_i2_r8_d1_tst %_i4_r8_d1_tst %_i8_r8_d1_tst \
	%_i2_r8_d0_tst %_i4_r8_d0_tst %_i8_r8_d0_tst \
	%_i2_r16_d1_tst %_i4_r16_d1_tst %_i8_r16_d1_tst \
	%_i2_r16_d0_tst %_i4_r16_d0_tst %_i8_r16_d0_tst
	./checktest --warning $(TESTLOG)
	@printf "\n$@ ends!\n\n" >> $(TESTLOG)
	@printf "\n$@ ends!\n\n"
	make clean

# No idea why, but if we put $(COMMONSRC), $(NEWUOASRC), and $(TESTSRC) as prerequisites for %_tst,
# then the Makefile does not work with GitHub Actions, although it works on a normal machine.
# In both cases, the `make` is `GNU Make 4.2.1`.
%_tst: src header_%
	@printf "\n$@ starts.\n\n"
	@printf "\n$@ starts.\n\n" >> $(TESTLOG)
	$(FC) $(FFLAGS) -o $@ $(COMMONSRC) $(NEWUOASRC) $(TESTSRC) >> $(TESTLOG) 2>&1
	./$@ >> $(TESTLOG) 2>&1
	./checktest --error $(TESTLOG)
	@printf "\n$@ ends!\n\n" >> $(TESTLOG)
	@printf "\n$@ ends!\n\n"

src: ../newuoa/* ../common/* \
	$(NEWUOA)/* $(COMMON)/*
	@cp -r ../newuoa ./
	@cp -r ../common ./

# For REAL128:
# Flang and nvfortran do not support REAL128;
# AOCC Flang complains about a symbol lookup error: undefined symbol: "fort_rnumq_i8";
# ifx encounters an internal error when compiling with REAL128.
dtest ftest vtest xtest: QPAVLB = 0
atest gtest ntest itest stest: QPAVLB = 1
header_%: IK = $(shell expr 8 \* $$(echo $@ | sed "s/.*_i//" | sed "s/_.*//"))
header_%: RP = $(shell expr 8 \* $$(echo $@ | sed "s/.*_r//" | sed "s/_.*//"))
header_%: DBG = $(shell echo $@ | sed "s/.*_d//" | sed "s/_.*//")
header_%: $(HEADERS)
	@sed -i '0,/^#define __QP_AVAILABLE__ [0-9]*/s//#define __QP_AVAILABLE__ $(QPAVLB)/' $(HEADERS)
	@sed -i '0,/^#define __INTEGER_KIND__ [0-9]*/s//#define __INTEGER_KIND__ $(IK)/' $(HEADERS)
	@sed -i '0,/^#define __REAL_PRECISION__ [0-9]*/s//#define __REAL_PRECISION__ $(RP)/' $(HEADERS)
	@sed -i '0,/^#define __DEBUGGING__ [0-9]*/s//#define __DEBUGGING__ $(DBG)/' $(HEADERS)


####################################################################################################
# Cleaning up.
clean:
	rm -f *.o *.mod *.dbg *.cmdx *.cmod *.ilm *.stb
	rm -f *_tst *_tst_cobyla *_tst_uobyqa *_tst_newuoa *_tst_bobyqa *_tst_lincoa
	mkdir -p log
	@if [ `ls -1 *.log 2>/dev/null | wc -l ` -gt 0 ]; then mv *.log log > /dev/null 2>&1 ; fi
