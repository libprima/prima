# This Makefile intends to test the modern-Fortran version of Powell's solvers on as many compilers
# as possible.
#
# Coded by Zaikun ZHANG (www.zhangzk.net).
#
# Started: September 2021
#
# Last Modified: September 13, 2021


#$(VERBOSE).SILENT:

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
TESTSRC != cat ./ffiles.txt
COMMONSRC != sed "s/^/\.\/common\//" ../common/ffiles.txt
NEWUOASRC != sed "s/^/\.\/newuoa\//" ../newuoa/ffiles.txt
# Tests.
TESTS = ftest
# Log file.
$(foreach TST, $(TESTS), $(eval $(TST): TESTLOG = $(TST).log))

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

# LLVM Flang
FFORT = flang -std=f$(FSTD) -Mstandard -Wall -Wextra
f_i2_r4_d1_tst: FC = $(FFORT)


####################################################################################################
# Making a compiler-specific test
%test: %_i2_r4_d1_tst
	./checktest --warning $(TESTLOG)
	@printf "\n$@ ends!\n\n" >> $(TESTLOG)
	@printf "\n$@ ends!\n\n"
	@make clean

%_tst: src header_% $(COMMONSRC) $(NEWUOASRC) $(TESTSRC)
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
ftest: QPAVLB = 0
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
