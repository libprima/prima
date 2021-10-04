################################################################
# This Makefile tests the following compilers on empty arrays.
#
# af95: Absoft 64-bit Pro Fortran 21.0.0
# gfortran: GNU Fortran (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0
# ifort: ifort (IFORT) 2021.2.0 20210228
# nagfor: NAG Fortran Compiler Release 7.0(Yurakucho) Build 7036
# pgfortran: pgfortran (aka nvfortran) 21.3-0 LLVM 64-bit x86-64
# sunf95: Oracle Developer Studio 12.6
# g95: G95 (GCC 4.0.3 (g95 0.94!) Jan 17 2013)
#
# Tested on Ubuntu 20.04 with Linux 5.4.0-77-generic x86_64

.PHONY: test clean

test:
	make -s gtest
	make -s itest
	make -s ntest
	make -s ptest
	make -s stest
	make -s 9test
	make -s atest

gtest: FC = gfortran -Wall -Wextra -fcheck=all
itest: FC = ifort -warn all -check all
ntest: FC = nagfor -C
ptest: FC = pgfortran -C -Mbounds
stest: FC = sunf95 -w3 -xcheck=%all -C
9test: FC = g95 -Wall -Wextra -fbounds-check
atest: FC = af95 -no-pie -et -Rb

%test: testempty.f90
	$(FC) -g -O0 -o $@ $<
	./$@

clean:
	rm -f *.o *.mod *.dbg *test
################################################################
