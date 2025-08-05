This directory contains tests for the modernized and improved Fortran implementation of Powell's
derivative-free optimization solvers.

The major tests are:

"$COMPILER"test_i"$INTEGER_KIND"."$SOLVER",

where

COMPILER is g, i, n, f, v, d, or x, which stands for gfortran, ifort, nagfor, flang, nvfortran,
aocc flang, or ifx, respectively;
INTEGER_KIND is 2, 4, or 8, which stands for 16-bit, 32-bit, or 64-bit integers, respectively;
SOLVER is cobyla, uobyqa, newuoa, bobyqa, or lincoa.

For example, try

make clean && make gtest_i2.cobyla

N.B.: The scripts in this directory are written for testing and developing purposes. They do NOT
intend to illustrate the usage of the package. For such illustrations, see the "examples" directory.

Coded by Zaikun ZHANG (www.zhangzk.net).

Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).

Started in October 2021.
