This directory contains simple examples that illustrate how to use the modernized Fortran
implementation of Powell's derivative-free optimization solvers.

N.B.:

1. See the Makefiles for how to compile the examples.

2. The first example in every folder, example_1, uses the same objective function:
f(x1, x2) = (x1 - 5)**2 + (x2 - 4)**2.
The unconstrained minimizer is obviously (5, 4) with an optimal value of 0. This example uses a
trivial objective function in order to let the user focus on understanding the PRIMA API and usage.

3. The examples assume that the macros in ../common/ppf.h are set to their default values. In
particular, PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0 (default integer).

4. In the Makefiles, we impose Fortran 2018 standard in the compilation. It is our intention to be
compliant with Fortran 2008 and above.

5. As of March 2024, the examples run successfully with the following compilers on Ubuntu 22.04.
- AMD AOCC Flang 4.1.0
- Arm Fortran Compiler 22.1
- Classic Flang 17.0.2
- GNU gfortran 13.1.0
- Intel ifx 2024.0.2
- Intel ifort 2021.11.1
- NAG Fortran Compiler Release 7.1 (Hanzomon) Build 7143
- NVIDIA nvfortran 24.1
The following discontinued compilers are not supported: Absoft af95, g95, Oracle sunf95.

Coded by Zaikun ZHANG (www.zhangzk.net).

Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).

Started in July 2021.
