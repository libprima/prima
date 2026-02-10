This directory contains simple examples that illustrate how to use the modernized Fortran
implementation of Powell's derivative-free optimization solvers.

N.B.:

0. In production, if the dimension of the problem is big (e.g., more than 100), then compilers
should be instructed to compile PRIMA with the automatic arrays allocated on the heap. Otherwise,
those arrays may be allocated on the stack, which may lead to stack overflow. Since PRIMA is
designed to solve problems with expensive function evaluations, we do not worry about the
performance of heap arrays.

The compiler flags for heap arrays are as follows:
- AMD AOCC Flang: -fno-stack-arrays
- AMD AOMP Flang: -fno-stack-arrays -mmlir -fdynamic-heap-array
- Arm Fortran Compiler: -fno-stack-arrays -mmlir -fdynamic-heap-array
- LLVM Flang: -fno-stack-arrays -mmlir -fdynamic-heap-array
- GNU gfortran: -fno-stack-arrays
- Intel ifx: -heap-arrays
- Intel ifort: -heap-arrays
- NVIDIA nvfortran: -Mnostack_arrays
- NAG Fortran Compiler: Not needed. According to the NAG support, "the behaviour of nagfor is for
  small fixed-size arrays to go on the stack, and for variable-size arrays and fixed-size arrays to
  go on the heap".

If ever a segmentation fault occurs, check whether the above flags are used in the compilation.

N.B.: For LLVM Flang and relatives, `-mmlir -fdynamic-heap-array` is needed as of LLVM 21.1.8. See
https://github.com/llvm/llvm-project/issues/88344
https://github.com/zequipe/flang_heap_arrays

1. See the Makefiles for how to compile the examples.

2. The first example in every folder, example_1, uses the same objective function:
f(x1, x2) = (x1 - 5)**2 + (x2 - 4)**2.
The unconstrained minimizer is obviously (5, 4) with an optimal value of 0. This example uses a
trivial objective function in order to let the user focus on understanding the PRIMA API and usage.

3. The examples assume that the macros in ../common/ppf.h are set to their default values. In
particular, PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0 (default integer).

4. In the Makefiles, we impose Fortran 2018 standard in the compilation. It is our intention to be
compliant with Fortran 2008 and above.

5. As of January 2026, the examples run successfully with the following compilers on Ubuntu 24.04.
- AMD AOCC Flang 5.1
- AMD AOMP Flang 22.0
- Arm Fortran Compiler 23.10
- LLVM Flang 21.1
- GNU gfortran 14.2
- Intel ifx 2025.3
- Intel ifort 2021.11.1
- NVIDIA nvfortran 26.1
- NAG Fortran Compiler Release 7.2(Shin-Urayasu) Build 7231
The following discontinued compilers are not supported: Absoft af95, g95, Oracle sunf95.

Coded by Zaikun ZHANG (www.zhangzk.net).

Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).

Started in July 2021.
