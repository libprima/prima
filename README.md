<h2 align="center">PRIMA: Reference Implementation for Powell's Methods with Modernization and Amelioration</h2>
<p align="center">Dedicated to the late Professor <b><a href="https://www.zhangzk.net/powell.html">M. J. D. Powell</a></b> FRS (1936--2015)</p>

- [What](#what)
- [Why](#why)
- [How](#how)
- [Current status](#current-status)
    - [Modern Fortran](#modern-fortran)
    - [C](#c)
    - [Python](#python)
    - [MATLAB](#matlab)
    - [Julia](#julia)
    - [Other languages](#other-languages)
- [Bug fixes](#bug-fixes)
- [Improvements](#improvements)
- [Who was Powell?](#who-was-powell)
- [A "fun" fact](#a-fun-fact)
- [Acknowledgment](#acknowledgment)
- [Citing PRIMA](#citing-prima)
- [Charityware](#charityware)
- [Contact](#contact)
- [Mirrors](#mirrors)
    - [Gitee](https://gitee.com/libprima/prima)
    - [GitHub](https://github.com/libprima/prima)
    - [GitLab](https://gitlab.com/libprima/prima)
- [Star history](#star-history)


### What

PRIMA is a package for **solving general nonlinear optimization problems without using derivatives**.
It provides the reference implementation for Powell's renowned derivative-free optimization methods, i.e., COBYLA, UOBYQA, NEWUOA, BOBYQA, and LINCOA.
The "P" in the name stands for [**P**owell](https://www.zhangzk.net/powell.html),
and "RIMA" is an acronym for "**R**eference **I**mplementation with **M**odernization and **A**melioration".

PRIMA is part of a research project funded by the
[Hong Kong Research Grants Council](https://www.ugc.edu.hk/eng/rgc) and
the [Department of Applied Mathematics](https://www.polyu.edu.hk/ama) (AMA) at the
[Hong Kong Polytechnic University](https://www.polyu.edu.hk) (PolyU).
The current version is ready to be used [in Fortran](#modern-fortran), [in C](#c),
[in Python](https://github.com/libprima/prima#python),
[in MATLAB](https://github.com/libprima/prima/blob/main/README_mat.md),
and [in Julia](https://juliahub.com/ui/Packages/General/PRIMA).

PRIMA was initiated by [Zaikun Zhang](https://www.zhangzk.net) in July 2020, based on
the [PDFO](https://www.pdfo.net) package by [Tom M. Ragonneau](https://tomragonneau.com/) and Zaikun Zhang.

See [Zaikun Zhang's talk](https://raw.githubusercontent.com/ztalks/20230825-iciam23/main/20230825-iciam.pdf)
on PRIMA at [The 10th International Congress on Industrial and Applied Mathematics](https://iciam2023.org/) for more information.


### Why

Professor Powell carefully implemented his derivative-free optimization methods into publicly available solvers,
which are genuine masterpieces. They are widely used by engineers and scientists. For instance,
see Section 1 of [a recent paper on Powell's solvers](https://arxiv.org/pdf/2302.13246.pdf)
as well as the Google searches of [COBYLA](https://www.google.com/search?q=cobyla)
and [BOBYQA](https://www.google.com/search?q=bobyqa).

However, Professor Powell's implementation was done in [Fortran 77](./fortran/original).
The code is nontrivial to understand or maintain, let alone extend.
For many practitioners, this has become an obstacle to exploiting these solvers in their
applications. Even worse, it has hindered researchers from exploring the wealth left by Professor
Powell. By all means, it is
[necessary to make the solvers available in languages other than Fortran](https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-23-23992)
promptly, first wrapping Powell's code, which is the objective of [PDFO](https://www.pdfo.net),
and then providing native and modernized implementations, which is the mission of PRIMA.

Before he passed, Professor Powell had asked me and
[Professor Nick Gould](https://www.numerical.rl.ac.uk/people/nimg) to maintain his solvers.
This is an honorable mission. To make the solvers more accessible, I started PRIMA.
It is a project similar to the translation, interpretation, and annotation of Euclid’s
*Elements*. It will make Powell's solvers easily understandable to everyone, not only the experts.
Few people remember [who translated *Elements*](https://en.wikipedia.org/wiki/Euclid%27s_Elements#Translations),
but it is a job that must be done.

PRIMA aims to provide the reference implementation of Powell's methods in modern languages,
including [**modern** Fortran](https://fortran-lang.org) (F2008 or newer), C/C++, Python, MATLAB,
Julia, and R. It will be a **faithful** implementation, in the sense that the code will be
mathematically equivalent to Powell’s, **except for** the
[bug fixes](#bug-fixes) and [improvements](#improvements) made intentionally.

The focus is to implement these methods in a **structured** and **modularized** way so that they
are **understandable**, **maintainable**, **extendable**, **fault-tolerant**, and **future-proof**.
The code will **have no GOTO** (of course)
and will **use matrix-vector procedures instead of loops** whenever possible.
In doing so, PRIMA codes the algorithms **in a way that we would present them on a blackboard**.
Such an implementation will enable us to get a deeper understanding of Powell's methods and
pave the way for new developments based on them.

There do exist "translations" of Powell's Fortran 77 code in other languages. For example,
[NLopt](https://github.com/stevengj/nlopt) contains a C version of COBYLA, NEWUOA, and BOBYQA,
but the C code in NLopt is translated from the Fortran 77 code straightforwardly, if
not automatically by [f2c](https://netlib.org/f2c/f2c.pdf), and hence inherits the style, structure,
and probably [bugs](#bug-fixes) of the original Fortran 77 implementation.
Note, however, that
[Py-BOBYQA](https://numericalalgorithmsgroup.github.io/pybobyqa/) is a **true translation** of BOBYQA
to Python, with significant improvements.


### How

The mission of PRIMA is nontrivial due to the delicacy of Powell's algorithms and the unique style
of his code. To ensure the faithfulness of PRIMA,
the **modern** Fortran version was started by refactoring Powell's code into the free form via a small
[MATLAB tool](./matlab/setup_tools/freeform.m).
However, such refactored code is far from what is desired, because it inherits completely
the structure and style of Powell's code except for the layout. Significant modifications are needed
to reorganize (indeed, to **rewrite**) the code. To maintain the faithfulness and quality of the
reference implementation, extensive tests are conducted after each and every tiny modification,
using the [CUTEst](https://github.com/ralna/CUTEst) problems via [MatCUTEst](https://github.com/matcutest/matcutest).
The tests do not only verify the faithfulness of the implementation but also check that **the solvers
behave properly even if they are invoked with improper inputs or [encounter failures of function
evaluations](https://github.com/libprima/prima/blob/main/matlab/tests/private/tough.m)**.
[**Stress tests**](#stress-tests) are also conducted
periodically to verify that the solvers work correctly without running into errors when applied to
**excessively large problems**.

The tests are **automated** by
[GitHub Actions](https://docs.github.com/en/actions). As of August 2023, more than
45,000 "workflows" have been successfully run by GitHub Actions. Normally, each workflow consists of \~ 5
([sometimes more than 200](https://github.com/primalib/prima/actions/runs/5763631681))
**randomized** tests,
each test taking from tens of minutes to several hours (the maximum
is 6 hours, after which the test will be canceled automatically). In other words,
PRIMA has been verified by more than 200,000 hours (or **more than 20 years**) of randomized tests.
**Code must be battle-tested before becoming software.**

Since each GitHub Team account can only run at most 60 GitHub Actions workflows concurrently, I have
to distribute this large amount of tests to multiple Team accounts as follows.

- [Tests](https://github.com/libprima/prima/actions) at [libprima/prima](https://github.com/libprima/prima)

    - [![CMake build](https://github.com/libprima/prima/actions/workflows/cmake.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake.yml)
    - [![CMake build, nagfor](https://github.com/libprima/prima/actions/workflows/cmake_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake_nagfor.yml)
    - [![CMake build, macOS ARM64](https://github.com/libprima/prima/actions/workflows/cmake_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake_mac.yml)
    - [![CMake build, macOS ARM64, nagfor](https://github.com/libprima/prima/actions/workflows/cmake_mac_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/cmake_mac_nagfor.yml)
    - [![Test nagfor, macOS ARM64](https://github.com/libprima/prima/actions/workflows/test_nagfor_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/test_nagfor_mac.yml)
    - [![Test matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/test_matlab_mac.yml)
    - [![Stress test, matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/stress_test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/stress_test_matlab_mac.yml)
    - [![Parallel test, matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/parallel_test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/parallel_test_matlab_mac.yml)
    - [![Recursive test, matlab, macOS ARM64](https://github.com/libprima/prima/actions/workflows/recursive_test_matlab_mac.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/recursive_test_matlab_mac.yml)
    - [![Test nagfor](https://github.com/libprima/prima/actions/workflows/test_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/test_nagfor.yml)
    - [![Build Python wheels](https://github.com/libprima/prima/actions/workflows/build_python.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/build_python.yml)
    - [![Compile MEX](https://github.com/libprima/prima/actions/workflows/compile_mex.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/compile_mex.yml)
    - [![Lint the Fortran code and the MEX gateways with nagfor](https://github.com/libprima/prima/actions/workflows/lint_nagfor.yml/badge.svg)](https://github.com/libprima/prima/actions/workflows/lint_nagfor.yml)

- <a name="verification"></a> [Tests](https://github.com/libsprima/prima/actions) at [libsprima/prima](https://github.com/libsprima/prima)
    - [![Verification, small](https://github.com/libsprima/prima/actions/workflows/verify_small.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/verify_small.yml)
    - [![Verification, big](https://github.com/libsprima/prima/actions/workflows/verify_big.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/verify_big.yml)
    - [![Verification, large](https://github.com/libsprima/prima/actions/workflows/verify_large.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/verify_large.yml)
    - [![Plot performance profiles for lincoa](https://github.com/libsprima/prima/actions/workflows/profile_lincoa_small.yml/badge.svg)](https://github.com/libsprima/prima/actions/workflows/profile_lincoa_small.yml)

- <a name="profiling"></a> [Tests](https://github.com/primapack/prima/actions) at [primapack/prima](https://github.com/primapack/prima)

    - [![Plot performance profiles for cobyla, small](https://github.com/primapack/prima/actions/workflows/profile_cobyla_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_cobyla_small.yml)
    - [![Plot performance profiles for uobyqa, small](https://github.com/primapack/prima/actions/workflows/profile_uobyqa_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_uobyqa_small.yml)
    - [![Plot performance profiles for newuoa, small](https://github.com/primapack/prima/actions/workflows/profile_newuoa_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_newuoa_small.yml)
    - [![Plot performance profiles for bobyqa, small](https://github.com/primapack/prima/actions/workflows/profile_bobyqa_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_bobyqa_small.yml)
    - [![Plot performance profiles for PRIMA, small](https://github.com/primapack/prima/actions/workflows/profile_prima_small.yml/badge.svg)](https://github.com/primapack/prima/actions/workflows/profile_prima_small.yml)

- [Tests](https://github.com/fortlab/prima/actions) at [fortlab/prima](https://github.com/fortlab/prima)

    - [![Plot performance profiles for all problems, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_all_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_all_sq.yml)
    - [![Plot performance profiles for cobyla, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_cobyla_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_cobyla_small_sq.yml)
    - [![Plot performance profiles for uobyqa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_uobyqa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_uobyqa_small_sq.yml)
    - [![Plot performance profiles for newuoa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_newuoa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_newuoa_small_sq.yml)
    - [![Plot performance profiles for bobyqa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_bobyqa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_bobyqa_small_sq.yml)
    - [![Plot performance profiles for lincoa, small, single and quadruple](https://github.com/fortlab/prima/actions/workflows/profile_lincoa_small_sq.yml/badge.svg)](https://github.com/fortlab/prima/actions/workflows/profile_lincoa_small_sq.yml)

- [Tests](https://github.com/primalib/prima/actions) at [primalib/prima](https://github.com/primalib/prima)

    - [![Verification, archiva](https://github.com/primalib/prima/actions/workflows/verify_archiva.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/verify_archiva.yml)
    - [![Plot performance profiles for all problems](https://github.com/primalib/prima/actions/workflows/profile_all.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_all.yml)
    - [![Plot performance profiles, single](https://github.com/primalib/prima/actions/workflows/profile_single.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_single.yml)
    - [![Plot performance profiles, quadruple](https://github.com/primalib/prima/actions/workflows/profile_quadruple.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_quadruple.yml)
    - [![Plot performance profiles, int16 and int64](https://github.com/primalib/prima/actions/workflows/profile_integer_kind.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_integer_kind.yml)
    - [![Plot performance profiles, infnan](https://github.com/primalib/prima/actions/workflows/profile_infnan.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_infnan.yml)
    - [![Plot performance profiles, various compiler options](https://github.com/primalib/prima/actions/workflows/profile_compiler_options.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_compiler_options.yml)
    - [![Plot performance profiles, various npt](https://github.com/primalib/prima/actions/workflows/profile_npt.yml/badge.svg)](https://github.com/primalib/prima/actions/workflows/profile_npt.yml)

- [Tests](https://github.com/opt4ai/prima/actions) at [opt4ai/prima](https://github.com/opt4ai/prima)

    - [![Test MATLAB](https://github.com/opt4ai/prima/actions/workflows/test_matlab.yml/badge.svg)](https://github.com/opt4ai/prima/actions/workflows/test_matlab.yml)
    - [![Test MATLAB, Linux](https://github.com/opt4ai/prima/actions/workflows/test_matlab_linux.yml/badge.svg)](https://github.com/opt4ai/prima/actions/workflows/test_matlab_linux.yml)

- [Tests](https://github.com/opt2ai/prima/actions) at [opt2ai/prima](https://github.com/opt2ai/prima)

    - [![Test MATLAB, macOS Intel](https://github.com/opt2ai/prima/actions/workflows/test_matlab_mac_intel.yml/badge.svg)](https://github.com/opt2ai/prima/actions/workflows/test_matlab_mac_intel.yml)
    - [![Test MATLAB, Windows](https://github.com/opt2ai/prima/actions/workflows/test_matlab_windows.yml/badge.svg)](https://github.com/opt2ai/prima/actions/workflows/test_matlab_windows.yml)

- [Tests](https://github.com/sprimalib/prima/actions) at [sprimalib/prima](https://github.com/sprimalib/prima)<a name="stress-tests"></a>

    - <a name="stress"></a>[![Stress test on large problems, MATLAB](https://github.com/sprimalib/prima/actions/workflows/stress_test_matlab.yml/badge.svg)](https://github.com/sprimalib/prima/actions/workflows/stress_test_matlab.yml)
    - [![Stress test on large problems, Fortran](https://github.com/sprimalib/prima/actions/workflows/stress_test_fortran.yml/badge.svg)](https://github.com/sprimalib/prima/actions/workflows/stress_test_fortran.yml)

- [Tests](https://github.com/zequipe/prima/actions) at [zequipe/prima](https://github.com/zequipe/prima)
    - [![Parallel test, MATLAB](https://github.com/zequipe/prima/actions/workflows/parallel_test_matlab.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/parallel_test_matlab.yml)
    - [![Recursive test, MATLAB](https://github.com/zequipe/prima/actions/workflows/recursive_test_matlab.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/recursive_test_matlab.yml)
    - [![Test Flang](https://github.com/zequipe/prima/actions/workflows/test_flang.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_flang.yml)
    - [![Test Flang in AMD AOCC](https://github.com/zequipe/prima/actions/workflows/test_aflang.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_aflang.yml)
    - [![Test nvfortran](https://github.com/zequipe/prima/actions/workflows/test_nvfortran.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_nvfortran.yml)
    - [![Test Oracle sunf95](https://github.com/zequipe/prima/actions/workflows/test_sunf95.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_sunf95.yml)
    - [![Test g95](https://github.com/zequipe/prima/actions/workflows/test_g95.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/test_g95.yml)
    - [![Test RESCUE and IDZ, classical](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_classical.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_classical.yml)
    - [![Test RESCUE and IDZ, modernized](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_modernized.yml/badge.svg)](https://github.com/zequipe/prima/actions/workflows/profile_rescue_idz_modernized.yml)

- [Tests](https://github.com/equipez/prima/actions) at [equipez/prima](https://github.com/equipez/prima)
    - [![Test gfortran on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_gfortran_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_gfortran_kunpeng.yml)
    - [![Test Flang on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_flang_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_flang_kunpeng.yml)
    - [![Test nvfortran on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_nvfortran_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_nvfortran_kunpeng.yml)
    - [![Test armflang on Kunpeng](https://github.com/equipez/prima/actions/workflows/test_armflang_kunpeng.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_armflang_kunpeng.yml)
    - [![CMake build on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/cmake_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/cmake_pi.yml)
    - [![Test gfortran on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_gfortran_pi64.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_gfortran_pi64.yml)
    - [![Test Flang on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_flang_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_flang_pi.yml)
    - [![Test armflang on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_armflang_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_armflang_pi.yml)
    - [![Test nvfortran on Raspberry Pi](https://github.com/equipez/prima/actions/workflows/test_nvfortran_pi.yml/badge.svg)](https://github.com/equipez/prima/actions/workflows/test_nvfortran_pi.yml)

- [Tests](https://github.com/s-prima/prima/actions) at [s-prima/prima](https://github.com/s-prima/prima)

    - [![Lint the Fortran code and the MEX gateways on GitHub hosted runners](https://github.com/s-prima/prima/actions/workflows/lint_hosted.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/lint_hosted.yml)
    - [![Test gfortran](https://github.com/s-prima/prima/actions/workflows/test_gfortran.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/test_gfortran.yml)
    - [![Test ifort](https://github.com/s-prima/prima/actions/workflows/test_ifort.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/test_ifort.yml)
    - [![Test ifx](https://github.com/s-prima/prima/actions/workflows/test_ifx.yml/badge.svg)](https://github.com/s-prima/prima/actions/workflows/test_ifx.yml)


### Current status

#### Modern Fortran

After almost **three** years of intensive coding, the [modern Fortran version](./fortran) of
PRIMA was finished by December 2022.
It can be compiled using CMake as follows.<a name="cmake"></a>
```bash
git clone --depth 1 https://github.com/libprima/prima.git
cd prima
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build --target install
```
This should create the `primaf` library for Fortran usage, located in the `install/lib/` directory
to be used with the module files in `install/include/prima/mod/`.
In case CMake fails to find your Fortran compiler,
you can indicate it by specifying `-DCMAKE_Fortran_COMPILER=/path/to/your/Fortran/compiler`.
Similarly, set `-DCMAKE_C_COMPILER=/path/to/your/C/compiler` for your C compiler if needed.

Examples on how to use the library from an external code are available in [`fortran/examples/`](https://github.com/libprima/prima/tree/main/fortran/examples).
Below is an illustration with COBYLA.
```bash
cd fortran/examples/cobyla
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install -DPRIMA_DIR=$PWD/../../../install/lib/cmake/prima/
cmake --build build --target install
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/../../../install/lib ./install/bin/cobyla_example_1
```

#### C

A C binding to the Fortran library is available in the [`c/` folder](https://github.com/libprima/prima/tree/main/c).
In the same way as the Fortran library, it can be [compiled using CMake](#cmake),
which should also create the `primac` library for C compilation, located in `install/lib/` to be used with the `prima.h` header in `install/include/prima/`.

Examples on how to use the library from an external code are available in [`c/examples/`](https://github.com/libprima/prima/tree/main/c/examples).
Below is an illustration with COBYLA.
```bash
cd c/examples/cobyla
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install -DPRIMA_DIR=$PWD/../../../install/lib/cmake/prima/
cmake --build build --target install
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/../../../install/lib ./install/bin/cobyla_example
```

#### Python

- An [interface](./python) is provided for [using the **modern** Fortran implementation in Python](./python/examples/rosenbrock.py).
- The inclusion of PRIMA into SciPy is [under discussion](https://github.com/scipy/scipy/issues/18118). It will replace the [buggy](#bug-fixes) and unmaintained Fortran 77 version of [COBYLA underlying `scipy.optimize.minimize`](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyla.html#optimize-minimize-cobyla), and make the other four solvers available to all SciPy users.
- A [native Python implementation of PRIMA](https://github.com/libprima/prima/pull/37) is under development.

#### MATLAB

- An [interface](./matlab/interfaces/prima.m) is provided for [using the **modern** Fortran implementation in MATLAB](./README_mat.md).
- <a name="newuoa_mat"></a>A [pure MATLAB version of NEWUOA](./matlab/interfaces/+newuoa_mat/) is implemented. It was
  generated straightforwardly (indeed, **automatically**) from an earlier version of the
  **modern** Fortran code (with the help of Mr. Galann Pennec). The other four solvers will be
  implemented in MATLAB similarly.

#### Julia
- A [Julia interface](https://juliahub.com/ui/Packages/General/PRIMA) is provided
by [`PRIMA.jl`](https://github.com/libprima/prima.jl).
It is registered in the General Registry of Julia as
[`PRIMA`](https://github.com/JuliaRegistries/General/tree/master/P/PRIMA).

#### Other languages

- Interfaces for using the modern Fortran implementation in other languages will be available later.
- Given the **modern** Fortran version, **native implementations** in other languages
become **much easier**, because we now have a structured and modularized implementation as a reference.
My team will implement the methods in other languages in this way.
This is the main motivation for developing the **modern** Fortran version first &mdash;
to provide a modernized reference implementation for the development in other languages.

### Bug fixes

PRIMA has fixed some **serious** issues in the **original Fortran 77 implementation** of Powell's methods.
Note that all of them are problems in the Fortran 77 code rather than flaws in the algorithms.

<!---[NLopt.jl](https://github.com/JuliaOpt/NLopt.jl), -->
The examples given below are bugs or requests sent to [SciPy](https://github.com/scipy/scipy),
[NLopt](https://github.com/stevengj/nlopt),
[nloptr](https://github.com/astamm/nloptr),
[OpenTURNS](https://github.com/openturns/openturns),
etc., which are reputable packages that wrap/interface the **original Fortran 77 implementation**
of Powell's solver. Inevitably, they suffer from the bugs in the Fortran 77 code.

- The Fortran 77 solvers may get **stuck** in infinite loops.

     - [optimize: COBYLA hangs / infinite loop #8998](https://github.com/scipy/scipy/issues/8998)
     - [BUG: Scipy.optimize / COBYLA hangs on some CPUs #15527](https://github.com/scipy/scipy/issues/15527)

	 - [COBYLA freezes (though maxeval and maxtime are given) #370](https://github.com/stevengj/nlopt/issues/370)

	 - [COBYLA hangs #118](https://github.com/stevengj/nlopt/issues/118)

	 - [NEWUOA_BOUND stuck in infinite loop inside MMA #117](https://github.com/stevengj/nlopt/issues/117)

     - [Cobyla freezes in 0T1.16rc1 #1651](https://github.com/openturns/openturns/issues/1651)

     - [Optimization freezes #25](https://github.com/astamm/nloptr/issues/25)

     - [BOBYQA gets stuck in infinite loop. #7](https://github.com/cureos/csnumerics/issues/7)

     - [Cobyla turns into infinite loop and never finishes #8](https://github.com/cureos/csnumerics/issues/8)

     - [Algorithm turns into infinite loop and never finishes #3](https://github.com/xypron/jcobyla/issues/3)

     - [The Fortran 77 version of UOBYQA encounters infinite cyclings very often if PRIMA_REAL_PRECISION is 32](https://github.com/libprima/prima/issues/98)

- The Fortran 77 solvers may **crash** with [segmentation faults](https://en.wikipedia.org/wiki/Segmentation_fault)
  due to uninitialized variables that are used as indices.

     - [Fix all uninitialized variable warnings #134](https://github.com/stevengj/nlopt/issues/134)

	 - [BOBYQA uninitialised variables in rare cases #133](https://github.com/stevengj/nlopt/issues/133)

	 - [Use of uninitialized variable in BOBYQA altmov #36](https://github.com/stevengj/nlopt/issues/36)

- Fortran 77 COBYLA may **not return the best point** that is evaluated; sometimes, the returned point can have a
large constraint violation even though the starting point is feasible.

	 - [nlopt COBYLA optimizer gives unexpected output #182](https://github.com/stevengj/nlopt/issues/182)

	 - [Last Result Returned Not Optimized Result #110](https://github.com/stevengj/nlopt/issues/110)

	 - [COBYLA returns last evaluated function which might not be minimum #57](https://github.com/stevengj/nlopt/issues/57)

     - [Successful termination when constraints violated #1](https://github.com/cureos/jcobyla/issues/1)

<!---
- Thread-safety
    - [scipy.optimize.minimize(method='COBYLA') not threadsafe #9658](https://github.com/scipy/scipy/issues/9658)

    - [BUG: Make cobyla threadsafe #3](https://github.com/sturlamolden/scipy/pull/3)
-->


### Improvements

Thanks to the improvements introduced into the new implementation, PRIMA outperforms Powell's
original code in terms of the **number of function evaluations**, which is the standard performance
indicator in derivative-free optimization.
Below are the [performance profiles](https://arxiv.org/pdf/cs/0102001.pdf)
of the PRIMA solvers compared with Powell's implementation, the convergence tolerance being $\tau = 10^{-6}$.
Roughly speaking, performance profiles plot the percentage of test problems solved against the budget,
which is measured relative to the cost of the most efficient solver in the comparison.
A **higher** curve indicates a **better** solver.
See [Benchmarking Derivative-Free Optimization Algorithms](https://www.mcs.anl.gov/~wild/dfo/benchmarking)
([J. J. Moré](https://www.anl.gov/profile/jorge-j-more) and [S. M. Wild](https://www.anl.gov/profile/stefan-m-wild))
for more information.


- NEWUOA on unconstrained CUTEst problems of at most 200 variables
<img src="./benchmark/latest/prima_newuoa.png" style="width:26em;"/>

- BOBYQA on bound-constrained CUTEst problems of at most 200 variables
<img src="./benchmark/latest/prima_bobyqa.png" style="width:26em;"/>

- LINCOA on linearly constrained CUTEst problems of at most 200 variables and 20000 constraints
<img src="./benchmark/latest/prima_lincoa.png" style="width:26em;"/>

- COBYLA on nonlinearly constrained CUTEst problems of at most 100 variables and 10000 constraints
<img src="./benchmark/latest/prima_cobyla.png" style="width:26em;"/>

- UOBYQA on unconstrained CUTEst problems of at most 100 variables
<img src="./benchmark/latest/prima_uobyqa.png" style="width:26em;"/>


### Who was Powell?

[Michael James David Powell FRS](https://en.wikipedia.org/wiki/Michael_J._D._Powell) was
["a British numerical analyst who was among the pioneers of computational mathematics"](https://royalsocietypublishing.org/doi/full/10.1098/rsbm.2017.0023).
He was the inventor/early contributor of
[quasi-Newton method](https://en.wikipedia.org/wiki/Quasi-Newton_method),
[trust region method](https://en.wikipedia.org/wiki/Trust_region),
[augmented Lagrangian method](https://en.wikipedia.org/wiki/Augmented_Lagrangian_method),
and [SQP method](https://en.wikipedia.org/wiki/Sequential_quadratic_programming).
Each of them is a pillar of modern numerical optimization. He also made significant contributions
to [approximation theory and methods](https://www.cambridge.org/highereducation/books/approximation-theory-and-methods/66FD8CD6F18FE1ED499A8CA9A05F2A5A#overview).

Among numerous honors, Powell was one of the two recipients of the first
[Dantzig Prize](https://en.wikipedia.org/wiki/Dantzig_Prize)
from the Mathematical Programming Society (MOS) and Society for Industrial and Applied Mathematics (SIAM).
This is considered the highest award in optimization.


### A "fun" fact

In the past years, while working on PRIMA, I have spotted a dozen of [bugs in reputable Fortran compilers](https://github.com/zequipe/test_compiler)
and three [bugs in MATLAB](https://github.com/zequipe/test_matlab). Each of them represents days of **bitter** debugging, which finally led to the conclusion
that it was not a problem in my code but a flaw in the Fortran compilers or in MATLAB. From a very unusual angle, this reflects how intensive
the coding has been.

The bitterness behind this "fun" fact is exactly why I work on PRIMA: I hope that all
the frustrations that I have experienced will not happen to any user of Powell's methods anymore.
I hope I am the last one in the world to decode a maze of 244 GOTOs in 7939 lines of Fortran 77 code &mdash;
I did this for three years and I do not want anyone else to do it again.


### Acknowledgment

PRIMA is dedicated to the memory of the late [Professor Powell](https://www.zhangzk.net/powell.html) with gratitude for his inspiration and
for the wealth he left to us.

I am grateful to [Professor Ya-xiang Yuan](http://lsec.cc.ac.cn/~yyx) for his everlasting encouragement and support.

The development of PRIMA would have been a mission impossible without the groundwork laid by the [PDFO](https://www.pdfo.net)
package of [Tom M. Ragonneau](https://tomragonneau.com/) and Zaikun Zhang. PDFO is Chapter 3 of
Ragonneau's [thesis](https://tomragonneau.com/documents/thesis.pdf) co-supervised by Zaikun Zhang and Professor [Xiaojun Chen](https://www.polyu.edu.hk/ama/staff/xjchen/ChenXJ.htm), with financial support from
the [Hong Kong Ph.D. Fellowship Scheme](https://cerg1.ugc.edu.hk/hkpfs/index.html) (ref. PF18-24698).

PRIMA is a long-term project, which would not have been sustainable without the continued funds from the
[Hong Kong Research Grants Council](https://www.ugc.edu.hk/eng/rgc/) (ref. PolyU 253012/17P, PolyU 153054/20P,
PolyU 153066/21P, and PolyU 153086/23P) and [The Hong Kong Polytechnic University](https://www.polyu.edu.hk/) (PolyU),
in particular the [Department of Applied Mathematics](https://www.polyu.edu.hk/ama) (AMA).


### Citing PRIMA

PRIMA has taken me significant energy and time. I will be delighted if it is useful to you. All I need is a citation / acknowledgment.
Note that PRIMA contains [bug fixes](#bug-fixes) and [improvements](#improvements) that do not exist in Powell's original
implementation of the solvers. Results produced by PRIMA may not be reproducible using the original solvers.

If you use PRIMA, please cite it as follows. The citation will be pointed to my paper on PRIMA when I finish it.

[1] Z. Zhang, PRIMA: Reference Implementation for Powell's Methods with Modernization and Amelioration, available at http://www.libprima.net, [DOI: 10.5281/zenodo.8052654](https://doi.org/10.5281/zenodo.8052654), 2023

```bibtex
@misc{Zhang_2023,
    title        = {{PRIMA: Reference Implementation for Powell's Methods with Modernization and Amelioration}},
    author       = {Zhang, Z.},
    howpublished = {available at http://www.libprima.net, DOI: 10.5281/zenodo.8052654},
    year         = {2023}
}
```

In addition, Powell’s methods can be cited as follows.

[2] M. J. D. Powell, A direct search optimization method that models the
objective and constraint functions by linear interpolation, In Advances
in *Optimization and Numerical Analysis*, *eds.* S. Gomez and J. P. Hennart,
pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994

[3] M. J. D. Powell, UOBYQA: unconstrained optimization by quadratic
approximation, *Math. Program.*, 92(B):555--582, 2002

[4] M. J. D. Powell, The NEWUOA software for unconstrained optimization
without derivatives, In *Large-Scale Nonlinear Optimization*, *eds.* G. Di Pillo
and M. Roma, pages 255--297, Springer, New York, US, 2006

[5] M. J. D. Powell, The BOBYQA algorithm for bound constrained
optimization without derivatives, Technical Report DAMTP 2009/NA06,
Department of Applied Mathematics and Theoretical Physics, Cambridge
University, Cambridge, UK, 2009

[6] T. M. Ragonneau and Z. Zhang,
[PDFO: a cross-platform package for Powell's derivative-free optimization solvers](https://arxiv.org/pdf/2302.13246.pdf),
arXiv:2302.13246, 2023


**Remarks**

- LINCOA seeks the least value of a nonlinear function subject to
linear inequality constraints without using derivatives of the objective
function. Powell did not publish a paper to introduce the algorithm.

- [The paper [6]](https://arxiv.org/pdf/2302.13246.pdf) introduces [the PDFO package](https://www.pdfo.net)
rather than PRIMA. Nevertheless, it provides a good introduction to Powell's methods.


### Charityware

PRIMA is [charityware](https://en.wikipedia.org/wiki/Careware), distributed for free under its
[license](https://github.com/libprima/prima/blob/main/LICENCE.txt).
If you appreciate it, you may consider making a donation to a charity that you trust
(in addition to [citing \& acknowledging PRIMA](https://github.com/libprima/prima#citing-prima)).
This is only a suggestion, not an obligation.

The inspiration comes from [Vim](https://www.vim.org/), with which Zaikun Zhang typed all his PRIMA code.


### Contact

In case of problems, [open a GitHub issue](https://github.com/libprima/prima/issues) or [contact
Zaikun Zhang](https://www.zhangzk.net).


### Mirrors

- Gitee: [https://gitee.com/libprima/prima](https://gitee.com/libprima/prima)

- GitHub: [https://github.com/libprima/prima](https://github.com/libprima/prima)

- GitLab: [https://gitlab.com/libprima/prima](https://gitlab.com/libprima/prima)


### <a href="https://star-history.com/#libprima/prima&Date">Star history</a>

[stardev](https://stardev.io/) ranking: [28 among 30,975](https://stardev.io/top/repos/fortran?developer=libprima&repo=prima) Fortran repos as of Jan. 2024.

<img src="https://api.star-history.com/svg?repos=libprima/prima&type=Date">


<p align="center"><strong>Thank you for your support.</strong></p>
