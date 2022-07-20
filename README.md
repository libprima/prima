This is the reference implementation of Powell's derivative-free optimization solvers,
namely COBYLA, UOBYQA, NEWUOA, BOBYQA, and LINCOA.

This package is part of a research project funded by the
[Hong Kong Research Grants Council](https://www.ugc.edu.hk/eng/rgc) and
the [Hong Kong Polytechnic University](https://www.polyu.edu.hk) (PolyU).
It is still under intensive development, and there is no release yet. If you want to use the
above-mentioned solvers, see the [website](https://www.pdfo.net)
and [repository](https://github.com/pdfo/pdfo) of PDFO instead.

The goal is to implement these solvers in modern languages --- first [**modern** Fortran](https://fortran-lang.org)
(F2003 or later), and then MATLAB, Python, and probably Julia and R. It will be a faithful implementation, in the
sense that the new code will be mathematically equivalent to Powell’s, except for the bug fixes and
improvements that we make intentionally.

The focus is to implement Powell’s solvers in a modularized and structured way so that they are
**readable**, **maintainable**, and **extendable**. The new code will have no GOTO (of course) and will use
matrix-vector procedures instead of loops whenever possible.

This is not a trivial mission due to the delicacy of Powell's algorithms and the unique style of his code.
We started The Fortran code by refactoring Powell's code into the free form via a small
[MATLAB tool](https://github.com/zequipe/pdfo_ref/blob/master/matlab/setup_tools/freeform.m) written
by ourselves. However, such refactored code is far from what we want, because it inherits
completely the structure and style of Powell's code except for the format. Extensive modifications
are needed to reorganize (indeed, to rewrite) the code.
To maintain the faithfulness and quality of our implementation, intensive tests are conducted
after every tiny modification, the test problems coming from the [CUTEst set](https://github.com/ralna/CUTEst).
The tests are automated with the help of
[GitHub Actions](https://en.wikipedia.org/wiki/Explorative_strategies). As of July 2022, more than 20,000
"workflows" have been run by GitHub Actions
(see https://github.com/zequipe/gitpersonal/actions and https://github.com/zequipe/pdfo_ref/actions).
Normally, each workflow consists of \~ 5 **randomized** tests
that are conducted in parallel, each test taking from tens of minutes to several hours (the maximum is
6 hours, after which the workflow will be canceled automatically). In other words, our
implementation has been tested by  \~ $10^5$ hours (or \~ $10$ years) of randomized computation.


Dedicated to late Professor [M. J. D. Powell](https://www.zhangzk.net/powell.html) FRS (1936--2015).
