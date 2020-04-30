***********************************************************************
Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
            and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
            Department of Applied Mathematics,
            The Hong Kong Polytechnic University

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

We look forward to your feedback! Thank you very much!

March 2020, Hong Kong
***********************************************************************

This is the README file for the Python version of PDFO on Windows.
See https://www.pdfo.net for more information.

For the moment, F2PY does not work well with Python 2 on Windows. If you want
to use the package on Windows, make sure to do everything in any version of
Python 3.


0. Prerequisites

To use the Python version of PDFO, you need Python, NumPy, F2PY, C++ Build
Tools and ifort, which can be installed in the following way.

0.1. Install Python (version 3.0 or later) according to https://www.python.org .

0.2. We recommend to install the latest version of SciPy. See
https://www.scipy.org/install.html .

0.3. Install C++ Build Tools according to
https://visualstudio.microsoft.com/visual-cpp-build-tools/ . Your installation
has to include the dependencies "C++ Build Tools".

0.4. Install Intel Fortran compiler according to
https://software.intel.com/en-us/fortran-compilers . Your installation has to
include Python librairies (true by default).


1. Installation

1.1. Recommanded installation via PyPI

PDFO can be installed in the following way. In an Intel shell environment
(either 32 of 64bits, depending on your Python version), execute the following
command:

python -m pip install pdfo

If this command runs successfully, PDFO is installed. You can now test the
package by executing the following command in any shell:

python -m unittest pdfo.testpdfo

1.2. Manual installation (only if 1.1. has not been executing)

Alternatively, you can download the source files at https://www.pdfo.net .

1.2.1. Decompress the source code package of PDFO if you have not done so. You
will obtain a folder containing setup.py. Place this folder at the location
where you want PDFO to be installed.

1.2.2. In an Intel shell environment (either 32 of 64bits, depending on your
Python version), change your directory to the above-mentioned folder, and
execute the following command:

python -m pip install .

If this command runs successfully, PDFO is installed. You can now test the
package by executing the following command in any shell:

python -m unittest pdfo.testpdfo


2. Usage

2.1. PDFO provides the following Python functions:
pdfo, uobyqa, newuoa, bobyqa, lincoa, cobyla.

2.2. The "pdfo" function can automatically identify the type of your problem
and the call one of Powell's solvers. The other five functions call the solver
indicated by their names. It is highly recommended to use "pdfo" instead of
"uobyqa", "newuoa", etc.

2.3. The "pdfo" function is designed to be compatible with the "minimize"
function available in scipy.optimize. You can call "pdfo" in exactly the same
way as calling "minimize", without the derivative arguments (PDFO does not use
derivatives).

2.4. For detailed syntax of these functions, use the standard "help" command of
Python. For example,

>>> from pdfo import pdfo
>>> help(pdfo)

will tell you how to use "pdfo".


3. References

[1] M. J. D. Powell, A direct search optimization method that models the
objective and constraint functions by linear interpolation, In Advances
in Optimization and Numerical Analysis, eds. S. Gomez and J. P. Hennart,
pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994

[2] M. J. D. Powell, UOBYQA: unconstrained optimization by quadratic
approximation, Math. Program., 92(B):555--582, 2002

[3] M. J. D. Powell, The NEWUOA software for unconstrained optimization
without derivatives, In Large-Scale Nonlinear Optimization, eds. G. Di Pillo
and M. Roma, pages 255--297, Springer, New York, US, 2006

[4] M. J. D. Powell, The BOBYQA algorithm for bound constrained
optimization without derivatives, Technical Report DAMTP 2009/NA06,
Department of Applied Mathematics and Theoretical Physics, Cambridge
University, Cambridge, UK, 2009

Remark: LINCOA seeks the least value of a nonlinear function subject to
linear inequality constraints without using derivatives of the objective
function. Powell did not publish a paper to introduce the algorithm.