***********************************************************************
Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
            and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
            Department of Applied Mathematics,
            The Hong Kong Polytechnic University

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

We look forward to your feedback! Thank you very much!

March 2020, Hong Kong
***********************************************************************

This is the README file for the Python version of PDFO on Linux or Mac.
See https://www.pdfo.net for more information.


0. Prerequisites

To use the Python version of PDFO, you need Python, NumPy, F2PY, and 
gfortran, which can be installed in the following way.

0.1. Install Python (version 2.7 or later) according to https://www.python.org .

0.2. Install NumPy (version 1.10.0 or later) for your Python. NumPy provides 
F2PY. We recommend to install the latest version of SciPy. Then NumPy will be 
installed by default. See https://www.scipy.org/install.html .

0.3. Install gfortran using your package manager, e.g., apt on Debian/Ubuntu,
yum on Fedora/RHEL/CentOS, and Homebrew on Mac.


1. Installation

PDFO can be installed by the setup.py script in the following way.

1.1. Decompress the source code package of PDFO if you have not done so. You
will obtain a folder containing setup.py. Place this folder at the location
where you want PDFO to be installed.

1.2. In a command shell, change your directory to the above-mentioned folder,
and execute the following command:

python setup.py

If this command runs successfully, PDFO is installed.


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
