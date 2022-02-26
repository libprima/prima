***********************************************************************
Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
            and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
            Department of Applied Mathematics,
            The Hong Kong Polytechnic University

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

We look forward to your feedback! Thank you very much!

June 2020, Hong Kong
***********************************************************************

This is the README file for the Python version of PDFO on Windows.
See https://www.pdfo.net for more information.


0. Prerequisites

To use the Python version of PDFO on Windows, you need Python 3 (only 3.7 and
earlier, and PDFO does not support Python 2 on Windows), NumPy, F2PY, Microsoft
Visual Studio, and Intel Fortran compiler (ifort), which can be installed in
the following way.

0.1. Install Python (either version 3.0 and onwards, 3.7 and earlier) according
to https://www.python.org . Note that the Python installed from the Microsoft
Store may not work. In addition, the version of Python has to be compatible with
the Intel Fortran compiler; otherwise, Python may complain of missing dll.

0.2. Install the latest version of SciPy. See https://www.scipy.org/install.html .
Then NumPy will be installed by default. NumPy provides F2PY.

0.3. Install Microsoft Visual Studio according to 
https://visualstudio.microsoft.com . Make sure to include "C++ Build Tools"
in the Workloads pane, and "Microsoft Windows SDK" in the Individual Components
pane.

0.4. Install the Intel Fortran compiler (ifort) according to
https://software.intel.com/en-us/fortran-compilers .
Your installation has to include the Intel Distribution for Python (included
by default). Note that ifort may not support the latest version of Python.
As of June 15, 2020, the most recent release of ifort supports Python 3.7
according to the Key Specifications of the Intel Distribution for Python.


1. Installation

PDFO can be installed via PyPI. 

Install pip in your system ( https://pip.pypa.io/en/stable/installing/).
In the Intel command shell (either 32 or 64-bit, depending on your Python
version), execute

python -m pip install pdfo

Make sure to run the above command in the correct Intel Command shell,
as it is important to guarantee the compatibility between the compiler
and your Python. If your Python launcher is not python, adapt the
command accordingly. If this command runs successfully, PDFO is
installed. You may verify the installation by

python -m unittest pdfo.testpdfo

Note that an Ubuntu terminal is available in Windows 10:
https://ubuntu.com/tutorials/tutorial-ubuntu-on-windows
Within such a terminal, you can use the Linux version of PDFO following
README_py_unix.txt. The prerequisites are much easier to fulfill and no
commercial software such as ifort is needed.


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

3. Uninstall

PDFO can be uninstalled by executing the following command in a command shell:

python -m pip uninstall pdfo


4. References

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
