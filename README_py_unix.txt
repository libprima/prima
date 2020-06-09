***********************************************************************
Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
            and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
            Department of Applied Mathematics,
            The Hong Kong Polytechnic University

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

We look forward to your feedback! Thank you very much!

June 2020, Hong Kong
***********************************************************************

This is the README file for the Python version of PDFO on Linux or Mac.
See https://www.pdfo.net for more information.


0. Prerequisites

To use the Python version of PDFO on Linux or Mac, you need Python (including
the headers), NumPy, F2PY, and gfortran, which can be installed in the following way.

0.1. Install Python (version 2.7 or above) according to https://www.python.org .
To include the headers on Linux, python3-dev or python-dev should be installed
depending on the Python you want to use; on Mac, the headers are included by
default if you install Python using Homebrew.

0.2. Install the latest version of SciPy. See https://www.scipy.org/install.html .
SciPy includes NumPy, which provides F2PY.

0.3. Install gfortran using your package manager, e.g., apt on Debian/Ubuntu,
yum on Fedora/RHEL/CentOS, and Homebrew on Mac.

0.4. On Mac, it is necessary to install Xcode. See https://developer.apple.com/xcode/ .
Xcode provides the libraries needed by F2PY. 


1. Installation

Note: You may need to replace 'python' with 'python3' in all the commands 
below if you use Python 3. 

1.1. Method 1: Installation via PyPI (recommended)

Install pip in your system ( https://pip.pypa.io/en/stable/installing ). 
Then install PDFO by the following command in a command shell (e.g., bash, zsh):

python -m pip install pdfo

If this command runs successfully, PDFO is installed. You can then test the 
installation by executing

python -m unittest pdfo.testpdfo

1.2. Method 2: Manual installation 

Download the source code package from https://www.pdfo.net . Decompressing 
the package, you will obtain a folder. In a command shell (e.g., bash, zsh),
change your directory to this folder, and then install PDFO by running

python -m pip install ./

which needs pip ( https://pip.pypa.io/en/stable/installing ); alternatively,
you can install PDFO using setup.py in the above-mentioned folder by running

python setup.py install --user

The installation can be tested by

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
