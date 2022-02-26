***********************************************************************
Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
            and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
            Department of Applied Mathematics,
            The Hong Kong Polytechnic University

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

We look forward to your feedback! Thank you very much!

August 2021, Hong Kong
***********************************************************************

This is the README file for the Python version of PDFO.
See https://www.pdfo.net for more information.


1. Installation

1.1. Recommended installation
To use the Python version of PDFO on Linux, Mac, or Windows, you need Python
(version 3.7 or above).

It is highly recommended to install PDFO via PyPI.

Install pip in your system ( https://pip.pypa.io/en/stable/installing/ ) if
you Python version does not incude it. Then execute

python -m pip install pdfo

in a command shell (e.g., the terminal for Linux and macOS, or the Command
Shell for Windows). If your Python launcher is not python, adapt the command
accordingly. If this command runs successfully, PDFO is installed. You may
verify the installation by

python -m unittest pdfo.testpdfo

If you are an Anaconda user, PDFO is also available through the conda installer
( https://anaconda.org/conda-forge/pdfo ). However, it is not managed by us.

1.2. Alternative installation (using source distribution)
Alternatively, although deeply discouraged, PDFO can be installed from the
source code. It requires you to install additional Python headers, a Fortran
compiler (e.g., gfortran), and F2PY (provided by NumPy). Download and
decompress the source code package; you will obtain a folder containing
setup.py; in a command shell, change your directory to this folder; then
install PDFO by executing

python -m pip install .


2. Usage

2.1. PDFO provides the following Python functions:
pdfo, uobyqa, newuoa, bobyqa, lincoa, cobyla.

2.2. The "pdfo" function can automatically identify the type of your problem
and call one of Powell's solvers. The other five functions call the solver
indicated by their names. It is highly recommended to use "pdfo" instead of
"uobyqa", "newuoa", etc.

2.3. The "pdfo" function is designed to be compatible with the "minimize"
function available in scipy.optimize. You can call "pdfo" in exactly the same
way as calling "minimize", without the derivative arguments (PDFO does not use
derivatives).

2.4. For thedetailed syntax of these functions, use the standard "help" command
of Python. For example,

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
