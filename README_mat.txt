This is the README file for using PRIMA under MATLAB.

0. Prerequisites

PRIMA supports MATLAB R2018a and later releases. To use PRIMA, you need first
set up the MEX of your MATLAB so that it can compile Fortran (N.B.: the setup
of MEX is a pure MATLAB usage problem and it has nothing to do with PRIMA).

0.1. To see whether your MEX is ready, run the following code in MATLAB:

mex('-setup', '-v', 'FORTRAN'); mex('-v', fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.F'));

If this completes successfully, then your MEX is ready. Otherwise, it is not.

0.2. To configure MEX for compiling Fortran, you may refer to

https://github.com/equipez/setup_mex

which is a package providing scripts that attempt to facilitate setting up MEX.
In case it does not work, then please check the official documentation of MEX at
https://www.mathworks.com/help/matlab/ref/mex.html .
It will require you to install a supported Fortran compiler on your system.
See https://www.mathworks.com/support/requirements/previous-releases.html .
Note that MathWorks (rather than PRIMA) is quite rigid concerning the version
of your compiler, which has to be compatible with the release of your MATLAB;
the latest compiler is NOT necessarily supported by your MATLAB. On Windows,
in addition to the Fortran compiler, MathWorks needs you to install the
Microsoft Visual Studio with the "Desktop development with C++" workload and the
Microsoft Windows SDK. Follow the official documentation of MathWorks closely.


1. Installation

Download and decompress the source code package of PRIMA. You will obtain
a folder containing setup.m. Place this folder at the location where you
want PRIMA to be installed. In MATLAB, change the directory to this folder,
and execute the following command:

setup

If this command runs successfully, PRIMA is installed. You may execute the
following command in MATLAB to verify the installation:

testprima


2. Usage

2.1. PRIMA provides a MATLAB function `prima`, which can solve general
constrained or unconstrained optimization problems without using derivatives.

2.2. The `prima` function can automatically identify the type of your problem
and then call one of Powell's solvers, namely COBYLA, UOBYQA, NEWUOA, BOBYQA,
and LINCOA. The user can also specify the solver by setting the `solver` field
of the options passed to `prima`.

2.3. The `prima` function is designed to be compatible with the `fmincon`
function available in the Optimization Toolbox of MATLAB. You can call `prima`
in exactly the same way as calling `fmincon`.

2.4. For detailed syntax of `prima`, use the standard `help` command of MATLAB:

help prima

will tell you how to use `prima`.


3. Uninstall

PRIMA can be uninstalled using the setup.m script by executing the following
command in MATLAB:

setup uninstall


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
