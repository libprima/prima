***********************************************************************
Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
            and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
            Department of Applied Mathematics,
            The Hong Kong Polytechnic University

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

We look forward to your feedback! Thank you very much!

June 2020, Hong Kong
***********************************************************************

This is the README file for the MATLAB version of PDFO.
See https://www.pdfo.net for more information.


0. Prerequisites

PDFO supports MATLAB R2014a and later releases. To use PDFO, you need first
configure the MEX of your MATLAB so that it can compile Fortran. 

0.1. To see whether your MEX is ready, run the following code in MATLAB:

mex('-setup', '-v', 'FORTRAN'); mex('-v', fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.F'));

If this completes successfully, then your MEX is ready. Otherwise, it is not.

0.2. To configure MEX for compiling Fortran, see
https://www.mathworks.com/help/matlab/ref/mex.html .
It will require you to install a supported Fortran compiler on your system.
See https://www.mathworks.com/support/requirements/previous-releases.html .
Note that MathWorks (rather than PDFO) is quite rigid concerning the version
of your compiler, which has to be compatible with the release of your MATLAB;
the latest compiler is NOT necessarily supported by your MATLAB. On Windows, 
in addition to the Fortran compiler, MathWorks needs you to install the 
Microsoft Visual Studio and the Microsoft Windows SDK. Follow the official 
documentation of MathWorks closely. 


1. Installation

Download and decompress the source code package of PDFO. You will obtain
a folder containing setup.m. Place this folder at the location where you
want PDFO to be installed. In MATLAB, change the directory to this folder,
and execute the following command:

setup

If this command runs successfully, PDFO is installed. You may execute the
following command in MATLAB to verify the installation:

testpdfo


2. Usage 

2.1. PDFO provides the following MATLAB functions:
pdfo, uobyqa, newuoa, bobyqa, lincoa, cobyla.

2.2. The "pdfo" function can automatically identify the type of your problem
and then call one of Powell's solvers. The other five functions call the solver
indicated by their names. It is highly recommended to use "pdfo" instead of
"uobyqa", "newuoa", etc. 

2.3. The "pdfo" function is designed to be compatible with the "fmincon"
function available in the Optimization Toolbox of MATLAB. You can call "pdfo"
in exactly the same way as calling "fmincon". In addition, "pdfo" can be
called in some flexible ways that are not supported by "fmincon".

2.4. For detailed syntax of these functions, use the standard "help" command
of MATLAB. For example,

help pdfo

will tell you how to use "pdfo".


3. Uninstall

PDFO can be uninstalled using the setup.m script by executing the following
command in MATLAB:

setup unistall


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
