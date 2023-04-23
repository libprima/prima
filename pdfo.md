# PDFO for MATLAB: Powell's Derivative-Free Optimization solvers

Dedicated to the late Professor [M. J. D. Powell](https://www.zhangzk.net/powell.html)
FRS (1936&ndash;2015).

We look forward to your feedback! Thank you very much!

This is the README for the MATLAB version of PDFO.
See https://www.pdfo.net for more information.

## Prerequisites

PDFO supports MATLAB R2014a and later releases. To use PDFO, you need first
configure the [MEX](https://www.mathworks.com/help/matlab/ref/mex.html) of your
MATLAB so that it can compile Fortran.

To see whether your MEX is ready, run the following code in MATLAB:

```matlab
mex('-setup', '-v', 'FORTRAN'); mex('-v', fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.F'));
```

If this completes successfully, then your MEX is ready. Otherwise, it is not, and
you may try the [`setup_mex` package](https://github.com/equipez/setup_mex) at
```bash
https://github.com/equipez/setup_mex
```
It will help you to set MEX up on Windows or macOS (the setup of MEX is trivial on Linux).
In case `setup_mex` does not work, you need to consult a local MATLAB expert or the technical support of
MathWorks about "[how to set up MEX](https://www.mathworks.com/help/matlab/ref/mex.html)", which is
**not** part of PDFO.


## Installation

Download and decompress the [source code package](https://www.pdfo.net/docs.html#download).
You will obtain a folder containing `setup.m`. Place this folder at the location
where you  want PDFO to be installed. In MATLAB, change the directory to this
folder, and execute the following command:

```matlab
setup
```

If this command runs successfully, PDFO is installed. You may execute the
following command in MATLAB to verify the installation:

```matlab
testpdfo
```


## Usage

PDFO provides the following MATLAB functions:
`pdfo`, `uobyqa`, `newuoa`, `bobyqa`, `lincoa`, `cobyla`.

The `pdfo` function can automatically identify the type of your problem
and then call one of [Powell's](https://www.zhangzk.net/powell.html) solvers.
The other five functions call the solver indicated by their names. It is highly
recommended using `pdfo` instead of `uobyqa`, `newuoa`, etc.

The `pdfo` function is designed to be compatible with the `fmincon`
function available in the [Optimization Toolbox](https://www.mathworks.com/products/optimization.html)
of MATLAB. You can call `pdfo` in exactly the same way as calling `fmincon`. In
addition, `pdfo` can be  called in some flexible ways that are not supported by
`fmincon`.

For detailed syntax of these functions, use the standard `help` command
of MATLAB. For example,

```matlab
help pdfo
```

will tell you how to use `pdfo`.

## Uninstall

PDFO can be uninstalled using the setup.m script by executing the following
command in MATLAB:

```matlab
setup uninstall
```

## References

[1] M. J. D. Powell, A direct search optimization method that models the
objective and constraint functions by linear interpolation, In Advances
in Optimization and Numerical Analysis, eds. S. Gomez and J. P. Hennart,
pages 51&ndash;67, Springer Verlag, Dordrecht, Netherlands, 1994

[2] M. J. D. Powell, UOBYQA: unconstrained optimization by quadratic
approximation, Math. Program., 92(B):555&ndash;582, 2002

[3] M. J. D. Powell, The NEWUOA software for unconstrained optimization
without derivatives, In Large-Scale Nonlinear Optimization, eds. G. Di Pillo
and M. Roma, pages 255&ndash;297, Springer, New York, US, 2006

[4] M. J. D. Powell, The BOBYQA algorithm for bound constrained
optimization without derivatives, Technical Report DAMTP 2009/NA06,
Department of Applied Mathematics and Theoretical Physics, Cambridge
University, Cambridge, UK, 2009

*Remark:* LINCOA seeks the least value of a nonlinear function subject to
linear inequality constraints without using derivatives of the objective
function. Powell did not publish a paper to introduce the algorithm.
