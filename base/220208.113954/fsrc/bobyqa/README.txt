This is a modernized and improved Fortran implementation of the BOBYQA algorithm described in

M. J. D. Powell, The BOBYQA algorithm for bound constrained optimization without derivatives,
Technical Report DAMTP 2009/NA06, Department of Applied Mathematics and Theoretical Physics,
Cambridge University, Cambridge, UK, 2009

Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper.

Due to the dependency among the modules, the files must be compiled in the order indicated in
ffiles.txt. Before that, one must compile the Fotran files in the directory "common" according
to the ffiles.txt therein.

See the directory "examples" for an illustration about how to call the solver.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

Started in February 2022.
