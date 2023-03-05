This is the reference implementation of Powell's UOBYQA algorithm described in

M. J. D. Powell, UOBYQA: unconstrained optimization by quadratic approximation, Math. Program.,
92(B):555--582, 2002

Coded by Zaikun ZHANG (www.zhangzk.net) based the above-mentioned paper and Powell's code, with
modernization, bug fixes, and improvements.

Due to the dependency among the modules, the files must be compiled in the order indicated in
ffiles.txt. Before that, one must compile the Fotran files in the directory "common" according
to the ffiles.txt therein.

See the directory "examples" for an illustration about how to call the solver.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

Started in February 2022.
