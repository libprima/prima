This is the reference implementation of Powell's COBYLA algorithm described in

M. J. D. Powell, A direct search optimization method that models the objective and constraint
functions by linear interpolation, In Advances in Optimization and Numerical Analysis, eds. S. Gomez
and J. P. Hennart, pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994

Coded by Zaikun ZHANG (www.zhangzk.net) based on the above-mentioned paper and Powell's code, with
modernization and improvements.

Due to the dependency among the modules, the files must be compiled in the order indicated in
ffiles.txt. Before that, one must compile the Fotran files in the directory "common" according
to the ffiles.txt therein.

See the directory "examples" for an illustration about how to call the solver.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

Started in July 2021.
