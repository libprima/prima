This is a modern Fortran implementation of the COBYLA algorithm described in

M. J. D. Powell, A direct search optimization method that models the objective and constraint
functions by linear interpolation, In Advances in Optimization and Numerical Analysis, eds. S. Gomez
and J. P. Hennart, pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994

Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper.

Started in July 2021.

Due to the dependency among the modules, the files must be compiled in the order indicated in
ffiles.txt. Before that, one must compile the Fotran files in the directory "common" according
to the ffiles.txt therein.

See the directory "examples" for an illustration about how to call COBYLA.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
