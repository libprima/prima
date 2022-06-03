This is the reference implementation of Powell's LINCOA algorithm.

Coded by Zaikun ZHANG (www.zhangzk.net) based on the paper

M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
Math. Program. Comput., 7:237--267, 2015

and Powell's code, with modernization and improvements.

Powell did not publish a paper to introduce the algorithm. The above-mentioned paper does not
describe LINCOA but discusses how to solve linearly-constrained trust-region subproblems.

Due to the dependency among the modules, the files must be compiled in the order indicated in
ffiles.txt. Before that, one must compile the Fotran files in the directory "common" according
to the ffiles.txt therein.

See the directory "examples" for an illustration about how to call the solver.

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

Started in February 2022.
