The code was sent by Professor Powell to Zaikun Zhang on December 14th, 2013.  
The file "email.txt" is the original email. For more information on LINCOA, 
you might contact Professor Powell (mjdp@cam.ac.uk).

December 15th, 2013                   Zaikun Zhang (www.zhangzk.net) 


Below are the remarks from Professor Powell.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

     The Fortran version of LINCOA is attached. Its purpose is to seek the
least value of a function F of several variables subject to general linear
inequality constraints on the variables, when derivatives of F are not
available. The name LINCOA denotes LINearly Constrained Optimization
Algorithm. F is specified by the user through a subroutine called CALFUN.
The algorithm is intended to change the variables to values that are close
to a local constrained minimum of F. The user, however, should assume
responsibility for finding out if the calculations are adequate. It may be
helpful to employ several starting points in the space of the variables and
to try different values of the parameters NPT and RHOEND. I intend to write
a paper that explains briefly the main features of the software.

     LINCOA is not suitable for very large numbers of variables because no
attention is given to any sparsity. A few calculations with 1000 variables,
however, have been run successfully overnight, and the performance of LINCOA
is satisfactory usually for small numbers of variables. Several calculations
of the objective function may be required at points that do not satisfy the
linear constraints, especially if an equality constraint is expressed as
two inequalities.

     The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the "PtsinTet" problem, in order to provide
an example for testing. Then LINCOA and its six auxiliary routines, namely
LINCOB, GETACT, PRELIM, QMSTEP, TRSTEP and UPDATE, are given. Finally, the
output from the author's computer for the PtsinTet problem is listed.

     In addition to providing CALFUN, the linear constraints, and an initial
vector of variables, the user has to set the values of RHOBEG, RHOEND and
NPT. After scaling the individual variables, so that the magnitudes of their
expected changes are similar, RHOBEG is the initial steplength for changes
to the variables, a reasonable choice being the mesh size of a coarse grid
search. Further, RHOEND should be suitable for a search on a very fine grid.
Typically, the final vector of variables is within distance 10*RHOEND of
a local minimum. The parameter NPT specifies the number of interpolation
conditions on each quadratic model, the value NPT=2*N+1 being recommended
for a start, where N is the number of variables.

     The way of calling LINCOA should be clear from the PtsinTet example
and from the comments near the beginning of SUBROUTINE LINCOA. There are no
restrictions on or charges for the use of the software. I hope that the time
and effort I have spent on developing the package will be helpful to much
research and to many applications.

December 6th, 2013                    M.J.D. Powell (mjdp@cam.ac.uk)
