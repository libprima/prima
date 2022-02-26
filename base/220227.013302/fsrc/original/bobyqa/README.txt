The code was sent by Professor Powell to Zaikun Zhang on December 16th, 2013.  
The file "email.txt" is the original email. For more information on BOBYQA, 
you might contact Professor Powell (mjdp@cam.ac.uk).

December 16th, 2013                   Zaikun Zhang (www.zhangzk.net) 


Below are the remarks from Professor Powell.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

     The Fortran version of BOBYQA is attached. Its purpose is to seek
the least value of a function F of several variables, when derivatives
are not available, where F is specified by the user through a subroutine
called CALFUN. The name BOBYQA denotes Bound Approximation BY Quadratic
Approximation, the constraints being lower and upper bounds on every
variable, which can be set to huge values for unconstrained variables.
The algorithm is intended to change the variables to values that are close
to a local minimum of F. The user, however, should assume responsibility for
finding out if the calculations are satisfactory, by considering carefully
the values of F that occur. Details of the method of BOBYQA are given in
the report "The BOBYQA algorithm for bound constrained optimization without
derivatives", which can be reached from the "damtp.cam.ac.uk" home-page on
the web, by clicking on "Research at DAMTP", then on "Numerical Analysis"
and then on "Reports", the number of the report being 2009/NA06.

     The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the "Invdist2" problem, in order to provide
an example for testing. Then BOBYQA and its six auxiliary routines, namely
BOBYQB, ALTMOV, PRELIM, RESCUE, TRSBOX and UPDATE, are given. Finally, the
computed output that the author obtained for the Invdist2 problems is listed.

     In addition to providing CALFUN, an initial vector of variables and
the lower and upper bounds, the user has to set the values of the parameters
RHOBEG, RHOEND and NPT. After scaling the individual variables if necessary,
so that the magnitudes of their expected changes are similar, RHOBEG is the
initial steplength for changes to the variables, a reasonable choice being
the mesh size of a coarse grid search. Further, RHOEND should be suitable for
a search on a very fine grid. Typically, the software calculates a vector
of variables that is within distance 10*RHOEND of a local minimum. Another
consideration is that every trial vector of variables is forced to satisfy
the lower and upper bounds, but there has to be room to make a search in all
directions. Therefore an error return occurs if the difference between the
bounds on any variable is less than 2*RHOBEG. The parameter NPT specifies
the number of interpolation conditions on each quadratic model, the value
NPT=2*N+1 being recommended for a start, where N is the number of variables.
It is often worthwhile to try other choices too, but much larger values tend
to be inefficient, because the amount of routine work of each iteration is
of magnitude NPT**2, and because the achievement of adequate accuracy in some
matrix calculations becomes more difficult. Some excellent numerical results
have been found in the case NPT=N+6 even with more than 100 variables.

     The way of calling BOBYQA should be clear from the Invdist2 examples
and from the comments near the beginning of SUBROUTINE BOBYQA. There are no
restrictions on or charges for the use of the software. I hope that the time
and effort I have spent on developing the package will be helpful to much
research and to many applications.

January 5th, 2009                    M.J.D. Powell (mjdp@cam.ac.uk)
