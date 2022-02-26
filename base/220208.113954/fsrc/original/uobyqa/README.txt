The code was sent by Professor Powell to Zaikun Zhang on April 6th, 2015.  
The file "email.txt" is the original email. For more information on
UOBYQA, you might contact Professor Powell (mjdp@cam.ac.uk).

April 6th, 2015                   Zaikun Zhang (www.zhangzk.net) 


Below are the remarks from Professor Powell.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

     The Fortran version of UOBYQA, written by M.J.D. Powell, is attached.
Its purpose is to seek the least value of a function F of several variables,
when derivatives are not available, where F is specified by the user through
a subroutine called CALFUN. The algorithm is intended to change the variables
to values that are close to a local minimum of F. The user, however, should
assume responsibility for finding out if the calculations are satisfactory,
by giving careful attention to the values of F that occur. The details of
the method are described in "UOBYQA: unconstrained optimization by quadratic
approximation" by M.J.D. Powell, Mathematical Programming Series B, Volume
92, pages 555-582 (2002).

     The attachments in sequence are a suitable Makefile, followed by a main
program and a CALFUN routine for the Chebyquad problems, in order to provide
an example for testing. Then UOBYQA and its three auxiliary routines, namely
UOBYQB, TRSTEP and LAGMAX, are given. Finally, the computed output that the
author obtained for the Chebyquad problems is listed.

     The way of calling UOBYQA should be clear from the given example and
from the comments of that subroutine. It is hoped that the software will
be helpful to much future research and to many applications. There are no
restrictions on or charges for its use. If you wish to refer to it, please
cite the paper that is mentioned above.
