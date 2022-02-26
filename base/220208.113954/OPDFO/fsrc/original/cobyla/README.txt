The code was sent by Professor Powell to Zaikun Zhang on December 16th, 2013.  
The file "email.txt" is the original email. The makefile was by Zhang. 
For more information on COBYLA, you might contact Professor Powell (mjdp@cam.ac.uk).

December 16th, 2013                   Zaikun Zhang (www.zhangzk.net) 


Below are the remarks from Professor Powell.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

                                 COBYLA
                                 ~~~~~~

     Here is a single-precision Fortran implementation of the algorithm for
constrained optimization that is the subject of the report I have written on
"A direct search optimization method that models the objective and constraint
functions by linear interpolation". This report has the number DAMTP 1992/NA5,
University of Cambridge, and it has been published in the proceedings of the
conference on Numerical Analysis and Optimization that was held in Oaxaca,
Mexico in January, 1992, which is the book "Advances in Optimization and
Numerical Analysis" (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer
Academic Publishers (1994).

     The instructions for using the Fortran code are given in the comments of
SUBROUTINE COBYLA, which is the interface between the user and the main
calculation that is done by SUBROUTINE COBYLB. There is a need for a linear
programming problem to be solved subject to a Euclidean norm trust region
constraint. Therefore SUBROUTINE TRSTLP is provided too, but you may have some
software that you prefer to use instead. These 3 subroutines are separated by
lines of hyphens below. Further, there follows the main program, the CALCFC
subroutine and the output that are appropriate to the numerical examples that
are discussed in the last section of DAMTP 1992/NA5. Please note, however,
that some cosmetic restructuring of the software has caused the given output
to differ slightly from Table 1 of the report.

     There are no restrictions on the use of the software, nor do I offer any
guarantees of success. Indeed, at the time of writing this note I had applied
it only to test problems that have up to 10 variables.

                        Mike Powell (May 7th, 1992).
