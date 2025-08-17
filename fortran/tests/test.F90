program test
!--------------------------------------------------------------------------------------------------!
! This program tests the modernized version of Powell's solvers.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Sun 17 Aug 2025 07:48:05 AM CST
!--------------------------------------------------------------------------------------------------!

#if !defined PRIMA_TESTDIM
#define PRIMA_TESTDIM 'small'
#endif

use, non_intrinsic :: datetime_mod, only : year, week
use, non_intrinsic :: test_solver_mod, only : test_solver
implicit none

integer :: seed

! SEED is the random seed for the tests.
! The default value of SEED below is altered weekly to test the solvers as much as possible.
seed = 100 * modulo(year(), 100) + week()
#if defined PRIMA_TESTSEED
seed = PRIMA_TESTSEED
#endif

print *, 'The random seed is', seed

! PRIMA_TESTDIM is the dimension of the test problem. It can be 'small', 'big', or 'large'.
! When it is 'small', then `test_solver` also accepts `mindim`, `maxdim`, and `dimstride`
! to specify the dimensions of the test problems; when it is 'big' or 'large', then dimension of
! the test is hard coded in `test_*.f90` for each solver.
print *, 'The test dimension is ', PRIMA_TESTDIM

call test_solver(randseed=seed, testdim=PRIMA_TESTDIM)

end program test
