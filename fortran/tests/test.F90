program test
!--------------------------------------------------------------------------------------------------!
! This program tests the modernized version of Powell's solvers.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Saturday, June 24, 2023 PM11:54:28
!--------------------------------------------------------------------------------------------------!

#if !defined PRIMA_TESTDIM
#define PRIMA_TESTDIM small
#endif

use, non_intrinsic :: datetime_mod, only : year, week
use, non_intrinsic :: test_solver_mod, only : test_solver
implicit none

integer :: yw

! YW is the random seed for the tests. It is altered weekly to test the solvers as much as possible.
yw = 100 * modulo(year(), 100) + week()
print *, 'The random seed is', yw

call test_solver(randseed=yw, testdim=PRIMA_TESTDIM)

end program test
