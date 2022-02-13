program test
!--------------------------------------------------------------------------------------------------!
! This program tests the modernized version of Powell's solvers.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Sunday, February 13, 2022 PM05:03:36
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: datetime_mod, only : year, week
use, non_intrinsic :: test_solver_mod, only : test_solver
implicit none

integer :: yw

! YW is the random seed for the tests. It is altered weekly to test the solvers as much as possible.
yw = 100 * modulo(year(), 100) + week()
print *, 'The random seed is', yw

call test_solver(randseed=yw)

end program test
