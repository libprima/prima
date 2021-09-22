program test
!--------------------------------------------------------------------------------------------------!
! This program tests the modernized version of Powell's solvers.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Wednesday, September 22, 2021 PM12:36:32
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: test_newuoa_mod, only : test_newuoa
implicit none

call test_newuoa()

end program test
