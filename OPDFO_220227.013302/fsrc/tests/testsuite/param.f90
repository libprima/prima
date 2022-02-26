module param_mod
!--------------------------------------------------------------------------------------------------!
! This module defines some default values for the tests.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Sunday, December 05, 2021 PM07:12:34
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, TENTH
implicit none
private
public :: DIMSTRIDE_DFT, MINDIM_DFT, MAXDIM_DFT, NRAND_DFT, NOISE_LEVEL_DFT, NOISE_TYPE_DFT

! Use an odd stride so that both odd and even dimensional problems will be tested.
integer(IK), parameter :: DIMSTRIDE_DFT = 3
! Testing univariate problems can help us to uncover some bugs that can only occur in extreme cases.
integer(IK), parameter :: MINDIM_DFT = 1
integer(IK), parameter :: MAXDIM_DFT = 20
integer(IK), parameter :: NRAND_DFT = 5
real(RP), parameter :: NOISE_LEVEL_DFT = TENTH
character(len=*), parameter :: NOISE_TYPE_DFT = 'gaussian'

end module param_mod
