module param_mod
!--------------------------------------------------------------------------------------------------!
! This module defines some default values for the tests.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Saturday, October 09, 2021 PM12:37:30
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, TENTH
implicit none
private
public :: DIMSTRIDE_DFT, MINDIM_DFT, MAXDIM_DFT, NRAND_DFT, NOISE, F_NOISE, X_NOISE

integer(IK), parameter :: DIMSTRIDE_DFT = 1
integer(IK), parameter :: MAXDIM_DFT = 20
! Testing univariate problems can help us to uncover some bugs that can only occur in extreme cases.
integer(IK), parameter :: MINDIM_DFT = 1
integer(IK), parameter :: NRAND_DFT = 5
real(RP), parameter :: F_NOISE = 0.2_RP
real(RP), parameter :: NOISE = TENTH
real(RP), parameter :: X_NOISE = TENTH

end module param_mod
