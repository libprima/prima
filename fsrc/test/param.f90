module param_mod
!--------------------------------------------------------------------------------------------------!
! This module defines some default values for the tests.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Thursday, September 30, 2021 AM12:04:09
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, TENTH
implicit none
private
public :: DIMSTRIDE_DFT, MINDIM_DFT, MAXDIM_DFT, NRAND_DFT, NOISE, F_NOISE, X_NOISE

integer(IK), parameter :: DIMSTRIDE_DFT = 1
integer(IK), parameter :: MAXDIM_DFT = 20
integer(IK), parameter :: MINDIM_DFT = 1
integer(IK), parameter :: NRAND_DFT = 5
real(RP), parameter :: F_NOISE = 0.2_RP
real(RP), parameter :: NOISE = TENTH
real(RP), parameter :: X_NOISE = TENTH

end module param_mod
