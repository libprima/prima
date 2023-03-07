module redrho_mod
!--------------------------------------------------------------------------------------------------!
! This module provides a function that calculates RHO when it needs to be reduced.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Wednesday, March 08, 2023 AM01:10:18
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: redrho


contains


function redrho(rho_in, rhoend) result(rho)
!--------------------------------------------------------------------------------------------------!
! This function calculates RHO when it needs to be reduced.
! The scheme is shared by UOBYQA, NEWUOA, BOBYQA, LINCOA. For COBYLA, Powell's code reduces RHO by
! `RHO = HALF * RHO; IF (RHO <= 1.5_RP * RHOEND) RHO = RHOEND`, as specified in (11) of the COBYLA
! paper. However, this scheme seems to work better, especially after we introduce DELTA.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: rho_in
real(RP), intent(in) :: rhoend

! Outputs
real(RP) :: rho

real(RP) :: rho_ratio
character(len=*), parameter :: srname = 'REDRHO'

! Preconditions
if (DEBUGGING) then
    call assert(rho_in > rhoend .and. rhoend > 0, 'RHO_IN > RHOEND > 0', srname)
end if

!====================!
! Calculation starts !
!====================!

rho_ratio = rho_in / rhoend

if (rho_ratio > 250.0_RP) then
    rho = TENTH * rho_in
else if (rho_ratio <= 16.0_RP) then
    rho = rhoend
else
    rho = sqrt(rho_ratio) * rhoend  !rho = sqrt(rho_in * rhoend)
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(rho_in > rho .and. rho >= rhoend, 'RHO_IN > RHO >= RHOEND', srname)
end if
end function redrho


end module redrho_mod
