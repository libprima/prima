module resolution_mod
!--------------------------------------------------------------------------------------------------!
! This module provides procedures to enhance the resolution of the solver by reducing RHO and DELTA,
! and updating CPEN, if applicable.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Sunday, October 10, 2021 AM03:18:58
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: resenhance

interface resenhance
    module procedure resenhance_unc, resenhance_nlc
end interface


contains


subroutine resenhance_unc(rhoend, delta, rho)
!--------------------------------------------------------------------------------------------------!
! This subroutine enhances the resolution of the solver in the unconstrained case.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO, HALF, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(inout) :: delta
real(RP), intent(inout) :: rho

! Local variables
real(RP) :: delta_old
real(RP) :: rho_old
real(RP) :: rho_ratio
character(len=*), parameter :: srname = 'RESENHANCE_UNC'

! Preconditions
if (DEBUGGING) then
    call assert(rhoend >= ZERO, 'RHOEND >= 0', srname)
    call assert(rho > rhoend, 'RHO > RHOEND', srname)
    call assert(delta >= rho, 'DELTA >= RHO', srname)
end if

!====================!
! Calculation starts !
!====================!

! Record the values of RHO and DELTA for checking.
rho_old = rho
delta_old = delta

delta = HALF * rho
rho_ratio = rho / rhoend
if (rho_ratio <= 16.0_RP) then
    rho = rhoend
elseif (rho_ratio <= 250.0_RP) then
    rho = sqrt(rho_ratio) * rhoend
else
    rho = TENTH * rho
end if
delta = max(delta, rho)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(rho < rho_old .and. rho >= rhoend, 'RHOEND <= RHO < RHO_OLD', srname)
    call assert(delta < delta_old .and. delta >= rho, 'RHO <= DELTA < DELTA_OLD', srname)
end if
end subroutine resenhance_unc


subroutine resenhance_nlc(conmat, fval, rhoend, cpen, rho)
!--------------------------------------------------------------------------------------------------!
! This subroutine enhances the resolution of the solver in the nonlinearly constrained case.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: conmat(:, :)
real(RP), intent(in) :: fval(:)
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(inout) :: cpen
real(RP), intent(inout) :: rho

! Local variables
real(RP) :: cmax(size(conmat, 1))
real(RP) :: cmin(size(conmat, 1))
real(RP) :: cpen_old
real(RP) :: denom
real(RP) :: rho_old
character(len=*), parameter :: srname = 'RESENHANCE_NLC'

! Preconditions
if (DEBUGGING) then
    call assert(size(fval) >= 1, 'SIZE(FVAL) >= 1', srname)
    call assert(size(conmat, 2) == size(fval), 'SIZE(CONMAT, 2) == SIZE(FVAL)', srname)
    call assert(rho > rhoend, 'RHO > RHOEND', srname)
    call assert(cpen >= ZERO, 'CPEN >= 0', srname)
end if

!====================!
! Calculation starts !
!====================!

! Record the values of RHO and CPEN for checking.
rho_old = rho
cpen_old = cpen

! See equation (11) in Section 3 of the COBYLA paper for the update of RHO.
rho = HALF * rho
if (rho <= 1.5E0_RP * rhoend) then
    rho = rhoend
end if

! See equations (12)--(13) in Section 3 of the COBYLA paper for the update of CPEN.
! If the original CPEN = 0, then the updated CPEN is also 0.
cmin = minval(conmat, dim=2)
cmax = maxval(conmat, dim=2)
if (any(cmin < HALF * cmax)) then
    denom = minval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
    cpen = min(cpen, (maxval(fval) - minval(fval)) / denom)
else
    cpen = ZERO
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(rho < rho_old .and. rho >= rhoend, 'RHEND <= RHO < RHO_OLD', srname)
    call assert(cpen <= cpen_old .and. cpen >= ZERO, '0 <= CPEN <= CPEN_OLD', srname)
end if
end subroutine resenhance_nlc

end module resolution_mod
