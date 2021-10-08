module resolution_mod

implicit none
private
public :: resenhance


interface resenhance
    module procedure resenhance_unc, resenhance_nlc
end interface


contains


subroutine resenhance_unc(rhoend, delta, rho)

use, non_intrinsic :: consts_mod, only : RP, HALF, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Input
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(inout) :: delta
real(RP), intent(inout) :: rho

! Local variables
real(RP) :: rho_ratio
character(len=*), parameter :: srname = 'RESENHANCE_UNC'


if (DEBUGGING) then
    if (rho <= rhoend) then
        call errstop(srname, 'RHO <= RHOEND')
    end if
end if

delta = HALF * rho
rho_ratio = rho / rhoend
if (rho_ratio <= 16.0_RP) then
    rho = rhoend
else if (rho_ratio <= 250.0_RP) then
    rho = sqrt(rho_ratio) * rhoend
else
    rho = TENTH * rho
end if
delta = max(delta, rho)

end subroutine resenhance_unc


subroutine resenhance_nlc(conmat, fval, rhoend, cpen, rho)

use, non_intrinsic :: consts_mod, only : RP, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop
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
real(RP) :: denom
character(len=*), parameter :: srname = 'RESENHANCE_NLC'


if (DEBUGGING) then
    if (size(fval) < 1) then
        call errstop(srname, 'SIZE(FVAL) < 1')
    end if
    if (size(conmat, 2) /= size(fval)) then
        call errstop(srname, 'SIZE(CONMAT, 2) /= SIZE(FVAL)')
    end if
    if (rho <= rhoend) then
        call errstop(srname, 'RHO <= RHOEND')
    end if
end if

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

end subroutine resenhance_nlc

end module resolution_mod
