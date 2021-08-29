module resolution_mod

implicit none
private
public :: resenhance


contains


subroutine resenhance(conmat, fval, rhoend, cpen, rho)

use consts_mod, only : RP, ZERO, HALF, DEBUGGING
use debug_mod, only : errstop, verisize
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
character(len=*), parameter :: srname = 'RESENHANCE'


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

end subroutine resenhance

end module resolution_mod
