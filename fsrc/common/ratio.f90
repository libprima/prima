module ratio_mod

implicit none
private
public :: redrat


contains


function redrat(ared, pred, rshrink) result(ratio)

use consts_mod, only : RP, ZERO, ONE, HALF, HUGENUM, DEBUGGING
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
implicit none

! Inputs
real(RP) :: ared
real(RP) :: pred
real(RP) :: rshrink  ! When RATIO <= RSHRINK, DELTA will be shrunk.

! Output
real(RP) :: ratio

! Local variables
character(len=*), parameter :: srname = 'REDRAT'

if (DEBUGGING) then
    if (rshrink <= ZERO) then
        call errstop(srname, 'The ratio threshold for shrinking the trust-region radius is not positive')
    end if
end if

if (is_nan(pred) .or. pred <= ZERO) then
    ! The trust-region subproblem solver fails in this rare case. Instead of terminating as Powell's
    ! original code does, we set RATIO as follows so that the solver may continue to progress.
    if (ared > ZERO) then
        ! The trial point will be accepted, but the trust-region radius will be shrunk.
        ratio = HALF * rshrink
    else
        ! Signify a bad trust-region step, so that the solver will check whether to take the
        ! geometry step or reduce rho.
        ratio = -HUGENUM
    end if
elseif (is_posinf (pred) .and. is_posinf(ared)) then
    ratio = ONE
else
    ratio = ared / pred
end if

end function redrat


end module ratio_mod
