module checkexit_mod

implicit none
private
public :: checkexit

interface checkexit
    module procedure checkexit_unc, checkexit_con
end interface checkexit


contains


function checkexit_unc(maxfun, nf, f, ftarget, x) result(info)

! Generic modules
use consts_mod, only : RP, IK, DEBUGGING
use info_mod, only : INFO_DFT, NAN_X, FTARGET_ACHIEVED, MAXFUN_REACHED
use infnan_mod, only : is_nan
use debug_mod, only : errstop

implicit none

! Inputs
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: x(:)

! Output
integer(IK) :: info

character(len=*), parameter :: srname = 'CHECKEXIT_UNC'

if (DEBUGGING) then
    if (NAN_X == INFO_DFT .or. FTARGET_ACHIEVED == INFO_DFT .or. MAXFUN_REACHED == INFO_DFT) then
        call errstop(srname, 'Invalid info codes.')
    end if
end if

info = INFO_DFT  ! Default info.

if (any(is_nan(x))) then
    info = NAN_X
end if
! N.B.: With the moderated extreme barrier, F or CSTRV cannot be Inf/NaN.
!if (is_nan(f) .or. is_posinf(f)) then
!    info = NAN_INF_F
!end if
if (f <= ftarget) then
    info = FTARGET_ACHIEVED
end if
if (nf >= maxfun) then
    info = MAXFUN_REACHED
end if

end function checkexit_unc


function checkexit_con(maxfun, nf, cstrv, ctol, f, ftarget, x) result(info)

! Generic modules
use consts_mod, only : RP, IK, DEBUGGING
use info_mod, only : INFO_DFT, NAN_X, FTARGET_ACHIEVED, MAXFUN_REACHED
use infnan_mod, only : is_nan
use debug_mod, only : errstop

implicit none

! Inputs
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: nf
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: ctol
real(RP), intent(in) :: f
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: x(:)

! Output
integer(IK) :: info

character(len=*), parameter :: srname = 'CHECKEXIT_CON'

if (DEBUGGING) then
    if (NAN_X == INFO_DFT .or. FTARGET_ACHIEVED == INFO_DFT .or. MAXFUN_REACHED == INFO_DFT) then
        call errstop(srname, 'Invalid info codes.')
    end if
end if

info = INFO_DFT  ! Default info.

if (any(is_nan(x))) then
    info = NAN_X
end if
! N.B.: With the moderated extreme barrier, F or CSTRV cannot be Inf/NaN.
!if (is_nan(f) .or. is_posinf(f) .or. is_nan(cstrv) .or. is_posinf(cstrv)) then
!    info = NAN_INF_F
!end if
if (f <= ftarget .and. cstrv <= ctol) then
    info = FTARGET_ACHIEVED
end if
if (nf >= maxfun) then
    info = MAXFUN_REACHED
end if

end function checkexit_con

end module checkexit_mod
