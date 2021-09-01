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
use info_mod, only : INFO_DFT, NAN_X, NaN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED
use infnan_mod, only : is_nan, is_posinf
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
    ! NAN_X should never happen if the initial X does not contain NaN and the subroutines generating
    ! trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
    if (any(is_nan(x))) then
        call errstop(srname, 'NaN occurs in X')
    end if
    ! With the moderated extreme barrier, F cannot be Inf/NaN.
    if (is_nan(f) .or. is_posinf(f)) then
        call errstop(srname, 'NaN occurs in F')
    end if
end if

info = INFO_DFT  ! Default info.

! This should never happen unless there is a bug.
if (any(is_nan(x))) then
    info = NAN_X
end if
! With the moderated extreme barrier, this should never happen unless there is a bug.
if (is_nan(f) .or. is_posinf(f)) then
    info = NAN_INF_F
end if
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
use info_mod, only : INFO_DFT, NAN_X, NAN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED
use infnan_mod, only : is_nan, is_posinf
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
    ! NAN_X should never happen if the initial X does not contain NaN and the subroutines generating
    ! trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
    if (any(is_nan(x))) then
        call errstop(srname, 'NaN occurs in X')
    end if
    ! With the moderated extreme barrier, F or CSTRV cannot be Inf/NaN.
    if (is_nan(f) .or. is_posinf(f) .or. is_nan(cstrv) .or. is_posinf(cstrv)) then
        call errstop(srname, 'NaN occurs in F or CSTRV')
    end if
end if

info = INFO_DFT  ! Default info.

! This should never happen unless there is a bug.
if (any(is_nan(x))) then
    info = NAN_X
end if
! With the moderated extreme barrier, this should never happen unless there is a bug.
if (is_nan(f) .or. is_posinf(f) .or. is_nan(cstrv) .or. is_posinf(cstrv)) then
    info = NAN_INF_F
end if
if (f <= ftarget .and. cstrv <= ctol) then
    info = FTARGET_ACHIEVED
end if
if (nf >= maxfun) then
    info = MAXFUN_REACHED
end if

end function checkexit_con

end module checkexit_mod
