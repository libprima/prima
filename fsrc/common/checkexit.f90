module checkexit_mod
!--------------------------------------------------------------------------------------------------!
! This module checks whether to exit the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Friday, November 19, 2021 PM03:20:24
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: checkexit

interface checkexit
    module procedure checkexit_unc, checkexit_con
end interface checkexit


contains


function checkexit_unc(maxfun, nf, f, ftarget, x) result(info)
!--------------------------------------------------------------------------------------------------!
! This module checks whether to exit the solver in the unconstrained case.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_inf
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED

implicit none

! Inputs
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: x(:)

! Outputs
integer(IK) :: info

! Local variables
character(len=*), parameter :: srname = 'CHECKEXIT_UNC'

! Preconditions
if (DEBUGGING) then
    call assert(.not. any([NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED] == INFO_DFT), &
        & 'NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, and MAXFUN_REACHED differ from INFO_DFT', srname)
    ! X does not contain NaN if the initial X does not contain NaN and the subroutines generating
    ! trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! With the moderated extreme barrier, F cannot be NaN/+Inf.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

info = INFO_DFT  ! Default info, indicating that the solver should not exit.

! Although X should not contain NaN unless there is a bug, we write the following for security.
! X can be Inf, as finite + finite can be Inf numerically.
if (any(is_nan(x) .or. is_inf(x))) then
    info = NAN_INF_X
end if

! Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
if (is_nan(f) .or. is_posinf(f)) then
    info = NAN_INF_F
end if

if (f <= ftarget) then
    info = FTARGET_ACHIEVED
end if

if (nf >= maxfun) then
    info = MAXFUN_REACHED
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(any([INFO_DFT, NAN_INF_X, FTARGET_ACHIEVED, MAXFUN_REACHED] == info), &
        & 'INFO is NAN_INF_X, FTARGET_ACHIEVED, MAXFUN_REACHED, or INFO_DFT', srname)
end if

end function checkexit_unc


function checkexit_con(maxfun, nf, cstrv, ctol, f, ftarget, x) result(info)
!--------------------------------------------------------------------------------------------------!
! This module checks whether to exit the solver in the constrained case.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_inf
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED

implicit none

! Inputs
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: nf
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: ctol
real(RP), intent(in) :: f
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: x(:)

! Outputs
integer(IK) :: info

! Local variables
character(len=*), parameter :: srname = 'CHECKEXIT_CON'

! Preconditions
if (DEBUGGING) then
    call assert(.not. any([NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED] == INFO_DFT), &
        & 'NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, and MAXFUN_REACHED differ from INFO_DFT', srname)
    ! X does not contain NaN if the initial X does not contain NaN and the subroutines generating
    ! trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! With the moderated extreme barrier, F or CSTRV cannot be NaN/+Inf.
    call assert(.not. (is_nan(f) .or. is_posinf(f) .or. is_nan(cstrv) .or. is_posinf(cstrv)), &
        & 'F or CSTRV is not NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

info = INFO_DFT   ! Default info, indicating that the solver should not exit.

! Although X should not contain NaN unless there is a bug, we write the following for security.
! X can be Inf, as finite + finite can be Inf numerically.
if (any(is_nan(x) .or. is_inf(x))) then
    info = NAN_INF_X
end if

! Although NAN_INF_F should not happen unless there is a bug, we include the following for security.
if (is_nan(f) .or. is_posinf(f) .or. is_nan(cstrv) .or. is_posinf(cstrv)) then
    info = NAN_INF_F
end if

if (f <= ftarget .and. cstrv <= ctol) then
    info = FTARGET_ACHIEVED
end if

if (nf >= maxfun) then
    info = MAXFUN_REACHED
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(any([INFO_DFT, NAN_INF_F, FTARGET_ACHIEVED, MAXFUN_REACHED] == info), &
        & 'INFO is NAN_INF_X, FTARGET_ACHIEVED, MAXFUN_REACHED, or INFO_DFT', srname)
end if

end function checkexit_con


end module checkexit_mod
