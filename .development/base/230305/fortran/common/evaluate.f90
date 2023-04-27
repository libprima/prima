module evaluate_mod
!--------------------------------------------------------------------------------------------------!
! This is a module evaluating the objective/constraint function with Nan/Inf handling.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: August 2021
!
! Last Modified: Tuesday, December 13, 2022 PM12:54:25
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: moderatex
public :: moderatef
public :: moderatec
public :: evaluate

interface evaluate
    module procedure evaluatef, evaluatefc
end interface evaluate


contains


function moderatex(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function moderates a decision variable. It replaces NaN by 0 and Inf/-Inf by REALMAX/-REALMAX.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO, REALMAX
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : trueloc
implicit none

! Inputs
real(RP), intent(in) :: x(:)
! Outputs
real(RP) :: y(size(x))

y = x
y(trueloc(is_nan(x))) = ZERO
y = max(-REALMAX, min(REALMAX, y))
end function moderatex


pure elemental function moderatef(f) result(y)
!--------------------------------------------------------------------------------------------------!
! This function moderates the function value of a MINIMIZATION problem. It replaces NaN and any
! value above FUNCMAX by FUNCMAX.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, FUNCMAX
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: f
! Outputs
real(RP) :: y

y = f
if (is_nan(f)) then
    y = FUNCMAX
end if
y = min(FUNCMAX, y)
! We may moderate huge negative function values, but we decide not to.
!y = max(-FUNCMAX, min(FUNCMAX, y))
end function moderatef


function moderatec(c) result(y)
!--------------------------------------------------------------------------------------------------!
! This function moderates the constraint value, the constraint demanding this value to be NONNEGATIVE.
! It replaces NaN and any value below -CONSTRMAX by -CONSTRMAX, and any value above CONSTRMAX by CONSTRMAX.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, CONSTRMAX
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : trueloc
implicit none

! Inputs
real(RP), intent(in) :: c(:)
! Outputs
real(RP) :: y(size(c))

y = c
y(trueloc(is_nan(c))) = -CONSTRMAX
y = max(-CONSTRMAX, min(CONSTRMAX, y))
end function moderatec


subroutine evaluatef(calfun, x, f)
!--------------------------------------------------------------------------------------------------!
! This function evaluates CALFUN at X, setting F to the objective function value. Nan/Inf are
! handled by a moderated extreme barrier.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: pintrf_mod, only : OBJ
implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(in) :: x(:)

! Output
real(RP), intent(out) :: f

! Local variables
character(len=*), parameter :: srname = 'EVALUATEF'

! Preconditions
if (DEBUGGING) then
    ! X should not contain NaN if the initial X does not contain NaN and the subroutines generating
    ! trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
end if

!====================!
! Calculation starts !
!====================!

if (any(is_nan(x))) then
    ! Although this should not happen unless there is a bug, we include this case for robustness.
    f = sum(x)  ! Set F to NaN
else
    call calfun(moderatex(x), f)  ! Evaluate F; We moderate X before doing so.

    ! Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
    ! finite value. This is naive. Better approaches surely exist.
    f = moderatef(f)

    ! We may moderate huge negative values of F (NOT an extreme barrier), but we decide not to.
    !f = max(-FUNCMAX, f)
end if


!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    ! With X not containing NaN, and with the moderated extreme barrier, F cannot be NaN/+Inf.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
end if

end subroutine evaluatef


subroutine evaluatefc(calcfc, x, f, constr, cstrv)
!--------------------------------------------------------------------------------------------------!
! This function evaluates CALCFC at X, setting F to the objective function value, CONSTR to the
! constraint value, and CSTRV to the constraint violation. Nan/Inf are handled by a moderated
! extreme barrier.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: pintrf_mod, only : OBJCON
implicit none

! Inputs
procedure(OBJCON) :: calcfc ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: cstrv

! Local variables
character(len=*), parameter :: srname = 'EVALUATEFC'

! Preconditions
if (DEBUGGING) then
    ! X should not contain NaN if the initial X does not contain NaN and the subroutines generating
    ! trust-region/geometry steps work properly so that they never produce a step containing NaN/Inf.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
end if

!====================!
! Calculation starts !
!====================!

if (any(is_nan(x))) then
    ! Although this should not happen unless there is a bug, we include this case for robustness.
    ! Set F, CONSTR, and CSTRV to NaN.
    f = sum(x)
    constr = f
    cstrv = f
else
    call calcfc(moderatex(x), f, constr)  ! Evaluate F and CONSTR; We moderate X before doing so.

    ! Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
    ! finite value. This is naive, and better approaches surely exist.
    f = moderatef(f)
    constr = moderatec(constr)
    ! We may moderate huge negative values of F (NOT an extreme barrier), but we decide not to.
    !f = max(-FUNCMAX, f)

    ! Evaluate the constraint violation for constraints CONSTR(X) >= 0.
    cstrv = maxval([-constr, ZERO])
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    ! With X not containing NaN, and with the moderated extreme barrier, F cannot be NaN/+Inf, and
    ! CONSTR cannot be NaN/-Inf.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    call assert(.not. any(is_nan(constr) .or. is_neginf(constr)), &
        & 'CONSTR does not containt NaN/-Inf', srname)
    call assert(.not. (cstrv < 0 .or. is_nan(cstrv) .or. is_posinf(cstrv)), &
        & 'CSTRV is nonnegative and not NaN/+Inf', srname)
end if

end subroutine evaluatefc


end module evaluate_mod
