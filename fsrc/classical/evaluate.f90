module evaluate_mod
!--------------------------------------------------------------------------------------------------!
! This is a module evaluating the objective/constraint function with Nan/Inf handling.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: August 2021
!
! Last Modified: Friday, January 21, 2022 AM04:06:13
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
private
public :: evalf
public :: evalfc
public :: xhist, fhist, chist, conhist
public :: nf
public :: fc_x0_provided, x0, f_x0, constr_x0

real(RP), allocatable :: xhist(:, :), fhist(:), chist(:), conhist(:, :)
integer(IK) :: nf

! N.B.: FC_X0_PROVIDED, X0, F_X0 and CONSTR_X0 are only used in nonlinear constrained problems,
! where the user may provide the function & constraint values of the starting point X0.
logical :: fc_x0_provided
real(RP), allocatable :: x0(:)
real(RP) :: f_x0
real(RP), allocatable :: constr_x0(:)


contains


subroutine evalf(calfun, x, f)
!--------------------------------------------------------------------------------------------------!
! This function evaluates CALFUN at X, setting F to the objective function value. Nan/Inf are
! handled by a moderated extreme barrier.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, HUGEFUN
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: pintrf_mod, only : OBJ
implicit none

! Inputs
procedure(OBJ) :: calfun
real(RP), intent(in) :: x(:)

! Output
real(RP), intent(out) :: f

! Local variables
character(len=*), parameter :: srname = 'EVALF'
integer(IK) :: khist
integer(IK) :: maxfhist
integer(IK) :: maxxhist

if (any(is_nan(x))) then
    ! Although this should not happen unless there is a bug, we include this case for security.
    f = sum(x)  ! Set F to NaN
else
    call calfun(x, f)  ! Evaluate F.

    ! Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
    ! finite value. This is naive. Better approaches surely exist.
    if (f > HUGEFUN .or. is_nan(f)) then
        f = HUGEFUN
    end if

    !! We may moderate huge negative values of F (NOT an extreme barrier), but we decide not to.
    !!f = max(-HUGEFUN, f)
end if

! Update NF and the history.
nf = nf + int(1, kind(nf))
maxxhist = int(size(xhist, 2), kind(maxxhist))
if (maxxhist >= 1) then
    khist = modulo(nf - 1_IK, maxxhist) + 1_IK
    xhist(:, khist) = x
end if
maxfhist = int(size(fhist), kind(maxfhist))
if (maxfhist >= 1) then
    khist = modulo(nf - 1_IK, maxfhist) + 1_IK
    fhist(khist) = f
end if

end subroutine evalf


subroutine evalfc(calcfc, x, f, constr, cstrv)
!--------------------------------------------------------------------------------------------------!
! This function evaluates CALCFC at X, setting F to the objective function value, CONSTR to the
! constraint value, and CSTRV to the constraint violation. Nan/Inf are handled by a moderated
! extreme barrier.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGEFUN, HUGECON
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : norm
use, non_intrinsic :: pintrf_mod, only : OBJCON
implicit none

! Inputs
procedure(OBJCON) :: calcfc
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: cstrv

! Local variables
character(len=*), parameter :: srname = 'EVALFC'
integer(IK) :: khist
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxxhist
logical :: evaluated

if (any(is_nan(x))) then
    ! Although this should not happen unless there is a bug, we include this case for security.
    ! Set F, CONSTR, and CSTRV to NaN.
    f = sum(x)
    constr = f
    cstrv = f
else
    evaluated = .false.
    if (fc_x0_provided .and. allocated(x0) .and. allocated(constr_x0)) then
        if (norm(x - x0) <= 0) then
            f = f_x0
            constr = constr_x0
            evaluated = .true.
        end if
    end if
    if (.not. evaluated) then
        call calcfc(x, f, constr)  ! Evaluate F and CONSTR.
    end if

    ! Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
    ! finite value. This is naive, and better approaches surely exist.
    if (f > HUGEFUN .or. is_nan(f)) then
        f = HUGEFUN
    end if
    where (constr < -HUGECON .or. is_nan(constr))
        ! The constraint is CONSTR(X) >= 0, so NaN should be replaced with a large negative value.
        constr = -HUGECON  ! MATLAB code: constr(constr < -HUGECON | isnan(constr)) = -HUGECON
    end where

    ! Moderate huge positive values of CONSTR, or they may lead to Inf/NaN in subsequent calculations.
    ! This is NOT an extreme barrier.
    constr = min(HUGECON, constr)
    !! We may moderate F similarly, but we decide not to.
    !!f = max(-HUGEFUN, f)

    ! Evaluate the constraint violation for constraints CONSTR(X) >= 0.
    cstrv = maxval([-constr, ZERO])
end if

! Update NF and the history.
nf = nf + int(1, kind(nf))
maxxhist = int(size(xhist, 2), kind(maxxhist))
if (maxxhist >= 1) then
    khist = modulo(nf - 1_IK, maxxhist) + 1_IK
    xhist(:, khist) = x
end if
maxfhist = int(size(fhist), kind(maxfhist))
if (maxfhist >= 1) then
    khist = modulo(nf - 1_IK, maxfhist) + 1_IK
    fhist(khist) = f
end if
maxchist = int(size(chist), kind(maxchist))
if (maxchist >= 1) then
    khist = modulo(nf - 1_IK, maxchist) + 1_IK
    chist(khist) = cstrv
end if
maxconhist = int(size(conhist, 2), kind(maxconhist))
if (maxconhist >= 1) then
    khist = modulo(nf - 1_IK, maxconhist) + 1_IK
    conhist(:, khist) = constr
end if

end subroutine evalfc


end module evaluate_mod
