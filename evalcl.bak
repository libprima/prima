module evalcl_mod
!--------------------------------------------------------------------------------------------------!
! This is a module evaluating the objective/constraint functions for the CLASSICAL MODE of Powell's
! solvers. It does the following in addition to the function evaluations.
! 1. Handle Nan/Inf with moderated extreme barriers.
! 2. Record the number of function evaluations (NF) and the history of evaluations.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: August 2021
!
! Last Modified: Wednesday, February 09, 2022 AM12:50:03
!--------------------------------------------------------------------------------------------------!

! N.B.: The module is NOT used by the modernized version of Powell's solvers. It is used only when we
! invoke Powell's Fortran 77 code, i.e., when we call the solvers in the CLASSICAL MODE. We provide
! only limited support for the classical mode. Regarding the current module, the limitation is
! reflected in the following aspects.
! 1. The implementation of this module is NOT thread safe due to the module variables NF, XHIST etc.
! 2. The preconditions and postconditions are only enforced minimally compared to the modern version.
! 3. The argument NF in RANGEHIST mask the module variables with the same name, which is not good
! practice. Similar things can be said about XHIST etc.
! 4. The definition of EVALUATE relies on a procedure of the same name from the module evaluate_mod.
! Even though we rename the latter as EVAL, some may find it not perfect practice.

use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
private
public :: evaluate
public :: rangehist
public :: nf
public :: xhist, fhist, chist, conhist
public :: fc_x0_provided, x0, f_x0, constr_x0

interface evaluate
    module procedure evaluatef, evaluatefc
end interface evaluate

integer(IK) :: nf
real(RP), allocatable :: xhist(:, :), fhist(:), chist(:), conhist(:, :)

! N.B.: FC_X0_PROVIDED, X0, F_X0 and CONSTR_X0 are only used in nonlinear constrained problems,
! where the user may provide the function & constraint values of the starting point X0.
logical :: fc_x0_provided
real(RP), allocatable :: x0(:)
real(RP) :: f_x0
real(RP), allocatable :: constr_x0(:)


contains


subroutine evaluatef(calfun, x, f)
!--------------------------------------------------------------------------------------------------!
! This function evaluates CALFUN at X, setting F to the objective function value. Nan/Inf are
! handled by a moderated extreme barrier. In addition, the function updates NF and the history.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: evaluate_mod, only : eval => evaluate
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: pintrf_mod, only : OBJ
implicit none

! Inputs
procedure(OBJ) :: calfun
real(RP), intent(in) :: x(:)

! Output
real(RP), intent(out) :: f

! Local variables
integer(IK) :: khist
integer(IK) :: maxfhist
integer(IK) :: maxxhist

if (any(is_nan(x))) then   ! Set F to NaN
    f = sum(x)
else
    call eval(calfun, x, f)
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

end subroutine evaluatef


subroutine evaluatefc(calcfc, x, f, constr, cstrv)
!--------------------------------------------------------------------------------------------------!
! This function evaluates CALCFC at X, setting F to the objective function value, CONSTR to the
! constraint value, and CSTRV to the constraint violation. Nan/Inf are handled by a moderated
! extreme barrier. In addition, the function updates NF and the history.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
use, non_intrinsic :: evaluate_mod, only : moderatef, moderatec, eval => evaluate
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
integer(IK) :: khist
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxxhist
logical :: evaluated

if (any(is_nan(x))) then  ! Set F, CONSTR, and CSTRV to NaN.
    f = sum(x)
    constr = f
    cstrv = f
else
    evaluated = .false.
    if (fc_x0_provided .and. allocated(x0) .and. allocated(constr_x0)) then
        if (norm(x - x0) <= 0) then
            f = moderatef(f_x0)
            constr = moderatec(constr_x0)
            cstrv = maxval([-constr, ZERO])
            evaluated = .true.
        end if
    end if
    if (.not. evaluated) then
        call eval(calcfc, x, f, constr, cstrv)
    end if
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

end subroutine evaluatefc


subroutine rangehist(nf, xhist, fhist, chist, conhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine arranges FHIST, XHIST, CHIST, and CONHIST in the chronological order.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: nf

! In-outputs
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)
real(RP), intent(inout), optional :: chist(:)
real(RP), intent(inout), optional :: conhist(:, :)

! Local variables
integer(IK) :: khist
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
character(len=*), parameter :: srname = 'RANGEHIST'

! Sizes
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
if (present(chist)) then
    maxchist = int(size(chist), kind(maxchist))
else
    maxchist = 0_IK
end if
if (present(conhist)) then
    maxconhist = int(size(conhist, 2), kind(maxconhist))
else
    maxconhist = 0_IK
end if
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)

! Preconditions
if (DEBUGGING) then
    ! Check the sizes of XHIST, FHIST, CHIST, CONHIST.
    call assert(maxxhist * (maxxhist - maxhist) == 0, 'SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    call assert(maxconhist * (maxconhist - maxhist) == 0, 'SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
end if

!====================!
! Calculation starts !
!====================!

! The ranging should be done only if 0 < MAXXHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxxhist > 0 .and. maxxhist < nf) then
    ! We could replace MODULO(NF - 1_IK, MAXXHIST) + 1_IK) with MODULO(NF - 1_IK, MAXHIST) + 1_IK)
    ! based on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    khist = modulo(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! N.B.:
    ! 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
    ! elements are used to construct the array.
    ! 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
    ! 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
end if
! The ranging should be done only if 0 < MAXFHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxfhist > 0 .and. maxfhist < nf) then
    khist = modulo(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if
! The ranging should be done only if 0 < MAXCONHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxconhist > 0 .and. maxconhist < nf) then
    khist = modulo(nf - 1_IK, maxconhist) + 1_IK
    conhist = reshape([conhist(:, khist + 1:maxconhist), conhist(:, 1:khist)], shape(conhist))
end if
! The ranging should be done only if 0 < MAXCHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxchist > 0 .and. maxchist < nf) then
    khist = modulo(nf - 1_IK, maxchist) + 1_IK
    chist = [chist(khist + 1:maxchist), chist(1:khist)]
end if

end subroutine rangehist


end module evalcl_mod
