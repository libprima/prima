module selectx_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines that ensure the returned X is optimal among all the calculated
! points in the sense that no other point achieves both lower function value and lower constraint
! violation at the same time. The module is needed only in the constrained case.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Saturday, March 04, 2023 PM01:37:33
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: savefilt, selectx, isbetter

interface isbetter
    module procedure isbetter10, isbetter01
end interface isbetter

contains


subroutine savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, constr, confilt)
!--------------------------------------------------------------------------------------------------!
! This subroutine saves X, F, and CSTRV in XFILT, FFILT, and CFILT (and CONSTR in CONFILT if they
! are present), unless a vector in XFILT(:, 1:NFILT) is better than X. If X is better than some
! vectors in XFILT(:, 1:NFILT), then these vectors will be removed. If X is not better than any of X
! FILT(:, 1:NFILT) but NFILT=MAXFILT, then we remove a column from XFILT according to the merit
! function PHI = FFILT + CWEIGHT * MAX(CFILT - CTOL, ZERO).
! N.B.:
! 1. Only XFILT(:, 1:NFILT) and FFILT(:, 1:NFILT) etc contains valid information, while
! XFILT(:, NFILT+1:MAXFILT) and FFILT(:, NFILT+1:MAXFILT) etc are not initialized yet.
! 2. We decide whether an X is better than another by the ISBETTER function.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: linalg_mod, only : trueloc

implicit none

! Inputs
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: ctol
real(RP), intent(in) :: cweight
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)  ! N
real(RP), intent(in), optional :: constr(:)  ! M

! In-outputs
integer(IK), intent(inout) :: nfilt
real(RP), intent(inout) :: cfilt(:)  ! MAXFILT
real(RP), intent(inout) :: ffilt(:)  ! MAXFILT
real(RP), intent(inout) :: xfilt(:, :) ! (N, MAXFILT)
real(RP), intent(inout), optional :: confilt(:, :)  ! (M, MAXFILT)

! Local variables
character(len=*), parameter :: srname = 'SAVEFILT'
integer(IK) :: index_to_keep(size(ffilt))
integer(IK) :: m
integer(IK) :: maxfilt
integer(IK) :: n
integer(IK) :: kworst
logical :: keep(nfilt)
real(RP) :: cfilt_shifted(size(ffilt))
real(RP) :: cref
real(RP) :: fref
real(RP) :: phi(size(ffilt))
real(RP) :: phimax

! Sizes
if (present(constr)) then
    m = int(size(constr), kind(m))
else
    m = 0
end if
n = int(size(x), kind(n))
maxfilt = int(size(ffilt), kind(maxfilt))

! Preconditions
if (DEBUGGING) then
    ! Check the size of X.
    call assert(n >= 1, 'N >= 1', srname)
    ! Check NFILT
    call assert(nfilt >= 0 .and. nfilt <= maxfilt, '0 <= NFILT <= MAXFILT', srname)
    ! Check the sizes of XFILT, FFILT, CFILT.
    call assert(maxfilt >= 1, 'MAXFILT >= 1', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == maxfilt, 'SIZE(XFILT) == [N, MAXFILT]', srname)
    call assert(size(cfilt) == maxfilt, 'SIZE(CFILT) == MAXFILT', srname)
    ! Check the values of XFILT, FFILT, CFILT.
    call assert(.not. any(is_nan(xfilt(:, 1:nfilt))), 'XFILT does not contain NaN', srname)
    call assert(.not. any(is_nan(ffilt(1:nfilt)) .or. is_posinf(ffilt(1:nfilt))), &
        & 'FFILT does not contain NaN/+Inf', srname)
    call assert(.not. any(cfilt(1:nfilt) < 0 .or. is_nan(cfilt(1:nfilt)) .or. is_posinf(cfilt(1:nfilt))), &
        & 'CFILT does not contain nonnegative values of NaN/+Inf', srname)
    ! Check the values of X, F, CSTRV.
    ! X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! F cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    ! CSTRV cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (cstrv < 0 .or. is_nan(cstrv) .or. is_posinf(cstrv)), 'CSTRV is nonnegative and not NaN/+Inf', srname)
    ! Check CONSTR and CONFILT.
    call assert(present(constr) .eqv. present(confilt), 'CONSTR and CONFILT are both present or both absent', srname)
    if (present(constr)) then
        ! CONSTR cannot contain NaN/-Inf due to the moderated extreme barrier.
        call assert(.not. any(is_nan(constr) .or. is_neginf(constr)), 'CONSTR does not contain NaN/-Inf', srname)
        call assert(size(confilt, 1) == m .and. size(confilt, 2) == maxfilt, 'SIZE(CONFILT) == [M, MAXFILT]', srname)
        call assert(.not. any(is_nan(confilt(:, 1:nfilt)) .or. is_neginf(confilt(:, 1:nfilt))), &
            & 'CONFILT does not contain NaN/-Inf', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! Return immediately if any column of XFILT is better than X.
if (any(isbetter(ffilt(1:nfilt), cfilt(1:nfilt), f, cstrv, ctol))) then
    return
end if

! Decide which columns of XFILT to keep.
keep = (.not. isbetter(f, cstrv, ffilt(1:nfilt), cfilt(1:nfilt), ctol))

! If NFILT == MAXFILT and X is not better than any column of XFILT, then we remove the worst column
! of XFILT according to the merit function PHI = FFILT + CWEIGHT * MAX(CFILT - CTOL, ZERO).
if (count(keep) == maxfilt) then  ! In this case, NFILT = SIZE(KEEP) = COUNT(KEEP) = MAXFILT > 0.
    cfilt_shifted = max(cfilt - ctol, ZERO)
    if (cweight <= 0) then
        phi = ffilt
    elseif (is_posinf(cweight)) then
        phi = cfilt_shifted
        ! We should not use CFILT here; if MAX(CFILT_SHIFTED) is attained at multiple indices, then
        ! we will check FFILT to exhaust the remaining degree of freedom.
    else
        phi = max(ffilt, -REALMAX) + cweight * cfilt_shifted
        ! MAX(FFILT, -REALMAX) makes sure that PHI will not contain NaN (unless there is a bug).
    end if
    ! We select X to maximize PHI. In case there are multiple maximizers, we take the one with the
    ! largest CSTRV_SHIFTED; if there are more than one choices, we take the one with the largest F;
    ! if there are several candidates, we take the one with the largest CSTRV; if the last comparison
    ! still leads to more than one possibilities, then they are equally bad and we choose the first.
    ! N.B.:
    ! 1. This process is the opposite of selecting KOPT in SELECTX.
    ! 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 == CSTRV_SHIFTED_2 do
    ! not ensure that F_1 == F_2!
    phimax = maxval(phi)
    cref = maxval(cfilt_shifted, mask=(phi >= phimax))
    fref = maxval(ffilt, mask=(cfilt_shifted >= cref))
    kworst = int(maxloc(cfilt, mask=(ffilt >= fref), dim=1), kind(kworst))
    !!MATLAB: cmax = max(cfilt(ffilt >= fref)); kworst = find(ffilt >= fref & ~(cfilt < cmax), 1,'first');
    if (kworst < 1 .or. kworst > size(keep)) then  ! For security. Should not happen.
        kworst = 1
    end if
    keep(kworst) = .false.
end if

nfilt = int(count(keep), kind(nfilt))
index_to_keep(1:nfilt) = trueloc(keep)
xfilt(:, 1:nfilt) = xfilt(:, index_to_keep(1:nfilt))
ffilt(1:nfilt) = ffilt(index_to_keep(1:nfilt))
cfilt(1:nfilt) = cfilt(index_to_keep(1:nfilt))
if (present(confilt) .and. present(constr)) then
    confilt(:, 1:nfilt) = confilt(:, index_to_keep(1:nfilt))
end if

nfilt = nfilt + 1_IK
xfilt(:, nfilt) = x
ffilt(nfilt) = f
cfilt(nfilt) = cstrv
if (present(confilt) .and. present(constr)) then
    confilt(:, nfilt) = constr
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    ! Check NFILT and the sizes of XFILT, FFILT, CFILT.
    call assert(nfilt >= 1 .and. nfilt <= maxfilt, '1 <= NFILT <= MAXFILT', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == maxfilt, 'SIZE(XFILT) == [N, MAXFILT]', srname)
    call assert(size(ffilt) == maxfilt, 'SIZE(FFILT) = MAXFILT', srname)
    call assert(size(cfilt) == maxfilt, 'SIZE(CFILT) = MAXFILT', srname)
    ! Check the values of XFILT, FFILT, CFILT.
    call assert(.not. any(is_nan(xfilt(:, 1:nfilt))), 'XFILT does not contain NaN', srname)
    call assert(.not. any(is_nan(ffilt(1:nfilt)) .or. is_posinf(ffilt(1:nfilt))), &
        & 'FFILT does not contain NaN/+Inf', srname)
    call assert(.not. any(cfilt(1:nfilt) < 0 .or. is_nan(cfilt(1:nfilt)) .or. is_posinf(cfilt(1:nfilt))), &
        & 'CFILT does not contain nonnegative values of NaN/+Inf', srname)
    ! Check that no point in the filter is better than X, and X is better than no point.
    call assert(.not. any(isbetter(ffilt(1:nfilt), cfilt(1:nfilt), f, cstrv, ctol)), &
        & 'No point in the filter is better than X', srname)
    call assert(.not. any(isbetter(f, cstrv, ffilt(1:nfilt), cfilt(1:nfilt), ctol)), &
        & 'X is better than no point in the filter', srname)
    ! Check CONFILT.
    if (present(confilt)) then
        call assert(size(confilt, 1) == m .and. size(confilt, 2) == maxfilt, 'SIZE(CONFILT) == [M, MAXFILT]', srname)
        call assert(.not. any(is_nan(confilt(:, 1:nfilt)) .or. is_neginf(confilt(:, 1:nfilt))), &
            & 'CONFILT does not contain NaN/-Inf', srname)
    end if
end if

end subroutine savefilt


function selectx(fhist, chist, cweight, ctol) result(kopt)
!--------------------------------------------------------------------------------------------------!
! This subroutine selects X according to the FHIST and CHIST, which represents (a part of) history
! of F and CSTRV. Normally, FHIST and CHIST are not the full history but only a filter, e.g., FFILT
! and CFILT generated by SAVEFILT. However, we name them as FHIST and CHIST because the [F, CSTRV]
! in a filter should not dominate each other, but this subroutine does NOT assume such a property.
! N.B.: CTOL is the tolerance of constraint violation (CSTRV). A point is considered feasible if
! its constraint violation is at most CTOL. Note that CTOL is absolute, not relative.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : IK, RP, EPS, REALMAX, FUNCMAX, CONSTRMAX, ZERO, TWO, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
real(RP), intent(in) :: cweight
real(RP), intent(in) :: chist(:)
real(RP), intent(in) :: ctol
real(RP), intent(in) :: fhist(:)

! Outputs
integer(IK) :: kopt

! Local variables
character(len=*), parameter :: srname = 'SELECTX'
integer(IK) :: nhist
real(RP) :: chist_shifted(size(fhist))
real(RP) :: cmin
real(RP) :: cref
real(RP) :: fref
real(RP) :: phi(size(fhist))
real(RP) :: phimin

! Sizes
nhist = int(size(fhist), IK)

! Preconditions
if (DEBUGGING) then
    call assert(nhist >= 1, 'SIZE(FHIST) >= 1', srname)
    call assert(size(chist) == nhist, 'SIZE(FHIST) == SIZE(CHIST)', srname)
    call assert(.not. any(is_nan(fhist) .or. is_posinf(fhist)), &
        & 'FHIST does not contain NaN/+Inf', srname)
    call assert(.not. any(chist < 0 .or. is_nan(chist) .or. is_posinf(chist)), &
        & 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
    call assert(cweight >= 0, 'CWEIGHT >= 0', srname)
    call assert(ctol >= 0, 'CTOL >= 0', srname)
end if

!====================!
! Calculation starts !
!====================!

! We select X among the points with F < FREF and CSTRV < CREF.
! Do NOT use F <= FREF, because F == FREF (FUNCMAX or REALMAX) may mean F == INF in practice!
if (any(fhist < FUNCMAX .and. chist < CONSTRMAX)) then
    fref = FUNCMAX
    cref = CONSTRMAX
elseif (any(fhist < REALMAX .and. chist < CONSTRMAX)) then
    fref = REALMAX
    cref = CONSTRMAX
elseif (any(fhist < FUNCMAX .and. chist < REALMAX)) then
    fref = FUNCMAX
    cref = REALMAX
else
    fref = REALMAX
    cref = REALMAX
end if

if (.not. any(fhist < fref .and. chist < cref)) then
    kopt = nhist
else
    ! Shift the constraint violations by CTOL, so that CSTRV <= CTOL is regarded as no violation.
    chist_shifted = max(chist - ctol, ZERO)
    ! CMIN is the minimal shifted constraint violation attained in the history.
    cmin = minval(chist_shifted, mask=(fhist < fref))
    ! We consider only the points whose shifted constraint violations are at most the CREF below.
    ! N.B.: Without taking MAX(EPS, .), CREF would be 0 if CMIN = 0. In that case, asking for
    ! CSTRV_SHIFTED < CREF would be WRONG!
    cref = max(EPS, TWO * cmin)
    ! We use the following PHI as our merit function to select X.
    if (cweight <= 0) then
        phi = fhist
    elseif (is_posinf(cweight)) then
        phi = chist_shifted
        ! We should not use CHIST here; if MIN(CHIST_SHIFTED) is attained at multiple indices, then
        ! we will check FHIST to exhaust the remaining degree of freedom.
    else
        phi = max(fhist, -REALMAX) + cweight * chist_shifted
        ! MAX(FHIST, -REALMAX) makes sure that PHI will not contain NaN (unless there is a bug).
    end if
    ! We select X to minimize PHI subject to F < FREF and CSTRV_SHIFTED <= CREF (see the comments
    ! above for the reason of taking "<" and "<=" in these two constraints). In case there are
    ! multiple minimizers, we take the one with the least CSTRV_SHIFTED; if there are more than one
    ! choices, we take the one with the least F; if there are several candidates, we take the one
    ! with the least CSTRV; if the last comparison still leads to more than one possibilities, then
    ! they are equally good and we choose the first.
    ! N.B.:
    ! 1. This process is the opposite of selecting KWORST in SAVEFILT.
    ! 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 == CSTRV_SHIFTED_2 do
    ! not ensure that F_1 == F_2!
    phimin = minval(phi, mask=(fhist < fref .and. chist_shifted <= cref))
    cref = minval(chist_shifted, mask=(fhist < fref .and. phi <= phimin))
    fref = minval(fhist, mask=(chist_shifted <= cref))
    kopt = int(minloc(chist, mask=(fhist <= fref), dim=1), kind(kopt))
    !!MATLAB: cmin = min(chist(fhist <= fref)); kopt = find(fhist <= fref & ~(chist > cmin), 1,'first');
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(kopt >= 1 .and. kopt <= nhist, '1 <= KOPT <= SIZE(FHIST)', srname)
    call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), fhist(kopt), chist(kopt), ctol)), &
        & 'No point in the history is better than X', srname)
end if

end function selectx


function isbetter00(f1, c1, f2, c2, ctol) result(is_better)
!--------------------------------------------------------------------------------------------------!
! This function compares whether FC1 = (F1, C1) is (strictly) better than FC2 = (F2, C2), which
! basically means that (F1 < F2 and C1 <= C2) or (F1 <= F2 and C1 < C2).
! It takes care of the cases where some of these values are NaN or Inf, even though some cases
! should never happen due to the moderated extreme barrier.
! At return, BETTER = TRUE if and only if (F1, C1) is better than (F2, C2).
! Here, C means constraint violation, which is a nonnegative number.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, TEN, EPS, CONSTRMAX, REALMAX, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: f1
real(RP), intent(in) :: c1
real(RP), intent(in) :: f2
real(RP), intent(in) :: c2
real(RP), intent(in) :: ctol

! Outputs
logical :: is_better

! Local variables
character(len=*), parameter :: srname = 'ISBETTER'
real(RP) :: cref

! Preconditions
if (DEBUGGING) then
    call assert(.not. any(is_nan([f1, c1]) .or. is_posinf([f2, c2])), 'FC1 does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan([f2, c2]) .or. is_posinf([f2, c2])), 'FC2 does not contain NaN/+Inf', srname)
    call assert(c1 >= 0 .and. c2 >= 0, 'C1 >= 0, C2 >= 0', srname)
    call assert(ctol >= 0, 'CTOL >= 0', srname)
end if

!====================!
! Calculation starts !
!====================!

is_better = .false.
! Even though NaN/+Inf should not occur in FC1 or FC2 due to the moderated extreme barrier, for
! security and robustness, the code below does not make this assumption.
is_better = is_better .or. (any(is_nan([f1, c1])) .and. .not. any(is_nan([f2, c2])))
is_better = is_better .or. (f1 < f2 .and. c1 <= c2)
is_better = is_better .or. (f1 <= f2 .and. c1 < c2)
! If C1 <= CTOL and C2 is significantly larger/worse than CTOL, i.e., C2 > MAX(CTOL,CREF),
! then FC1 is better than FC2 as long as F1 < REALMAX. Normally CREF >= CTOL so MAX(CTOL, CREF)
! is indeed CREF. However, this may not be true if CTOL > 1E-1*CONSTRMAX.
cref = TEN * max(EPS, min(ctol, 1.0E-2_RP * CONSTRMAX))  ! The MIN avoids overflow.
is_better = is_better .or. (f1 < REALMAX .and. c1 <= ctol .and. (c2 > max(ctol, cref) .or. is_nan(c2)))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(.not. (is_better .and. f1 >= f2 .and. c1 >= c2), &
        & '[F1, C1] >= [F2, C2] and IS_BETTER cannot be both true', srname)
    call assert(is_better .or. .not. (f1 <= f2 .and. c1 < c2), &
        & 'if [F1, C1] <= [F2, C2] but not equal, then IS_BETTER must be true', srname)
    call assert(is_better .or. .not. (f1 < f2 .and. c1 <= c2), &
        & 'if [F1, C1] <= [F2, C2] but not equal, then IS_BETTER must be true', srname)
end if

end function isbetter00


function isbetter10(f1, c1, f2, c2, ctol) result(is_better)
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
real(RP), intent(in) :: f1(:)
real(RP), intent(in) :: c1(:)
real(RP), intent(in) :: f2
real(RP), intent(in) :: c2
real(RP), intent(in) :: ctol

! Outputs
logical, allocatable :: is_better(:)

! Local variables
character(len=*), parameter :: srname = 'ISBETTER10'
integer(IK) :: i
integer(IK) :: nfc

! Sizes
nfc = int(size(f1), kind(nfc))

! Preconditions
if (DEBUGGING) then
    call assert(nfc >= 0, 'NFC >= 0', srname)
    call assert(size(f1) == size(c1), 'SIZE(F1) == SIZE(C1)', srname)
    call assert(.not. any(is_nan(f1) .or. is_posinf(f1)), 'F1 does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(c1) .or. is_posinf(c1)), 'C1 does not contain NaN/+Inf', srname)
    call assert(.not. (is_nan(f2) .or. is_posinf(f2)), 'F2 is not NaN/+Inf', srname)
    call assert(.not. (is_nan(c2) .or. is_posinf(c2)), 'C2 is not NaN/+Inf', srname)
    call assert(all(c1 >= 0) .and. c2 >= 0, 'C1 >= 0, C2 >= 0', srname)
    call assert(ctol >= 0, 'CTOL >= 0', srname)
end if

!====================!
! Calculation starts !
!====================!

call safealloc(is_better, nfc)
is_better = [(isbetter00(f1(i), c1(i), f2, c2, ctol), i=1, nfc)]

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(is_better) == size(f1), 'SIZE(IS_BETTER) == SIZE(F1)', srname)
end if

end function isbetter10


function isbetter01(f1, c1, f2, c2, ctol) result(is_better)
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
real(RP), intent(in) :: f1
real(RP), intent(in) :: c1
real(RP), intent(in) :: f2(:)
real(RP), intent(in) :: c2(:)
real(RP), intent(in) :: ctol

! Outputs
logical, allocatable :: is_better(:)

! Local variables
character(len=*), parameter :: srname = 'ISBETTER01'
integer(IK) :: i
integer(IK) :: nfc

! Sizes
nfc = int(size(f2), kind(nfc))

! Preconditions
if (DEBUGGING) then
    call assert(nfc >= 0, 'NFC >= 0', srname)
    call assert(.not. (is_nan(f1) .or. is_posinf(f1)), 'F1 is not NaN/+Inf', srname)
    call assert(.not. (is_nan(c1) .or. is_posinf(c1)), 'C1 is not NaN/+Inf', srname)
    call assert(size(f2) == size(c2), 'SIZE(F2) == SIZE(C2)', srname)
    call assert(.not. any(is_nan(f2) .or. is_posinf(f2)), 'F2 does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(c2) .or. is_posinf(c2)), 'C2 does not contain NaN/+Inf', srname)
    call assert(c1 >= 0 .and. all(c2 >= 0), 'C1 >= 0, C2 >= 0', srname)
    call assert(ctol >= 0, 'CTOL >= 0', srname)
end if

!====================!
! Calculation starts !
!====================!

call safealloc(is_better, nfc)
is_better = [(isbetter00(f1, c1, f2(i), c2(i), ctol), i=1, nfc)]

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(is_better) == size(f2), 'SIZE(IS_BETTER) == SIZE(F2)', srname)
end if

end function isbetter01


end module selectx_mod
