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
! Last Modified: Monday, March 14, 2022 PM12:38:26
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: savefilt, selectx, isbetter


contains


subroutine savefilt(constr, cstrv, ctol, cweight, f, x, nfilt, cfilt, confilt, ffilt, xfilt)
!--------------------------------------------------------------------------------------------------!
! This subroutine saves X, F, CONSTR, and CSTRV in XFILT, FFILT, CONFILT, and CFILT, unless a vector
! in XFILT(:, 1:NFILT) is better than X. If X is better than some vectors in XFILT(:, 1:NFILT), then
! these vectors will be removed. If X is not better than any of XFILT(:, 1:NFILT) but NFILT=MAXFILT,
! then we remove a column from XFILT according to PHI = FFILT + CWEIGHT * MAX(CFILT - CTOL, ZERO).
! N.B.:
! 1. Only XFILT(:, 1:NFILT) and FFILT(:, 1:NFILT) etc contains valid information, while
! XFILT(:, NFILT+1:MAXFILT) and FFILT(:, NFILT+1:MAXFILT) etc are not initialized yet.
! 2. We decide whether an X is better than another by the ISBETTER function.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: linalg_mod, only : trueloc

implicit none

! Inputs
real(RP), intent(in) :: constr(:)  ! M
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: ctol
real(RP), intent(in) :: cweight
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)  ! N

! In-outputs
integer(IK), intent(inout) :: nfilt
real(RP), intent(inout) :: cfilt(:)  ! MAXFILT
real(RP), intent(inout) :: confilt(:, :)  ! (M, MAXFILT)
real(RP), intent(inout) :: ffilt(:)  ! MAXFILT
real(RP), intent(inout) :: xfilt(:, :) ! (N, MAXFILT)

! Local variables
character(len=*), parameter :: srname = 'SAVEFILT'
integer(IK) :: i
integer(IK) :: index_to_keep(size(ffilt))
integer(IK) :: m
integer(IK) :: maxfilt
integer(IK) :: n
integer(IK) :: kworst
logical :: better(nfilt)
logical :: keep(nfilt)
real(RP) :: cfilt_shifted(size(ffilt))
real(RP) :: cref
real(RP) :: fref
real(RP) :: phi(size(ffilt))
real(RP) :: phimax

! Sizes
m = int(size(constr), kind(m))
n = int(size(x), kind(n))
maxfilt = int(size(ffilt), kind(maxfilt))

! Preconditions
if (DEBUGGING) then
    ! Check the size of X.
    call assert(n >= 1, 'N >= 1', srname)
    ! Check NFILT
    call assert(nfilt >= 0 .and. nfilt <= maxfilt, '0 <= NFILT <= MAXFILT', srname)
    ! Check the sizes of XFILT, FFILT, CONFILT, CFILT.
    call assert(maxfilt >= 1, 'MAXFILT >= 1', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == maxfilt, 'SIZE(XFILT) == [N, MAXFILT]', srname)
    call assert(size(confilt, 1) == m .and. size(confilt, 2) == maxfilt, 'SIZE(CONFILT) == [M, MAXFILT]', srname)
    call assert(size(cfilt) == maxfilt, 'SIZE(CFILT) == MAXFILT', srname)
    ! Check the values of XFILT, FFILT, CONFILT, CFILT.
    call assert(.not. any(is_nan(xfilt(:, 1:nfilt))), 'XFILT does not contain NaN', srname)
    call assert(.not. any(is_nan(ffilt(1:nfilt)) .or. is_posinf(ffilt(1:nfilt))), &
        & 'FFILT does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(confilt(:, 1:nfilt)) .or. is_neginf(confilt(:, 1:nfilt))), &
        & 'CONFILT does not contain NaN/-Inf', srname)
    call assert(.not. any(cfilt(1:nfilt) < 0 .or. is_nan(cfilt(1:nfilt)) .or. is_posinf(cfilt(1:nfilt))), &
        & 'CFILT does not contain nonnegative values of NaN/+Inf', srname)
    ! Check the values of X, F, CONSTR, CSTRV.
    ! X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! F cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    ! CONSTR cannot contain NaN/-Inf due to the moderated extreme barrier.
    call assert(.not. any(is_nan(constr) .or. is_neginf(constr)), 'CONSTR does not contain NaN/-Inf', srname)
    ! CSTRV cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (cstrv < 0 .or. is_nan(cstrv) .or. is_posinf(cstrv)), 'CSTRV is nonnegative and not NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! Return immediately if any column of XFILT is better than X.
! BETTER is defined by an array constructor with an implied do loop.
better = [(isbetter([ffilt(i), cfilt(i)], [f, cstrv], ctol), i=1, nfilt)]
if (any(better)) then
    return
end if

! Decide which columns of XFILT to keep. We use again an array constructor with an implied do loop.
keep = [(.not. isbetter([f, cstrv], [ffilt(i), cfilt(i)], ctol), i=1, nfilt)]

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
        phi = max(ffilt, -HUGENUM) + cweight * cfilt_shifted
        ! MAX(FFILT, -HUGENUM) makes sure that PHI will not contain NaN (unless there is a bug).
    end if
    ! We select X to maximize PHI. In case there are multiple maximizers, we take the one with the
    ! largest CSTRV_SHIFTED; if there are more than one choices, we take the one with the largest F;
    ! if there are several candidates, we take the one with the largest CSTRV; if the last comparison
    ! still leads to more than one possibilities, then they are equally bad and we choose the first.
    ! N.B.:
    ! 1. This process is the opposite of selecting KOPT in SELECTX.
    ! 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 == CSTRV_SHIFTED_2 do
    ! not ensure that F_1 == F_2!!!
    phimax = maxval(phi)
    cref = maxval(cfilt_shifted, mask=(phi >= phimax))
    fref = maxval(ffilt, mask=(cfilt_shifted >= cref))
    kworst = int(maxloc(cfilt, mask=(ffilt >= fref), dim=1), kind(kworst))
    if (kworst < 1 .or. kworst > size(keep)) then  ! For security. Should not happen.
        kworst = 1_IK
    end if
    keep(kworst) = .false.
end if

nfilt = int(count(keep), kind(nfilt))
index_to_keep(1:nfilt) = trueloc(keep)
xfilt(:, 1:nfilt) = xfilt(:, index_to_keep(1:nfilt))
ffilt(1:nfilt) = ffilt(index_to_keep(1:nfilt))
confilt(:, 1:nfilt) = confilt(:, index_to_keep(1:nfilt))
cfilt(1:nfilt) = cfilt(index_to_keep(1:nfilt))

nfilt = nfilt + 1_IK
xfilt(:, nfilt) = x
ffilt(nfilt) = f
confilt(:, nfilt) = constr
cfilt(nfilt) = cstrv

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    ! Check NFILT and the sizes of XFILT, FFILT, CONFILT, CFILT.
    call assert(nfilt >= 1 .and. nfilt <= maxfilt, '1 <= NFILT <= MAXFILT', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == maxfilt, 'SIZE(XFILT) == [N, MAXFILT]', srname)
    call assert(size(ffilt) == maxfilt, 'SIZE(FFILT) = MAXFILT', srname)
    call assert(size(confilt, 1) == m .and. size(confilt, 2) == maxfilt, 'SIZE(CONFILT) == [M, MAXFILT]', srname)
    call assert(size(cfilt) == maxfilt, 'SIZE(CFILT) = MAXFILT', srname)
    ! Check the values of XFILT, FFILT, CONFILT, CFILT.
    call assert(.not. any(is_nan(xfilt(:, 1:nfilt))), 'XFILT does not contain NaN', srname)
    call assert(.not. any(is_nan(ffilt(1:nfilt)) .or. is_posinf(ffilt(1:nfilt))), &
        & 'FFILT does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(confilt(:, 1:nfilt)) .or. is_neginf(confilt(:, 1:nfilt))), &
        & 'CONFILT does not contain NaN/-Inf', srname)
    call assert(.not. any(cfilt(1:nfilt) < 0 .or. is_nan(cfilt(1:nfilt)) .or. is_posinf(cfilt(1:nfilt))), &
        & 'CFILT does not contain nonnegative values of NaN/+Inf', srname)
    ! Check that no point in the filter is better than X, and X is better than no point.
    call assert(.not. any([(isbetter([ffilt(i), cfilt(i)], [f, cstrv], ctol), i=1, nfilt)]), &
        & 'No point in the filter is better than X', srname)
    call assert(.not. any([(isbetter([f, cstrv], [ffilt(i), cfilt(i)], ctol), i=1, nfilt)]), &
        & 'X is better than no point in the filter', srname)
end if

end subroutine savefilt


function selectx(fhist, chist, cweight, ctol) result(kopt)
!--------------------------------------------------------------------------------------------------!
! This subroutine selects X according to the FHIST and CHIST, which represents (a piece of) history
! of F and CSTRV. Normally, FHIST and CHIST are not the full history but only a filter, e.g., FFILT
! and CFILT generated by SAVEFILT. However, we name them as FHIST and CHIST because the [F, CSTRV]
! in a filter should not dominate each other, but this subroutine does NOT assume such a structure.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : IK, RP, HUGENUM, HUGEFUN, HUGECON, ZERO, TWO, DEBUGGING
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
integer(IK) :: i
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
!--------------------------------------------------------------------------------------------------!
!! Do NOT use F <= FREF, because F == FREF (HUGEFUN or HUGENUM) may mean F == INF in practice !!
!--------------------------------------------------------------------------------------------------!
if (any(fhist < HUGEFUN .and. chist < HUGECON)) then
    fref = HUGEFUN
    cref = HUGECON
elseif (any(fhist < HUGENUM .and. chist < HUGECON)) then
    fref = HUGENUM
    cref = HUGECON
elseif (any(fhist < HUGEFUN .and. chist < HUGENUM)) then
    fref = HUGEFUN
    cref = HUGENUM
else
    fref = HUGENUM
    cref = HUGENUM
end if

if (.not. any(fhist < fref .and. chist < cref)) then
    kopt = nhist
else
    ! Shift the constraint violations by CTOL, so that CSTRV <= CTOL is regarded as no violation.
    chist_shifted = max(chist - ctol, ZERO)
    ! CMIN is the minimal shifted constraint violation attained in the history.
    cmin = minval(chist_shifted, mask=(fhist < fref))
    ! We consider only the points whose shifted constraint violations are at most the CREF below.
    cref = TWO * cmin  ! CREF = 0 if CMIN = 0; thus asking for CSTRV_SHIFTED < CREF is WRONG!
    ! We use the following PHI as our merit function to select X.
    if (cweight <= 0) then
        phi = fhist
    elseif (is_posinf(cweight)) then
        phi = chist_shifted
        ! We should not use CHIST here; if MIN(CHIST_SHIFTED) is attained at multiple indices, then
        ! we will check FHIST to exhaust the remaining degree of freedom.
    else
        phi = max(fhist, -HUGENUM) + cweight * chist_shifted
        ! MAX(FHIST, -HUGENUM) makes sure that PHI will not contain NaN (unless there is a bug).
    end if
    ! We select X to minimize PHI subject to F < FREF and CSTRV_SHIFTED <= CREF (seee the comments
    ! above for the reason of taking "<" and "<=" in these two constraints). In case there are
    ! multiple minimizers, we take the one with the least CSTRV_SHIFTED; if there are more than one
    ! choices, we take the one with the least F; if there are several candidates, we take the one
    ! with the least CSTRV; if the last comparison still leads to more than one possibilities, then
    ! they are equally good and we choose the first.
    ! N.B.:
    ! 1. This process is the opposite of selecting KWORST in SAVEFILT.
    ! 2. In finite-precision arithmetic, PHI_1 == PHI_2 and CSTRV_SHIFTED_1 == CSTRV_SHIFTED_2 do
    ! not ensure that F_1 == F_2!!!
    phimin = minval(phi, mask=(fhist < fref .and. chist_shifted <= cref))
    cref = minval(chist_shifted, mask=(fhist < fref .and. phi <= phimin))
    fref = minval(fhist, mask=(chist_shifted <= cref))
    kopt = int(minloc(chist, mask=(fhist <= fref), dim=1), kind(kopt))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(kopt >= 1 .and. kopt <= nhist, '1 <= KOPT <= SIZE(FHIST)', srname)
    call assert(.not. any([(isbetter([fhist(i), chist(i)], [fhist(kopt), chist(kopt)], ctol), &
        & i=1, nhist)]), 'No point in the history is better than X', srname)
end if

end function selectx


function isbetter(fc1, fc2, ctol) result(is_better)
!--------------------------------------------------------------------------------------------------!
! This function compares whether FC1 = (F1, CSTRV1) is (strictly) better than FC2 = (F2, CSTRV2),
! meaning that (F1 < F2 and CSTRV1 <= CSTRV2) or (F1 <= F2 and CSTRV1 < CSTRV2).
! It takes care of the cases where some of these values are NaN or Inf, even though some cases
! should never happen due to the moderated extreme barrier.
! At return, BETTER = TRUE if and only if (F1, CSTRV1) is better than (F2, CSTRV2).
! Here, CSTRV means constraint violation, which is a nonnegative number.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, TEN, EPS, HUGECON, HUGENUM, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: fc1(:)
real(RP), intent(in) :: fc2(:)
real(RP), intent(in) :: ctol

! Outputs
logical :: is_better

! Local variables
character(len=*), parameter :: srname = 'ISBETTER'
real(RP) :: cref

! Preconditions
if (DEBUGGING) then
    call assert(size(fc1) == 2 .and. size(fc2) == 2, 'SIZE(FC1) == 2 == SIZE(FC2)', srname)
    call assert(.not. any(is_nan(fc1) .or. is_posinf(fc1)), 'FC1 does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(fc2) .or. is_posinf(fc2)), 'FC2 does not contain NaN/+Inf', srname)
    call assert(fc1(2) >= 0 .and. fc2(2) >= 0, 'FC1(2) >= 0, FC(2) >= 0', srname)
    call assert(ctol >= 0, 'CTOL >= 0', srname)
end if

!====================!
! Calculation starts !
!====================!

is_better = .false.
! Even though NaN/+Inf should not occur in FC1 or FC2 due to the moderated extreme barrier, for
! security and robustness, the code below does not make this assumption.
is_better = is_better .or. (any(is_nan(fc2)) .and. .not. any(is_nan(fc2)))
is_better = is_better .or. (fc1(1) < fc2(1) .and. fc1(2) <= fc2(2))
is_better = is_better .or. (fc1(1) <= fc2(1) .and. fc1(2) < fc2(2))
! If FC1(2) <= CTOL and FC2(2) is significantly larger/worse than CTOL, i.e., FC(2) > MAX(CTOL,CREF),
! then FC1 is better than FC2 as long as FC(1) < HUGENUM. Normally CREF >= CTOL so MAX(CTOL, CREF)
! is indeed CREF. However, this may not be true if CTOL > 1E-1*HUGECON.
cref = TEN * max(EPS, min(ctol, 1.0E-2_RP * HUGECON))  ! The MIN avoids overflow.
is_better = is_better .or. (fc1(1) < HUGENUM .and. fc1(2) <= ctol .and. &
    & (fc2(2) > max(ctol, cref) .or. is_nan(fc2(2))))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(.not. (is_better .and. fc1(1) >= fc2(1) .and. fc1(2) >= fc2(2)), &
        & 'FC1 >= FC2 and IS_BETTER cannot be both true', srname)
    call assert(is_better .or. .not. (fc1(1) <= fc2(1) .and. fc1(2) < fc2(2)), &
        & 'if FC1 <= FC2 but not equal, then IS_BETTER must be true', srname)
    call assert(is_better .or. .not. (fc1(1) < fc2(1) .and. fc1(2) <= fc2(2)), &
        & 'if FC1 <= FC2 but not equal, then IS_BETTER must be true', srname)
end if

end function isbetter


end module selectx_mod
