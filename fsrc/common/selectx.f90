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
! Last Modified: Thursday, December 02, 2021 PM02:05:42
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: savefilt, selectx, isbetter


contains


subroutine savefilt(constr, cstrv, ctol, f, x, nfilt, cfilt, confilt, ffilt, xfilt)
!--------------------------------------------------------------------------------------------------!
! This subroutine saves X, F, CONSTR, and CSTRV in XFILT, FFILT, CONFILT, and CFILT, unless a vector
! in XFILT(:, 1:NFILT) is better than X. If X is better than some vectors in XFILT(:, 1:NFILT), then
! these vectors will be removed. If X is not better than any of XFILT(:, 1:NFILT) but NFILT=NFILTMAX,
! then we remove XFILT(:,1), which is the oldest vector in XFILT(:, 1:NFILT).
! N.B.:
! 1. Only XFILT(:, 1:NFILT) and FFILT(:, 1:NFILT) etc contains valid information, while
! XFILT(:, NFILT+1:NFILTMAX) and FFILT(:, NFILT+1:NFILTMAX) etc are not initialized yet.
! 2. We decide whether an X is better than another by the ISBETTER function.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: memory_mod, only : safealloc

implicit none

! Inputs
real(RP), intent(in) :: constr(:)  ! M
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: ctol
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)  ! N

! In-outputs
integer(IK), intent(inout) :: nfilt
real(RP), intent(inout) :: cfilt(:)  ! NFILTMAX
real(RP), intent(inout) :: confilt(:, :)  ! (M, NFILTMAX)
real(RP), intent(inout) :: ffilt(:)  ! NFILTMAX
real(RP), intent(inout) :: xfilt(:, :) ! (N, NFILTMAX)

! Local variables
character(len=*), parameter :: srname = 'SAVEFILT'
integer(IK) :: i
integer(IK) :: m
integer(IK) :: n
integer(IK) :: nfiltmax
integer(IK), allocatable :: index_to_keep(:)
logical :: better(nfilt)
logical :: keep(nfilt)

! Sizes
m = int(size(constr), kind(m))
n = int(size(x), kind(n))
nfiltmax = int(size(ffilt), kind(nfiltmax))

! Preconditions
if (DEBUGGING) then
    ! Check the size of X.
    call assert(n >= 1, 'N >= 1', srname)
    ! Check NFILT
    call assert(nfilt <= nfiltmax, 'NFILT <= NFILTMAX', srname)
    ! Check the sizes of XFILT, FFILT, CONFILT, CFILT.
    call assert(nfiltmax >= 1, 'NFILTMAX >= 1', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == nfiltmax, 'SIZE(XFILT) == [N, NFILTMAX]', srname)
    call assert(size(confilt, 1) == m .and. size(confilt, 2) == nfiltmax, 'SIZE(CONFILT) == [M, NFILTMAX]', srname)
    call assert(size(cfilt) == nfiltmax, 'SIZE(CFILT) == NFILTMAX', srname)
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
! BETTER is defined by the array constructor with an implied do loop.
better = [(isbetter([ffilt(i), cfilt(i)], [f, cstrv], ctol), i=1, nfilt)]
if (any(better)) then
    return
end if

! Decide which columns of XFILT to keep. We use again the array constructor with an implied do loop.
keep = [(.not. isbetter([f, cstrv], [ffilt(i), cfilt(i)], ctol), i=1, nfilt)]
! If X is not better than any column of XFILT, then we remove the first (oldest) column of XFILT.
if (count(keep) == nfiltmax) then
    keep(1) = .false.
end if

!--------------------------------------------------!
!----The SAFEALLOC line is removable in F2003.-----!
call safealloc(index_to_keep, nfilt)
!--------------------------------------------------!
index_to_keep = pack([(int(i, IK), i=1, nfilt)], mask=keep)
nfilt = int(count(keep) + 1, kind(nfilt))
xfilt(:, 1:nfilt - 1) = xfilt(:, index_to_keep)
ffilt(1:nfilt - 1) = ffilt(index_to_keep)
confilt(:, 1:nfilt - 1) = confilt(:, index_to_keep)
cfilt(1:nfilt - 1) = cfilt(index_to_keep)
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer to deallocate
! them immediately when they finish their jobs.
deallocate (index_to_keep)
xfilt(:, nfilt) = x
ffilt(nfilt) = f
confilt(:, nfilt) = constr
cfilt(nfilt) = cstrv

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    ! Check NFILT and the sizes of XFILT, FFILT, CONFILT, CFILT
    call assert(nfilt >= 0 .and. nfilt <= nfiltmax, '0 <= NFILT <= NFILTMAX', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == nfiltmax, 'SIZE(XFILT) == [N, NFILTMAX]', srname)
    call assert(size(ffilt) == nfiltmax, 'SIZE(FFILT) = NFILTMAX', srname)
    call assert(size(confilt, 1) == m .and. size(confilt, 2) == nfiltmax, 'SIZE(CONFILT) == [M, NFILTMAX]', srname)
    call assert(size(cfilt) == nfiltmax, 'SIZE(CFILT) = NFILTMAX', srname)
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


function selectx(fhist, chist, cpen, ctol) result(kopt)
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
real(RP), intent(in) :: cpen
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
    call assert(cpen >= 0, 'CPEN >= 0', srname)
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
    ! We select X among the points whose shifted constraint violations are at most CREF.
    cref = TWO * cmin  ! CREF = 0 if CMIN = 0; thus asking for CSTRV_SHIFTED < CREF is unreasonable.
    ! We use the following PHI as our merit function to select X.
    phi = fhist + cpen * chist_shifted
    ! We select X to minimize PHI subject to F < FREF and CSTRV_SHIFTED <= CREF. In case there are
    ! multiple minimizers, we take the one with the least CSTRV; if there are still more than one
    ! choices, we take the one with the least F; if the last comparison still leads to more than
    ! one possibilities, then they are all equally good, and we choose the first one of them.
    ! N.B.: In finite-precision arithmetic, PHI_1 = PHI_2 and CSTRV_SHIFTED_1 = CSTRV_SHIFTED_2 do
    ! not ensure that F_1 = F_2!!!
    phimin = minval(phi, mask=(fhist < fref .and. chist_shifted <= cref))
    cref = minval(chist, mask=(fhist < fref .and. phi <= phimin))
    kopt = int(minloc(fhist, mask=(chist <= cref), dim=1), kind(kopt))
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

use, non_intrinsic :: consts_mod, only : RP, TEN, HUGENUM, DEBUGGING
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
is_better = is_better .or. (.not. any(is_nan(fc1)) .and. any(is_nan(fc2)))
is_better = is_better .or. (fc1(1) < fc2(1) .and. fc1(2) <= fc2(2))
is_better = is_better .or. (fc1(1) <= fc2(1) .and. fc1(2) < fc2(2))
cref = TEN * max(ctol, epsilon(ctol))
is_better = is_better .or. (fc1(1) < HUGENUM .and. fc1(2) <= ctol .and. &
    & ((fc2(2) > cref) .or. is_nan(fc2(2))))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(.not. (fc1(1) >= fc2(1) .and. fc1(2) >= fc2(2) .and. is_better), &
        & 'FC1 >= FC2 and IS_BETTER cannot be both true', srname)
    call assert((.not. (fc1(1) <= fc2(1) .and. fc1(2) < fc2(2))) .or. is_better, &
        & 'if FC1 <= FC2 but not equal, then IS_BETTER must be true', srname)
    call assert((.not. (fc1(1) < fc2(1) .and. fc1(2) <= fc2(2))) .or. is_better, &
        & 'if FC1 <= FC2 but not equal, then IS_BETTER must be true', srname)
end if

end function isbetter


end module selectx_mod
