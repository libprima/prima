module selectx_mod

implicit none
private
public :: selectx, isbetter

contains


function selectx(fhist, chist, cpen, ctol) result(kopt)
! This subroutine selects X according to the history of F and CSTRV. Normally, FHIST and CHIST are
! not the complete history but only a filter.

use consts_mod, only : IK, RP, HUGENUM, HUGEFUN, HUGECON, ZERO, TWO, DEBUGGING
use debug_mod, only : errstop, verisize

implicit none

! Inputs
real(RP), intent(in) :: cpen
real(RP), intent(in) :: chist(:)
real(RP), intent(in) :: ctol
real(RP), intent(in) :: fhist(:)

! Output
integer(IK) :: kopt

! Local variables
real(RP) :: cmin
real(RP) :: cref
real(RP) :: chist_shifted(size(fhist))
real(RP) :: fref
real(RP) :: phi(size(fhist))
real(RP) :: phimin
character(len=*), parameter :: srname = 'SELECTX'


! Get the verify the sizes.
if (DEBUGGING) then
    if (size(fhist) == 0) then
        call errstop(srname, 'SIZE(FHIST) == 0')
    end if
    call verisize(chist, size(fhist))
end if

! We select X among the points with F < FREF and CSTRV < CREF.
! Do NOT use F <= FREF, because F = FREF (HUGEFUN or HUGENUM) may mean F = INF in practice.
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
    kopt = size(fhist)
else
    ! Shift the constraint violations by CTOL, so that CSTRV <= CTOL is regarded as no violation.
    chist_shifted = max(chist - ctol, ZERO)
    ! CMIN is the minimal shifted constraint violation attained in the history.
    cmin = minval(chist_shifted, mask=(fhist < fref))
    ! We select X among the points whose shifted constraint violations are at most CREF.
    cref = TWO * cmin  ! CREF = ZERO if CMIN = ZERO; asking for CSTRV_SHIFTED < CREF is unreasonable.
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

end function selectx


function isbetter(fc1, fc2, ctol) result(is_better)
! This function compares whether FC1 = (F1, CSTRV1) is (strictly) better than FC2 = (F2, CSTRV2),
! meaning that (F1 < F2 and CSTRV1 <= CSTRV2) or (F1 <= F2 and CSTRV1 < CSTRV2).
! It takes care of the cases where some of these values are NaN or Inf.
! At return, BETTER = TRUE if and only if (F1, CSTRV1) is better than (F2, CSTRV2).
! Here, CSTRV means constraint violation, which is a nonnegative number.

! Generic modules
use consts_mod, only : RP, TEN, HUGENUM, DEBUGGING
use infnan_mod, only : is_nan
use debug_mod, only : errstop
implicit none

! Inputs
real(RP), intent(in) :: fc1(:)
real(RP), intent(in) :: fc2(:)
real(RP), intent(in) :: ctol

! Output
real(RP) :: cref
logical :: is_better

! Local variables
character(len=*), parameter :: srname = 'ISBETTER'


! Verify the sizes
if (DEBUGGING .and. (size(fc1) /= 2 .or. size(fc2) /= 2)) then
    call errstop(srname, 'SIZE(FC1) or SIZE(FC2) is not 2')
end if

is_better = .false.
is_better = is_better .or. (.not. any(is_nan(fc1)) .and. any(is_nan(fc2)))
is_better = is_better .or. (fc1(1) < fc2(1) .and. fc1(2) <= fc2(2))
is_better = is_better .or. (fc1(1) <= fc2(1) .and. fc1(2) < fc2(2))
cref = TEN * max(ctol, epsilon(ctol))
is_better = is_better .or. (fc1(1) < HUGENUM .and. fc1(2) <= ctol .and. ((fc2(2) > cref) .or. is_nan(fc2(2))))

end function isbetter


end module selectx_mod
