module hist_mod

implicit none
private
public :: savehist, selectx

contains


subroutine savehist(x, f, con, cstrv, xsav, fsav, consav, csav, nsav, ctol)   !?????
! This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless a vector in XSAV(:, 1:NSAV) is
! better than XDROP. If XDROP is better than some vectors in XSAV(:, 1:NSAV), then these vectors
! will be removed. If XDROP is not better than any of XSAV(:, 1:NSAV) but NSAV=NSAVMAX, then we
! remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NSAV).
!
! When COBYLA calls this subroutine, XDROP is a vector to be "dropped", and  DATDROP contains its
! function/constraint information (constraint value in the first M entries, DATDROP(M+1) = F(XDROP),
! and DATDROP(M+2) = CSTRV(X)). XSAV and DATSAV save at most NSAVMAX vectors "dropped" by COBYLB
! and their function/constraint information. Only XSAV(:, 1:NSAV) and DATSAV(:, 1:NSAV) contains
! such vectors, while XSAV(:, NSAV+1:NSAVMAX) and DATSAV(:, NSAV+1:NSAVMAX) are not initialized yet.
!
! Note: We decide whether X is better than the function/constraint of Y according to the ISBETTER
! function with CPEN = -ONE. Due to the implementation of ISBETTER,
! X is better than Y with CPEN < 0
! ==> X is better than Y with any CPEN >= 0,
! ==> X is better than Y regardless of CPEN.

use consts_mod, only : RP, IK, ONE, DEBUGGING, SRNLEN
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop, verisize
use output_mod, only : retmssg, rhomssg, fmssg
use memory_mod, only : safealloc
use lina_mod, only : calquad, inprod

implicit none

! Inputs
real(RP), intent(IN) :: ctol
real(RP), intent(IN) :: con(:)  ! M
real(RP), intent(IN) :: cstrv
real(RP), intent(IN) :: f
real(RP), intent(IN) :: x(:)  ! N

! In-outputs
integer(IK), intent(INOUT) :: nsav
real(RP), intent(INOUT) :: consav(:, :)  ! (M, NSAVMAX)
real(RP), intent(INOUT) :: csav(:)  ! NSAVMAX
real(RP), intent(INOUT) :: fsav(:)  ! NSAVMAX
real(RP), intent(INOUT) :: xsav(:, :) ! (N, NSAVMAX)

! Local variables
integer(IK) :: m
integer(IK) :: n
integer(IK) :: nsavmax
integer(IK) :: i
integer(IK), allocatable :: index_to_keep(:)
logical :: better(nsav)
logical :: keep(nsav)
character(len=SRNLEN), parameter :: srname = 'SAVEHIST'

m = size(con)
n = size(x)
nsavmax = size(fsav)

if (DEBUGGING) then
    if (n <= 0) then
        call errstop(srname, 'SIZE(X) is invalid')
    end if
    call verisize(xsav, n, nsavmax)
    call verisize(consav, m, nsavmax)
    call verisize(csav, nsavmax)
end if

if (nsavmax <= 0) then
    return  ! Do nothing if NSAVMAX = 0
end if

! Return immediately if any column of XSAV is better than XDROP.
! BETTER is defined by the array constructor with an implicit do loop.
better = [(isbetter([fsav(i), csav(i)], [f, cstrv], ctol), i=1, nsav)]
if (any(better)) then
    return
end if

! Decide which columns of XSAV to keep. We use again the array constructor with an implicit do loop.
keep = [(.not. isbetter([f, cstrv], [fsav(i), csav(i)], ctol), i=1, nsav)]
! If X is not better than any column of XSAV, then we remove the first (oldest) column of XSAV.
if (count(keep) == nsavmax) then
    keep(1) = .false.
end if

!--------------------------------------------------!
!----The SAFEALLOC line is removable in F2003.-----!
call safealloc(index_to_keep, nsav)
!--------------------------------------------------!
index_to_keep = pack([(i, i=1, nsav)], mask=keep)
nsav = count(keep) + 1
xsav(:, 1:nsav - 1) = xsav(:, index_to_keep)
fsav(1:nsav - 1) = fsav(index_to_keep)
consav(:, 1:nsav - 1) = consav(:, index_to_keep)
csav(1:nsav - 1) = csav(index_to_keep)
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet
! we prefer to deallocate them immediately when they finish their jobs.
deallocate (index_to_keep)
xsav(:, nsav) = x
fsav(nsav) = f
consav(:, nsav) = con
csav(nsav) = cstrv

end subroutine savehist


function selectx(fhist, chist, cpen, ctol) result(kopt)

use consts_mod, only : IK, RP, HUGENUM, HUGEFUN, HUGECON, ZERO, TWO, DEBUGGING, SRNLEN
use infnan_mod, only : is_nan, is_posinf
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
character(len=SRNLEN), parameter :: srname = 'SELECTX'


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
use consts_mod, only : RP, TEN, HUGENUM, DEBUGGING, SRNLEN
use infnan_mod, only : is_nan
use debug_mod, only : errstop
implicit none

! Inputs
real(RP), intent(IN) :: fc1(:)
real(RP), intent(IN) :: fc2(:)
real(RP), intent(IN) :: ctol

! Output
real(RP) :: cref
logical :: is_better

! Local variables
character(len=SRNLEN), parameter :: srname = 'ISBETTER'


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


end module hist_mod
