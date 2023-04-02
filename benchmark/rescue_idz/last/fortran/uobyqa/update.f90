module update_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Sunday, March 05, 2023 PM09:58:57
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(knew, d, f, moderr, kopt, fval, pl, pq, xpt)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates PL, PQ, XPT, KOPT, and FVAL when XPT(:, KNEW) becomes XNEW.
! See Section 4 of the UOBYQA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : outprod, norm
use, non_intrinsic :: powalg_mod, only : calvlag

! Inputs
integer(IK), intent(in) :: knew
real(RP), intent(in) :: d(:)  ! D(N)
real(RP), intent(in) :: f
real(RP), intent(in) :: moderr

! In-outputs
integer(IK), intent(inout) :: kopt
real(RP), intent(inout) :: fval(:)   ! FVAL(NPT)
real(RP), intent(inout) :: pl(:, :)  ! PL(NPT-1, NPT)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT-1)
real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)

! Local variables
character(len=*), parameter :: srname = 'UPDATE'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: plnew(size(pl, 1))
real(RP) :: vlag(size(xpt, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(npt == (n + 1) * (n + 2) / 2, 'NPT = (N+1)(N+2)/2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew >= 1 .or. f >= fval(kopt), 'KNEW >= 1 unless X is not improved', srname)
    call assert(knew /= kopt .or. f < fval(kopt), 'KNEW /= KOPT unless X is improved', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN or +Inf', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(pl, 1) == npt - 1 .and. size(pl, 2) == npt, 'SIZE(PL) == [NPT-1, NPT]', srname)
    call assert(size(pq) == npt - 1, 'SIZE(PQ) == NPT-1', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do essentially nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    return
end if

! Update the Lagrange functions.
vlag = calvlag(pl, d, xpt(:, kopt), kopt)
pl(:, knew) = pl(:, knew) / vlag(knew)
plnew = pl(:, knew)
pl = pl - outprod(plnew, vlag)
pl(:, knew) = plnew

! Update the quadratic model.
pq = pq + moderr * plnew

! Replace the interpolation point that has index KNEW by the point XNEW.
xpt(:, knew) = xpt(:, kopt) + d
fval(knew) = f

! KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW) = FVAL(KOPT) and KNEW < KOPT, then
! MINLOC(FVAL) = KNEW /= KOPT. Do not change KOPT in this case.
if (f < fval(kopt)) then
    kopt = knew
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt .and. all(is_finite(xpt)), &
        & 'SIZE(XPT) == [N, NPT], XPT is finite', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(pl, 1) == npt - 1 .and. size(pl, 2) == npt, 'SIZE(PL) == [NPT-1, NPT]', srname)
    call assert(size(pq) == npt - 1, 'SIZE(PQ) == NPT-1', srname)
end if

end subroutine update


end module update_mod
