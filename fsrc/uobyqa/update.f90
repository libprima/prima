module update_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the update of PL, PQ, XPT, KOPT, XOPT, and FOPT
! when XPT(:, KNEW) becomes XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Monday, November 28, 2022 PM05:48:15
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(knew, d, f, moderr, kopt, fopt, pl, pq, xopt, xpt)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates PL, PQ, XPT, KOPT, XOPT, and FOPT when XPT(:, KNEW) becomes XNEW = XOPT+D.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : outprod
use, non_intrinsic :: powalg_mod, only : calvlag

! Inputs
integer(IK), intent(in) :: knew
real(RP), intent(in) :: d(:)  ! D(N)
real(RP), intent(in) :: f
real(RP), intent(in) :: moderr

! In-outputs
integer(IK), intent(inout) :: kopt
real(RP), intent(inout) :: fopt
real(RP), intent(inout) :: pl(:, :)  ! PL(NPT-1, NPT)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT-1)
real(RP), intent(inout) :: xopt(:)  ! XOPT(N)
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
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

! Update the Lagrange functions and the quadratic model.
! It can happen that VLAG(KNEW) = 0 due to rounding.
vlag = calvlag(pl, d, xopt, kopt)
pl(:, knew) = pl(:, knew) / vlag(knew)
plnew = pl(:, knew)
pl = pl - outprod(plnew, vlag)
pl(:, knew) = plnew
pq = pq + moderr * plnew

! Replace the interpolation point that has index KNEW by the point XNEW.
! Update KOPT if F is the least calculated value of the objective function.
xpt(:, knew) = xopt + d
if (f < fopt) then
    kopt = knew
    fopt = f
    xopt = xpt(:, kopt)
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt .and. all(is_finite(xpt)), &
        & 'SIZE(XPT) == [N, NPT], XPT is finite', srname)
end if

end subroutine update


end module update_mod
