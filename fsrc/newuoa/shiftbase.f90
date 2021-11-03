module shiftbase_mod
!--------------------------------------------------------------------------------------------------!
! This module contanis a subroutine that shifts the base point from XBASE to XBASE + XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Wednesday, November 03, 2021 PM11:48:28
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: shiftbase


contains


subroutine shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt)
!--------------------------------------------------------------------------------------------------!
! SHIFTBASE shifts the base point for XBASE to XBASE + XOPT and updates GQ, HQ, and BMAT
! accordingly. PQ and ZMAT remain the same after the shifting. See Section 7 of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : r1update, r2update, inprod, matprod, issymmetric, hessmul

implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: pq(:)   ! PQ(NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! In-outputs
real(RP), intent(inout) :: bmat(:, :)   ! BMAT(N, NPT + N)
real(RP), intent(inout) :: gq(:)    ! GQ(N)
real(RP), intent(inout) :: hq(:, :) ! HQ(N, N)
real(RP), intent(inout) :: xbase(:) ! XBASE(N)
real(RP), intent(inout) :: xopt(:)  ! XOPT(N)
real(RP), intent(inout) :: xpt(:, :)    ! XPT(N, NPT)

! Local variables
character(len=*), parameter :: srname = 'SHIFTBASE'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bmatk(size(bmat, 1))
real(RP) :: qxoptq
real(RP) :: sumz(size(zmat, 2))
real(RP) :: t
real(RP) :: vlag(size(xopt))
real(RP) :: w1(size(pq))
real(RP) :: w2(size(xopt))
real(RP) :: xoptsq
real(RP) :: xpq(size(xopt))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetrize', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(gq) == n, 'SIZE(GQ) == N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), &
        & 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT)==N, XOPT is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

xoptsq = inprod(xopt, xopt)
qxoptq = QUART * xoptsq

! Update GQ.
gq = hessmul(hq, pq, xpt, xopt) + gq

! Update HQ. See (7.14) of the NEWUOA paper.
w1 = matprod(xopt, xpt) - HALF * xoptsq
! W1 equals MATPROD(XPT, XOPT) after XPT is updated TEMPORARILY as follows.
xpt = xpt - HALF * spread(xopt, dim=2, ncopies=npt)  ! TEMPORARY
xpq = matprod(xpt, pq)
call r2update(hq, ONE, xopt, xpq)  ! Implement R2UPDATE properly so that HQ is symmetric.

! Update BMAT. See (7.11)--(7.12) of the NEWUOA paper.
! First, make the changes to BMAT that do not depend on ZMAT.
do k = 1, npt
    bmatk = bmat(:, k)
    w2 = w1(k) * xpt(:, k) + qxoptq * xopt
    ! Implement R2UPDATE properly so that BMAT(:, NPT+1:NPT+N) is symmetric.
    call r2update(bmat(:, npt + 1:npt + n), ONE, bmatk, w2)
end do
! Then the revisions of BMAT that depend on ZMAT are calculated.
sumz = sum(zmat, dim=1)
do k = 1, int(npt - n - 1, kind(k))
    vlag = qxoptq * sumz(k) * xopt + matprod(xpt, w1 * zmat(:, k))
    if (k <= idz - 1) then
        t = -ONE
    else
        t = ONE
    end if
    call r1update(bmat(:, 1:npt), t, vlag, zmat(:, k))
    ! Implement R1UPDATE properly so that BMAT(:, NPT+1:NPT+N) is symmetric.
    call r1update(bmat(:, npt + 1:npt + n), t, vlag)
end do

! The following instructions complete the shift of XBASE. Recall the we have already subtracted
! HALF*XOPT from XPT. Therefore, overall, the new XPT is XPT - XOPT.
xpt = xpt - HALF * spread(xopt, dim=2, ncopies=npt)
xbase = xbase + xopt
xopt = ZERO

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(issymmetric(hq), 'HQ is symmetric', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
end if

end subroutine shiftbase


end module shiftbase_mod
