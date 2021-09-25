! SHIFTBASE_MOD is a module contaning a subroutine that shifts the base
! point from XBASE to XBASE + XPT.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.
!
! Last Modified: Friday, July 23, 2021 PM12:02:24

module shiftbase_mod

implicit none
private
public :: shiftbase


contains


subroutine shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt)
! SHIFTBASE shifts the base point for XBASE to XBASE + XOPT and updates
! GQ, HQ, and BMAT accordingly. PQ and ZMAT remain the same after the
! shifting. See Section 7 of the NEWUOA paper.

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, HALF, QUART, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
use lina_mod, only : Ax_plus_y, r1update, r2update, inprod, matprod, hessmul

implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: pq(:)  ! PQ(NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! In-outputs
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: gq(:)  ! GQ(N)
real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(inout) :: xbase(:)  ! XBASE(N)
real(RP), intent(inout) :: xopt(:)  ! XOPT(N)
real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)

! Local variables
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bmatk(size(bmat, 1))
real(RP) :: qxoptq
real(RP) :: sumz(size(zmat, 2))
real(RP) :: vlag(size(xopt))
real(RP) :: w1(size(pq))
real(RP) :: w2(size(xopt))
real(RP) :: xoptsq
real(RP) :: xpq(size(xopt))
character(len=SRNLEN), parameter :: srname = 'SHIFTBASE'


! Get and verify the sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

if (DEBUGGING) then
    if (n == 0 .or. npt < n + 2) then
        call errstop(srname, 'SIZE(XPT) is invalid')
    end if
    call verisize(pq, npt)
    call verisize(zmat, npt, int(npt - n - 1, kind(n)))
    call verisize(bmat, n, npt + n)
    call verisize(gq, n)
    call verisize(hq, n, n)
    call verisize(xopt, n)
    call verisize(xbase, n)
end if


xoptsq = inprod(xopt, xopt)
qxoptq = QUART * xoptsq

!---------------------------------------------------------------------------!
!--------------! gq = hessmul(hq, pq, xpt, xopt) + gq !---------------------!
!--------------! OR (not numerically equivalent)!---------------------------!
gq = matprod(hq, xopt) + (matprod(xpt, pq * matprod(xopt, xpt)) + gq) !-!
!gq = Ax_plus_y(xpt, pq * matprod(xopt, xpt), gq)
!gq = Ax_plus_y(hq, xopt, gq)
!---------------------------------------------------------------------------!

w1 = matprod(xopt, xpt) - HALF * xoptsq
! W1 equals MATPROD(XPT, XOPT) after XPT is updated as follows.
xpt = xpt - HALF * spread(xopt, dim=2, ncopies=npt)
!do k = 1, npt
!    xpt(:, k) = xpt(:, k) - HALF*xopt
!end do

! Update HQ. It has to be done after the above revision to XPT!!!
xpq = matprod(xpt, pq)

!----------------------------------------------------------------!
! Implement R2UPDATE properly so that it ensures HQ is symmetric.
call r2update(hq, ONE, xopt, xpq)
!----------------------------------------------------------------!


! Make the changes to BMAT that do not depend on ZMAT.
do k = 1, npt
    bmatk = bmat(:, k)
    w2 = w1(k) * xpt(:, k) + qxoptq * xopt
    ! Implement R2UPDATE properly so that it ensures
    ! bmat(:, npt+1:npt+n) is symmetric.
    call r2update(bmat(:, npt + 1:npt + n), ONE, bmatk, w2)
end do

! Then the revisions of BMAT that depend on ZMAT are calculated.
sumz = sum(zmat, dim=1)
do k = 1, int(idz - 1, kind(k))
!----------------------------------------------------------------------!
    vlag = qxoptq * sumz(k) * xopt + matprod(xpt, w1 * zmat(:, k)) !--!
!    vlag = Ax_plus_y(xpt, w1 * zmat(:, k), qxoptq * sumz(k) * xopt)
!----------------------------------------------------------------------!
    call r1update(bmat(:, 1:npt), -ONE, vlag, zmat(:, k))
    ! Implement R1UPDATE properly so that it ensures
    ! bmat(:, npt+1:npt+n) is symmetric.
    call r1update(bmat(:, npt + 1:npt + n), -ONE, vlag)
end do
do k = idz, int(npt - n - 1, kind(k))
!----------------------------------------------------------------------!
    vlag = qxoptq * sumz(k) * xopt + matprod(xpt, w1 * zmat(:, k)) !--!
!    vlag = Ax_plus_y(xpt, w1 * zmat(:, k), qxoptq * sumz(k) * xopt)
!----------------------------------------------------------------------!
    call r1update(bmat(:, 1:npt), ONE, vlag, zmat(:, k))
    ! Implement R1UPDATE properly so that it ensures
    ! bmat(:, npt+1:npt+n) is symmetric.
    call r1update(bmat(:, npt + 1:npt + n), ONE, vlag)
end do

!----------------------------------------------------------------!
! The following instructions complete the shift of XBASE.
! Recall the we have already subtracted HALF*XOPT from XPT.
! Therefore, overall, the new XPT is XPT - XOPT.
xpt = xpt - HALF * spread(xopt, dim=2, ncopies=npt)
!do k = 1, npt
!    xpt(:, k) = xpt(:, k) - HALF*xopt
!end do

xbase = xbase + xopt
xopt = ZERO

end subroutine shiftbase


end module shiftbase_mod
