module shiftbase_mod
!--------------------------------------------------------------------------------------------------!
! This module contanis a subroutine that shifts the base point from XBASE to XBASE + XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Tuesday, April 05, 2022 AM08:34:31
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
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: BY(N, N), SXPT(NPT), V(N), VXOPT(N, N), XPTXAV(N, NPT), YMAT(N, NPT),
! YZMAT(N, NPT-N-1), YZMAT_C(N, NPT-N-1).
! Size of local arrays: REAL(RP)*(4*N*NPT+NPT-N), LARGE!!!
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod, outprod, issymmetric, hess_mul

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
real(RP) :: by(size(xopt), size(xopt))
real(RP) :: qxoptq
real(RP) :: sxpt(size(pq))
real(RP) :: v(size(xopt))
real(RP) :: vxopt(size(xopt), size(xopt))
real(RP) :: xoptsq
real(RP) :: xptxav(size(xpt, 1), size(xpt, 2))
real(RP) :: ymat(size(xopt), size(pq))
real(RP) :: yzmat(size(xopt), size(zmat, 2))
real(RP) :: yzmat_c(size(xopt), size(zmat, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

xoptsq = inprod(xopt, xopt)

xptxav = xpt - HALF * spread(xopt, dim=2, ncopies=npt)
!!MATLAB: xptxav = xpt - xopt/2  % xopt should be a column!! Implicit expansion
sxpt = matprod(xopt, xptxav)
!sxpt = matprod(xopt, xpt) - HALF * xoptsq

! Update BMAT. See (7.11)--(7.12) of the NEWUOA paper and the elaborations around.

! First, make the changes to BMAT that do not depend on ZMAT.
qxoptq = QUART * xoptsq
do k = 1, npt
    ymat(:, k) = sxpt(k) * xptxav(:, k) + qxoptq * xopt  ! Should it be called VLAG or W2?
end do
!!MATLAB: ymat = xptxav .* sxpt + qxoptq * xopt  % sxpt should be a row, xopt should be a column
by = matprod(bmat(:, 1:npt), transpose(ymat))  ! BMAT(:, 1:NPT) is not updated yet.
bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + (by + transpose(by))
!call symmetrize(bmat(:, npt + 1:npt + n))  ! Do this if the update above does not ensure symmetry

! Then the revisions of BMAT that depend on ZMAT are calculated.
yzmat = matprod(ymat, zmat)
yzmat_c(:, 1:idz - 1) = -yzmat(:, 1:idz - 1)
yzmat_c(:, idz:npt - n - 1) = yzmat(:, idz:npt - n - 1)
bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + matprod(yzmat, transpose(yzmat_c))
bmat(:, 1:npt) = bmat(:, 1:npt) + matprod(yzmat_c, transpose(zmat))

! Update the quadratic model. Only GQ and HQ need revision. For HQ, see (7.14) of the NEWUOA paper.
gq = hess_mul(hq, pq, xpt, xopt) + gq
v = matprod(xptxav, pq)  ! Vector V in (7.14) of the NEWUOA paper
vxopt = outprod(v, xopt)  !!MATLAB: vxopt = v * xopt';  % v and xopt should be both columns
hq = hq + (vxopt + transpose(vxopt))
!call symmetrize(hq)  ! Do this if the update above does not ensure symmetry

! The following instructions complete the shift of XBASE.
xpt = xpt - spread(xopt, dim=2, ncopies=npt)
!!MATLAB: xpt = xpt - xopt  % xopt should be a column!! Implicit expansion
xbase = xbase + xopt
xopt = ZERO

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

end subroutine shiftbase


end module shiftbase_mod
