module shiftbase_mod
!--------------------------------------------------------------------------------------------------!
! This module contanis a subroutine that shifts the base point from XBASE to XBASE + XPT. It is
! used in NEWUOA, BOBYQA, and LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Tuesday, March 21, 2023 AM10:38:48
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: shiftbase

interface shiftbase
    module procedure shiftbase_lfqint, shiftbase_qint
end interface shiftbase

contains


subroutine shiftbase_lfqint(kopt, xbase, xpt, zmat, bmat, pq, hq, idz)
!--------------------------------------------------------------------------------------------------!
! This subroutine shifts the base point from XBASE to XBASE + XOPT and updates BMAT and HQ
! accordingly. PQ and ZMAT remain the same after the shifting. See Section 7 of the NEWUOA paper.
! N.B.:
! 1. In Powell's implementation of NEWUOA, the quadratic model is represented by [GQ, PQ, HQ], where
! GQ is the gradient of the quadratic model at XBASE. In that case, GQ should be updated by
! GQ = GQ + HESS_MUL(XOPT, XPT, PQ, HQ) in this subroutine, where PQ is the un-updated version.
! However, Powell implemented BOBYQA and LINCOA without GQ but with GOPT, which is the gradient at
! XBASE + XOPT. Note that GOPT remains unchanged when XBASE is shifted. In our implementation,
! NEWUOA also uses GOPT instead of GQ.
! 2. [IDZ, ZMAT] provides the factorization of Omega in (3.17) of the NEWUOA paper; in specific,
! Omega = sum_{i=1}^{NPT-N-1} s_i*ZMAT(:,i)*ZMAT(:,i)^T, s_i = -1 if i < IDZ, and si = 1 if i >= IDZ.
! In precise arithmetic, IDZ should be always 1; to cope with rounding errors, NEWUOA and LINCOA
! allow IDZ = -1 (see (4.18)--(4.20) of the NEWUOA paper); in BOBYQA, IDZ is always 1, and the
! rounding errors are handled by the RESCUE subroutine (Sec. 5 of the BOBYQA paper).
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: BY(N, N), SXPT(NPT), V(N), VXOPT(N, N), XPTXAV(N, NPT), YMAT(N, NPT),
! YZMAT(N, NPT-N-1), YZMAT_C(N, NPT-N-1).
! Size of local arrays: REAL(RP)*(4*N*NPT+NPT-N), LARGE!
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod, outprod, issymmetric

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: pq(:)   ! PQ(NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
integer(IK), intent(in), optional :: idz  ! Absent in BOBYQA, being equivalent to IDZ = 1

! In-outputs
real(RP), intent(inout) :: bmat(:, :)   ! BMAT(N, NPT + N)
real(RP), intent(inout) :: hq(:, :) ! HQ(N, N)
real(RP), intent(inout) :: xbase(:) ! XBASE(N)
real(RP), intent(inout) :: xpt(:, :)    ! XPT(N, NPT)

! Local variables
character(len=*), parameter :: srname = 'SHIFTBASE_LFQINT'
integer(IK) :: idz_loc
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bymat(size(xbase), size(xbase))
!real(RP) :: htol
real(RP) :: qxoptq
real(RP) :: sxpt(size(xpt, 2))
real(RP) :: v(size(xbase))
real(RP) :: vxopt(size(xbase), size(xbase))
real(RP) :: xopt(size(xbase))
real(RP) :: xoptsq
real(RP) :: xptxav(size(xpt, 1), size(xpt, 2))
real(RP) :: ymat(size(xpt, 1), size(xpt, 2))
real(RP) :: yzmat(size(xbase), size(zmat, 2))
real(RP) :: yzmat_c(size(xbase), size(zmat, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Read IDZ, which is absent from BOBYQA, being equivalent to IDZ = 1.
idz_loc = 1
if (present(idz)) then
    idz_loc = idz
end if

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(idz_loc >= 1 .and. idz_loc <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    ! The following test cannot be passed.
    !htol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E10_RP * EPS)) ! Tolerance for error in H
    !call assert(errh(idz_loc, bmat, zmat, xpt) <= htol, 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
end if

!====================!
! Calculation starts !
!====================!

! Read XOPT.
xopt = xpt(:, kopt)
xoptsq = inprod(xopt, xopt)

! Update BMAT. See (7.11)--(7.12) of the NEWUOA paper and the elaborations around.
! XPTXAV corresponds to XPT - XAV in the NEWUOA paper, with XAV = (X0 + XOPT)/2.
xptxav = xpt - HALF * spread(xopt, dim=2, ncopies=npt)
!!MATLAB: xptxav = xpt - xopt/2  % xopt should be a column! Implicit expansion
!sxpt = matprod(xopt, xptxav)
sxpt = matprod(xopt, xpt) - HALF * xoptsq  ! This one seems to work better numerically.

! First, make the changes to BMAT that do not depend on ZMAT.
qxoptq = QUART * xoptsq
do k = 1, npt
    ymat(:, k) = sxpt(k) * xptxav(:, k) + qxoptq * xopt
end do
!!MATLAB: ymat = xptxav .* sxpt + qxoptq * xopt  % sxpt should be a row, xopt should be a column
!ymat(:, kopt) = HALF * xoptsq * xopt ! This makes no difference according to a test on 20220406
bymat = matprod(bmat(:, 1:npt), transpose(ymat))  ! BMAT(:, 1:NPT) is not updated yet.
bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + (bymat + transpose(bymat))
! Then the revisions of BMAT that depend on ZMAT are calculated.
yzmat = matprod(ymat, zmat)
yzmat_c = yzmat
yzmat_c(:, 1:idz_loc - 1) = -yzmat(:, 1:idz_loc - 1)  ! IDZ_LOC is usually small. So this assignment is cheap.
bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + matprod(yzmat, transpose(yzmat_c))
bmat(:, 1:npt) = bmat(:, 1:npt) + matprod(yzmat_c, transpose(zmat))

! Update the quadratic model. Note that PQ remains unchanged. For HQ, see (7.14) of the NEWUOA paper.
!v = matprod(xptxav, pq)  ! Vector V in (7.14) of the NEWUOA paper
v = matprod(xpt, pq) - HALF * sum(pq) * xopt ! This one seems to work better numerically.
vxopt = outprod(v, xopt)  !!MATLAB: vxopt = v * xopt';  % v and xopt should be both columns
hq = (vxopt + transpose(vxopt)) + hq !call r2update(hq, ONE, xopt, v)
!call symmetrize(hq)  ! Do this if the update above does not ensure symmetry.

! The following instructions complete the shift of XBASE.
xbase = xbase + xopt
xpt = xpt - spread(xopt, dim=2, ncopies=npt)
xpt(:, kopt) = ZERO
!!MATLAB: xpt = xpt - xopt; xpt(:, kopt) = 0;  % xopt should be a column! Implicit expansion

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(abs(xpt(:, kopt)) <= 0), 'XPT(:, KOPT) == 0', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    ! The following test cannot be passed.
    !call assert(errh(idz_loc, bmat, zmat, xpt) <= htol, 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
end if

end subroutine shiftbase_lfqint


subroutine shiftbase_qint(kopt, pl, pq, xbase, xpt)
!--------------------------------------------------------------------------------------------------!
! This subroutine shifts the base point from XBASE to XBASE + XOPT, and make the corresponding
! changes to the gradients of the Lagrange functions and the quadratic model. See the discussion
! below (40) of the UOBYQA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : smat_mul_vec

implicit none

! Inputs
integer(IK), intent(in) :: kopt

! In-outputs
real(RP), intent(inout) :: xbase(:)  ! XBASE(N)
real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(inout) :: pl(:, :)  ! PL(NPT-1, NPT)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT-1)

! Local variables
character(len=*), parameter :: srname = 'SHIFTBASE_QINT'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: xopt(size(xbase))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(npt == (n + 1) * (n + 2) / 2, 'NPT = (N+1)(N+2)/2', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(pl, 1) == npt - 1 .and. size(pl, 2) == npt, 'SIZE(PL) == [NPT-1, NPT]', srname)
    call assert(size(pq) == npt - 1, 'SIZE(PQ) == NPT-1', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt .and. all(is_finite(xpt)), &
        & 'SIZE(XPT) == [N, NPT], XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

! Shift the base point from XBASE to XBASE + XOPT.
xopt = xpt(:, kopt)
xbase = xbase + xopt
xpt = xpt - spread(xopt, dim=2, ncopies=npt)
xpt(:, kopt) = ZERO

! Update the gradient of the model
pq(1:n) = pq(1:n) + smat_mul_vec(pq(n + 1:npt - 1), xopt)

! Upda the gradient of the Lagrange functions.
do k = 1, npt
    pl(1:n, k) = pl(1:n, k) + smat_mul_vec(pl(n + 1:npt - 1, k), xopt)
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(pl, 1) == npt - 1 .and. size(pl, 2) == npt, 'SIZE(PL) == [NPT-1, NPT]', srname)
    call assert(size(pq) == npt - 1, 'SIZE(PQ) == NPT-1', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt .and. all(is_finite(xpt)), &
        & 'SIZE(XPT) == [N, NPT], XPT is finite', srname)
end if

end subroutine shiftbase_qint


end module shiftbase_mod
