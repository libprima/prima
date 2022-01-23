module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines for initializing FVAL, XBASE, XPT, GQ, HQ, PQ, IDZ, ZMAT, BMAT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Thursday, January 06, 2022 PM08:47:00
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initxf, initq, inith


contains


subroutine initxf(calfun, iprint, maxfun, ftarget, rhobeg, x0, ij, kopt, nf, fhist, fval, xbase, &
    & xhist, xpt, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine does the initialization regarding the interpolation points & their function values
!
! N.B.:
! 1. Remark on IJ:
! If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [NPT-2*N-1, 2].
! IJ contains integers between 2 and 2*N + 1. For each K > 2*N + 1, XPT(:, K) is
! XPT(:, IJ(K, 1)) + XPT(:, IJ(K, 2)). Let I = IJ(K, 1) - 1 if such a number is <= N + 1; otherwise,
! let IJ = IJ(K, 1) - N - 1; define J using IJ(K, 2) similarly. Then all the entries of XPT(:, K)
! are zero except that the I and J entries are RHOBEG or -RHOBEG. Indeed, XPT(I, K) is RHOBEG if
! IJ(K, 1) <= N + 1 and -RHOBEG otherwise; XPT(J, K) is similar. Consequently, the Hessian of the
! quadratic model will get a possibly nonzero (I, J) entry.
! 2. At return,
! INFO = INFO_DFT: initialization finishes normally
! INFO = FTARGET_ACHIEVED: return because F <= FTARGET
! INFO = NAN_INF_X: return because X contains NaN
! INFO = NAN_INF_F: return because F is either NaN or +Inf
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : eye, sort
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)   ! X0(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: ij(:, :)    ! IJ(MAX(0_IK, NPT-2*N-1_IK), 2)
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: fval(:)    ! FVAL(NPT)
real(RP), intent(out) :: fhist(:)   ! FHIST(MAXFHIST)
real(RP), intent(out) :: xbase(:)   ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)    ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)

! Local variables
character(len=*), parameter :: solver = 'NEWUOA'
character(len=*), parameter :: srname = 'INITXF'
integer(IK) :: k
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: subinfo
logical :: evaluated(size(fval))
real(RP) :: f
real(RP) :: x(size(x0))

! Sizes
n = int(size(x0), kind(n))
npt = int(size(fval), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(all(is_finite(x0)), 'X0 is finite', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
end if

!====================!
! Calculation starts !
!====================!

info = INFO_DFT  ! Default info.

! Initialize XBASE to X0.
xbase = x0

! Begin the initialization procedure. The coordinates of the displacement of the next initial
! interpolation point from XBASE are set in XPT(:, :).

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously. However, the loop here
! is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
evaluated = .false.
! Initialize FVAL to HUGENUM. Otherwise, compilers may complain that FVAL is not (completely)
! initialized if the initialization aborts due to abnormality (see CHECKEXIT).
fval = HUGENUM

! Set XPT, FVAL, KOPT, and XOPT.

! Set XPT.
xpt(:, 1) = ZERO
xpt(:, 2:n + 1) = rhobeg * eye(n)
xpt(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1_IK)  ! XPT(:, 2*N+2 : NPT) = ZERO if it is nonempty.

! Set FVAL(1 : 2*N + 1) by evaluating F. Totally parallelizable except for FMSG.
do k = 1, min(npt, int(2 * n + 1, kind(npt)))
    x = xpt(:, k) + xbase
    call evaluate(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call fmsg(solver, iprint, k, f, x)
    ! Save X and F into the history.
    call savehist(k, f, x, fhist, xhist)
    ! Check whether to exit.
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do

! Set IJ.
! In general, when NPT = (N+1)*(N+2)/2, we can initialize IJ(2*N + 2 : NPT, :) to ANY permutation
! of {(I, J) : 1 <= J < I <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first
! NPT - (2*N + 1) elements of such a permutation. Powell took the following permutation.
ij(:, 1) = int(([(k, k=2_IK * n + 2_IK, npt)] - n - 2) / n, IK) ! [(K, K=1, NPT)] = [1, 2, ..., NPT]
ij(:, 2) = int([(k, k=2_IK * n + 2_IK, npt)] - (ij(:, 1) + 1) * n - 1, IK)
ij(:, 1) = modulo(ij(:, 1) + ij(:, 2) - 1_IK, n) + 1_IK  ! MODULO(K-1,N) + 1 = K-N for K in [N+1,2N]
! The next line ensures IJ(:, 1) > IJ(:, 2).
ij = sort(ij, 2, 'descend')
! Increment IJ by 1. This 1 comes from the fact that XPT(:, 1) corresponds to the base point XBASE.
ij = ij + 1_IK
! Further revise IJ according to FVAL(2 : 2*N + 1).
! N.B.:
! 1. For each K below, the following lines revises IJ(K, :) as follows:
! change IJ(K, 1) to IJ(K, 1) + N if FVAL(IJ(K, 1) + N) < FVAL(IJ(K, 1));
! change IJ(K, 2) to IJ(K, 2) + N if FVAL(IJ(K, 2) + N) < FVAL(IJ(K, 2)).
! 2. The idea of this revision is as follows: Let [I, J] = IJ(K, :) with the IJ BEFORE the revision;
! XPT(:, K) is the sum of either {XPT(:, I+1) or XPT(:, I+N+1)} + {XPT(:, J+1) or XPT(:, J+N+1)},
! each choice being made in favor of the point that has a lower function value, with the hope that
! such a choice will more likely render an XPT(:, K) with a low function value.
! 3. This revision is OPTIONAL. Due to this revision, the definition of XPT(:, 2*N + 2 : NPT) relies
! on FVAL(2 : 2*N + 1), and it is the sole origin of the such dependency. If we remove the revision
! IJ, then the evaluations of FVAL(1 : NPT) can be merged, and they are totally PARALLELIZABLE; this
! can be beneficial if the function evaluations are expensive, which is likely the case.
! 4. MATLAB can index a vector using a 2D array of indices (not in Fortran), thus the MATLAB code is
! % IJ(FVAL(IJ + N) < FVAL(IJ)) = IJ + N;
where (fval(ij(:, 1) + n) < fval(ij(:, 1)))
    ij(:, 1) = ij(:, 1) + n
end where
where (fval(ij(:, 2) + n) < fval(ij(:, 2)))
    ij(:, 2) = ij(:, 2) + n
end where

! Set XPT(:, 2*N + 2 : NPT). It depends on IJ and hence on FVAL(2 : 2*N + 1).
xpt(:, 2 * n + 2:npt) = xpt(:, ij(:, 1)) + xpt(:, ij(:, 2))

! Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except for FMSG.
if (info == INFO_DFT) then
    do k = int(2 * n + 2, kind(k)), npt
        x = xpt(:, k) + xbase
        call evaluate(calfun, x, f)
        call fmsg(solver, iprint, k, f, x)
        evaluated(k) = .true.
        fval(k) = f
        ! Save X and F into the history.
        call savehist(k, f, x, fhist, xhist)
        ! Check whether to exit.
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end do
end if

! Set NF, KOPT
nf = int(count(evaluated), kind(nf))
kopt = int(minloc(fval, dim=1, mask=evaluated), kind(kopt))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(all(ij >= 2 .and. ij <= 2 * n + 1), '1 <= IJ <= N', srname)
    call assert(all(modulo(ij(:, 1) - 2_IK, n) > modulo(ij(:, 2) - 2_IK, n)), &
        & 'MODULO(IJ(:, 1) - 2, N) > MODULO(IJ(:, 2) - 2, N)', srname)
    call assert(nf <= npt, 'NF <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= nf, '1 <= KOPT <= NF', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(evaluated .and. (is_nan(fval) .or. is_posinf(fval))), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
end if

end subroutine initxf


subroutine initq(ij, fval, xpt, gq, hq, pq, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes the quadratic model, which is represented by (GQ, HQ, PQ) so that its
! gradient at XBASE is GQ; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL
use, non_intrinsic :: linalg_mod, only : issymmetric

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :)     ! IJ(NPT, 2)
real(RP), intent(in) :: fval(:)     ! FVAL(NPT)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK), intent(out) :: info
real(RP), intent(out) :: gq(:)  ! GQ(N)
real(RP), intent(out) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)

! Local variables
character(len=*), parameter :: srname = 'INITQ'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: fbeg
real(RP) :: rhobeg
real(RP) :: xi
real(RP) :: xj

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(all(ij >= 2 .and. ij <= 2 * n + 1), '1 <= IJ <= N', srname)
    call assert(all(modulo(ij(:, 1) - 2_IK, n) > modulo(ij(:, 2) - 2_IK, n)), &
        & 'MODULO(IJ(:, 1) - 2, N) > MODULO(IJ(:, 2) - 2, N)', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) = [N, N]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

gq = ZERO
hq = ZERO
pq = ZERO  ! We will not update PQ. It is ZERO at return.

rhobeg = maxval(abs(xpt(:, 2)))  ! Read RHOBEG from XPT.
fbeg = fval(1)

! Set GQ by forward difference.
gq(1:n) = (fval(2:n + 1) - fbeg) / rhobeg
! If possible, revise GQ to central difference.
k = min(int(npt - n - 1, kind(n)), n)
gq(1:k) = HALF * (gq(1:k) + (fbeg - fval(n + 2:n + 1 + k)) / rhobeg)

! Set the diagonal of HQ by 2nd-order central finite difference.
do k = 1, min(int(npt - n - 1, kind(n)), n)
    hq(k, k) = ((fval(k + 1) - fbeg) / rhobeg - (fbeg - fval(k + n + 1)) / rhobeg) / rhobeg
end do
! When NPT > 2*N + 1, set the off-diagonal entries of HQ.
do k = 1, npt - 2_IK * n - 1_IK
    ! With the I, J, XI, and XJ defined below, we have
    ! FVAL(K+2*n+1) = F(XBASE + XI*e_I + XJ*e_J),
    ! FVAL(IJ(K, 1)) = F(XBASE + XI*e_I),
    ! FVAL(IJ(K, 2)) = F(XBASE + XJ*e_J).
    ! Thus the HQ(I,J) defined below approximates frac{partial^2}{partial X_I partial X_J} F(XBASE).
    i = modulo(ij(k, 1) - 2_IK, n) + 1_IK
    j = modulo(ij(k, 2) - 2_IK, n) + 1_IK
    xi = xpt(i, k + 2 * n + 1)
    xj = xpt(j, k + 2 * n + 1)
    hq(i, j) = (fbeg - fval(ij(k, 1)) - fval(ij(k, 2)) + fval(k + 2 * n + 1)) / (xi * xj)
    hq(j, i) = hq(i, j)
end do



if (any(is_nan(gq)) .or. any(is_nan(hq)) .or. any(is_nan(pq))) then
    info = NAN_MODEL
else
    info = INFO_DFT
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
end if

end subroutine initq


subroutine inith(ij, xpt, idz, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes IDZ, BMAT, and ZMAT.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL
use, non_intrinsic :: linalg_mod, only : issymmetric, eye

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :) ! IJ(NPT, 2)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK), intent(out) :: idz
integer(IK), intent(out) :: info
real(RP), intent(out) :: bmat(:, :) ! BMAT(N, NPT + N)
real(RP), intent(out) :: zmat(:, :) ! ZMAT(NPT, NPT - N - 1)

! Local variables
character(len=*), parameter :: srname = 'INITH'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: rhobeg
real(RP) :: rhosq

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(all(ij >= 2 .and. ij <= 2 * n + 1), '1 <= IJ <= N', srname)
    call assert(all(modulo(ij(:, 1) - 2_IK, n) > modulo(ij(:, 2) - 2_IK, n)), &
        & 'MODULO(IJ(:, 1) - 2, N) > MODULO(IJ(:, 2) - 2, N)', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! Set IDZ = 1. It will not be changed in the following.
idz = 1_IK

! Some values to be used for setting BMAT and ZMAT.
rhobeg = maxval(abs(xpt(:, 2)))  ! Read RHOBEG from XPT.
rhosq = rhobeg**2

! Set BMAT.
bmat = ZERO
if (npt <= 2 * n + 1) then
    ! Set BMAT(1 : NPT-N-1, :)
    bmat(1:npt - n - 1, 2:npt - n) = (HALF / rhobeg) * eye(npt - n - 1_IK)
    bmat(1:npt - n - 1, n + 2:npt) = -(HALF / rhobeg) * eye(npt - n - 1_IK)
    ! Set BMAT(NPT-N : N, :)
    bmat(npt - n:n, 1) = -ONE / rhobeg
    bmat(npt - n:n, npt - n + 1:n + 1) = (ONE / rhobeg) * eye(2_IK * n - npt + 1_IK)
    bmat(npt - n:n, 2 * npt - n:npt + n) = -(HALF * rhosq) * eye(2_IK * n - npt + 1_IK)
else
    bmat(:, 2:n + 1) = (HALF / rhobeg) * eye(n)
    bmat(:, n + 2:2 * n + 1) = -(HALF / rhobeg) * eye(n)
end if

! Set ZMAT.
zmat = ZERO
if (npt <= 2 * n + 1) then
    zmat(1, :) = -sqrt(TWO) / rhosq
    zmat(2:npt - n, :) = (sqrt(HALF) / rhosq) * eye(npt - n - 1_IK)
    zmat(n + 2:npt, :) = (sqrt(HALF) / rhosq) * eye(npt - n - 1_IK)
else
    ! Set ZMAT(:, 1:N).
    zmat(1, 1:n) = -sqrt(TWO) / rhosq
    zmat(2:n + 1, 1:n) = (sqrt(HALF) / rhosq) * eye(n)
    zmat(n + 2:2 * n + 1, 1:n) = (sqrt(HALF) / rhosq) * eye(n)
    ! Set ZMAT(:, N+1 : NPT-N-1).
    zmat(1, n + 1:npt - n - 1) = ONE / rhosq
    zmat(2 * n + 2:npt, n + 1:npt - n - 1) = (ONE / rhosq) * eye(npt - 2_IK * n - 1_IK)
    do k = 1, npt - 2_IK * n - 1_IK
        zmat(ij(k, :), k + n) = -ONE / rhosq
    end do
end if

if (any(is_nan(bmat)) .or. any(is_nan(zmat))) then
    info = NAN_MODEL
else
    info = INFO_DFT
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

end subroutine inith


end module initialize_mod
