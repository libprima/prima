module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines for initializing FVAL, XBASE, XPT, GQ, HQ, PQ, IDZ, ZMAT, BMAT.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Last Modified: Saturday, September 11, 2021 PM01:47:38
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
! When K > 2*N + 1, all the entries of XPT(:, K) will be zero except that the IJ(1, K) and IJ(2, K)
! entries will be set to RHOBEG or -RHOBEG. Consequently, the Hessian of the quadratic model will
! get a possibly nonzero (IJ(1, K), IJ(2, K)) entry. Indeed, IJ(:, 1 : 2*N + 1) is never used.
! 2. At return,
! INFO = INFO_DFT: initialization finishes normally
! INFO = FTARGET_ACHIEVED: return because F <= FTARGET
! INFO = NAN_X: return because X contains NaN
! INFO = NAN_INF_F: return because F is either NaN or +Inf
!--------------------------------------------------------------------------------------------------!

! Generic modules
use checkexit_mod, only : checkexit
use consts_mod, only : RP, IK, ZERO, DEBUGGING
use debug_mod, only : assert
use evaluate_mod, only : evalf
use history_mod, only : savehist
use infnan_mod, only : is_finite
use info_mod, only : INFO_DFT
use linalg_mod, only : eye
use output_mod, only : fmssg
use pintrf_mod, only : FUN

implicit none

! Inputs
procedure(FUN) :: calfun
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)   ! X0(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: ij(:, :)    ! IJ(2, NPT)
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
integer(IK) :: i
integer(IK) :: itemp
integer(IK) :: j
integer(IK) :: k
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: npt_revised
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
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXXHIST', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXXHIST', srname)
    call assert(size(ij, 1) == 2 .and. size(ij, 2) == npt, 'SIZE(IJ) == [2, NPT]', srname)
    call assert(rhobeg > ZERO, 'RHOBEG > 0', srname)
    call assert(all(is_finite(x0)), 'X0 is finite', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
end if

!====================!
! Calculation starts !
!====================!

info = INFO_DFT  ! Default info.

! We set IJ = 1 in case the initialization aborts due to abnormality. If we do not do this, IJ will
! be undefined if the initialization aborts.
ij = 1_IK

! Initialize XBASE to X0.
xbase = x0

! Begin the initialization procedure. The coordinates of the displacement of the next initial
! interpolation point from XBASE are set in XPT(:, :).

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously. However, the loop here
! is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
evaluated = .false.

! NPT_REVISED = NPT unless it turns out necessary to return abnormally (NaN/Inf, or F <= FTARGET).
npt_revised = npt


! Set XPT, FVAL, KOPT, and XOPT.

! Set XPT(:, 1) = ZERO
xpt(:, 1) = ZERO
! Set XPT(:, 2 : N + 1).
xpt(:, 2:n + 1) = rhobeg * eye(n)
! Set XPT(:, N+2 : NPT). XPT(:, 2*N+2 : NPT) is initialized to ZERO if it is nonempty.
xpt(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1)

! Set FVAL(1 : 2*N + 1) by evaluating F. Totally parallelizable except for FMSSG.
do k = 1, min(npt, int(2 * n + 1, kind(npt)))
    x = xpt(:, k) + xbase
    call evalf(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call fmssg(iprint, k, f, x, solver)
    ! Save X and F into the history.
    call savehist(k, f, x, fhist, xhist)
    ! Check whether to exit.
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        npt_revised = 0_IK
        exit
    end if
end do

! Set XPT(:, 2*N + 2 : NPT). It depends on FVAL(2 : 2*N + 1).
do k = int(2 * n + 2, kind(k)), npt_revised
    ! Decide IJ(:, K).  In general, when NPT = (N+1)*(N+2)/2, we can set IJ(:, 2*N + 2 : NPT) to
    ! ANY permutation of {(I, J) : 1 <= J < I <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the
    ! first NPT - (2*N + 1) elements of such a permutation. Powell took the following permutation:
    itemp = int((k - n - 2) / n, kind(itemp))
    j = int(k - (itemp + 1) * n - 1, kind(j))
    i = j + itemp
    if (i > n) then
        itemp = j
        j = i - n
        i = itemp
    end if
    ij(:, k) = [i, j]  ! MATLAB code: IJ(:, K) = [I; J];

    ! The following lines set XPT(;, K) to
    ! XPT(:, I + 1) or XPT(:, I + N + 1)
    ! +
    ! XPT(:, J + 1) or XPT(:, J + N + 1),
    ! depending on the values of FVAL(I + 1), FVAL(I + N + 1), FVAL(J + 1), and FVAL(J + N + 1).
    ! It is the only places where the definition of XPT(:, 2*N + 2 : NPT) relies on F(2 : 2*N + 1).
    ! If we set XPT(:, K) to XPT(:, I + 1) + XPT(:, J + 1) regardless of FVAL, then the evaluations
    ! of FVAL(1 : NPT) can be merged, and they are totally parallelizable; this can be beneficial if
    ! the function evaluations are expensive, which is likely the case.
    if (fval(i + 1) <= fval(i + n + 1)) then
        xpt(i, k) = xpt(i, i + 1)
    else
        xpt(i, k) = xpt(i, i + n + 1)
    end if
    if (fval(j + 1) <= fval(j + n + 1)) then
        xpt(j, k) = xpt(j, j + 1)
    else
        xpt(j, k) = xpt(j, j + n + 1)
    end if
end do

! Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except for FMSSG.
do k = int(2 * n + 2, kind(k)), npt_revised
    x = xpt(:, k) + xbase
    call evalf(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call fmssg(iprint, k, f, x, solver)
    ! Save X and F into the history.
    call savehist(k, f, x, fhist, xhist)
    ! Check whether to exit.
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do

! Set NF, KOPT
nf = int(count(evaluated), kind(nf))
kopt = int(minloc(fval, dim=1, mask=evaluated), kind(kopt))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(all(ij(:, 2 * n + 2:npt) >= 1 .and. ij(:, 2 * n + 2:npt) <= n), '1<=IJ<=N', srname)
    call assert(all(ij(1, 2 * n + 2:npt_revised) > ij(2, 2 * n + 2:npt_revised)), &
        & 'IJ(1, 2N+2:END) > IJ(2, 2N+2:END)', srname)
    call assert(nf <= npt, 'NF <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= nf, '1 <= KOPT <= NF', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

end subroutine initxf


subroutine initq(ij, fval, xpt, gq, hq, pq, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes the quadratic model, which is represented by (GQ, HQ, PQ) so that its
! gradient at XBASE is GQ; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use consts_mod, only : RP, IK, ZERO, HALF, DEBUGGING
use debug_mod, only : assert
use infnan_mod, only : is_nan, is_finite
use info_mod, only : INFO_DFT, NAN_MODEL
use linalg_mod, only : issymmetric

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :)     ! IJ(2, NPT)
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
real(RP) :: fi
real(RP) :: fj
real(RP) :: rhobeg
real(RP) :: xi
real(RP) :: xj

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(size(fval) == npt, 'SIZE(FVAL) == NPT', srname)
    call assert(size(ij, 1) == 2 .and. size(ij, 2) == npt, 'SIZE(IJ) == [2, NPT]', srname)
    call assert(all(ij(:, 2 * n + 2:npt) >= 1 .and. ij(:, 2 * n + 2:npt) <= n), '1<=IJ<=N', srname)
    call assert(size(gq) == n, 'SIZE(GQ) == N', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) == [n, n]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
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
do k = int(2 * n + 2, kind(k)), npt
    ! I, J, XI, and XJ will be used below.
    i = ij(1, k)
    j = ij(2, k)
    xi = xpt(i, k)
    xj = xpt(j, k)
    if (xi * xpt(i, i + 1) > ZERO) then
        fi = fval(i + 1)
    else
        fi = fval(i + n + 1)
    end if
    if (xj * xpt(j, j + 1) > ZERO) then
        fj = fval(j + 1)
    else
        fj = fval(j + n + 1)
    end if
    ! With the XI, XJ, FI, and FJ found above, we have
    ! FVAL(K) = F(XBASE + XI + XJ), FI = F(XBASE + XI), FJ = F(XBASE + XJ).
    ! Thus the HQ(I,J) defined below approximates frac{partial^2}{partial X_I partial X_J} F(XBASE).
    hq(i, j) = (fbeg - fi - fj + fval(k)) / (xi * xj)
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
    call assert(issymmetric(hq), 'HQ is symmetric', srname)
end if

end subroutine initq


subroutine inith(ij, xpt, idz, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes IDZ, BMAT, and ZMAT.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, HALF, DEBUGGING
use debug_mod, only : assert
use infnan_mod, only : is_nan, is_finite
use info_mod, only : INFO_DFT, NAN_MODEL
use linalg_mod, only : issymmetric, eye

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :)     ! IJ(2, NPT)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK), intent(out) :: idz
integer(IK), intent(out) :: info
real(RP), intent(out) :: bmat(:, :)     ! BMAT(N, NPT + N)
real(RP), intent(out) :: zmat(:, :)     ! ZMAT(NPT, NPT - N - 1)

! Local variables
character(len=*), parameter :: srname = 'INITH'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: recip
real(RP) :: reciq
real(RP) :: rhobeg

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(size(ij, 1) == 2 .and. size(ij, 2) == npt, 'SIZE(IJ) == [2, NPT]', srname)
    call assert(all(ij(:, 2 * n + 2:npt) >= 1 .and. ij(:, 2 * n + 2:npt) <= n), '1<=IJ<=N', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! Set IDZ = 1. It will not be changed in the following.
idz = 1_IK

! Some values to be used for setting BMAT and ZMAT.
rhobeg = maxval(abs(xpt(:, 2)))  ! Read RHOBEG from XPT.

! Set BMAT.
bmat = ZERO
if (npt <= 2 * n + 1) then
    bmat(1:npt - n - 1, 2:npt - n) = (HALF / rhobeg) * eye(npt - n - 1)
    bmat(1:npt - n - 1, n + 2:npt) = -(HALF / rhobeg) * eye(npt - n - 1)
    bmat(npt - n:n, 1) = -ONE / rhobeg
    bmat(npt - n:n, npt - n + 1:n + 1) = (ONE / rhobeg) * eye(2 * n - npt + 1)
    bmat(npt - n:n, 2 * npt - n:npt + n) = -(HALF * rhobeg**2) * eye(2 * n - npt + 1)
else
    bmat(:, 2:n + 1) = (HALF / rhobeg) * eye(n)
    bmat(:, n + 2:2 * n + 1) = -(HALF / rhobeg) * eye(n)
end if

! Set ZMAT.
recip = ONE / rhobeg**2
reciq = sqrt(HALF) / rhobeg**2
zmat = ZERO
if (npt <= 2 * n + 1) then
    zmat(1, :) = -reciq - reciq
    zmat(2:npt - n, :) = reciq * eye(npt - n - 1)
    zmat(n + 2:npt, :) = reciq * eye(npt - n - 1)
else
    zmat(1, 1:n) = -reciq - reciq
    zmat(2:n + 1, 1:n) = reciq * eye(n)
    zmat(n + 2:2 * n + 1, 1:n) = reciq * eye(n)
    zmat(1, n + 1:npt - n - 1) = recip
    zmat(2 * n + 2:npt, n + 1:npt - n - 1) = recip * eye(npt - 2 * n - 1)
    do k = int(2 * n + 2, kind(k)), int(npt, kind(k))
        where (xpt(ij(:, k), k) > ZERO)
            zmat(ij(:, k) + 1, k - n - 1) = -recip
        elsewhere
            zmat(ij(:, k) + n + 1, k - n - 1) = -recip
        end where
        ! MATLAB code:
        ! L = (XPT(IJ(:, K), K) > ZERO);
        ! ZMAT(IJ(L, K) + 1, K - N - 1) = -RECIP;
        ! ZMAT(IJ(~L, K) + N + 1, K - N - 1) = -RECIP;
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
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT - N', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
end if

end subroutine inith


end module initialize_mod
