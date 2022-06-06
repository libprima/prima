module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Tuesday, June 07, 2022 AM12:04:42
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initialize


contains


subroutine initialize(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x0, &
    & kopt, nf, bmat, f, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, zmat, info)

use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : matprod, trueloc, issymmetric, diag, eye, sort
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: xl(:)  ! XL(N)
real(RP), intent(in) :: xu(:)  ! XU(N)

! In-outputs
real(RP), intent(inout) :: x0(:)  ! X(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: fval(:)  ! FVAL(NPT)
real(RP), intent(out) :: gopt(:)  ! GOPT(N)
real(RP), intent(out) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)
real(RP), intent(out) :: sl(:)  ! SL(N)
real(RP), intent(out) :: su(:)  ! SU(N)
real(RP), intent(out) :: xbase(:)  ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(out) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)


! Local variables
character(len=*), parameter :: solver = 'BOBYQA'
character(len=*), parameter :: srname = 'INITIIALIZE'
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: subinfo
real(RP) :: x(size(xpt, 1)), xa(min(size(xpt, 1), size(xpt, 2) - size(xpt, 1) - 1)), xb(size(xa))
real(RP) :: fbase, rhosq, stepa, stepb, xi, xj
integer(IK) :: ij(max(0, size(xpt, 2) - 2 * size(xpt, 1) - 1), 2), itemp, k, i, j, ndiag
logical :: evaluated(size(xpt, 2))

! Sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(size(fval) == npt, 'SIZE(FVAL) == NPT', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(size(x0) == n .and. all(is_finite(x0)), 'SIZE(X0) == N, X0 is finite', srname)
    call assert(all(x0 >= xl .and. (x0 <= xl .or. x0 >= xl + rhobeg)), 'X0 == XL or X0 >= XL + RHOBEG', solver)
    call assert(all(x0 <= xu .and. (x0 >= xu .or. x0 <= xu - rhobeg)), 'X0 == XU or X0 >= XU - RHOBEG', solver)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) == [N, N]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
end if

!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!

! Initialize INFO to the default value. At return, an INFO different from this value will indicate
! an abnormal return.
info = INFO_DFT

! The lower and upper bounds on moves from X0 are set now, in order to provide useful and exact
! information about components of X0 that become within distance RHOBEG from their bounds.
sl = xl - x0
su = xu - x0
! After the preprocessing subroutine PREPROC, SL <= 0 and the nonzero entries of SL should be less
! than -RHOBEG, while SU >= 0 and the nonzeros of SU should be larger than RHOBEG. However, this may
! not be true due to rounding. The following lines revise SL and SU to ensure it. X0 is also revised
! accordingly. In precise arithmetic, the "revisions" do not change SL, SU, or X0.
where (sl < 0)
    sl = min(sl, -rhobeg)
elsewhere
    x0 = xl
    sl = ZERO
    su = xu - xl
end where
where (su > 0)
    su = max(su, rhobeg)
elsewhere
    x0 = xu
    sl = xl - xu
    su = ZERO
end where
!!MATLAB code for revising X, SL, and SU:
!!sl(sl < 0) = min(sl(sl < 0), -rhobeg);
!!x0(sl >= 0) = xl(sl >= 0);
!!sl(sl >= 0) = 0;
!!su(sl >= 0) = xu(sl >= 0) - xl(sl >= 0);
!!su(su > 0) = max(su(su > 0), rhobeg);
!!x0(su <= 0) = xu(su <= 0);
!!sl(su <= 0) = xl(su <= 0) - xu(su <= 0);
!!su(su <= 0) = 0;

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously.
evaluated = .false.
! Initialize FVAL to HUGENUM. Otherwise, compilers may complain that FVAL is not (completely)
! initialized if the initialization aborts due to abnormality (see CHECKEXIT).
fval = HUGENUM

! Set XBASE to the initial vector of variables.
xbase = x0
hq = ZERO
pq = ZERO
bmat = ZERO
zmat = ZERO

xpt = ZERO
! Set XPT(:, 2 : N+1)
do k = 1, n
    xpt(k, k + 1) = rhobeg
    if (su(k) <= 0) then  ! SU(K) == 0
        xpt(k, k + 1) = -rhobeg
    end if
end do
! Set XPT(:, N+2 : MIN(2*N + 1, NPT)).
do k = 1, min(npt - n - 1_IK, n)
    xpt(k, k + n + 1) = -rhobeg
    if (sl(k) >= 0) then  ! SL(K) == 0
        xpt(k, k + n + 1) = min(TWO * rhobeg, su(k))
    end if
    if (su(k) <= 0) then  ! SU(K) == 0
        xpt(k, k + n + 1) = max(-TWO * rhobeg, sl(k))
    end if
end do

! Set FVAL(1 : MIN(2*N + 1, NPT)) by evaluating F. Totally parallelizable except for FMSG.
do k = 1, min(npt, int(2 * n + 1, kind(npt)))
    x = min(max(xl, xbase + xpt(:, k)), xu)
    x(trueloc(xpt(:, k) <= sl)) = xl(trueloc(xpt(:, k) <= sl))
    x(trueloc(xpt(:, k) >= su)) = xu(trueloc(xpt(:, k) >= su))
    call evaluate(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call savehist(k, x, xhist, f, fhist)
    call fmsg(solver, iprint, k, f, x)
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do

! For the K between 2 and N + 1, switch XPT(:, K) and XPT(:, K+1) if XPT(K-1, K) and XPT(K-1, K+N)
! have different signs and FVAL(K) <= FVAL(K+N). This is OPTIONAL. It provides a bias towards
! potentially lower values of F when defining XPT(:, 2*N + 2 : NPT). We may drop the requirement on
! the signs, but Powell's code has such a requirement.
do k = 2, min(npt - n, int(n + 1, kind(npt)))
    if (xpt(k - 1, k) * xpt(k - 1, k + n) < 0 .and. fval(k + n) < fval(k)) then
        fval([k, k + n]) = fval([k + n, k])
        xpt(:, [k, k + n]) = xpt(:, [k + n, k])
        ! Indeed, only XPT(K-1, [K, K+N]) needs switching, as the other entries are zero.
    end if
end do

! Set IJ.
! In general, when NPT = (N+1)*(N+2)/2, we can initialize IJ(1 : NPT - (2*N+1), :) to ANY permutation
! of {(I, J) : 1 <= J < I <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first
! NPT - (2*N + 1) elements of such a permutation. Powell took the following one.
ij(:, 1) = int([(k, k=n, npt - n - 2_IK)] / n, IK)
!!MATLAB: ij(:, 1) = floor((n : npt - n - 2) / n);
ij(:, 2) = int([(k, k=n, npt - n - 2_IK)] - n * ij(:, 1) + 1_IK, IK)
!!MATLAB: ij(:, 2) = (n : npt-n-2) - n*ij(:, 1) + 1
ij(:, 1) = modulo(ij(:, 1) + ij(:, 2) - 1_IK, n) + 1_IK  ! MODULO(K-1,N) + 1 = K-N for K in [N+1,2N]
! The next line ensures IJ(:, 1) > IJ(:, 2).
ij = sort(ij, 2, 'descend')
! Increment IJ by 1. This 1 comes from the fact that XPT(:, 1) corresponds to the base point XBASE.
ij = ij + 1_IK

! Set XPT(:, 2*N + 2 : NPT). Indeed, XPT(:, K) has only two nonzeros for each K >= 2*N+2.
xpt(:, 2 * n + 2:npt) = xpt(:, ij(:, 1)) + xpt(:, ij(:, 2))

! Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except for FMSG.
if (info == INFO_DFT) then
    do k = int(2 * n + 2, kind(k)), npt
        x = min(max(xl, xbase + xpt(:, k)), xu)
        x(trueloc(xpt(:, k) <= sl)) = xl(trueloc(xpt(:, k) <= sl))
        x(trueloc(xpt(:, k) >= su)) = xu(trueloc(xpt(:, k) >= su))
        call evaluate(calfun, x, f)
        evaluated(k) = .true.
        fval(k) = f
        call savehist(k, x, xhist, f, fhist)
        call fmsg(solver, iprint, k, f, x)
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end do
end if

nf = count(evaluated)
kopt = int(minloc(fval, mask=evaluated, dim=1), kind(kopt))

! Set the nonzero initial elements of BMAT and the quadratic model in the cases when NF is at most
! 2*N+1. If NF exceeds N+1, then the positions of the NF-th and (NF-N)-th interpolation points may
! be switched, in order that the function value at the first of them contributes to the off-diagonal
! second derivative terms of the initial quadratic model.
if (all(evaluated)) then
    hq = ZERO
    pq = ZERO
    fbase = fval(1)
    gopt = (fval(2:n + 1) - fbase) / diag(xpt(:, 2:n + 1))
    ! The interpolation conditions decide the first NDIAG diagonal 2nd derivatives of the initial quadratic model.
    ndiag = min(n, npt - n - 1_IK)
    xa = diag(xpt(:, 2:ndiag + 1))
    xb = diag(xpt(:, n + 2:n + ndiag + 1))
    gopt(1:ndiag) = (((fval(2:ndiag + 1) - fbase) / xa) * xb - ((fval(n + 2:n + ndiag + 1) - fbase) / xb) * xa) / (xb - xa)
    do k = 1, ndiag
        hq(k, k) = TWO * ((fval(k + 1) - fbase) / xa(k) - (fval(n + k + 1) - fbase) / xb(k)) / (xa(k) - xb(k))
    end do

    ! Set the off-diagonal second derivatives of the initial quadratic model.
    do k = 1, npt - 2_IK * n - 1_IK
        i = ij(k, 1) - 1
        j = ij(k, 2) - 1
        xi = xpt(i, k + 2 * n + 1)
        xj = xpt(j, k + 2 * n + 1)
        hq(i, j) = (fbase - fval(ij(k, 1)) - fval(ij(k, 2)) + fval(k + 2 * n + 1)) / (xi * xj)
        hq(j, i) = hq(i, j)
    end do
    if (kopt /= 1) then
        gopt = gopt + matprod(hq, xpt(:, kopt))
    end if
end if

if (all(evaluated)) then
    rhosq = rhobeg * rhobeg
    ! The interpolation set decides the first NDIAG diagonal 2nd derivatives of the Lagrange polynomials.
    ndiag = min(n, npt - n - 1_IK)
    xa = diag(xpt(:, 2:ndiag + 1))
    xb = diag(xpt(:, n + 2:n + ndiag + 1))

    bmat = ZERO
    ! Set BMAT(1 : NDIAG, :)
    bmat(1:ndiag, 1) = -(xa + xb) / (xa * xb)
    do k = 1, ndiag
        bmat(k, k + n + 1) = -HALF / xpt(k, k + 1)
        bmat(k, k + 1) = -bmat(k, 1) - bmat(k, k + n + 1)
    end do
    ! Set BMAT(NDIAG+1 : N, :)
    bmat(npt - n:n, 1) = -ONE / diag(xpt(npt - n:n, npt - n + 1:n + 1))
    do k = npt - n, n
        bmat(k, k + 1) = ONE / xpt(k, k + 1)
        bmat(k, npt + k) = -HALF * rhosq
    end do

    zmat = ZERO
    ! Set ZMAT(:, 1 : NDIAG)
    zmat(1, 1:ndiag) = sqrt(TWO) / (xa * xb)
    zmat(n + 2:ndiag + n + 1, 1:ndiag) = (sqrt(HALF) / rhosq) * eye(ndiag)
    do k = 1, ndiag
        zmat(k + 1, k) = -zmat(1, k) - zmat(k + n + 1, k)
    end do
    ! Set ZMAT(:, NDIAG+1 : NPT-N-1)
    zmat(1, n + 1:npt - n - 1) = ONE / rhosq
    zmat(2 * n + 2:npt, n + 1:npt - n - 1) = (ONE / rhosq) * eye(npt - 2_IK * n - 1_IK)
    do k = 1, npt - 2_IK * n - 1_IK
        zmat(ij(k, :), k + n) = -ONE / rhosq
    end do
end if

! Postconditions
if (DEBUGGING) then
    call assert(nf <= npt, 'NF <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(all(xbase >= xl .and. xbase <= xu), 'XL <= XBASE <= XU', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt)) .and. &
        & all(xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    call assert(size(fval) == npt .and. .not. any(evaluated .and. (is_nan(fval) .or. is_posinf(fval))), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(.not. any(fval < fval(kopt) .and. evaluated), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(size(sl) == n .and. all(sl <= 0), 'SIZE(SL) == N, SL <= 0', srname)
    call assert(size(su) == n .and. all(su >= 0), 'SIZE(SU) == N, SU >= 0', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
end if

end subroutine initialize


end module initialize_mod
