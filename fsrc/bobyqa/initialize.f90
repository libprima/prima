module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of BOBYQA, described in Section 2 of the BOBYQA paper.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, June 13, 2022 PM05:39:46
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initxf, initq, inith


contains


subroutine initxf(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x0, ij, kopt, nf, fhist, fval, &
    & sl, su, xbase, xhist, xpt, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine does the initialization about the interpolation points & their function values.
!
! N.B.:
! 1. Remark on IJ:
! If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [NPT-2*N-1, 2].
! IJ contains integers between 1 and N. For each K > 2*N + 1, XPT(:, K) is
! XPT(:, IJ(K, 1) + 1) + XPT(:, IJ(K, 2) + 1). The 1 in IJ + 1 comes from the fact that XPT(:, 1)
! corresponds to the base point XBASE. Let I = IJ(K, 1) and J = IJ(K, 2). Then all the
! entries of XPT(:, K) are zero except for the I and J entries. Consequently, the Hessian of the
! quadratic model will get a possibly nonzero (I, J) entry.
! 2. At return,
! INFO = INFO_DFT: initialization finishes normally
! INFO = FTARGET_ACHIEVED: return because F <= FTARGET
! INFO = NAN_INF_X: return because X contains NaN
! INFO = NAN_INF_F: return because F is either NaN or +Inf
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : trueloc
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
integer(IK), intent(out) :: ij(:, :)    ! IJ(MAX(0_IK, NPT-2*N-1_IK), 2)
integer(IK), intent(out) :: info
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: fval(:)  ! FVAL(NPT)
real(RP), intent(out) :: sl(:)  ! SL(N)
real(RP), intent(out) :: su(:)  ! SU(N)
real(RP), intent(out) :: xbase(:)  ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)

! Local variables
character(len=*), parameter :: solver = 'BOBYQA'
character(len=*), parameter :: srname = 'INITIIALIZE'
integer(IK) :: k
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: subinfo
logical :: evaluated(size(xpt, 2))
real(RP) :: f
real(RP) :: x(size(xpt, 1))

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
    call assert(all(x0 >= xl .and. (x0 <= xl .or. x0 >= xl + rhobeg)), 'X0 == XL or X0 >= XL + RHOBEG', srname)
    call assert(all(x0 <= xu .and. (x0 >= xu .or. x0 <= xu - rhobeg)), 'X0 == XU or X0 >= XU - RHOBEG', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize INFO to the default value. At return, an INFO different from this value will indicate
! an abnormal return.
info = INFO_DFT

! SL and SU are the lower and upper bounds on feasible moves from X0.
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

! Initialize XBASE to X0.
xbase = x0

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously.
evaluated = .false.
! Initialize FVAL to HUGENUM. Otherwise, compilers may complain that FVAL is not (completely)
! initialized if the initialization aborts due to abnormality (see CHECKEXIT).
fval = HUGENUM

! Set XPT(:, 2 : N+1)
xpt = ZERO
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
! have different signs and FVAL(K) <= FVAL(K+N). This provides a bias towards potentially lower
! values of F when defining XPT(:, 2*N + 2 : NPT). We may drop the requirement on the signs, but
! Powell's code has such a requirement.
! N.B.:
! 1. The switching is OPTIONAL. It we remove it, then the evaluations of FVAL(1 : NPT) can be
! merged, and they are totally PARALLELIZABLE; this can be beneficial if the function evaluations
! are expensive, which is likely the case.
! 2. The initialization of NEWUOA revises IJ (see below) instead of XPT and FVAL. Mathematically,
! it is equivalent; practically, after XPT is revised, the initialization of the quadratic model and
! the Lagrange polynomials needs revision as well, as, e.g., XPT(:, 2:N) is not RHOBEG*EYE(N) anymore.
do k = 2, min(npt - n, int(n + 1, kind(npt)))
    if (xpt(k - 1, k) * xpt(k - 1, k + n) < 0 .and. fval(k + n) < fval(k)) then
        fval([k, k + n]) = fval([k + n, k])
        xpt(:, [k, k + n]) = xpt(:, [k + n, k])
        ! Indeed, only XPT(K-1, [K, K+N]) needs switching, as the other entries are zero.
    end if
end do

! Set IJ.
! In general, when NPT = (N+1)*(N+2)/2, we can set IJ(1 : NPT - (2*N+1), :) to ANY permutation
! of {{I, J} : 1 <= I /= J <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first NPT - (2*N+1)
! elements of such a permutation. The following IJ is defined according to Powell's code. See also
! Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
! N.B.: We do not distinguish between {I, J} and {J, I}, which represent the same set. If we want to
! ensure an order, e.g., IJ(:, 1) > IJ(:, 2) (so that the (IJ(K, 1), IJ(K, 2)) position is in the
! lower triangular part of a matrix), then we can sort IJ, e.g., by IJ = SORT(IJ, 2, 'DESCEND').
ij(:, 1) = int([(k, k=n, npt - n - 2_IK)] / n, IK)
ij(:, 2) = int([(k, k=n, npt - n - 2_IK)] - n * ij(:, 1) + 1_IK, IK)
ij(:, 1) = modulo(ij(:, 1) + ij(:, 2) - 1_IK, n) + 1_IK  ! MODULO(K-1,N) + 1 = K-N for K in [N+1,2N]
!!MATLAB: (N.B.: Fortran MODULO == MATLAB `mod`, Fortran MOD == MATLAB `rem`)
!!ij(:, 1) = floor((n : npt-n-2) / n);
!!ij(:, 2) = (n : npt-n-2) - n*ij(:, 1) + 1;
!!ij(:, 1) = mod(ij(:, 1) + ij(:, 2) - 1, n) + 1;  % mod(k-1,n) + 1 = k-n for k in [n+1,2n]

! Set XPT(:, 2*N + 2 : NPT). It depends on XPT(:, 1 : 2*N + 1) and hence on FVAL(1: 2*N + 1).
! Indeed, XPT(:, K) has only two nonzeros for each K >= 2*N+2.
! N.B.: The 1 in IJ + 1 comes from the fact that XPT(:, 1) corresponds to XBASE.
xpt(:, 2 * n + 2:npt) = xpt(:, ij(:, 1) + 1) + xpt(:, ij(:, 2) + 1)

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

! Set NF, KOPT
nf = int(count(evaluated), kind(nf))
kopt = int(minloc(fval, mask=evaluated, dim=1), kind(kopt))
!!MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nf <= npt, 'NF <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(all(ij >= 1 .and. ij <= n), '1 <= IJ <= N', srname)
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
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
end if

end subroutine initxf


subroutine initq(ij, fval, xpt, gopt, hq, pq, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes the quadratic model, which is represented by (GOPT, HQ, PQ) so that
! its gradient at XOPT is GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'. Here,
! XOPT = XBASE + XPT(:, KOPT).
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL
use, non_intrinsic :: linalg_mod, only : issymmetric, diag, matprod

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :)     ! IJ(NPT, 2)
real(RP), intent(in) :: fval(:)     ! FVAL(NPT)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK), intent(out), optional :: info
real(RP), intent(out) :: gopt(:)  ! GOPT(N)
real(RP), intent(out) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)

! Local variables
character(len=*), parameter :: srname = 'INITQ'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: k
integer(IK) :: kopt
integer(IK) :: n
integer(IK) :: ndiag
integer(IK) :: npt
real(RP) :: fbase
real(RP) :: xa(min(size(xpt, 1), size(xpt, 2) - size(xpt, 1) - 1))
real(RP) :: xb(size(xa))
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
    call assert(all(ij >= 1 .and. ij <= n), '1 <= IJ <= N', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) = N', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) = [N, N]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

fbase = fval(1)  ! FBASE is the function value at XBASE.

! Set GOPT by the forward difference.
gopt = (fval(2:n + 1) - fbase) / diag(xpt(:, 2:n + 1))

! The interpolation conditions decide GOPT(1:NDIAG) and the first NDIAG diagonal 2nd derivatives of
! the initial quadratic model by a quadratic interpolation on three points.
ndiag = min(n, npt - n - 1_IK)
xa = diag(xpt(:, 2:ndiag + 1))
xb = diag(xpt(:, n + 2:n + ndiag + 1))

! Revise GOPT(1:NDIAG) to the value provided by the three-point interpolation.
gopt(1:ndiag) = (gopt(1:ndiag) * xb - ((fval(n + 2:n + ndiag + 1) - fbase) / xb) * xa) / (xb - xa)

! Set the diagonal of HQ by the three-point interpolation. If we do this before the revision of
! GQ(1:NDIAG), we can avoid the calculation of FVAL(K + 1) - FBASE) / RHOBEG. But we prefer to
! decouple the initialization of GQ and HQ. We are not concerned by this amount of flops.
hq = ZERO
do k = 1, ndiag
    hq(k, k) = TWO * ((fval(k + 1) - fbase) / xa(k) - (fval(n + k + 1) - fbase) / xb(k)) / (xa(k) - xb(k))
end do
!!MATLAB:
!!hdiag = 2*((fval(2 : ndiag+1) - fbase) / xa - (fval(n+2 : n+ndiag+1) - fbase) / xb) / (xa-xb)
!!hq(1:ndiag, 1:ndiag) = diag(hdiag)

! When NPT > 2*N + 1, set the off-diagonal entries of HQ.
do k = 1, npt - 2_IK * n - 1_IK
    i = ij(k, 1)
    j = ij(k, 2)
    xi = xpt(i, k + 2 * n + 1)
    xj = xpt(j, k + 2 * n + 1)
    ! N.B.: The 1 in I+1 and J+1 comes from the fact that XPT(:, 1) corresponds to XBASE.
    hq(i, j) = (fbase - fval(i + 1) - fval(j + 1) + fval(k + 2 * n + 1)) / (xi * xj)
    hq(j, i) = hq(i, j)
end do

kopt = int(minloc(fval, dim=1), kind(kopt))
if (kopt /= 1) then
    gopt = gopt + matprod(hq, xpt(:, kopt))
end if

pq = ZERO

if (present(info)) then
    if (is_nan(sum(abs(gopt)) + sum(abs(hq)))) then
        info = NAN_MODEL
    else
        info = INFO_DFT
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(gopt) == n, 'SIZE(GOPT) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
end if

end subroutine initq


subroutine inith(ij, xpt, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes BMAT and ZMAT.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): NONE
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL
use, non_intrinsic :: linalg_mod, only : issymmetric, diag

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :) ! IJ(NPT, 2)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK), intent(out), optional :: info
real(RP), intent(out) :: bmat(:, :) ! BMAT(N, NPT + N)
real(RP), intent(out) :: zmat(:, :) ! ZMAT(NPT, NPT - N - 1)

! Local variables
character(len=*), parameter :: srname = 'INITH'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: ndiag
integer(IK) :: npt
real(RP) :: rhobeg
real(RP) :: rhosq
real(RP) :: xa(min(size(xpt, 1), size(xpt, 2) - size(xpt, 1) - 1))
real(RP) :: xb(size(xa))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(all(ij >= 1 .and. ij <= n), '1 <= IJ <= N', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! Some values to be used for setting BMAT and ZMAT.
rhobeg = maxval(abs(xpt(:, 2)))  ! Read RHOBEG from XPT.
rhosq = rhobeg**2

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
do k = ndiag + 1_IK, n
    bmat(k, 1) = -ONE / xpt(k, k + 1)
    bmat(k, k + 1) = -bmat(k, 1)
    bmat(k, npt + k) = -HALF * rhosq
end do

zmat = ZERO
! Set ZMAT(:, 1 : NDIAG)
zmat(1, 1:ndiag) = sqrt(TWO) / (xa * xb)
do k = 1, ndiag
    zmat(k + 1, k) = -zmat(1, k) - sqrt(HALF) / rhosq
    zmat(k + n + 1, k) = sqrt(HALF) / rhosq
end do
! Set ZMAT(:, NDIAG+1 : NPT-N-1)
do k = ndiag + 1_IK, npt - n - 1_IK
    zmat(1, k) = ONE / rhosq
    zmat(k + n + 1, k) = ONE / rhosq
    zmat(ij(k - n, :) + 1, k) = -ONE / rhosq
end do

if (present(info)) then
    if (is_nan(sum(abs(bmat)) + sum(abs(zmat)))) then
        info = NAN_MODEL
    else
        info = INFO_DFT
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

end subroutine inith


end module initialize_mod
