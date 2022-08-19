module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, June 16, 2022 AM12:00:45
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initxf, inith


contains


subroutine initxf(calfun, iprint, maxfun, A_orig, amat, b_orig, ctol, ftarget, rhobeg, x0, b, &
    & ij, kopt, nf, chist, fhist, fval, xbase, xhist, xpt, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine does the initialization about the interpolation points & their function values.
!
! N.B.:
! 1. Remark on IJ:
! If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [2, NPT-2*N-1].
! IJ contains integers between 1 and N. For each K > 2*N + 1, XPT(:, K) is
! XPT(:, IJ(1, K) + 1) + XPT(:, IJ(2, K) + 1). The 1 in IJ + 1 comes from the fact that XPT(:, 1)
! corresponds to the base point XBASE. Let I = IJ(1, K) and J = IJ(2, K). Then all the entries of
! XPT(:, K) are zero except for the I and J entries. Consequently, the Hessian of the quadratic
! model will get a possibly nonzero (I, J) entry.
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
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : matprod, maximum, eye
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : setij

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: A_orig(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: b_orig(:)  ! B_ORIG(M)
real(RP), intent(in) :: ctol
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)  ! X0(N)

! In-outputs
real(RP), intent(inout) :: b(:)  ! B(M)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: ij(:, :)  ! IJ(2, MAX(0_IK, NPT-2*N-1))
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)  ! CHIST(MAXCHIST)
real(RP), intent(out) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out) :: fval(:)  ! FVAL(NPT)
real(RP), intent(out) :: xbase(:)  ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)

! Local variables
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'INITXF'
integer(IK) :: j
integer(IK) :: k
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: subinfo
logical :: evaluated(size(xpt, 2))
logical :: feasible(size(xpt, 2))
real(RP) :: constr(size(b))
real(RP) :: cstrv
real(RP) :: f
real(RP) :: mincv
real(RP) :: x(size(x0))

! Sizes.
m = int(size(b), kind(m))
n = int(size(x), kind(n))
npt = int(size(fval), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = int(max(maxxhist, maxfhist, maxchist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(size(A_orig, 1) == n .and. size(A_orig, 2) == m, 'SIZE(A_ORIG) == [N, M]', srname)
    call assert(size(b_orig) == m, 'SIZE(B_ORIG) == M', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize INFO to the default value. At return, an INFO different from this value will indicate
! an abnormal return.
info = INFO_DFT

! Initialize XBASE to X0.
xbase = x0

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously. However, the loop here
! is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
evaluated = .false.

! Initialize FVAL to HUGENUM. Otherwise, compilers may complain that FVAL is not (completely)
! initialized if the initialization aborts due to abnormality (see CHECKEXIT).
fval = HUGENUM

! Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT], but they may be altered
! later to make a constraint violation sufficiently large.
xpt(:, 1) = ZERO
xpt(:, 2:n + 1) = rhobeg * eye(n)
xpt(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1_IK)  ! XPT(:, 2*N+2 : NPT) = ZERO if it is nonempty.

! Set IJ.
! In general, when NPT = (N+1)*(N+2)/2, we can set IJ(:, 1 : NPT - (2*N+1)) to ANY permutation
! of {{I, J} : 1 <= I /= J <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first NPT - (2*N+1)
! elements of such a permutation. The following IJ is defined according to Powell's code. See also
! Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
ij = setij(n, npt)

! Set XPT(:, 2*N + 2 : NPT).
! Indeed, XPT(:, K) has only two nonzeros for each K >= 2*N + 2,
! N.B.: The 1 in IJ + 1 comes from the fact that XPT(:, 1) corresponds to XBASE.
xpt(:, 2 * n + 2:npt) = xpt(:, ij(1, :) + 1) + xpt(:, ij(2, :) + 1)

! Update the constraint right-hand sides to allow for the shift XBASE.
b = b - matprod(xbase, amat)

! Go through the initial points, shifting every infeasible point if necessary so that its constraint
! violation is at least 0.2*RHOBEG. This is OPTIONAL. According to a test on 20220614, it seems to
! improve the performance of LINCOA modestly.
mincv = 0.2_RP * rhobeg
do k = 2, npt  ! LINCOA always start with a feasible point. So we do this only for K >= 2.
    ! Internally, we use AMAT and B to evaluate the constraints.
    constr = matprod(xpt(:, k), amat) - b
    if (all(constr < mincv) .and. any(constr > 0)) then
        j = int(maxloc(constr, dim=1), IK)
        xpt(:, k) = xpt(:, k) + (mincv - constr(j)) * amat(:, j)
    end if
end do

! Set FVAL by evaluating F. Totally parallelizable except for FMSG.
do k = 1, npt
    x = xbase + xpt(:, k)
    call evaluate(calfun, x, f)
    ! For the output, we use A_ORIG and B_ORIG to evaluate the constraints.
    constr = matprod(x, A_orig) - b_orig
    cstrv = maximum([ZERO, constr])
    call fmsg(solver, iprint, k, f, x, cstrv, constr)
    call savehist(k, x, xhist, f, fhist, cstrv, chist)
    evaluated(k) = .true.
    feasible(k) = all(constr <= 0)
    fval(k) = f
    subinfo = checkexit(maxfun, k, cstrv, ctol, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do
feasible(1) = .true.  ! LINCOA always start with a feasible point.

nf = int(count(evaluated), kind(nf))
kopt = int(minloc(fval, mask=(evaluated .and. feasible), dim=1), kind(kopt))
!!MATLAB:
!!fopt = min(fval(evaluated & feasible));
!!kopt = find(evaluated & feasible & ~(fval > fopt), 1, 'first');

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(ij, 1) == 2 .and. size(ij, 2) == max(0_IK, npt - 2_IK * n - 1_IK), &
        & 'SIZE(IJ) == [2, NPT - 2*N - 1]', srname)
    call assert(all(ij >= 1 .and. ij <= 2 * n), '1 <= IJ <= 2*N', srname)
    call assert(all(ij(1, :) /= ij(2, :)), 'IJ(1, :) /= IJ(:, 2)', srname)
    call assert(nf <= npt, 'NF <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= nf, '1 <= KOPT <= NF', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(evaluated .and. (is_nan(fval) .or. is_posinf(fval))), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(.not. any(evaluated .and. feasible .and. fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
end if

end subroutine initxf


subroutine inith(ij, rhobeg, xpt, idz, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes [IDZ, BMAT, ZMAT] which represents the matrix H in (3.12) of the
! NEWUOA paper (see also (2.7) of the BOBYQA paper).
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL
use, non_intrinsic :: linalg_mod, only : eye, issymmetric
use, non_intrinsic :: powalg_mod, only : updateh

implicit none

! Inputs
integer(IK), intent(in) :: ij(:, :)    ! IJ(2, MAX(0, NPT-2*N-1))
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! Outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out) :: idz
real(RP), intent(out) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(out) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Local variables
character(len=*), parameter :: srname = 'INITH'
integer(IK) :: k
integer(IK) :: kbase
integer(IK) :: n
integer(IK) :: npt
real(RP) :: recip
real(RP) :: reciq
real(RP) :: rhosq
real(RP) :: xpt_ref(size(xpt, 1), size(xpt, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(ij, 1) == 2 .and. size(ij, 2) == max(0_IK, npt - 2_IK * n - 1_IK), &
        & 'SIZE(IJ) == [2, NPT - 2*N - 1]', srname)
    call assert(all(ij >= 1 .and. ij <= 2 * n), '1 <= IJ <= 2*N', srname)
    call assert(all(ij(1, :) /= ij(2, :)), 'IJ(1, :) /= IJ(2, :)', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

rhosq = rhobeg**2

! Set BMAT corresponding to the XPT_REF that will be defined later.
recip = ONE / rhobeg
reciq = HALF / rhobeg
bmat = ZERO
if (npt <= 2 * n + 1) then
    ! Set BMAT(1 : NPT-N-1, :)
    bmat(1:npt - n - 1, 2:npt - n) = reciq * eye(npt - n - 1_IK)
    bmat(1:npt - n - 1, n + 2:npt) = -reciq * eye(npt - n - 1_IK)
    ! Set BMAT(NPT-N : N, :)
    bmat(npt - n:n, 1) = -recip
    bmat(npt - n:n, npt - n + 1:n + 1) = recip * eye(2_IK * n - npt + 1_IK)
    bmat(npt - n:n, 2 * npt - n:npt + n) = -(HALF * rhosq) * eye(2_IK * n - npt + 1_IK)
else
    bmat(:, 2:n + 1) = reciq * eye(n)
    bmat(:, n + 2:2 * n + 1) = -reciq * eye(n)
end if

! Set ZMAT corresponding to the XPT_REF that will be defined later.
recip = ONE / rhosq
reciq = sqrt(HALF) / rhosq
zmat = ZERO
if (npt <= 2 * n + 1) then
    zmat(1, :) = -reciq - reciq
    zmat(2:npt - n, :) = reciq * eye(npt - n - 1_IK)
    zmat(n + 2:npt, :) = reciq * eye(npt - n - 1_IK)
else
    ! Set ZMAT(:, 1:N).
    zmat(1, 1:n) = -reciq - reciq
    zmat(2:n + 1, 1:n) = reciq * eye(n)
    zmat(n + 2:2 * n + 1, 1:n) = reciq * eye(n)
    ! Set ZMAT(:, N+1 : NPT-N-1).
    zmat(1, n + 1:npt - n - 1) = recip
    zmat(2 * n + 2:npt, n + 1:npt - n - 1) = recip * eye(npt - 2_IK * n - 1_IK)
    do k = 1, npt - 2_IK * n - 1_IK
        zmat(ij(:, k) + 1, k + n) = -recip
    end do
end if

! Set IDZ = 1.
idz = 1_IK

! Up to now, [BMAT, ZMAT, IDZ] corresponds to the XPT_REF defined below.
xpt_ref(:, 1) = ZERO
xpt_ref(:, 2:n + 1) = rhobeg * eye(n)
xpt_ref(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1_IK)
xpt_ref(:, 2 * n + 2:npt) = xpt_ref(:, ij(1, :) + 1) + xpt_ref(:, ij(2, :) + 1)

! Update [BMAT, ZMAT, IDZ] so that it corresponds to XPT.
kbase = 1_IK
do k = 1, npt
    if (all(abs(xpt(:, k) - xpt_ref(:, k)) <= 0)) then  ! XPT(:, K) == XPT_REF(:, K)
        cycle
    end if
    call updateh(k, kbase, idz, xpt(:, k), xpt_ref, bmat, zmat)
    xpt_ref(:, k) = xpt(:, k)
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
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

end subroutine inith


end module initialize_mod
