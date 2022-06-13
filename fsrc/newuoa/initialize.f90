module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of NEWUOA, described in Section 3 of the NEWUOA paper.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
!
! Started: July 2020
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Last Modified: Monday, June 13, 2022 PM04:32:49
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initxf, initq, inith


contains


subroutine initxf(calfun, iprint, maxfun, ftarget, rhobeg, x0, ij, kopt, nf, fhist, fval, xbase, &
    & xhist, xpt, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine does the initialization about the interpolation points & their function values.
!
! N.B.:
! 1. Remark on IJ:
! If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [NPT-2*N-1, 2].
! IJ contains integers between 2 and 2*N + 1. For each K > 2*N + 1, XPT(:, K) is
! XPT(:, IJ(K, 1)) + XPT(:, IJ(K, 2)). Let I = IJ(K, 1) - 1 if such a number is <= N + 1; otherwise,
! let IJ = IJ(K, 1) - N - 1; define J using IJ(K, 2) similarly. Then all the entries of XPT(:, K)
! are zero except that the I and J entries are RHOBEG or -RHOBEG. Indeed, XPT(I, K) is RHOBEG if
! IJ(K, 1) <= N + 1 and -RHOBEG otherwise; XPT(J, K) is similar. Consequently, the Hessian of the
! quadratic model will get a possibly nonzero (I, J) entry. In the code, IJ is defined according to
! Powell's original code as well as Section 3 of the NEWUOA paper and (2.4) of the BOBYQA paper.
! 2. At return,
! INFO = INFO_DFT: initialization finishes normally
! INFO = FTARGET_ACHIEVED: return because F <= FTARGET
! INFO = NAN_INF_X: return because X contains NaN
! INFO = NAN_INF_F: return because F is either NaN or +Inf
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! LOGICAL :: EVALUATED(NPT)
! REAL(RP) :: X(N)
! Size of local arrays: LOGICAL*NPT + REAL(RP)*N
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : eye
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)   ! X0(N)

! Outputs
integer(IK), intent(out) :: ij(:, :)    ! IJ(MAX(0_IK, NPT-2*N-1_IK), 2)
integer(IK), intent(out) :: info
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: fhist(:)   ! FHIST(MAXFHIST)
real(RP), intent(out) :: fval(:)    ! FVAL(NPT)
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
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
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
    call assert(size(fval) == npt, 'SIZE(FVAL) == NPT', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(size(x0) == n .and. all(is_finite(x0)), 'SIZE(X0) == N, X0 is finite', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
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

! Initialize XPT(:, 1: MIN(2*N + 1, NPT)).
xpt(:, 1) = ZERO
xpt(:, 2:n + 1) = rhobeg * eye(n)
! After the following line, XPT(:, 2*N+2 : NPT) = ZERO if it is nonempty. It will be revised later
! according to FVAL(2 : 2*N + 1).
xpt(:, n + 2:npt) = -rhobeg * eye(n, npt - n - 1_IK)

! Set FVAL(1 : min(2*N + 1, NPT)) by evaluating F. Totally parallelizable except for FMSG.
do k = 1, min(npt, 2_IK * n + 1_IK)
    x = xpt(:, k) + xbase
    call evaluate(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call fmsg(solver, iprint, k, f, x)
    ! Save X and F into the history.
    call savehist(k, x, xhist, f, fhist)
    ! Check whether to exit.
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do

! Set IJ.
! In general, when NPT = (N+1)*(N+2)/2, we can initialize IJ(1 : NPT - (2*N+1), :) to ANY permutation
! of {{I, J} : 1 <= J /= I <= N}; when NPT < (N+1)*(N+2)/2, we can set it to the first NPT - (2*N+1)
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
! such a choice will more likely render an XPT(:, K) with a lower function value.
! 3. This revision is OPTIONAL. Due to this revision, the definition of XPT(:, 2*N + 2 : NPT) relies
! on FVAL(2 : 2*N + 1), and it is the sole origin of the such dependency. If we remove the revision
! IJ, then the evaluations of FVAL(1 : NPT) can be merged, and they are totally PARALLELIZABLE; this
! can be beneficial if the function evaluations are expensive, which is likely the case.
! 4. MATLAB (but not Fortran) can index a vector using a 2D array of indices, thus the MATLAB code is
!!MATLAB: ij(fval(ij + n) < fval(ij)) = ij(fval(ij +n) < fval(ij)) + n;
where (fval(ij(:, 1) + n) < fval(ij(:, 1))) ij(:, 1) = ij(:, 1) + n
where (fval(ij(:, 2) + n) < fval(ij(:, 2))) ij(:, 2) = ij(:, 2) + n

! Set XPT(:, 2*N + 2 : NPT). It depends on IJ and hence on FVAL(2 : 2*N + 1). Indeed, XPT(:, K) has
! only two nonzeros for each K >= 2*N+2.
xpt(:, 2 * n + 2:npt) = xpt(:, ij(:, 1)) + xpt(:, ij(:, 2))

! Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except for FMSG.
if (info == INFO_DFT) then
    do k = 2_IK * n + 2_IK, npt
        x = xpt(:, k) + xbase
        call evaluate(calfun, x, f)
        call fmsg(solver, iprint, k, f, x)
        evaluated(k) = .true.
        fval(k) = f
        ! Save X and F into the history.
        call savehist(k, x, xhist, f, fhist)
        ! Check whether to exit.
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end do
end if

! Set NF, KOPT
nf = int(count(evaluated), kind(nf))  !!MATLAB: nf = sum(evaluated);
kopt = int(minloc(fval, mask=evaluated, dim=1), kind(kopt))
!!MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(ij, 1) == max(0_IK, npt - 2_IK * n - 1_IK) .and. size(ij, 2) == 2, &
        & 'SIZE(IJ) == [NPT - 2*N - 1, 2]', srname)
    call assert(all(ij >= 2 .and. ij <= 2 * n + 1), '2 <= IJ <= 2*N + 1', srname)
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
! List of local arrays (including function-output arrays; likely to be stored on the stack): NONE
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
integer(IK), intent(out), optional :: info
real(RP), intent(out) :: gq(:)  ! GQ(N)
real(RP), intent(out) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)

! Local variables
character(len=*), parameter :: srname = 'INITQ'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: k
integer(IK) :: n
integer(IK) :: ndiag
integer(IK) :: npt
real(RP) :: fbase
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
    call assert(all(ij >= 2 .and. ij <= 2 * n + 1), '2 <= IJ <= 2*N + 1', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) = [N, N]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

rhobeg = maxval(abs(xpt(:, 2)))  ! Read RHOBEG from XPT.
fbase = fval(1)  ! FBASE is the function value at XBASE.

! Set GQ by the forward difference.
gq(1:n) = (fval(2:n + 1) - fbase) / rhobeg

! The interpolation conditions decide GQ(1:NDIAG) and the first NDIAG diagonal 2nd derivatives of
! the initial quadratic model by a quadratic interpolation on three points, which is equivalent to
! the central finite difference.
ndiag = min(npt - n - 1_IK, n)

! Revise GQ(1:NDIAG) to the value provided by the central finite difference.
gq(1:ndiag) = HALF * (gq(1:ndiag) + (fbase - fval(n + 2:n + 1 + ndiag)) / rhobeg)

! Set the diagonal of HQ by the 2nd-order central finite difference. If we do this before the
! revision of GQ(1:NDIAG), we can avoid the calculation of FVAL(K + 1) - FBASE) / RHOBEG. But we
! prefer to decouple the initialization of GQ and HQ. We are not concerned by this amount of flops.
hq = ZERO
do k = 1, ndiag
    hq(k, k) = ((fval(k + 1) - fbase) / rhobeg - (fbase - fval(k + n + 1)) / rhobeg) / rhobeg
end do
!!MATLAB:
!!hdiag = ((fval(2 : ndiag+1) - fbase) / rhobeg - (fbase - fval(n+2 : n+ndiag+1)) / rhobeg) / rhobeg
!!hq(1:ndiag, 1:ndiag) = diag(hdiag)

! When NPT > 2*N + 1, set the off-diagonal entries of HQ.
do k = 1, npt - 2_IK * n - 1_IK
    ! With the I, J, XI, and XJ defined below, we have
    ! FVAL(K+2*n+1) = F(XBASE + XI*e_I + XJ*e_J),
    ! FVAL(IJ(K, 1)) = F(XBASE + XI*e_I),
    ! FVAL(IJ(K, 2)) = F(XBASE + XJ*e_J).
    ! Thus the HQ(I,J) defined below approximates frac{partial^2}{partial X_I partial X_J} F(XBASE).
    ! N.B.: Here, exchanging I and J will not lead to any change in precise arithmetic. Powell's
    ! code exchanges I and J if needed to ensure that  I > J. This is because Powell's code saves HQ
    ! as a 1D array that contains the lower triangular part of this symmetric matrix.
    i = modulo(ij(k, 1) - 2_IK, n) + 1_IK
    j = modulo(ij(k, 2) - 2_IK, n) + 1_IK
    xi = xpt(i, k + 2 * n + 1)
    xj = xpt(j, k + 2 * n + 1)
    !hq(i, j) = (fbase - fval(ij(k, 1)) - fval(ij(k, 2)) + fval(k + 2 * n + 1)) / (xi * xj)
    hq(i, j) = (fbase - (fval(ij(k, 1)) + fval(ij(k, 2))) + fval(k + 2 * n + 1)) / (xi * xj)
    hq(j, i) = hq(i, j)
end do

pq = ZERO

if (present(info)) then
    if (is_nan(sum(abs(gq)) + sum(abs(hq)))) then
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
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
end if

end subroutine initq


subroutine inith(ij, xpt, idz, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes IDZ, BMAT, and ZMAT.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): NONE
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, DEBUGGING
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
integer(IK), intent(out), optional :: info
real(RP), intent(out) :: bmat(:, :) ! BMAT(N, NPT + N)
real(RP), intent(out) :: zmat(:, :) ! ZMAT(NPT, NPT - N - 1)

! Local variables
character(len=*), parameter :: srname = 'INITH'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: recip
real(RP) :: reciq
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
    call assert(all(ij >= 2 .and. ij <= 2 * n + 1), '2 <= IJ <= 2*N + 1', srname)
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

! Set BMAT.
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

! Set ZMAT.
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
        zmat(ij(k, :), k + n) = -recip
    end do
end if

! Set IDZ = 1.
idz = 1_IK

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
