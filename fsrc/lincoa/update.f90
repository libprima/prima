module update_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the update of [BMAT, ZMAT, IDZ] (represents H in the
! NEWUOA paper; there is no LINCOA paper), [XPT, FVAL, KOPT, XOPT, FOPT], [GQ, HQ, PQ] (the
! quadratic model), and RESCON when XPT(:, KNEW) is updated to XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's LICOA code.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, November 21, 2022 PM12:43:59
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updatexf, updateq, tryqalt, updateres


contains


subroutine updatexf(knew, freduced, d, f, kopt, fval, xpt, fopt, xopt)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates [XPT, FVAL, KOPT, XOPT, FOPT] so that X(:, KNEW) is updated to XOPT + D.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): NONE
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : norm

implicit none

! Inputs
integer(IK), intent(in) :: knew
real(RP), intent(in) :: d(:)     ! D(N)
real(RP), intent(in) :: f

! In-outputs
integer(IK), intent(inout) :: kopt
logical, intent(in) :: freduced
real(RP), intent(inout) :: fval(:)  ! FVAL(NPT)
real(RP), intent(inout) :: xpt(:, :)! XPT(N, NPT)

! Outputs
real(RP), intent(out) :: fopt
real(RP), intent(out) :: xopt(:)    ! XOPT(N)

! Local variables
character(len=*), parameter :: srname = 'UPDATEXF'
integer(IK) :: n
integer(IK) :: npt

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew >= 1 .or. f >= fval(kopt), 'KNEW >= 1 unless F >= FVAL(KOPT)', srname)
    call assert(knew /= kopt .or. f < fval(kopt), 'KNEW /= KOPT unless F < FVAL(KOPT)', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN or +Inf', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(size(xopt) == n, 'SIZE(XOPT) == N', srname)
    ! N.B.: Do NOT test the value of FOPT or XOPT. Being INTENT(OUT), they are UNDEFINED up to here.
end if

!====================!
! Calculation starts !
!====================!

! Do essentially nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    ! We must set XOPT and FOPT. Otherwise, they are UNDEFINED because we declare them as INTENT(OUT).
    xopt = xpt(:, kopt)
    fopt = fval(kopt)
    return
end if

xpt(:, knew) = xpt(:, kopt) + d
fval(knew) = f

if (freduced) then
    kopt = knew
end if

! Even if KOPT remains unchanged, we still need to update XOPT and FOPT, because it may happen that
! KNEW = KOPT, so that XPT(:, KOPT) has been updated to XNEW = XOPT + D.
xopt = xpt(:, kopt)
fopt = fval(kopt)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(abs(f - fval(knew)) <= 0, 'F == FVAL(KNEW)', srname)
    call assert(abs(fopt - fval(kopt)) <= 0, 'FOPT == FVAL(KOPT)', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(norm(xopt - xpt(:, kopt)) <= 0, 'XOPT == XPT(:, KOPT)', srname)
end if

end subroutine updatexf


subroutine updateq(idz, knew, freduced, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
! where XOSAV is the upupdated XOPT, namedly the XOPT before UPDATEXF is called.
! See Section 4 of the NEWUOA paper (there is no LINCOA paper).
! N.B.:
! XNEW is encoded in [BMAT, ZMAT, IDZ] after UPDATEH being called, and it also equals XPT(:, KNEW)
! after UPDATEXF being called. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): PQINC
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_posinf, is_nan
use, non_intrinsic :: linalg_mod, only : r1update, issymmetric
use, non_intrinsic :: powalg_mod, only : quadinc, omega_col, hess_mul

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
logical, intent(in) :: freduced
real(RP), intent(in) :: bmat(:, :) ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:) ! D(:)
real(RP), intent(in) :: moderr
real(RP), intent(in) :: xdrop(:)  ! XDROP(N)
real(RP), intent(in) :: xosav(:)  ! XOSAV(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! In-outputs
real(RP), intent(inout) :: gopt(:)  ! GOPT(N)
real(RP), intent(inout) :: hq(:, :) ! HQ(N, N)
real(RP), intent(inout) :: pq(:)    ! PQ(NPT)

! Local variables
character(len=*), parameter :: srname = 'UPDATEQ'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: pqinc(size(pq))

! Sizes
n = int(size(gopt), kind(n))
npt = int(size(pq), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(size(xdrop) == n .and. all(is_finite(xdrop)), 'SIZE(XDROP) == N, XDROP is finite', srname)
    call assert(size(xosav) == n .and. all(is_finite(xosav)), 'SIZE(XOSAV) == N, XOSAV is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    return
end if

! The unupdated model corresponding to [GOPT, HQ, PQ] interpolates F at all points in XPT except for
! XNEW. The error is MODERR = [F(XNEW)-F(XOPT)] - [Q(XNEW)-Q(XOPT)].

! Absorb PQ(KNEW)*XDROP*XDROP^T into the explicit part of the Hessian.
! Implement R1UPDATE properly so that it ensures that HQ is symmetric.
call r1update(hq, pq(knew), xdrop)
pq(knew) = ZERO

! Update the implicit part of the Hessian.
pqinc = moderr * omega_col(idz, zmat, knew)
pq = pq + pqinc

! Update the gradient, which needs the updated XPT.
gopt = gopt + moderr * bmat(:, knew) + hess_mul(xosav, xpt, pqinc)

! Further update GOPT if FREDUCED is TRUE, as XOPT changes from XOSAV to XNEW = XOSAV + D.
if (freduced) then
    gopt = gopt + hess_mul(d, xpt, pq, hq)
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

end subroutine updateq


subroutine tryqalt(idz, fval, ratio, bmat, zmat, itest, gq, hq, pq)
!--------------------------------------------------------------------------------------------------!
! TRYQALT tests whether to replace Q by the alternative model, namely the model that minimizes
! the F-norm of the Hessian subject to the interpolation conditions. It does the replacement
! when certain criteria are satisfied (i.e., when ITEST = 3). See Section 8 of the NEWUOA paper.
! N.B.: Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: GALT(N)
! Size of local arrays: REAL(RP)*(N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : inprod, matprod, issymmetric
use, non_intrinsic :: powalg_mod, only : omega_mul

implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: fval(:)     ! FVAL(NPT)
real(RP), intent(in) :: ratio
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! In-output
integer(IK), intent(inout) :: itest
real(RP), intent(inout) :: gq(:)    ! GQ(N)
real(RP), intent(inout) :: hq(:, :) ! HQ(N, N)
real(RP), intent(inout) :: pq(:)    ! PQ(NPT)
! N.B.:
! GQ, HQ, and PQ should be INTENT(INOUT) instead of INTENT(OUT). According to the Fortran 2018
! standard, an INTENT(OUT) dummy argument becomes undefined on invocation of the procedure.
! Therefore, if the procedure does not define such an argument, its value becomes undefined,
! which is the case for HQ and PQ when ITEST < 3 at exit. In addition, the information in GQ is
! needed for definining ITEST, so it must be INTENT(INOUT).

! Local variables
character(len=*), parameter :: srname = 'TRYQALT'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: galt(size(gq))

! Debugging variables
!real(RP) :: intp_tol

! Sizes
n = int(size(gq), kind(n))
npt = int(size(pq), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    ! By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
    ! NaN, which should NOT happen due to the moderated extreme barrier.
    call assert(.not. is_nan(ratio), 'RATIO is not NaN', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    ! [GQ, HQ, PQ] cannot pass the following test if FVAL contains extremely large values.
    !intp_tol = max(1.0E-8_RP, min(1.0E-1_RP, 1.0E10_RP * real(size(pq), RP) * EPS))
    !call wassert(errquad(fval, xpt, gq, pq, hq) <= intp_tol, 'Q interpolates FVAL at XPT', srname)
end if

!====================!
! Calculation starts !
!====================!

! In the NEWUOA paper, Powell replaces Q with Q_alt when RATIO <= 0.01 and ||G_alt|| <= 0.1||GQ||
! hold for 3 consecutive times (eq(8.4)). But Powell's code compares ABS(RATIO) instead of RATIO
! with 0.01. Here we use RATIO, which is more efficient as observed in Zaikun ZHANG's PhD thesis
! (Section 3.3.2).
!if (abs(ratio) > 1.0e-2_RP) then
if (ratio > 1.0E-2_RP) then
    itest = 0_IK
else
    galt = matprod(bmat(:, 1:npt), fval)
    if (inprod(gq, gq) < 1.0E2_RP * inprod(galt, galt)) then
        itest = 0_IK
    else
        itest = itest + 1_IK
    end if
end if

! Replace Q with Q_alt when ITEST >= 3.
if (itest >= 3) then
    gq = galt
    hq = ZERO
    pq = omega_mul(idz, zmat, fval)
    itest = 0_IK
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    ! [GQ, HQ, PQ] cannot pass the following test if FVAL contains extremely large values.
    !call wassert(errquad(fval, xpt, gq, pq, hq) <= intp_tol, 'QALT interpolates FVAL at XPT', srname)
end if

end subroutine tryqalt


subroutine updateres(amat, b, delta, dnorm, xopt, rescon)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates RESCON when XOPT has been updated by a step D.
! RESCON holds information about the constraint residuals at the current trust region center XOPT.
! 1. If if B(J) - AMAT(:, J)^T*XOPT <= DELTA, then RESCON(J) = B(J) - AMAT(:, J)^T*XOPT. Note that
! RESCON >= 0 in this case, because the algorithm keeps XOPT to be feasible.
! 2. Otherwise, RESCON(J) is a negative value that B(J) - AMAT(:, J)^T*XOPT >= |RESCON(J)| >= DELTA.
! RESCON can be updated without calculating the constraints that are far from being active, so that
! we only need to evaluate the constraints that are nearly active.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, trueloc

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: b(:)  ! B(M)
real(RP), intent(in) :: delta
real(RP), intent(in) :: dnorm  ! Norm of D
real(RP), intent(in) :: xopt(:)  ! XOPT(N); the updated value of XOPT

! In-outputs
real(RP), intent(inout) :: rescon(:)  ! RESCON(M)

! Local variables
character(len=*), parameter :: srname = 'UPDATERES'
integer(IK) :: m
integer(IK) :: n
logical :: mask(size(b))
real(RP) :: ax(size(b))

! Sizes
m = int(size(b), kind(m))
n = int(size(xopt), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(dnorm > 0, 'DNORM > 0', srname)
    call assert(all(is_finite(xopt)), 'XOPT is finite', srname)
    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    ! Zaikun 20221115: The following cannot pass?! Is it due to the update of DELTA?
    !call assert(all((rescon >= 0 .and. rescon <= delta) .or. rescon <= -delta), &
    !    & '0 <= RESCON <= DELTA or RESCON <= -DELTA', srname)
end if

!====================!
! Calculation starts !
!====================!

mask = (abs(rescon) < dnorm + delta)
ax(trueloc(mask)) = matprod(xopt, amat(:, trueloc(mask)))
where (mask)
    rescon = max(b - ax, ZERO)
elsewhere
    rescon = min(-abs(rescon) + dnorm, -delta)
end where
rescon(trueloc(rescon >= delta)) = -rescon(trueloc(rescon >= delta))

!!MATLAB:
!!mask = (abs(rescon) < delta + dnorm);
!!rescon(mask) = max(b(mask) - (xopt'*amat(:, mask))', 0);
!!rescon(~mask) = max(rescon(~mask) - dnorm, delta);
!!rescon(rescon >= delta) = -rescon(rescon >= delta);

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    !call assert(all((rescon >= 0 .and. rescon <= delta) .or. rescon <= -delta), &
    !    & '0 <= RESCON <= DELTA or RESCON <= -DELTA', srname)
end if
end subroutine updateres


end module update_mod
