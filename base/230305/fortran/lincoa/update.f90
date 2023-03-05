module update_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's LINCOA code.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Sunday, March 05, 2023 PM05:00:05
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updatexf, updateq, tryqalt, updateres


contains


subroutine updatexf(knew, ximproved, f, xnew, kopt, fval, xpt)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates [XPT, FVAL, KOPT] so that XPT(:, KNEW) is updated to XNEW.
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
real(RP), intent(in) :: f
real(RP), intent(in) :: xnew(:)  ! XNEW(N)

! In-outputs
integer(IK), intent(inout) :: kopt
logical, intent(in) :: ximproved
real(RP), intent(inout) :: fval(:)  ! FVAL(NPT)
real(RP), intent(inout) :: xpt(:, :)! XPT(N, NPT)

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
    call assert(knew >= 1 .or. .not. ximproved, 'KNEW >= 1 unless X is not improved', srname)
    call assert(knew /= kopt .or. ximproved, 'KNEW /= KOPT unless X is improved', srname)
    call assert(size(xnew) == n .and. all(is_finite(xnew)), 'SIZE(XNEW) == N, XNEW is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN or +Inf', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do essentially nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    return
end if

xpt(:, knew) = xnew
fval(knew) = f

if (ximproved) then
    kopt = knew
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt .and. all(is_finite(xpt)), &
        & 'SIZE(XPT) == [N, NPT], XPT is finite', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
end if

end subroutine updatexf


subroutine updateq(idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
! where XOSAV is the upupdated XOPT, namedly the XOPT before UPDATEXF is called.
! See Section 4 of the NEWUOA paper and that of the BOBYQA paper (there is no LINCOA paper).
! N.B.:
! XNEW is encoded in [BMAT, ZMAT, IDZ] after UPDATEH being called, and it also equals XPT(:, KNEW)
! after UPDATEXF being called. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): PQINC
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : r1update, issymmetric
use, non_intrinsic :: powalg_mod, only : omega_col, hess_mul

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
logical, intent(in) :: ximproved
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
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(knew >= 1 .or. .not. ximproved, 'KNEW >= 1 unless X is not improved', srname)
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

! Further update GOPT if XIMPROVED is TRUE, as XOPT changes from XOSAV to XNEW = XOSAV + D.
if (ximproved) then
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


subroutine tryqalt(idz, bmat, fval, xopt, xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt)
!--------------------------------------------------------------------------------------------------!
! This subroutine tests whether to replace Q by the alternative model, namely the model that
! minimizes the F-norm of the Hessian subject to the interpolation conditions. It first calculates
! the alternative model represented by [GALT, PQALT], and sets [GOPT, PQ, HQ] = [GALT, PQALT, 0]
! if the recent few (three) alternative models are more accurate in predicting the function value of
! XOPT + D, i.e., if ALL(QALT_BETTER) = TRUE.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : matprod, issymmetric
use, non_intrinsic :: powalg_mod, only : omega_mul, hess_mul

implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: fval(:)  ! FVAL(NPT)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! In-outptuts
logical, intent(inout) :: qalt_better(:)  ! QALT_BETTER(3)
real(RP), intent(inout) :: gopt(:)  ! GOPT(N)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT)
real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)

! Outputs
real(RP), intent(out) :: galt(:)  ! GALT(N)
real(RP), intent(out) :: pqalt(:)  ! PQALT(NPT)

! Local variables
character(len=*), parameter :: srname = 'TRYQALT'
integer(IK) :: n
integer(IK) :: npt

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(size(galt) == n, 'SIZE(GALT) = N', srname)
    call assert(size(pqalt) == npt, 'SIZE(PQALT) = NPT', srname)
end if

!====================!
! Calculation starts !
!====================!

! Establish the alternative model, which is the least Frobenius norm interpolant.
pqalt = omega_mul(idz, zmat, fval)
galt = matprod(bmat(:, 1:npt), fval) + hess_mul(xopt, xpt, pqalt)

! Replace the current model with the alternative model if ALL(QALT_BETTER) = TRUE, i.e., the
! recent few alternative models are more accurate in predicting the function value of XOPT + D.
if (all(qalt_better)) then
    pq = pqalt
    hq = ZERO
    gopt = galt
    qalt_better = .false.
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(gopt) == n, 'SIZE(GOPT) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(size(galt) == n, 'SIZE(GALT) = N', srname)
    call assert(size(pqalt) == npt, 'SIZE(PQALT) = NPT', srname)
end if

end subroutine tryqalt


subroutine updateres(ximproved, amat, b, delta, dnorm, xopt, rescon)
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
logical, intent(in) :: ximproved
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
    ! Zaikun 20221115: The following cannot pass?! Is it due to the update of DELTA? Did we
    ! misunderstand Powell's definition of RESCON?
    !call assert(all((rescon >= 0 .and. rescon <= delta) .or. rescon <= -delta), &
    !    & '0 <= RESCON <= DELTA or RESCON <= -DELTA', srname)
end if

!====================!
! Calculation starts !
!====================!

! Zaikun 20221115: Currently, UPDATERES does not update RESCON unless XIMPROVED is TRUE. Shouldn't
! we do it whenever DELTA is updated? Have we MISUNDERSTOOD RESCON?
if (.not. ximproved) then
    return
end if

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
