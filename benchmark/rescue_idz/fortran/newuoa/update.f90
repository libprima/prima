module update_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Wednesday, March 08, 2023 PM11:07:53
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updatexf, updateq, tryqalt


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
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
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

! KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW) = FVAL(KOPT) and KNEW < KOPT, then
! MINLOC(FVAL) = KNEW /= KOPT. Do not change KOPT in this case.
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
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
end if

end subroutine updatexf


subroutine updateq(idz, knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
! where XOSAV is the upupdated XOPT, namedly the XOPT before UPDATEXF is called.
! See Section 4 of the NEWUOA paper and that of the BOBYQA paper (there is no LINCOA paper).
! N.B.:
! 1. XNEW is encoded in [BMAT, ZMAT, IDZ] after UPDATEH being called, and it also equals XPT(:, KNEW)
! after UPDATEXF being called.
! 2. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
! 3. In Powell's implementation of NEWUOA, the quadratic model is represented by [GQ, PQ, HQ], where
! GQ is the gradient of the quadratic model at XBASE. However, Powell implemented BOBYQA and LINCOA
! without GQ but with GOPT, which is the gradient at XBASE + XOPT. In our implementation, we also
! use GOPT instead of GQ.
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


subroutine tryqalt(idz, bmat, fval, ratio, xopt, xpt, zmat, itest, gopt, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine tests whether to replace Q by the alternative model, namely the model that
! minimizes the F-norm of the Hessian subject to the interpolation conditions. It does the
! replacement if certain criteria are met (i.e., when ITEST = 3). See the paragraph around (8.4) of
! the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TEN, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric, trueloc
use, non_intrinsic :: powalg_mod, only : hess_mul, omega_mul

implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: fval(:)     ! FVAL(NPT)
real(RP), intent(in) :: ratio
real(RP), intent(in) :: xopt(:)     ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)   ! XOPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! In-output
integer(IK), intent(inout) :: itest
real(RP), intent(inout) :: gopt(:)    ! GOPT(N)
real(RP), intent(inout) :: hq(:, :) ! HQ(N, N)
real(RP), intent(inout) :: pq(:)    ! PQ(NPT)
! N.B.:
! GOPT, HQ, and PQ should be INTENT(INOUT) instead of INTENT(OUT). According to the Fortran 2018
! standard, an INTENT(OUT) dummy argument becomes undefined on invocation of the procedure.
! Therefore, if the procedure does not define such an argument, its value becomes undefined,
! which is the case for HQ and PQ when ITEST < 3 at exit. In addition, the information in GOPT is
! needed for definining ITEST, so it must be INTENT(INOUT).

! Local variables
character(len=*), parameter :: srname = 'TRYQALT'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: galt(size(gopt))
real(RP) :: pqalt(size(pq))

! Sizes
n = int(size(gopt), kind(n))
npt = int(size(pq), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    ! By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
    ! NaN, which should NOT happen due to the moderated extreme barrier.
    call assert(.not. is_nan(ratio), 'RATIO is not NaN', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the parameters of the least Frobenius norm interpolant to the current data.
pqalt = omega_mul(idz, zmat, fval)
galt = matprod(bmat(:, 1:npt), fval) + hess_mul(xopt, xpt, pqalt)

! Test whether to replace the new quadratic model by the least Frobenius norm interpolant, making
! the replacement if the test is satisfied. In the sequel, TEN seems to work a bit better than 100.
! In addition, Powell checked the magnitude of ABS(RATIO) instead of RATIO.
! !if (abs(ratio) > 0.01 .or. inprod(gopt, gopt) < 1.0E2_RP * inprod(galt, galt)) then ! Powell's code
if (ratio > TENTH .or. inprod(gopt, gopt) < TEN * inprod(galt, galt)) then
    itest = 0
else
    itest = itest + 1_IK
end if
if (itest >= 3) then
    gopt = galt
    pq = pqalt
    hq = ZERO
    itest = 0
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

end subroutine tryqalt


end module update_mod
