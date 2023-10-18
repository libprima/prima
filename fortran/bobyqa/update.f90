module update_bobyqa_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the updates when XPT(:, KNEW) becomes XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, October 18, 2023 PM10:15:42
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updateh, updatexf, updateq, tryqalt


contains


subroutine updateh(knew, kopt, d, xpt, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates arrays BMAT and ZMAT in order to replace the interpolation point
! XPT(:, KNEW) by XNEW = XPT(:, KOPT) + D. See Section 4 of the BOBYQA paper. [BMAT, ZMAT] describes
! the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of the coefficient matrix of the
! KKT system for the least-Frobenius norm interpolation problem: ZMAT holds a factorization of the
! leading NPT*NPT submatrix OMEGA of H, the factorization being OMEGA = ZMAT*ZMAT^T; BMAT holds the
! last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th
! column of H are not stored as they are unnecessary for the calculation.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: infos_mod, only : INFO_DFT, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : planerot, matprod, outprod, symmetrize, issymmetric
use, non_intrinsic :: powalg_mod, only : calbeta, calvlag
use, non_intrinsic :: string_mod, only : num2str

implicit none

! Inputs
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: d(:)  ! D(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
integer(IK), intent(out), optional :: info

! Local variables
character(len=*), parameter :: srname = 'UPDATEH'
integer(IK) :: j
integer(IK) :: n
integer(IK) :: npt
real(RP) :: alpha
real(RP) :: beta
real(RP) :: denom
real(RP) :: grot(2, 2)
real(RP) :: hcol(size(bmat, 2))
real(RP) :: sqrtdn
real(RP) :: tau
real(RP) :: v1(size(bmat, 1))
real(RP) :: v2(size(bmat, 1))
real(RP) :: vlag(size(bmat, 2))
real(RP) :: ztest

! Sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)

    do j = 1, npt
        hcol(1:npt) = matprod(zmat, zmat(j, :))
        hcol(npt + 1:npt + n) = bmat(:, j)
        call assert(sum(abs(hcol)) > 0, 'Column '//num2str(j)//' of H is nonzero', srname)
    end do

    call assert(all(is_finite(xpt)), 'XPT is finite', srname)

    ! The following is too expensive to check.
    !tol = 1.0E-2_RP
    !call wassert(errh(bmat, zmat, xpt) <= tol .or. RP == kind(0.0), &
    !    & 'H = W^{-1} in (2.7) of the BOBYQA paper', srname)
end if

!====================!
! Calculation starts !
!====================!

if (present(info)) then
    info = INFO_DFT
end if

! Do anything if KNEW is 0. This can only happen sometimes after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    return
end if

! Put the KNEW-th column of the unupdated H (except for the (NPT+1)th entry) into HCOL. Powell's
! code does this after ZMAT is rotated below, and then HCOL(1:NPT) = ZMAT(KNEW, 1) * ZMAT(:, 1),
! which saves flops but also introduces rounding errors due to the rotation.
hcol(1:npt) = matprod(zmat, zmat(knew, :))
hcol(npt + 1:npt + n) = bmat(:, knew)

! Calculate VLAG and BETA and other parameters for (4.9) and (4.14) of the BOBYQA paper.
beta = calbeta(kopt, bmat, d, xpt, zmat)
vlag = calvlag(kopt, bmat, d, xpt, zmat)

! In theory, DENOM can also be calculated after ZMAT is rotated below. However, this worsened the
! performance of BOBYQA in a test on 20220413.
alpha = hcol(knew)
tau = vlag(knew)
denom = alpha * beta + tau**2

! After the following line, VLAG = H*w - e_KNEW in the NEWUOA paper (where t = KNEW).
vlag(knew) = vlag(knew) - ONE

! Quite rarely, due to rounding errors, VLAG or BETA may not be finite, or DENOM may not be
! positive. In such cases, [BMAT, ZMAT] would be destroyed by the update, and hence we would rather
! not update them at all. Or should we simply terminate the algorithm?
if (.not. (is_finite(sum(abs(hcol)) + sum(abs(vlag)) + abs(beta)) .and. denom > 0)) then
    if (present(info)) then
        info = DAMAGING_ROUNDING
    end if
    return
end if

! Update the matrix BMAT. It implements the last N rows of (4.9) in the BOBYQA paper.
v1 = (alpha * vlag(npt + 1:npt + n) - tau * hcol(npt + 1:npt + n)) / denom
v2 = (-beta * hcol(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) / denom
bmat = bmat + outprod(v1, vlag) + outprod(v2, hcol) !call r2update(bmat, ONE, v1, vlag, ONE, v2, hcol)
! Numerically, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be symmetric.
call symmetrize(bmat(:, npt + 1:npt + n))

! Apply Givens rotations to put zeros in the KNEW-th row of ZMAT. After this, ZMAT(KNEW, :) contains
! only one nonzero at ZMAT(KNEW, 1). Entries of ZMAT are treated as 0 if the moduli are at most ZTEST.
ztest = 1.0E-20_RP * maxval(abs(zmat))  ! This threshold is by Powell
do j = 2, npt - n - 1_IK
    if (abs(zmat(knew, j)) > ztest) then
        grot = planerot(zmat(knew, [1_IK, j]))
        zmat(:, [1_IK, j]) = matprod(zmat(:, [1_IK, j]), transpose(grot))
    end if
    zmat(knew, j) = ZERO
end do

! Complete the updating of ZMAT. See (4.14) of the BOBYQA paper.
sqrtdn = sqrt(denom)
zmat(:, 1) = (tau / sqrtdn) * zmat(:, 1) - (zmat(knew, 1) / sqrtdn) * vlag(1:npt)
! Zaikun 20231012: Either of the following two lines worsens the performance of BOBYQA when the
! objective function is evaluated with 5 or less correct significance digits. Strange.
! !zmat(:, 1) = (tau * zmat(:, 1) - zmat(knew, 1) * vlag(1:npt)) / sqrtdn
! !zmat(knew, 1) = zknew1 / sqrtdn  ! ZKNEW1 is the unupdated ZMAT(KNEW, 1)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)

    do j = 1, npt
        hcol(1:npt) = matprod(zmat, zmat(j, :))
        hcol(npt + 1:npt + n) = bmat(:, j)
        call assert(sum(abs(hcol)) > 0, 'Column '//num2str(j)//' of H is nonzero', srname)
    end do

    ! The following is too expensive to check.
    ! !if (n * npt <= 50) then
    ! !    xpt_test = xpt
    ! !    xpt_test(:, knew) = xpt(:, kopt) + d
    ! !    call assert(errh(bmat, zmat, xpt_test) <= tol .or. RP == kind(0.0), &
    ! !        & 'H = W^{-1} in (2.7) of the BOBYQA paper', srname)
    ! !end if
end if
end subroutine updateh


subroutine updatexf(knew, ximproved, f, xnew, kopt, fval, xpt)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates [XPT, FVAL, KOPT] so that XPT(:, KNEW) is updated to XNEW.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf

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


subroutine updateq(knew, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates GOPT, HQ, and PQ when XPT(:, KNEW) changes from XDROP to XNEW = XOSAV + D,
! where XOSAV is the unupdated XOPT, namely the XOPT before UPDATEXF is called.
! See Section 4 of the NEWUOA paper and that of the BOBYQA paper (there is no LINCOA paper).
! N.B.:
! XNEW is encoded in [BMAT, ZMAT] after UPDATEH being called, and it also equals XPT(:, KNEW)
! after UPDATEXF being called. Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, r1update, issymmetric
use, non_intrinsic :: powalg_mod, only : hess_mul

implicit none

! Inputs
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
pqinc = moderr * matprod(zmat, zmat(knew, :))  ! pqinc = moderr * omega_col(1_IK, zmat, knew)
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


subroutine tryqalt(bmat, fval, ratio, sl, su, xopt, xpt, zmat, itest, gopt, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine tests whether to replace Q by the alternative model, namely the model that
! minimizes the F-norm of the Hessian subject to the interpolation conditions. It does the
! replacement if certain criteria are met (i.e., when ITEST = 3). See the paragraph around (6.12) of
! the BOBYQA paper.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TEN, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric, trueloc
use, non_intrinsic :: powalg_mod, only : hess_mul

implicit none

! Inputs
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: fval(:)     ! FVAL(NPT)
real(RP), intent(in) :: ratio
real(RP), intent(in) :: sl(:)       ! SL(N)
real(RP), intent(in) :: su(:)       ! SU(N)
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
! needed for defining ITEST, so it must be INTENT(INOUT).

! Local variables
character(len=*), parameter :: srname = 'TRYQALT'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: galt(size(gopt))
real(RP) :: pgalt(size(gopt))
real(RP) :: pgopt(size(gopt))
real(RP) :: pqalt(size(pq))

! Debugging variables
!real(RP) :: intp_tol

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

! Calculate the norm square of the projected gradient.
pgopt = gopt
pgopt(trueloc(xopt >= su)) = max(ZERO, gopt(trueloc(xopt >= su)))
pgopt(trueloc(xopt <= sl)) = min(ZERO, gopt(trueloc(xopt <= sl)))

! Calculate the parameters of the least Frobenius norm interpolant to the current data.
pqalt = matprod(zmat, matprod(fval, zmat))
galt = matprod(bmat(:, 1:npt), fval) + hess_mul(xopt, xpt, pqalt)

! Calculate the norm square of the projected alternative gradient.
pgalt = galt
pgalt(trueloc(xopt >= su)) = max(ZERO, galt(trueloc(xopt >= su)))
pgalt(trueloc(xopt <= sl)) = min(ZERO, galt(trueloc(xopt <= sl)))

! Test whether to replace the new quadratic model by the least Frobenius norm interpolant,
! making the replacement if the test is satisfied.
! N.B.: In the following IF, Powell's condition does not check RATIO. The condition here (with RATIO
! > TENTH)is adopted and adapted from NEWUOA, and it seems to improve the performance.
! !if (inprod(pgopt, pgopt) < TEN * inprod(pgalt, pgalt)) then  ! Powell's code
if (ratio > TENTH .or. inprod(pgopt, pgopt) < TEN * inprod(pgalt, pgalt)) then
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


end module update_bobyqa_mod
