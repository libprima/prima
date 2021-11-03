module update_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the update of IDZ, BMAT, ZMAT, GQ, HQ, PQ, FVAL, XPT,
! KOPT, FOPT, and XOPT when XPT(:, KNEW) is replaced by XNEW = XOPT + D.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Wednesday, November 03, 2021 AM11:20:38
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updateh, updateq, updatexf, tryqalt


contains


subroutine updateh(knew, kopt, idz, d, xpt, bmat, zmat)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates arrays BMAT and ZMAT together with IDZ, in order to replace the
! interpolation point XPT(:, KNEW) by XNEW = XOPT + D. See Section 4 of the NEWUOA paper.
! [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
! the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
! ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
! ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT-N-1)=1. BMAT will hold the last N ROWs of H
! except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th are not saved as
! they are unnecessary for the calculation.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, planerot, r2update, symmetrize, issymmetric

! Solver-specific modules
use, non_intrinsic :: vlagbeta_mod, only : vlagbeta

implicit none

! Inputs
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: d(:)  ! D(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
integer(IK), intent(inout) :: idz
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Local variables
character(len=*), parameter :: srname = 'UPDATEH'
integer(IK) :: j
integer(IK) :: ja
integer(IK) :: jb
integer(IK) :: jl
integer(IK) :: n
integer(IK) :: npt
logical :: reduce_idz
real(RP) :: alpha
real(RP) :: beta
real(RP) :: denom
real(RP) :: grot(2, 2)
!real(RP) :: htol
real(RP) :: scala
real(RP) :: scalb
real(RP) :: sqrtdn
real(RP) :: tau
real(RP) :: tausq
real(RP) :: temp
real(RP) :: tempa
real(RP) :: tempb
real(RP) :: v1(size(bmat, 1))
real(RP) :: v2(size(bmat, 1))
real(RP) :: vlag(size(bmat, 2))
real(RP) :: w(size(bmat, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '0 <= KOPT <= NPT', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N and D is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    ! The following test cannot be passed.
    !htol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E10_RP * EPS)) ! Tolerance for H
    !call assert(errh(idz, bmat, zmat, xpt) <= htol, 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew == 0) then
    return
end if

! Calculate VLAG and BETA according to D.
! VLAG contains the components of the vector Hw of the updating formula (4.11) in the NEWUOA paper,
! and BETA holds the value of the parameter that has this name.
! N.B.: Powell's original comments mention that VLAG is "the vector THETA*WCHECK + e_b of the
! updating formula (6.11)", which does not match the published version of the NEWUOA paper.
call vlagbeta(idz, kopt, bmat, d, xpt, zmat, beta, vlag)

! Apply rotations to put zeros in the KNEW-th row of ZMAT. A 2x2 rotation will be multiplied to ZMAT
! from the right so that ZMAT(KNEW, [JL, J]) becomes [SQRT(ZMAT(KNEW, JL)^2 + ZMAT(KNEW, J)^2), 0].
! As in MATLAB, PLANEROT(X) returns a 2x2 Givens matrix G for X in R^2 so that Y = G*X has Y(2) = 0.

! In the loop, if 2 <= J < IDZ, then JL = 1; if IDZ < J <= NPT - N - 1, then JL = IDZ.
jl = 1_IK
do j = 2, int(npt - n - 1, kind(j))
    if (j == idz) then
        jl = idz
        cycle
    end if
    if (abs(zmat(knew, j)) > ZERO) then
        grot = planerot(zmat(knew, [jl, j]))  ! MATLAB code: GROT = PLANEROT(ZMAT(KNEW, [JL, J])')
        zmat(:, [jl, j]) = matprod(zmat(:, [jl, j]), transpose(grot))
        zmat(knew, j) = ZERO
    end if
end do
! The value of JL after the loop is important below. Its value is determined by the current (i.e.,
! unupdated) value of IDZ. IDZ is an integer in {1, ..., NPT-N} such that S_J = -1 for J < IDZ while
! S_J = 1 for J >= IDZ in the factorization of OMEGA. See (3.17) and (4.16) of the NEWUOA paper.
! The value of JL has two possibilities:
! 1. JL = 1 iff IDZ = 1 or IDZ = NPT - N.
! 1.1. IDZ = 1 means that
! OMEGA = sum_{J=1}^{NPT-N-1} ZMAT(:, J)*ZMAT(:, J)' ;
! 1.2. IDZ = NPT - N means that
! OMEGA = - sum_{J=1}^{NPT-N-1} ZMAT(:, J)*ZMAT(:, J)' ;
! 2. JL = IDZ > 1 iff 2 <= IDZ <= NPT - N - 1.

! Put the first NPT components of the KNEW-th column of HLAG into W, and calculate the parameters of
! the updating formula.
tempa = zmat(knew, 1)
if (idz >= 2) then
    tempa = -tempa
end if

w(1:npt) = tempa * zmat(:, 1)
if (jl > 1) then
    tempb = zmat(knew, jl)
    w(1:npt) = w(1:npt) + tempb * zmat(:, jl)
end if

alpha = w(knew)
tau = vlag(knew)
tausq = tau * tau
denom = alpha * beta + tausq
! After the following line, VLAG = Hw - e_t in the NEWUOA paper.
vlag(knew) = vlag(knew) - ONE
sqrtdn = sqrt(abs(denom))

reduce_idz = .false.
if (jl == 1) then
    ! Complete the updating of ZMAT when there is only one nonzero element in ZMAT(KNEW, :) after
    ! the rotation. This is the normal case, because IDZ is 1 in precise arithmetic.
    !---------------------------------------------------------------------------------------------!
    ! Up to now, TEMPA = ZMAT(KNEW, 1) if IDZ = 1 and TEMPA = -ZMAT(KNEW, 1) if IDZ >= 2. However,
    ! according to (4.18) of the NEWUOA paper, TEMPB should always be ZMAT(KNEW, 1)/sqrtdn
    ! regardless of IDZ. Therefore, the following definition of TEMPB is inconsistent with (4.18).
    ! This is probably a BUG. See also Lemma 4 and (5.13) of Powell's paper "On updating the inverse
    ! of a KKT matrix". However, the inconsistency is hardly observable in practice, because JL = 1
    ! implies IDZ = 1 in precise arithmetic.
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    !tempb = tempa/sqrtdn
    !tempa = tau/sqrtdn
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Here is the corrected version (only TEMPB is changed).
    tempa = tau / sqrtdn
    tempb = zmat(knew, 1) / sqrtdn
    !---------------------------------------------------------------------------------------------!

    zmat(:, 1) = tempa * zmat(:, 1) - tempb * vlag(1:npt)

    !---------------------------------------------------------------------------------------------!
    ! The following six lines by Powell are obviously problematic --- SQRTDN is always nonnegative.
    ! According to (4.18) of the NEWUOA paper, "SQRTDN < ZERO" and "SQRTDN >= ZERO" below should be
    ! both revised to "DENOM < ZERO". See also the corresponding part of the LINCOA code. Note that
    ! the NEWUOA paper uses SIGMA to denote DENOM. Check also Lemma 4 and (5.13) of Powell's paper
    ! "On updating the inverse of a KKT matrix". It seems that the BOBYQA code does not have this
    ! part --- it does not have IDZ at all (why?). Anyway, these lines are not invoked very often in
    ! practice, because IDZ should always be 1 in precise arithmetic.
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    !if (idz == 1 .and. sqrtdn < ZERO) then
    !    idz = 2
    !end if
    !if (idz >= 2 .and. sqrtdn >= ZERO) then
    !    reduce_idz = .true.
    !end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! This is the corrected version. It duplicates` the
    ! corresponding part of the LINCOA code.
    if (denom < ZERO) then
        if (idz == 1) then
            ! This is the first place (out of two) where IDZ is
            ! increased. Note that IDZ = 2 <= NPT-N after the update.
            idz = 2_IK
        else
            ! This is the first place (out of two) where IDZ is
            ! decreased (by 1). Since IDZ >= 2 in this case, we have
            ! IDZ >= 1 after the update.
            reduce_idz = .true.
        end if
    end if
    !---------------------------------------------------------------------------------------------!
else
    ! Complete the updating of ZMAT in the alternative case. ZMAT(KNEW, :) has two nonzeros after
    ! the rotations.
    ja = 1_IK
    if (beta >= ZERO) then
        ja = jl
    end if
    jb = int(jl + 1 - ja, kind(jb))
    temp = zmat(knew, jb) / denom
    tempa = temp * beta
    tempb = temp * tau
    temp = zmat(knew, ja)
    scala = ONE / sqrt(abs(beta) * temp**2 + tausq)
    scalb = scala * sqrtdn
    zmat(:, ja) = scala * (tau * zmat(:, ja) - temp * vlag(1:npt))
    zmat(:, jb) = scalb * (zmat(:, jb) - tempa * w(1:npt) - tempb * vlag(1:npt))
    ! If and only if DENOM <= 0, IDZ will be revised according to the sign of BETA.
    ! See (4.19)--(4.20) of the NEWUOA paper.
    if (denom <= ZERO) then
        if (beta < ZERO) then
            ! This is the second place (out of two) where IDZ is increased. Since
            ! JL = IDZ <= NPT-N-1 in this case, we have IDZ <= NPT-N after the update.
            idz = int(idz + 1, kind(idz))
        end if
        if (beta >= ZERO) then
            ! This is the second place (out of two) where IDZ is decreased (by 1). Since IDZ >= 2
            ! in this case, we have IDZ >= 1 after the update.
            reduce_idz = .true.
        end if
    end if
end if

! IDZ is reduced in the following case. Then exchange ZMAT(:, 1) and ZMAT(:, IDZ).
if (reduce_idz) then
    idz = int(idz - 1, kind(idz))
    if (idz > 1) then
        ! If a vector subscript has two or more elements with the same value, an array section with
        ! that vector subscript is not definable and shall not be defined or become undefined.
        zmat(:, [1_IK, idz]) = zmat(:, [idz, 1_IK])
    end if
end if

! Finally, update the matrix BMAT. It implements the last N rows of (4.11) in the NEWUOA paper.
w(npt + 1:npt + n) = bmat(:, knew)
v1 = (alpha * vlag(npt + 1:npt + n) - tau * w(npt + 1:npt + n)) / denom
v2 = (-beta * w(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) / denom
call r2update(bmat, ONE, v1, vlag, ONE, v2, w)
! In floating-point arithmetic, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be
! symmetric. Symmetrization needed.
call symmetrize(bmat(:, npt + 1:npt + n))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    !xpt_test = xpt
    !xpt_test(:, knew) = xpt(:, kopt) + d
    ! The following test cannot be passed.
    !call assert(errh(idz, bmat, zmat, xpt_test) <= htol, &
    !    & 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
end if

end subroutine updateh


subroutine updateq(idz, knew, kopt, bmat, d, f, fval, xpt, zmat, gq, hq, pq)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates GQ, HQ, and PQ when XPT(:, KNEW) is replaced by XNEW = XOPT + D. See
! Section 4 of the NEWUOA paper.
! N.B.: Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_posinf, is_nan
use, non_intrinsic :: linalg_mod, only : r1update, Ax_plus_y, issymmetric, calquad

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :) ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:) ! D(N)
real(RP), intent(in) :: f
real(RP), intent(in) :: fval(:) ! FVAL(NPT)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! In-outputs
real(RP), intent(inout) :: gq(:)    ! GQ(N)
real(RP), intent(inout) :: hq(:, :) ! HQ(N, N)
real(RP), intent(inout) :: pq(:)    ! PQ(NPT)

! Local variables
character(len=*), parameter :: srname = 'UPDATEQ'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: moderr
real(RP) :: modez(size(zmat, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(knew >= 1 .or. f >= fval(kopt), 'KNEW >= 1 unless F >= FOPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt .or. f < fval(kopt), 'KNEW /= KOPT unless F < FOPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N and D is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN or +Inf', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    ! [GQ, HQ, PQ] cannot pass the following test if FVAL contains extremely large values.
    !call assert(errquad(gq, hq, pq, xpt, fval) <= intp_tol, 'Q interpolates FVAL at XPT', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew == 0) then
    return
end if

! The unupdated model corresponding to [GQ, HQ, PQ] interpolates F at all points in XPT except for
! XNEW, which will become XPT(:, KNEW). The error is MODERR = [F(XNEW)-F(XOPT)] - [Q(XNEW)-Q(XOPT)].
! In the following, CALQUAD = Q(XOPT + D) - Q(XOPT) = Q(XNEW) - Q(XOPT).
moderr = f - fval(kopt) + calquad(d, gq, hq, pq, xpt(:, kopt), xpt)

! Absorb PQ(KNEW)*XPT(:, KNEW)*XPT(:, KNEW)^T into the explicit part of the Hessian.
! Implement R1UPDATE properly so that it ensures HQ is symmetric.
call r1update(hq, pq(knew), xpt(:, knew))
pq(knew) = ZERO

! Update the implicit part of the Hessian.
modez = moderr * zmat(knew, :)
modez(1:idz - 1) = -modez(1:idz - 1)
!----------------------------------------------------------------!
!pq = pq + matprod(zmat, modez) !---------------------------------!
pq = Ax_plus_y(zmat, modez, pq)
!----------------------------------------------------------------!

! Update the gradient.
gq = gq + moderr * bmat(:, knew)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(issymmetric(hq), 'HQ is symmetric', srname)
end if

end subroutine updateq


subroutine updatexf(knew, d, f, kopt, fval, xpt, fopt, xopt)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates XPT, FVAL, KOPT, XOPT, and FOPT so that X(:, KNEW) is replaced by XOPT+D.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
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
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(knew >= 1 .or. f >= fval(kopt), 'KNEW >= 1 unless F >= FOPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt .or. f < fval(kopt), 'KNEW /= KOPT unless F < FOPT', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N and D is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN or +Inf', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(xopt) == n, 'SIZE(XOPT) == N', srname)
    ! N.B.: Do NOT test the value of FOPT or XOPT. Being INTENT(OUT), they are UNDEFINED up to here.
end if

!====================!
! Calculation starts !
!====================!

! Do essentially nothing when KNEW is 0. This can only happen after a trust-region step.
if (knew == 0) then
    ! We must set FOPT and XOPT. Otherwise, they are UNDEFINED because we declare them as INTENT(OUT).
    fopt = fval(kopt)
    xopt = xpt(:, kopt)
    return
end if

xpt(:, knew) = xpt(:, kopt) + d
fval(knew) = f

! KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW) = FVAL(KOPT) and KNEW < KOPT, then
! MINLOC(FVAL) = KNEW /= KOPT. Do not change KOPT in this case.
if (fval(knew) < fval(kopt)) then
    kopt = knew
end if

! Even if FVAL(KNEW) >= FVAL(KOPT) and KOPT remains unchanged, we still need to update XOPT and FOPT,
! because it may happen that KNEW = KOPT, so that XPT(:, KOPT) has been updated to XNEW = XOPT + D.
xopt = xpt(:, kopt)
fopt = fval(kopt)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(abs(f - fval(knew)) <= ZERO, 'F == FVAL(KNEW)', srname)
    call assert(abs(fopt - fval(kopt)) <= ZERO, 'FOPT == FVAL(KOPT)', srname)
    call assert(norm(xopt - xpt(:, kopt)) <= ZERO, 'XOPT == XPT(:, KOPT)', srname)
    call assert(.not. any(fval < fopt), '.NOT. ANY(FVAL < FOPT)', srname)
end if

end subroutine updatexf


subroutine tryqalt(idz, fval, ratio, bmat, zmat, itest, gq, hq, pq)
!--------------------------------------------------------------------------------------------------!
! TRYQALT tests whether to replace Q by the alternative model, namely the model that minimizes
! the F-norm of the Hessian subject to the interpolation conditions. It does the replacement
! when certain criteria are satisfied (i.e., when ITEST = 3). See Section 8 of the NEWUOA paper.
! N.B.: Indeed, we only need BMAT(:, KNEW) instead of the entire matrix.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : inprod, matprod, issymmetric

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
real(RP) :: fz(size(zmat, 2))
real(RP) :: galt(size(gq))

! Sizes
n = int(size(gq), kind(n))
npt = int(size(pq), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(.not. is_nan(ratio), 'RATION is not NaN', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
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
        itest = int(itest + 1, kind(itest))
    end if
end if

! Replace Q with Q_alt when ITEST >= 3.
if (itest >= 3) then
    gq = galt
    hq = ZERO
    fz = matprod(fval, zmat)
    fz(1:idz - 1) = -fz(1:idz - 1)
    pq = matprod(zmat, fz)
    itest = 0_IK
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(issymmetric(hq), 'HQ is symmetric', srname)
end if

end subroutine tryqalt


end module update_mod
