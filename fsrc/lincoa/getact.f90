module getact_mod
!--------------------------------------------------------------------------------------------------!
! This module provides the GETACT subroutine of LINCOA, which picks the current active set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, March 30, 2022 AM12:14:40
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, psd)
!--------------------------------------------------------------------------------------------------!
!!!----------------------------------------------------------!!!
! !!! THE FOLLOWING DESCRIPTION NEEDS VERIFICATION !!!
! !!! Note that the set J gets updated within this subroutine,
! !!! which seems inconsistent with the description given here.
!!!----------------------------------------------------------!!!
!
! This subroutine solves a linearly constrained projected problem (LCPP)
!
! min ||D + G|| subject to AMAT(:, J)^T * D <= 0.
!
! The solution is PSD, which is a projected steepest descent direction PSD for a linearly
! constrained trust-region subproblem (LCTRS)
!
! min Q(X_k + D)  subject to ||D|| <= Delta_k and AMAT^T*(X_k + D) <= B,
!
! where X_k is in R^N, B is in R^M, and AMAT is in R^{NxM}.
!
! In (LCPP), J is the index set defined in (3.3) of Powell (2015) as
!
! J = {j : B_j - A_j^T*Y <= 0.2*Delta_k*||A_j||, 1 <= j <= M} with A_j = AMAT(:, j),
!
! i.e., the index set of constraints in (LCTRS) that are nearly active (as per Powell, j is in J if
! and only if the distance from Y to the boundary of the j-th constraint is at most 0.2*Delta_k).
! Here, Y is the point where G is taken, namely G = grad Q(Y). Y is not necessarily X_k, but an
! iterate of the algorithm (e.g., truncated conjugate gradient) that solves (LCTRS). In LINCOA,
! ||A_j|| is 1 as the gradients of the linear constraints are normalized before LINCOA starts.
!
! The subroutine solves (LCPP) by the active set method of Goldfarb-Idnani 1983. It does not only
! calculate PSD, but also identify the active set of (LCPP) at the solution PSD, and maintains a QR
! factorization of A corresponding to the active set. More specifically, IACT(1:NACT) is a set of
! indices such that the columns of AMAT(:, IACT(1:NACT)) constitute a basis of the active constraint
! gradients, ans QFAC*RFAC(:, 1:NACT) is the QR factorization of AMAT(:, IACT(1:NACT)) such that
!
! SIZE(QFAC) = [M, M], SIZE(RFAC, 1) = M, diag(RFAC(:, 1:NACT)) > 0.
!
! NACT, IACT, QFAC and RFAC across invocations of GETACT for warm starts.
!
! SNORM, RESNEW, RESACT, and G are the same as the terms with these names in SUBROUTINE TRSTEP.
! The elements of RESNEW and RESACT are also kept up to date.
!
! VLAM is the vector of Lagrange multipliers of the calculation.
!
! See Section 3 of Powell (2015) for more information.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, EPS, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, eye, istriu, isorth, norm, lsqr, trueloc

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: snorm

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)

! Outputs
real(RP), intent(out) :: psd(:)  ! PSD(N)

! Local variables
character(len=*), parameter :: srname = 'GETACT'
integer(IK) :: i
integer(IK) :: ic
integer(IK) :: iter
integer(IK) :: l
integer(IK) :: m
integer(IK) :: maxiter
integer(IK) :: n
logical :: mask(size(amat, 2))
real(RP) :: apsd(size(amat, 2))
real(RP) :: dd
real(RP) :: ddsav
real(RP) :: dnorm
real(RP) :: gg
real(RP) :: fracmult(size(g))
real(RP) :: psdsav(size(psd))
real(RP) :: tdel
real(RP) :: tol
real(RP) :: violmx
real(RP) :: vlam(size(g))
real(RP) :: vmu(size(g))
real(RP) :: vmult

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(g), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(psd) == n, 'SIZE(PSD) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

! Quick return when M = 0.
if (m <= 0) then
    nact = 0_IK
    qfac = eye(n)
    psd = -g
    return
end if

! Set some constants.
gg = inprod(g, g)
tdel = 0.2_RP * snorm

! Set the initial QFAC to the identity matrix in the case NACT = 0.
if (nact == 0) then
    qfac = eye(n)
end if

! Remove any constraints from the initial active set whose residuals exceed TDEL.
! Compilers may complain if we do not set VLAM. The value does not matter, as it will be overwritten.
vlam = ZERO
do ic = nact, 1, -1
    if (resact(ic) > tdel) then
        ! Delete constraint IACT(IC) from the active set, and set NACT = NACT - 1.
        call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
    end if
end do

! Remove any constraints from the initial active set whose Lagrange multipliers are nonnegative,
! and set the surviving multipliers.
! The following loop will run for at most NACT times, since each call of DEL_ACT reduces NACT by 1.
do while (nact > 0)
    vlam(1:nact) = lsqr(g, qfac(:, 1:nact), rfac(1:nact, 1:nact))
    if (any(vlam(1:nact) >= 0)) then
        ic = maxval(trueloc(vlam(1:nact) >= 0))  ! MATLAB: ic = max(find(vlam(1:nact) >= 0))
        call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
    else
        exit
    end if
end do

! Set the new search direction D. Terminate if the 2-norm of D is ZERO or does not decrease, or if
! NACT=N holds. The situation NACT=N occurs for sufficiently large SNORM if the origin is in the
! convex hull of the constraint gradients.
! Start with initialization of  PSD, PSDSAV, DD, and DDSAV.
psd = ZERO  ! Must be set, in case NACT = N at this point.
psdsav = ZERO
dd = ZERO
ddsav = TWO * gg  ! By Powell. This value is used at iteration 1 to test whether DD >= DDSAV. Why?

! What is the theoretical maximal number of iterations in the following procedure? Powell's code for
! this part is essentially a `DO WHILE (NACT < N) ... END DO` loop. We enforce the following maximal
! number of iterations, which is never reached in our tests (indeed, even 2*N cannot be reached).
maxiter = 2_IK * (m + n)
do iter = 1_IK, maxiter
    ! When NACT == N, exit with PSD = 0. Indeed, with a correctly implemented matrix product, the
    ! lines below this IF should render DD = 0 and trigger the exit. We do it explicitly for clarity.
    if (nact >= n) then  ! Indeed, NACT > N should never happen.
        psd = ZERO
        exit
    end if

    psd(nact + 1:n) = matprod(g, qfac(:, nact + 1:n))
    psd = -matprod(qfac(:, nact + 1:n), psd(nact + 1:n)) ! Projection of -G to range(QFAC(:,NACT+1:N))
    dd = inprod(psd, psd)
    dnorm = sqrt(dd)

    if (dd == ZERO) then
        exit
    end if

    if (dd >= ddsav) then
        psd = ZERO  ! Zaikun 20220329: Powell wrote this. Why? What about PSD = PSDSAV?
        exit
    end if

    !---------------------------------------------------------------------------------------!
    ! Powell's code does not handle the following pathological cases.
    if (inprod(psd, g) > 0 .or. .not. is_finite(sum(abs(psd)))) then
        psd = psdsav
        exit
    end if
    ! In our tests, tolerating the following cases seems to render better numerical results.
    !!if (dd > gg) then
    !!    psd = (psd / dnorm) * sqrt(gg)
    !!    exit
    !!end if
    !!if (inprod(psd, g) < -gg) then
    !!    exit
    !!end if
    !---------------------------------------------------------------------------------------!

    psdsav = psd
    ddsav = dd

    ! Pick the next integer L or terminate; a positive L is the index of the most violated constraint.
    apsd = matprod(psd, amat)
    mask = (resnew > 0 .and. resnew <= tdel .and. apsd > (dnorm / snorm) * resnew)
    if (any(mask)) then
        l = int(maxloc(apsd, mask=mask, dim=1), IK)
        violmx = apsd(l)
        ! MATLAB: apsd(mask) = -Inf; [violmx , l] = max(apsd);
    else
        exit
    end if

    ! Terminate if a positive value of VIOLMX may be due to computer rounding errors.
    ! N.B.: 0. Powell wrote VIOLMX < 0.01*DNORM instead of VIOLMX <= 0.01*DNORM.
    ! 1. Theoretically (but not numerically), APSD(IACT(1:NACT)) = 0 or empty.
    ! 2. CAUTION: the inf-norm of APSD(IACT(1:NACT)) is NOT always MAXVAL(ABS(APSD(IACT(1:NACT)))),
    ! as the latter returns -HUGE(APSD) instead of 0 when NACT = 0! In MATLAB, max([]) = []; in
    ! Python, R, and Julia, the maximum of an empty array raises errors/warnings (as of 20220318).
    if (violmx <= 0.01_RP * dnorm .and. violmx <= TEN * norm(apsd(iact(1:nact)), 'inf')) then
        exit
    end if

    ! Add constraint L to the active set. ADD_ACT sets NACT = NACT + 1 and VLAM(NACT) = 0.
    call add_act(l, amat(:, l), iact, nact, qfac, resact, resnew, rfac, vlam)

    ! Set the components of the vector VMU if VIOLMX is positive.
    ! N.B.: 1. In theory, NACT > 0 is not needed in the condition below, because VIOLMX must be 0
    ! when NACT is 0. We keep NACT > 0 for security: when NACT <= 0, RFAC(NACT, NACT) is invalid.
    ! 2. The loop will run for at most NACT <= N times: if VIOLMX > 0, then IC > 0, and hence
    ! VLAM(IC) = 0, which implies that DEL_ACT will be called to reduce NACT by 1.
    do while (violmx > 0 .and. nact > 0)
        !------------------------------------------------------------------------------------------!
        ! Zaikun 20220329: What is VMU exactly???
        ! VMU(1:NACT) = LSQR(QFAC(:, 1:NACT),RFAC(1:NACT, 1:NACT),QFAC(:, NACT)) / RFAC(NACT, NACT)?
        vmu(nact) = ONE / rfac(nact, nact)**2  ! We must ensure NACT > 0. In theory, VMU(NACT) > 0.
        do i = nact - 1, 1, -1
            vmu(i) = -inprod(rfac(i, i + 1:nact), vmu(i + 1:nact)) / rfac(i, i)
        end do
        !------------------------------------------------------------------------------------------!

        ! Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
        ! N.B.: 1. VLAM(1:NACT-1) < 0 and VLAM(NACT) <= 0 by the updates of VLAM. 2. VMU(NACT) > 0.
        ! 3. Only the places where VMU(1:NACT) < 0 is relevant below, if any.
        fracmult = HUGENUM
        where (vmu(1:nact) < 0 .and. vlam(1:nact) < 0)
            fracmult(1:nact) = vlam(1:nact) / vmu(1:nact)
        end where
        vmult = minval([violmx, fracmult(1:nact)])
        ic = 0_IK
        if (any(fracmult(1:nact) <= vmult)) then
            ic = maxval(trueloc(fracmult(1:nact) <= vmult))
            ! MATLAB: ic = max(find(fracmult(1:nact)] <= vmult))
        end if
        ! N.B.: 0. The definition of IC given above is mathematically equivalent to the following.
        !!IC = MAXVAL(TRUELOC([VIOLMX, FRACMULT(1:NACT)] <= VMULT)) - 1_IK, or:
        !!IC = INT(MINLOC([VIOLMX, FRACMULT(1:NACT)], DIM=1, BACK=.TRUE.), IK) - 1_IK
        ! However, such implementations are problematic in the unlikely case of VMULT = NaN: IC will
        ! be -Inf in the first and unspecified in the second. The MATLAB counterpart of the first
        ! implementation will render IC = [] as `find` (the MATLAB version of TRUELOC) returns [].
        ! 1. The BACK argument in MINLOC is available in F2008. Not supported by Absoft as of 2022.
        ! 2. A motivation for backward MINLOC is to save computation in DEL_ACT below (what else?).

        violmx = max(violmx - vmult, ZERO)
        vlam(1:nact) = vlam(1:nact) - vmult * vmu(1:nact)
        if (ic > 0 .and. ic <= nact) then  ! Powell's code: IF (IC>0). We check IC<=NACT for safety.
            vlam(ic) = ZERO
        end if

        ! Reduce the active set if necessary, so that all components of the new VLAM are negative,
        ! with resetting of the residuals of the constraints that become inactive.
        do ic = nact, 1, -1
            if (vlam(ic) >= 0) then  ! Powell's version: IF (.NOT. VLAM(IC) < 0) THEN
                ! Delete the constraint with index IACT(IC) from the active set; set NACT = NACT-1.
                call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
            end if
        end do
    end do  ! End of DO WHILE (VIOLMX > 0 .AND. NACT > 0)

    !----------------------------------------------------------------------------------------------!
    ! NACT can become 0 at this point iif VLAM(1:NACT) >= 0 before calling DEL_ACT, which is true
    ! if NACT = 1 when the WHILE loop starts. However, we have never observed a failure of the
    ! assertion below as of 20220329. Why?
    !-----------------------------------------!
    call assert(nact > 0, 'NACT > 0', srname) !
    !-----------------------------------------!
    if (nact == 0) then
        exit
    end if
    !----------------------------------------------------------------------------------------------!
end do  ! End of DO WHILE (NACT < N)

! NACT == 0 can happen here. The following lines improve the performance of LINCOA. Powell's code
! does not take care of this case explicitly.
if (nact == 0) then
    psd = -g
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    ! During the development, we want to get alerted if ITER reaches MAXITER.
    call assert(iter < maxiter, 'ITER < MAXITER', srname)

    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)! Can NACT be 0?
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(psd) == n, 'SIZE(PSD) == N', srname)
    ! In theory, |PSD|^2 <= GG and -GG <= PSD^T*G <= 0.
    ! N.B. 1. Do not use DD, which may not be up to date. 2. PSD^T*G can be NaN if G is huge.
    call assert(inprod(psd, psd) <= TWO * gg, '|PSD|^2 <= 2*GG', srname)
    call assert(.not. (inprod(psd, g) > 0 .or. inprod(psd, g) < -TWO * gg), '-2*GG <= PSD^T*G <= 0', srname)
end if

end subroutine getact


subroutine add_act(l, c, iact, nact, qfac, resact, resnew, rfac, vlam)
!--------------------------------------------------------------------------------------------------!
! This subroutine adds the constraint with index L to the active set as the (NACT+ )-th active
! constriant, updates IACT, QFAC, etc accordingly, and increments NACT to NACT+1. Here, C is the
! gradient of the new active constraint.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : qradd, istriu, isorth
implicit none

! Inputs
integer(IK), intent(in) :: l
real(RP), intent(in) :: c(:)  ! C(N)

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)
real(RP), intent(inout) :: vlam(:)  ! VLAM(N)

! Local variables (debugging only)
character(len=*), parameter :: srname = 'ADD_ACT'
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol

! Sizes
m = size(iact)
n = size(vlam)

! Preconditions
if (DEBUGGING) then
    call assert(m >= 1, 'M >= 1', srname)  ! Should not be called when M == 0.
    call assert(n >= 1, 'N >= 1', srname)
    call assert(nact >= 0 .and. nact <= min(m, n) - 1_IK, '0 <= NACT <= MIN(M, N)-1', srname)
    call assert(l >= 1 .and. l <= m, '1 <= L <= M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(.not. any(iact(1:nact) == l), 'L is not in IACT(1:NACT)', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

!====================!
! Calculation starts !
!====================!

! QRADD applies Givens rotations to the last (N-NACT) columns of QFAC so that the first (NACT+1)
! columns of QFAC are the ones required for the addition of the L-th constraint, and add the
! appropriate column to RFAC.
! N.B.: QRADD always augment NACT by 1. This is different from the strategy in COBYLA. Is it ensured
! that C cannot be linearly represented by the gradients of the existing active constraints?
call qradd(c, qfac, rfac, nact)  ! NACT is increased by 1!

! Indeed, it suffices to pass RFAC(:, 1:NACT+1) to QRADD as follows.
!!call qradd(c, qfac, rfac(:, 1:nact + 1), nact)  ! NACT is increased by 1!

! Update IACT, RESACT, RESNEW, and VLAM. N.B.: NACT has been increased by 1 in QRADD.
iact(nact) = l
resact(nact) = resnew(l)  ! RESACT(NACT) = RESNEW(IACT(NACT))
resnew(l) = ZERO  ! RESNEW(IACT(NACT)) = ZERO  ! Why not TINYCV? See DECACT.
vlam(nact) = ZERO

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact >= 1 .and. nact <= min(m, n), '1 <= NACT <= MIN(M, N)', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

end subroutine add_act


subroutine del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
!--------------------------------------------------------------------------------------------------!
! This subroutine deletes the constraint with index IACT(IC) from the active set, updates IACT,
! QFAC, etc accordingly, and reduces NACT to NACT-1.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, EPS, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : qrexc, isorth, istriu
implicit none

! Inputs
integer(IK), intent(in) :: ic

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)
real(RP), intent(inout) :: vlam(:)  ! VLAM(N)

! Local variables (debugging only)
character(len=*), parameter :: srname = 'DEL_ACT'
integer(IK) :: l
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol

! Sizes
m = size(iact)
n = size(vlam)

! Preconditions
! Preconditions
if (DEBUGGING) then
    call assert(m >= 1, 'M >= 1', srname)  ! Should not be called when M == 0.
    call assert(n >= 1, 'N >= 1', srname)
    call assert(nact >= 1 .and. nact <= min(m, n), '1 <= NACT <= MIN(M, N)', srname)
    call assert(ic >= 1 .and. ic <= nact, '1 <= IC <= NACT', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
    l = iact(ic)  ! For debugging only
end if

!====================!
! Calculation starts !
!====================!

! The following instructions rearrange the active constraints so that the new value of IACT(NACT) is
! the old value of IACT(IC). QREXC implements the updates of QFAC and RFAC by sequence of Givens
! rotations. Then NACT is reduced by one.

call qrexc(qfac, rfac(:, 1:nact), ic)  ! QREXC does nothing if IC == NACT.
! Indeed, it suffices to pass QFAC(:, 1:NACT) and RFAC(1:NACT, 1:NACT) to QREXC as follows. However,
! compilers may create a temporary copy of RFAC(1:NACT, 1:NACT), which is not contiguous in memory.
!!call qrexc(qfac(:, 1:nact), rfac(1:nact, 1:nact), ic)

iact(ic:nact) = [iact(ic + 1:nact), iact(ic)]
resact(ic:nact) = [resact(ic + 1:nact), resact(ic)]
resnew(iact(nact)) = max(resact(nact), TINYCV)
vlam(ic:nact) = [vlam(ic + 1:nact), vlam(ic)]
nact = nact - 1_IK

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact >= 0 .and. nact <= min(m, n) - 1, '1 <= NACT <= MIN(M, N)-1', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(.not. any(iact(1:nact) == l), 'L is not in IACT(1:NACT)', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

end subroutine del_act


end module getact_mod
