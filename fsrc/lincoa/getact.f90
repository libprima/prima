module getact_mod
!--------------------------------------------------------------------------------------------------!
! This module provides the GETACT subroutine of LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, June 23, 2022 PM10:24:06
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, delta, g, iact, nact, qfac, resact, resnew, rfac, psd)
!--------------------------------------------------------------------------------------------------!
!!!-------------------------------------------------------------!!!
!!! THE FOLLOWING DESCRIPTION NEEDS VERIFICATION!               !!!
!!! Note that the set J gets updated within this subroutine,    !!!
!!! which seems inconsistent with the description below.        !!!
!!! See the lines below "Pick the next integer L or terminate". !!!
!!!-------------------------------------------------------------!!!
!
! This subroutine solves a linearly constrained projected problem (LCPP)
!
! min ||D + G|| subject to AMAT(:, j)^T * D <= 0 for j in J.
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
! i.e., the index set of the nearly active constraints of (LCTRS) (as per Powell, j is in J if
! and only if the distance from Y to the boundary of the j-th constraint is at most 0.2*Delta_k).
! Here, Y is the point where G is taken, namely G = grad Q(Y). Y is not necessarily X_k, but an
! iterate of the algorithm (e.g., truncated conjugate gradient) that solves (LCTRS). In LINCOA,
! ||A_j|| is 1 as the gradients of the linear constraints are normalized before LINCOA starts.
!
! The subroutine solves (LCPP) by the active set method of Goldfarb-Idnani (1983). It does not only
! calculate PSD, but also identify the active set of (LCPP) at the solution PSD, namely
!
! I = {j in J : AMAT(:, j)^T*PSD = 0} (see 3.5 of Powell (2015)),
!
! and maintains a QR factorization of A corresponding to the active set. More specifically,
! IACT(1:NACT) is a set of indices such that the columns of AMAT(:, IACT(1:NACT)) constitute a basis
! of the active constraint gradients, ans QFAC*RFAC(:, 1:NACT) is the QR factorization of
! AMAT(:, IACT(1:NACT)) such that
!
! SIZE(QFAC) = [N, N], SIZE(RFAC, 1) = N, diag(RFAC(:, 1:NACT)) > 0.
!
! NACT, IACT, QFAC and RFAC are maintained up to date across invocations of GETACT for warm starts.
!
! SNORM, RESNEW, RESACT, and G are the same as the terms with these names in SUBROUTINE TRSTEP.
! The elements of RESNEW and RESACT are also kept up to date. See Section 3 of Powell (2015).
!
! VLAM is the vector of Lagrange multipliers of the calculation.
!
! See Section 3 of Powell (2015) for more information.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, EPS, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, eye, istriu, isorth, norm, lsqr, solve, trueloc

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: delta
real(RP), intent(in) :: g(:)  ! G(N)

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
integer(IK) :: icon
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
real(RP) :: frac(size(g))
real(RP) :: psdsav(size(psd))
real(RP) :: tdel
real(RP) :: tol
real(RP) :: v(size(g))
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
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
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
tdel = 0.2_RP * delta

! Set the initial QFAC to the identity matrix in the case NACT = 0.
if (nact == 0) then
    qfac = eye(n)
end if

! Remove any constraints from the initial active set whose residuals exceed TDEL.
! Compilers may complain if VLAM is not set. The value does not matter, as it will be overwritten.
vlam = ZERO
do icon = nact, 1, -1
    if (resact(icon) > tdel) then
        ! Delete constraint IACT(ICON) from the active set, and set NACT = NACT - 1.
        call delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam)
    end if
end do

! Remove any constraints from the initial active set whose Lagrange multipliers are nonnegative,
! and set the surviving multipliers.
! The following loop will run for at most NACT times, since each call of DEL_ACT reduces NACT by 1.
do while (nact > 0)
    vlam(1:nact) = lsqr(g, qfac(:, 1:nact), rfac(1:nact, 1:nact))
    if (.not. any(vlam(1:nact) >= 0)) then
        exit
    end if
    icon = maxval(trueloc(vlam(1:nact) >= 0))
    !!MATLAB: icon = max(find(vlam(1:nact) >= 0)); % OR: icon = find(vlam(1:nact) >= 0, 1, 'last')
    call delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam)
end do
! Zaikun 20220330: What if NACT = 0 at this point?

! Set the new search direction D. Terminate if the 2-norm of D is ZERO or does not decrease, or if
! NACT=N holds. The situation NACT=N occurs for sufficiently large SNORM if the origin is in the
! convex hull of the constraint gradients.
! Start with initialization of PSDSAV and DDSAV.
psdsav = ZERO  ! Must be set, in case the loop exits due to abnormality at iteration 1.
ddsav = TWO * gg  ! By Powell. This value is used at iteration 1 to test whether DD >= DDSAV. Why?

! What is the theoretical maximal number of iterations in the following procedure? Powell's code for
! this part is essentially a `DO WHILE (NACT < N) ... END DO` loop. We enforce the following maximal
! number of iterations, which is never reached in our tests (indeed, even 2*N cannot be reached).
! The iteration counter ITER never appears in the code of the iterations, as its purpose is merely
! to impose an upper bound on the number of iterations.
maxiter = 2_IK * (m + n)
do iter = 1_IK, maxiter
    ! When NACT == N, exit with PSD = 0. Indeed, with a correctly implemented matrix product, the
    ! lines below this IF should render DD = 0 and trigger an exit. We make it explicit for clarity.
    if (nact >= n) then  ! Indeed, NACT > N should never happen.
        psd = ZERO
        exit
    end if

    ! Set PSD to the projection of -G to range(QFAC(:,NACT+1:N))
    psd = -matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
    !!MATLAB: psd = -qfac(:, nact + 1:n) * (g' * qfac(:, nact + 1:n))';
    !----------------------------------------------------------------------------------------------!
    ! Zaikun: The schemes below work evidently worse than the one above in a test on 20220417. Why?
    !-------------------------------------------------------------------------!
    ! VERSION 1:
    !!psd = matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact))) - g
    !-------------------------------------------------------------------------!
    ! VERSION 2:
    !!if (2 * nact < n) then
    !!    psd = matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact))) - g
    !!else
    !!    psd = -matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
    !!end if
    !-------------------------------------------------------------------------!
    !----------------------------------------------------------------------------------------------!

    dd = inprod(psd, psd)
    dnorm = sqrt(dd)

    if (dd == ZERO) then
        exit
    end if

    if (dd >= ddsav) then
        psd = ZERO  ! Zaikun 20220329: Powell wrote this. Why?
        !psd = psdsav  ! This does not seem to improve the performance.
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
    !!    psd = (sqrt(gg) / dnorm) * psd
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
    mask = (resnew > 0 .and. resnew <= tdel .and. apsd > (dnorm / delta) * resnew)
    !----------------------------------------------------------------------------------------------!
    ! N.B.: the definition of L and VIOLMX can be simplified as follows, but we prefer explicitness.
    !L = INT(MAXLOC(APSD, MASK=MASK, DIM=1), IK) ! MAXLOC(...) = 0 if MASK is all FALSE.
    !VIOLMX = MAXVAL(APSD, MASK=MASK)  ! MAXVAL(...) = -HUGE(APSD) if MASK is all FALSE.
    if (any(mask)) then
        l = int(maxloc(apsd, mask=mask, dim=1), IK)
        violmx = apsd(l)
    else
        l = 0_IK
        violmx = -HUGENUM
    end if
    !!MATLAB: apsd(mask) = -Inf; [violmx , l] = max(apsd);
    ! N.B.: the value of L will differ from the Fortran version if MASK is all FALSE, but this is OK
    ! because VIOLMX will be -Inf, which will trigger the `exit` below. This is tricky. Be cautious!
    !----------------------------------------------------------------------------------------------!

    ! Terminate if VIOLMX <= 0 (when MASK contains only FALSE) or a positive value of VIOLMX may be
    ! due to computer rounding errors.
    ! N.B.: 0. `L <= 0` makes no difference in the condition below. Added for robustness.
    ! 1. Powell wrote VIOLMX < 0.01*DNORM instead of VIOLMX <= 0.01*DNORM.
    ! 2. Theoretically (but not numerically), APSD(IACT(1:NACT)) = 0 or empty.
    ! 3. CAUTION: the Inf-norm of APSD(IACT(1:NACT)) is NOT always MAXVAL(ABS(APSD(IACT(1:NACT)))),
    ! as the latter returns -HUGE(APSD) instead of 0 when NACT = 0! In MATLAB, max([]) = []; in
    ! Python, R, and Julia, the maximum of an empty array raises errors/warnings (as of 20220318).
    if (all(.not. mask) .or. violmx <= min(0.01_RP * dnorm, TEN * norm(apsd(iact(1:nact)), 'inf'))) then
        exit
    end if

    ! Add constraint L to the active set. ADD_ACT sets NACT = NACT + 1 and VLAM(NACT) = 0.
    call addact(l, amat(:, l), iact, nact, qfac, resact, resnew, rfac, vlam)

    ! Set the components of the vector VMU if VIOLMX is positive.
    ! N.B.: 1. In theory, NACT > 0 is not needed in the condition below, because VIOLMX must be 0
    ! when NACT is 0. We keep NACT > 0 for security: when NACT <= 0, RFAC(NACT, NACT) is invalid.
    ! 2. The loop will run for at most NACT <= N times: if VIOLMX > 0, then ICON > 0, and hence
    ! VLAM(ICON) = 0, which implies that DEL_ACT will be called to reduce NACT by 1.
    do while (violmx > 0 .and. nact > 0)
        v(1:nact - 1) = ZERO
        v(nact) = ONE / rfac(nact, nact) ! This is why we must ensure NACT > 0.
        ! Solve the linear system RFAC(1:NACT, 1:NACT) * VMU(1:NACT) = V(1:NACT) .
        vmu(1:nact) = solve(rfac(1:nact, 1:nact), v(1:nact)) ! VMU(NACT) = V(NACT)/RFAC(NACT,NACT)>0
        !!MATLAB: vmu(1:nact) = rfac(1:nact, 1:nact) \ v(1:nact);

        ! Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
        ! N.B.: 1. VLAM(1:NACT-1) < 0 and VLAM(NACT) <= 0 by the updates of VLAM. 2. VMU(NACT) > 0.
        ! 3. Only the places where VMU(1:NACT) < 0 is relevant below, if any.
        frac = HUGENUM
        where (vmu(1:nact) < 0 .and. vlam(1:nact) < 0) frac(1:nact) = vlam(1:nact) / vmu(1:nact)
        !!MATLAB: frac = vlam / vmu; frac(vmu >= 0 | vlam >= 0) = Inf;
        vmult = minval([violmx, frac(1:nact)])
        icon = maxval([0_IK, trueloc(frac(1:nact) <= vmult)])
        !!MATLAB: icon = max([0; find(frac(1:nact) <= vmult)]); % find(frac(1:nact)<=vmult) can be empty

        ! N.B.: 0. The definition of ICON given above is mathematically equivalent to the following.
        !!ICON = MAXVAL(TRUELOC([VIOLMX, FRACMULT(1:NACT)] <= VMULT)) - 1_IK, OR
        !!ICON = INT(MINLOC([VIOLMX, FRACMULT(1:NACT)], DIM=1, BACK=.TRUE.), IK) - 1_IK
        ! However, such implementations are problematic in the unlikely case of VMULT = NaN: ICON
        ! will be -Inf in the first and unspecified in the second. The MATLAB counterpart of the
        ! first implementation will render ICON = [] as `find` (the MATLAB version of TRUELOC)
        ! returns [].
        ! 1. The BACK argument in MINLOC is available in F2008. Not supported by Absoft as of 2022.
        ! 2. A motivation for backward MINLOC is to save computation in DEL_ACT below (what else?).

        violmx = max(violmx - vmult, ZERO)
        vlam(1:nact) = vlam(1:nact) - vmult * vmu(1:nact)
        if (icon > 0 .and. icon <= nact) then  ! Powell: IF (ICON>0). We check ICON<=NACT for safety.
            vlam(icon) = ZERO
        end if

        ! Reduce the active set if necessary, so that all components of the new VLAM are negative,
        ! with resetting of the residuals of the constraints that become inactive.
        do icon = nact, 1, -1
            if (vlam(icon) >= 0) then  ! Powell's version: IF (.NOT. VLAM(ICON) < 0) THEN
                ! Delete the constraint with index IACT(ICON) from the active set; set NACT = NACT-1.
                call delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam)
            end if
        end do
    end do  ! End of DO WHILE (VIOLMX > 0 .AND. NACT > 0)

    !----------------------------------------------------------------------------------------------!
    ! NACT can become 0 at this point iif VLAM(1:NACT) >= 0 before calling DEL_ACT, which is true
    ! if NACT happens to be 1 when the WHILE loop starts. However, we have never observed a failure
    ! of the assertion below as of 20220329. Why?
    !-----------------------------------------!
    call assert(nact > 0, 'NACT > 0', srname) !
    !-----------------------------------------!
    if (nact == 0) then
        exit
    end if
    !----------------------------------------------------------------------------------------------!
end do  ! End of DO WHILE (NACT < N)

! It is possible to have NACT == 0 here. The following lines improve the performance of LINCOA.
! Powell's code does not take care of this case explicitly.
if (nact == 0) then
    qfac = eye(n)
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
    ! PSD = -G when NACT == 0; G may contain Inf/NaN.
    call assert(all(is_finite(psd)) .or. nact == 0, 'PSD is finite unless NACT == 0', srname)
    ! In theory, |PSD|^2 <= GG and -GG <= PSD^T*G <= 0.
    ! N.B. 1. Do not use DD, which may not be up to date. 2. PSD^T*G can be NaN if G is huge.
    call assert(inprod(psd, psd) <= TWO * gg, '|PSD|^2 <= 2*GG', srname)
    call assert(.not. (inprod(psd, g) > 0 .or. inprod(psd, g) < -TWO * gg), '-2*GG <= PSD^T*G <= 0', srname)
end if

end subroutine getact


subroutine addact(l, c, iact, nact, qfac, resact, resnew, rfac, vlam)
!--------------------------------------------------------------------------------------------------!
! This subroutine adds the constraint with index L to the active set as the (NACT+ )-th active
! constriant, updates IACT, QFAC, etc accordingly, and increments NACT to NACT+1. Here, C is the
! gradient of the new active constraint.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : istriu, isorth
use, non_intrinsic :: powalg_mod, only : qradd
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
integer(IK) :: nsave
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
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
    nsave = nact  ! For debugging only
end if

!====================!
! Calculation starts !
!====================!

! QRADD applies Givens rotations to the last (N-NACT) columns of QFAC so that the first (NACT+1)
! columns of QFAC are the ones required for the addition of the L-th constraint, and add the
! appropriate column to RFAC.
! N.B.: QRADD always augment NACT by 1, which differs from the corresponding subroutine in COBYLA.
! Is it ensured that C cannot be represented by the gradients of the existing active constraints?
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
    call assert(nact == nsave + 1, 'NACT = NSAVE + 1', srname)
    call assert(nact >= 1 .and. nact <= min(m, n), '1 <= NACT <= MIN(M, N)', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

end subroutine addact


subroutine delact(icon, iact, nact, qfac, resact, resnew, rfac, vlam)
!--------------------------------------------------------------------------------------------------!
! This subroutine deletes the constraint with index IACT(ICON) from the active set, updates IACT,
! QFAC, etc accordingly, and reduces NACT to NACT-1.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, EPS, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : isorth, istriu
use, non_intrinsic :: powalg_mod, only : qrexc
implicit none

! Inputs
integer(IK), intent(in) :: icon

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
integer(IK) :: nsave
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
    call assert(icon >= 1 .and. icon <= nact, '1 <= ICON <= NACT', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
    nsave = nact  ! For debugging only
    l = iact(icon)  ! For debugging only
end if

!====================!
! Calculation starts !
!====================!

! The following instructions rearrange the active constraints so that the new value of IACT(NACT) is
! the old value of IACT(ICON). QREXC implements the updates of QFAC and RFAC by a sequence of Givens
! rotations. Then NACT is reduced by one.

call qrexc(qfac, rfac(:, 1:nact), icon)  ! QREXC does nothing if ICON == NACT.
! Indeed, it suffices to pass QFAC(:, 1:NACT) and RFAC(1:NACT, 1:NACT) to QREXC as follows. However,
! compilers may create a temporary copy of RFAC(1:NACT, 1:NACT), which is not contiguous in memory.
!!call qrexc(qfac(:, 1:nact), rfac(1:nact, 1:nact), icon)

iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
resact(icon:nact) = [resact(icon + 1:nact), resact(icon)]
resnew(iact(nact)) = max(resact(nact), TINYCV)
vlam(icon:nact) = [vlam(icon + 1:nact), vlam(icon)]
nact = nact - 1_IK

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact == nsave - 1, 'NACT = NSAVE - 1', srname)
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

end subroutine delact


end module getact_mod
