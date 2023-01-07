module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of LINCOA.
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
! Last Modified: Monday, December 26, 2022 PM05:27:51
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstep, trrad


contains


subroutine trstep(amat, delta, gq_in, hq_in, pq_in, rescon, xpt, iact, nact, qfac, rfac, ngetact, s)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, TWO, HALF, EPS, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm, solve, isorth, istriu, &
    & issymmetric, trueloc
use, non_intrinsic :: powalg_mod, only : hess_mul

! Solver-specific modules
use, non_intrinsic :: getact_mod, only : getact

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: delta
real(RP), intent(in) :: gq_in(:)  ! GQ(N)
real(RP), intent(in) :: hq_in(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq_in(:)  ! PQ(NPT)
real(RP), intent(in) :: rescon(:)  ! RESCON(M)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M); Will be updated in GETACT
integer(IK), intent(inout) :: nact  ! Will be updated in GETACT
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N); Will be updated in GETACT
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N); Will be updated in GETACT

! Outputs
integer(IK), intent(out) :: ngetact
real(RP), intent(out) :: s(:)  ! S(N)

! Local variables
logical :: newact
character(len=*), parameter :: srname = 'TRSTEP'
real(RP) :: d(size(gq_in))
real(RP) :: dproj(size(gq_in))
real(RP) :: psd(size(gq_in))
real(RP) :: hd(size(gq_in))
real(RP) :: pg(size(gq_in))
real(RP) :: resact(size(amat, 2))
real(RP) :: resnew(size(amat, 2))
real(RP) :: tol
real(RP) :: g(size(gq_in))
real(RP), parameter :: ctest = 0.01_RP  ! Convergence test parameter.
real(RP) :: frac(size(amat, 2)), ad(size(amat, 2)), restmp(size(amat, 2)), alpha, alphm, alpht, &
& beta, dd, dg, dhd, ds, gamma, reduct, resid, scaling, delsq, ss, temp, sold(size(s))
integer(IK) :: maxiter, iter, itercg, jsav
integer(IK) :: m
integer(IK) :: n
integer(IK) :: npt
real(RP) :: pq(size(pq_in)), gq(size(gq_in)), hq(size(hq_in, 1), size(hq_in, 2))

scaling = maxval(abs(gq_in))
if (scaling > 1.0E12) then
    gq = gq_in * max(TWO * tiny(scaling), ONE / scaling)
    pq = pq_in * max(TWO * tiny(scaling), ONE / scaling)
    hq = hq_in * max(TWO * tiny(scaling), ONE / scaling)
else
    gq = gq_in
    pq = pq_in
    hq = hq_in
end if

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(gq_in), kind(n))
npt = int(size(pq_in), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(size(hq_in, 1) == n .and. issymmetric(hq_in), 'HQ is n-by-n and symmetric', srname)

    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call wassert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
end if

!====================!
! Calculation starts !
!====================!

g = gq

! Return if G is not finite. Otherwise, GETACT will fail in the debugging mode.
if (.not. is_finite(sum(abs(g)))) then
    ngetact = 0
    s = ZERO
    return
end if

!
!     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
!       are the same as the terms with these names in LINCOB. If RESCON(J)
!       is negative, then |RESCON(J)| must be no less than the trust region
!       radius, so that the J-th constraint can be ignored.
!     S is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current S.
!     RESNEW, RESACT, D, PG are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       zero value of RESNEW(J) indicates that the J-th constraint is active,
!       and otherwise RESNEW(J) is set to the greater of TINYCV and the actual
!       residual of the J-th constraint for the current S. RESACT holds
!       the residuals of the active constraints, which may be positive.
!       D is the search direction of each line search.
!
!     Set some numbers for the conjugate gradient iterations.
!
delsq = delta * delta

! Set the initial elements of RESNEW, RESACT and S.

! 1. RESNEW(J) < 0 indicates that the J-th constraint does not restrict the CG steps of the current
! trust region calculation. In other words, RESCON >= DELTA.
! 2. RESNEW(J) = 0 indicates that J is an entry of IACT(1:NACT).
! 3. RESNEW(J) > 0 means that RESNEW(J) = max(B(J) - AMAT(:, J)^T*(XOPT+S), TINYCV), where S is the
! step up to now, calculated by a sequence of (truncated) CG iterations.
! N.B.: The order of the following lines is important, as the later ones override the earlier.
resnew = rescon
resnew(trueloc(rescon >= 0)) = max(TINYCV, rescon(trueloc(rescon >= 0)))
resnew(trueloc(rescon >= delta)) = -ONE
!!MATLAB:
!!resnew = rescon; resnew(rescon >= 0) = max(TINYCV, rescon(rescon >= 0)); resnew(rescon >= delta) = -1;
resnew(iact(1:nact)) = ZERO

! RESACT contains the constraint residuals of the constraints in IACT(1:NACT), namely the values
! of B(J) - AMAT(:, J)^T*(XOPT+S) for the J in IACT(1:NACT). Here, IACT(1:NACT) is a set of
! indicates such that the columns of AMAT(:, IACT(1:NACT)) form a basis of the constraint gradients
! in the "active set". For the definition of the "active set", see (3.5) of Powell (2015) and the
! comments at the beginning of the GETACT subroutine.
! N.B.: Between two calls of GETACT, S is updated in the orthogonal complement of the "active"
! gradients (i.e., null space of the "active" constraints). Therefore, RESACT remains unchanged.
! RESACT is changed right after GETACT is called if the first search direction D is not PSD
! but PSD + GAMMA * DPROJ.
resact(1:nact) = rescon(iact(1:nact))

s = ZERO
ss = ZERO
reduct = ZERO
ngetact = 0
newact = .true.

! ITERCG is the number of CG iterations corresponding to the current "active set" obtained by
! calling GETACT. These CG iterations are restricted in the orthogonal complement of the active
! gradients (i.e., null space of the active constraints).
! The following initial value of ITERCG is an artificial value that is not used. It is to entertain
! Fortran compilers (can it be be removed?).
itercg = -1_IK

maxiter = min(10000_IK, 10_IK * (m + n))  ! What is the theoretical upper bound of ITER?
do iter = 1, maxiter  ! Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
    if (newact) then
        ! GETACT picks the active set for the current S. It also sets PSD to the vector closest to
        ! -G that is orthogonal to the normals of the active constraints. PSD is scaled to have
        ! length 0.2*DELTA. Then a move of PSD from S is allowed by the linear constraints: PSD
        ! reduces the values of the nearly active constraints; it changes the inactive constraints
        ! by at most 0.2*DELTA, but the residuals of these constraints at no less than 0.2*DELTA.
        ngetact = ngetact + 1
        call getact(amat, delta, g, iact, nact, qfac, resact, resnew, rfac, psd)
        dd = inprod(psd, psd)
        !if (dd <= 0 .or. is_nan(dd)) then
        if (dd <= EPS * delsq .or. is_nan(dd)) then
            exit
        end if
        scaling = 0.2_RP * delta / sqrt(dd)
        psd = scaling * psd

        ! If the modulus of the residual of an "active constraint" is substantial (i.e., more than
        ! 1.0E-4*DELTA), then modify the searching direction PSD by a projection step to the
        ! boundaries of the "active constraint". This modified step will reduce the constraint
        ! residuals of the "active constraints" (see the update of RESACT below). The motivation is
        ! that the constraints in the "active set" are presumed to be active, and hence should have
        ! zero residuals (no constraint is violated, as the current method is feasible). According
        ! to a test on 20220821, this modification is important for the performance of LINCOA.
        ! N.B.:
        ! 1. The residual of the constraint A*X <= B is defined as B - A*X. It is not the constraint
        ! violation. Indeed, the constraint violations of the iterates are 0 in the current method.
        ! 2. We prefer `ANY(X > Y)` to `MAXVAL(X) > Y`, as Fortran standards do not specify
        ! MAXVAL(X) when X contains NaN, and MATLAB/Python/R/Julia behave differently in this
        ! respect. Moreover, MATLAB defines max(X) = [] if X == [], differing from mathematics
        ! and other languages.
        gamma = ZERO  ! The steplength of the projection step to be taken.
        if (any(resact(1:nact) > 1.0E-4_RP * delta)) then
            ! Set DPROJ to the shortest move (projection step) from S to the boundaries of the
            ! active constraints. We will use DPROJ to modify PSD.
            dproj = matprod(qfac(:, 1:nact), solve(transpose(rfac(1:nact, 1:nact)), resact(1:nact)))
            !!MATLAB: dproj = qfac(:, 1:nact) * (rfac(1:nact, 1:nact)' \ resact(1:nact))

            ! The vector DPROJ is also the shortest move from S + PSD to the boundaries of the
            ! active constraints (this is because PSD is parallel to the boundaries of the active
            ! constraints). Set GAMMA to the greatest steplength of this move that satisfies both
            ! the trust region bound and the linear constraints.
            ds = inprod(dproj, s + psd)
            dd = sum(dproj**2)
            resid = delsq - sum((s + psd)**2)
            !if (resid > 0) then
            if (resid > 0 .and. dd > EPS * delsq .and. .not. is_nan(ds)) then
                ! Set GAMMA to the greatest value so that S + PSD + GAMMA*DPROJ satisfies the trust
                ! region bound.
                !temp = sqrt(ds * ds + dd * resid)
                temp = maxval([abs(ds), sqrt(dd * resid), sqrt(ds * ds + dd * resid)])
                if (ds <= 0) then
                    gamma = (temp - ds) / dd
                else
                    gamma = resid / (temp + ds)
                end if
                if (gamma <= 0 .or. .not. is_finite(gamma)) gamma = ZERO
                ! Reduce GAMMA so that the move along DPROJ also satisfies the linear constraints.
                ad = -ONE
                ad(trueloc(resnew > 0)) = matprod(dproj, amat(:, trueloc(resnew > 0)))
                frac = ONE
                restmp(trueloc(ad > 0)) = resnew(trueloc(ad > 0)) - matprod(psd, amat(:, trueloc(ad > 0)))
                frac(trueloc(ad > 0)) = restmp(trueloc(ad > 0)) / ad(trueloc(ad > 0))
                gamma = minval([gamma, ONE, frac])  ! GAMMA = MINVAL([GAMMA, ONE, FRAC(TRUELOC(AD>0))])
            end if
        end if

        ! Set the next direction for seeking a reduction in the model function subject to the trust
        ! region bound and the linear constraints.
        ! Do NOT write D = PSD + GAMMA*DPROJ, as DPROJ may contain NaN/Inf, in which case GAMMA = 0.
        if (gamma > 0) then
            d = psd + gamma * dproj  ! Modified searching direction.
            itercg = -1_IK
        else
            d = psd  ! Original searching direction.
            itercg = 0_IK
        end if
    end if

    ! Set ALPHA to the steplength from S along D to the trust region boundary. Return if the first
    ! derivative term of this step is sufficiently small or if no further progress is possible.
    itercg = itercg + 1
    ! After the above line, ITERCG = 0 iff GETACT has been just called, and D is not PSD but a
    ! modified step.
    resid = delsq - ss
    if (resid <= 0) then
        exit
    end if
    dg = inprod(d, g)
    ds = inprod(d, s)
    dd = inprod(d, d)
    !if (dg >= 0 .or. is_nan(dg)) then
    !if (dg >= -EPS * sqrt(dd) * sqrt(sum(g**2)) .or. is_nan(dg) .or. dd <= EPS**2 .or. is_nan(dd)) then
    if (dd <= EPS * delsq .or. is_nan(ds)) then
        exit
    end if
    !temp = sqrt(ds * ds + dd * resid)
    temp = maxval([abs(ds), sqrt(dd * resid), sqrt(ds * ds + dd * resid)])
    if (ds <= 0) then
        alpha = (temp - ds) / dd
    else
        alpha = resid / (temp + ds)
    end if
    if (alpha <= 0 .or. .not. is_finite(alpha)) exit
    !if (-alpha * dg <= ctest * reduct) then
    if (-alpha * dg <= ctest * reduct .or. is_nan(alpha * dg)) then
        exit
    end if

    hd = hess_mul(d, xpt, pq, hq)

    ! Set DHD to the curvature of the model along D. Then reduce ALPHA if necessary to the value
    ! that minimizes the model.
    dhd = inprod(d, hd)
    alpht = alpha
    if (dg + alpha * dhd > 0) then
        alpha = -dg / dhd
    end if

    ! Make a further reduction in ALPHA if necessary to preserve feasibility.
    alphm = alpha
    ad = -ONE
    ad(trueloc(resnew > 0)) = matprod(d, amat(:, trueloc(resnew > 0)))
    frac = alpha
    frac(trueloc(ad > 0)) = resnew(trueloc(ad > 0)) / ad(trueloc(ad > 0))
    frac(trueloc(is_nan(frac))) = alpha
    jsav = 0_IK
    if (any(frac < alpha)) then
        jsav = int(minloc(frac, dim=1), IK)
        alpha = frac(jsav)
    end if
    !----------------------------------------------------------------------------------------------!
    ! Alternatively, JSAV and ALPHA can be calculated as below.
    !JSAV = INT(MINLOC([ALPHA, FRAC], DIM=1), IK) - 1_IK
    !ALPHA = MINVAL([ALPHA, FRAC])  ! This line cannot be exchanged with the last.
    ! We prefer our implementation as the code is more explicit; in addition, it is more flexible:
    ! we can change the condition ANY(FRAC < ALPHA) to ANY(FRAC < (1 - EPS) * ALPHA) or
    ! ANY(FRAC < (1 + EPS) * ALPHA), depending on whether we believe a false positive or a false
    ! negative of JSAV > 0 is more harmful.
    !----------------------------------------------------------------------------------------------!

    ! Post-process ALPHA according to some prior information.
    ! N.B.:
    ! 1. Since we set ALPHA=1 when ITERCG=0, the ALPHA calculated above is needed only if ITERCG>0.
    ! 2. Zaikun 20220821: In theory, shouldn't this post-processing change nothing? According to
    ! a test on 20220821, it does change ALPHA sometimes. Strange! Why?
    if (itercg == 0) then  ! Iff GETACT has been called, and D is not PSD but a modified step.
        ! By the definition of D, ALPHA = ONE is the largest ALPHA so that S + ALPHA*D satisfies the
        ! linear and trust region constraints.
        alpha = ONE
    elseif (itercg == 1 .and. gamma <= 0) then  ! Iff GETACT has been called, and D is not modified.
        ! Due to the scaling of PSD, S+D satisfies the linear and trust region constraints.
        alpha = max(alpha, ONE)
    else
        alpha = max(alpha, ZERO)
    end if

    ! Set ALPHA to the minimum between ALPHA and ALPHM, namely the steplength obtained by minimizing
    ! the quadratic model along D.
    alpha = min(alpha, alphm)

    ! Update S, G.
    sold = s
    s = s + alpha * d
    ss = sum(s**2)
    if (.not. is_finite(ss)) then
        s = sold
        exit
    end if
    g = g + alpha * hd
    if (.not. is_finite(sum(abs(g)))) then
        exit
    end if

    ! Update RESNEW.
    restmp = resnew - alpha * ad  ! Only RESTMP(TRUELOC(RESNEW > 0)) is needed.
    resnew(trueloc(resnew > 0)) = max(TINYCV, restmp(trueloc(resnew > 0)))
    !!MATLAB: mask = (resnew > 0); resnew(mask) = max(TINYCV, resnew(mask) - alpha * ad(mask));

    ! Update RESACT. This is done iff GETACT has been called, and D is not PSD but a modified step.
    !----------------------------------------------------------------------------------------------!
    ! Zaikun 20220821: There seems be a typo here. Powell's original code does not take ALPHA into
    ! account. Then RESACT seems to correspond to S + D, where D is defined as PSD + GAMMA*DPROJ
    ! during the modification procedure after GETACT is called. Without this modification, RESACT
    ! would remain unchanged because D = PSD, which is in the null space of the active constraints.
    ! The GAMMA*DPROJ component in the modified step D reduces RESACT by GAMMA*RESACT. However,
    ! since S is updated to S + ALPHA*D, shouldn't RESACT be reduced by ALPHA*GAMMA*RESACT?
    ! Note that Powell chose to update RESACT after ALPHA is calculated (instead of right after
    ! GAMMA is calculated), which might be an indication that he wanted to take ALPHA into account.
    ! In the following code, we try correcting this apparent typo, but it has little impact on the
    ! performance of LINCOA according to a test on 20220821.
    if (itercg == 0) then
        resact(1:nact) = (ONE - alpha * gamma) * resact(1:nact)
        !resact(1:nact) = (ONE - gamma) * resact(1:nact)  ! Powell's code.
    end if
    !----------------------------------------------------------------------------------------------!

    ! Update REDUCT, the reduction up to now.
    reduct = reduct - alpha * (dg + HALF * alpha * dhd)
    if (reduct <= 0 .or. is_nan(reduct)) then
        s = sold
        exit
    end if

    ! Test for termination.
    if (alpha >= alpht .or. -alphm * (dg + HALF * alphm * dhd) <= ctest * reduct) then
        exit
    end if

    ! Branch to a new loop if there is a new active constraint.
    ! When JSAV > 0, Powell's code branches back with NEWACT = .TRUE. only if |S| <= 0.8*DELTA, and
    ! it exits if |S| > 0.8*DELTA, as mentioned at the end of Section 3 of Powell 2015. The
    ! motivation seems to avoid small steps that changes the active set, because GETACT is expensive
    ! in flops. However, according to a test on 20220820, removing this condition (essentially
    ! replacing it with |S| < DELTA) improves the performance of LINCOA a bit. This may lead to
    ! small steps, but tiny steps will lead to tiny reductions and trigger an exit.
    newact = (jsav > 0)
    if (newact) then
        cycle
    end if

    ! If N-NACT CG iterations has been taken in the current null space (corresponding to the
    ! current "active set"), then, in theory, a stationary point in this subspace has been found.
    ! If the "active set" is the true active set, then a stationary point of the
    ! linearly-constrained trust region subproblem is found. So a termination is reasonable.
    ! However, the "active set" is not precisely the true active set, is it? See (3.5) of Powell
    ! (2015) and the comments at the beginning of the GETACT subroutine. Also, we should take into
    ! account the modification after GETACT is called.
    if (itercg >= n - nact) then  ! ITERCG > N - NACT is impossible.
        exit
    end if

    ! Calculate the next search direction, which is conjugate to the previous one if ITERCG /= NACT.
    ! N.B.: NACT < 0 is impossible unless GETACT is buggy; NACT = 0 can happen, particularly if
    ! there is no constraint. In theory, the code for the second case below covers the first as well.
    if (nact <= 0) then
        pg = g
    else
        pg = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
        !!MATLAB: pg = qfac(:, nact+1:n) * (g' * qfac(:, nact+1:n))';
    end if

    if (itercg == 0) then  ! Iff GETACT has been called, and D is not PSD but a modified step.
        beta = ZERO
    else
        beta = inprod(pg, hd) / dhd
    end if
    d = -pg + beta * d
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(s) == n .and. all(is_finite(s)), 'SIZE(S) == N, S is finite', srname)
    ! Due to rounding, it may happen that |S| > DELTA, but |S| > 2*DELTA is highly improbable.
    call assert(norm(s) <= TWO * delta, '|S| <= 2*DELTA', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
end if

end subroutine trstep
!--------------------------------------------------------------------------------------------------!
! Zaikun 20220417:
! For PG, the schemes below work evidently worse than the one above in a test on 20220417. Why?
!-----------------------------------------------------------------------!
! VERSION 1:
! !pg = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
!-----------------------------------------------------------------------!
! VERSION 2:
! !if (2 * nact < n) then
! !    pg = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
! !else
! !    pg = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
! !end if
!-----------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


function trrad(delta_in, dnorm, eta1, eta2, gamma1, gamma2, gamma3, ratio) result(delta)
!--------------------------------------------------------------------------------------------------!
! This function updates the trust region radius according to RATIO and DNORM.
!--------------------------------------------------------------------------------------------------!

! Generic module
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Input
real(RP), intent(in) :: delta_in   ! Current trust-region radius
real(RP), intent(in) :: dnorm   ! Norm of current trust-region step
real(RP), intent(in) :: eta1    ! Ratio threshold for contraction
real(RP), intent(in) :: eta2    ! Ratio threshold for expansion
real(RP), intent(in) :: gamma1  ! Contraction factor
real(RP), intent(in) :: gamma2  ! Expansion factor
real(RP), intent(in) :: gamma3  ! Expansion factor
real(RP), intent(in) :: ratio   ! Reduction ratio

! Outputs
real(RP) :: delta

! Local variables
character(len=*), parameter :: srname = 'TRRAD'

! Preconditions
if (DEBUGGING) then
    call assert(delta_in >= dnorm .and. dnorm > 0, 'DELTA_IN >= DNORM > 0', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma3 > 1 .and. gamma2 >= gamma3, &
        & '0 < GAMMA1 < 1 < GAMMA3 <= GAMMA2', srname)
    ! By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
    ! NaN, which should NOT happen due to the moderated extreme barrier.
    call assert(.not. is_nan(ratio), 'RATIO is not NaN', srname)
end if

!====================!
! Calculation starts !
!====================!

if (ratio <= eta1) then
    delta = gamma1 * dnorm  ! Powell's UOBYQA/NEWUOA.
    !delta = gamma1 * delta_in  ! Powell's COBYLA/LINCOA.
    !delta = min(gamma1 * delta_in, dnorm)  ! Powell's BOBYQA.
elseif (ratio <= eta2) then
    delta = max(gamma1 * delta_in, dnorm)   ! Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA
else
    delta = min(max(gamma1 * delta_in, gamma2 * dnorm), gamma3 * delta_in)
    !delta = max(gamma1 * delta_in, gamma2 * dnorm)  ! Powell's NEWUOA/BOBYQA.
    !delta = max(delta_in, 1.25_RP * dnorm, dnorm + rho)  ! Powell's UOBYQA
    !delta = max(delta_in, gamma2 * dnorm)  ! Modified version. Works well for UOBYQA.
end if

! For noisy problems, the following may work better.
! !if (ratio <= eta1) then
! !    delta = gamma1 * dnorm
! !elseif (ratio <= eta2) then  ! Ensure DELTA >= DELTA_IN
! !    delta = delta_in
! !else  ! Ensure DELTA > DELTA_IN with a constant factor
! !    delta = max(delta_in * (1.0_RP + gamma2) / 2.0_RP, gamma2 * dnorm)
! !end if
!

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(delta > 0, 'DELTA > 0', srname)
end if

end function trrad


end module trustregion_mod
