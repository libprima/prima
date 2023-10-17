! TODO: Implement GETMODEL to get the model of the objective function and constraints, i.e., g and A.
module cobylb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of COBYLA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
!
! N.B. (Zaikun 20220131): Powell's implementation of COBYLA uses RHO rather than DELTA as the
! trust-region radius, and RHO is never increased. DELTA does not exist in Powell's COBYLA code.
! Following the idea in Powell's other solvers (UOBYQA, ..., LINCOA), our code uses DELTA as the
! trust-region radius, while RHO works a lower bound of DELTA and indicates the current resolution
! of the algorithm. DELTA is updated in a classical way subject to DELTA >= RHO, whereas RHO is
! updated as in Powell's COBYLA code and is never increased. The new implementation improves the
! performance of COBYLA.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Monday, October 16, 2023 AM01:15:52
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobylb


contains


subroutine cobylb(calcfc, iprint, maxfilt, maxfun, amat, bvec, ctol, cweight, eta1, eta2, ftarget, &
    & gamma1, gamma2, rhobeg, rhoend, constr, f, x, nf, chist, conhist, cstrv, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the actual calculations of COBYLA.
!
! IPRINT, MAXFILT, MAXFUN, MAXHIST, CTOL, CWEIGHT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG,
! RHOEND, X, NF, F, XHIST, FHIST, CHIST, CONHIST, CSTRV and INFO are identical to the corresponding
! arguments in subroutine COBYLA.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, QUART, TENTH, EPS, REALMAX, &
    & DEBUGGING, MIN_MAXFILT
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: infos_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: message_mod, only : retmsg, rhomsg, fmsg
use, non_intrinsic :: pintrf_mod, only : OBJCON
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: selectx_mod, only : savefilt, selectx, isbetter

! Solver-specific modules
use, non_intrinsic :: geometry_cobyla_mod, only : assess_geo, setdrop_geo, setdrop_tr, geostep
use, non_intrinsic :: initialize_cobyla_mod, only : initxfc, initfilt
use, non_intrinsic :: trustregion_cobyla_mod, only : trstlp, trrad
use, non_intrinsic :: update_cobyla_mod, only : updatexfc, updatepole

implicit none

! Inputs
procedure(OBJCON) :: calcfc ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfilt
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: amat(:, :) ! AMAT(N, M_LCON)
real(RP), intent(in) :: bvec(:) ! BVEC(M_LCON)
real(RP), intent(in) :: ctol
real(RP), intent(in) :: cweight
real(RP), intent(in) :: eta1
real(RP), intent(in) :: eta2
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: gamma1
real(RP), intent(in) :: gamma2
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
! On entry, [X, F, CONSTR] = [X0, F(X0), CONSTR(X0)]
real(RP), intent(inout) :: constr(:)    ! CONSTR(M)
real(RP), intent(inout) :: f
real(RP), intent(inout) :: x(:)     ! X(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)   ! CHIST(MAXCHIST)
real(RP), intent(out) :: conhist(:, :)  ! CONHIST(M, MAXCONHIST)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: fhist(:)   ! FHIST(MAXFHIST)
real(RP), intent(out) :: xhist(:, :)    ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'COBYLA'
character(len=*), parameter :: srname = 'COBYLB'
integer(IK) :: jdrop_geo
integer(IK) :: jdrop_tr
integer(IK) :: kopt
integer(IK) :: m
integer(IK) :: m_lcon
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxtr
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: nfilt
integer(IK) :: nhist
integer(IK) :: subinfo
integer(IK) :: tr
logical :: bad_trstep
logical :: adequate_geo
logical :: evaluated(size(x) + 1)
logical :: improve_geo
logical :: reduce_rho
logical :: shortd
logical :: trfail
logical :: ximproved
real(RP) :: A(size(x), size(constr)) ! A contains the approximate gradient for the constraints
real(RP) :: actrem
real(RP) :: cfilt(min(max(maxfilt, 1_IK), maxfun))
real(RP) :: confilt(size(constr), size(cfilt))
real(RP) :: conmat(size(constr), size(x) + 1)
real(RP) :: cpen  ! Penalty parameter for constraint in merit function (PARMU in Powell's code)
real(RP) :: cval(size(x) + 1)
real(RP) :: d(size(x))
real(RP) :: delta
real(RP) :: dnorm
real(RP) :: ffilt(size(cfilt))
real(RP) :: fval(size(x) + 1)
real(RP) :: g(size(x))
real(RP) :: gamma3
real(RP) :: prerec  ! Predicted reduction in constraint violation
real(RP) :: preref  ! Predicted reduction in objective Function
real(RP) :: prerem  ! Predicted reduction in merit function
real(RP) :: ratio  ! Reduction ratio: ACTREM/PREREM
real(RP) :: rho
real(RP) :: sim(size(x), size(x) + 1)
real(RP) :: simi(size(x), size(x))
real(RP) :: xfilt(size(x), size(cfilt))
! CPENMIN is the minimum of the penalty parameter CPEN for the L-infinity constraint violation in
! the merit function. Note that CPENMIN = 0 in Powell's implementation, which allows CPEN to be 0.
! Here, we take CPENMIN > 0 so that CPEN is always positive. This avoids the situation where PREREM
! becomes 0 when PREREF = 0 = CPEN. It brings two advantages as follows.
! 1. If the trust-region subproblem solver works correctly and the trust-region center is not
! optimal for the subproblem, then PREREM > 0 is guaranteed. This is because, in theory, PREREC >= 0
! and MAX(PREREC, PREREF) > 0 , and the definition of CPEN in GETCPEN ensures that PREREM > 0.
! 2. There is no need to revise ACTREM and PREREM when CPEN = 0 and F = FVAL(N+1) as in lines
! 312--314 of Powell's cobylb.f code. Powell's code revises ACTREM to CVAL(N + 1) - CSTRV and PREREM
! to PREREC in this case, which is crucial for feasibility problems.
real(RP), parameter :: cpenmin = EPS
! FACTOR_ALPHA, FACTOR_BETA, FACTOR_GAMMA, and FACTOR_DELTA are four factors that COBYLB uses
! when managing the simplex. Note the following.
! 1. FACTOR_ALPHA < FACTOR_GAMMA < 1 < FACTOR_DELTA <= FACTOR_BETA.
! 2. FACTOR_DELTA has nothing to do with DELTA, which is the trust-region radius.
! 3. FACTOR_GAMMA has nothing to do with GAMMA1 and GAMMA2, which are the contracting/expanding
! factors for updating the trust-region radius DELTA.
real(RP), parameter :: factor_alpha = QUART  ! The factor alpha in the COBYLA paper
real(RP), parameter :: factor_beta = 2.1_RP  ! The factor beta in the COBYLA paper
real(RP), parameter :: factor_delta = 1.1_RP  ! The factor delta in the COBYLA paper
real(RP), parameter :: factor_gamma = HALF  ! The factor gamma in the COBYLA paper

! Sizes
m_lcon = int(size(bvec), kind(m_lcon))
m = int(size(constr), kind(m))
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxconhist = int(size(conhist, 2), kind(maxconhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = int(max(maxxhist, maxfhist, maxconhist, maxchist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(maxfun >= n + 2, 'MAXFUN >= N + 2', srname)
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(ctol >= 0, 'CTOL >= 0', srname)
    call assert(cweight >= 0, 'CWEIGHT >= 0', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == size(bvec), 'SIZE(AMAT) == [N, SIZE(BVEC)]', srname)
    call assert(maxfilt >= min(MIN_MAXFILT, maxfun) .and. maxfilt <= maxfun, &
        & 'MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(conhist, 1) == m .and. maxconhist * (maxconhist - maxhist) == 0, &
        & 'SIZE(CONHIST, 1) == M, SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    call assert(factor_alpha > 0 .and. factor_alpha < factor_gamma .and. factor_gamma < 1, &
        & '0 < FACTOR_ALPHA < FACTOR_GAMMA < 1', srname)
    call assert(factor_beta >= factor_delta .and. factor_delta > 1, 'FACTOR_BETA >= FACTOR_DELTA > 1', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize SIM, SIMI, FVAL, CONMAT, and CVAL, together with the history, NF, and EVALUATED.
! After the initialization, SIM(:, N+1) holds the vertex of the initial simplex with the smallest
! function value (regardless of the constraint violation), and SIM(:, 1:N) holds the displacements
! from the other vertices to SIM(:, N+1). FVAL, CONMAT, and CVAL hold the function values,
! constraint values, and constraint violations on the vertices in the order corresponding to SIM.
call initxfc(calcfc_internal, iprint, maxfun, constr, ctol, f, ftarget, rhobeg, x, nf, chist, conhist, &
   & conmat, cval, fhist, fval, sim, simi, xhist, evaluated, subinfo)

! Initialize the filter, including XFILT, FFILT, CONFILT, CFILT, and NFILT.
! N.B.: The filter is used only when selecting which iterate to return. It does not interfere with
! the iterations. COBYLA is NOT a filter method but a trust-region method based on an L-infinity
! merit function. Powell's implementation does not use a filter to select the iterate, possibly
! returning a suboptimal iterate.
call initfilt(conmat, ctol, cweight, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)

! Check whether to return due to abnormal cases that may occur during the initialization.
if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Return the best calculated values of the variables.
    ! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
    kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol)
    x = xfilt(:, kopt)
    f = ffilt(kopt)
    constr = confilt(:, kopt)
    cstrv = cfilt(kopt)
    ! Arrange CHIST, CONHIST, FHIST, and XHIST so that they are in the chronological order.
    call rangehist(nf, xhist, fhist, chist, conhist)
    ! Print a return message according to IPRINT.
    call retmsg(solver, info, iprint, nf, f, x, cstrv, constr)
    ! Postconditions
    if (DEBUGGING) then
        call assert(nf <= maxfun, 'NF <= MAXFUN', srname)
        call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
        call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
        call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
        call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
        ! The last calculated X can be Inf (finite + finite can be Inf numerically).
        call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
        call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
            & 'FHIST does not contain NaN/+Inf', srname)
        call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxconhist, &
            & 'SIZE(CONHIST) == [M, MAXCONHIST]', srname)
        call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
            & is_posinf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/+Inf', srname)
        call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
        call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) &
            & .or. is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain negative values or NaN/+Inf', srname)
        nhist = minval([nf, maxfhist, maxchist])
        call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f, cstrv, ctol)), &
            & 'No point in the history is better than X', srname)
    end if
    return
end if

! Set some more initial values.
! We must initialize ACTREM and PREREM. Otherwise, when SHORTD = TRUE, compilers may raise a
! run-time error that they are undefined. But their values will not be used: when SHORTD = FALSE,
! they will be overwritten; when SHORTD = TRUE, the values are used only in BAD_TRSTEP, which is
! TRUE regardless of ACTREM or PREREM. Similar for PREREC, PREREF, PREREM, RATIO, and JDROP_TR.
! No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
! Our initialization of CPEN differs from Powell's in two ways. First, we use the ratio defined in
! (13) of Powell's COBYLA paper to initialize CPEN. Second, we impose CPEN >= CPENMIN > 0. Powell's
! code simply initializes CPEN to 0.
rho = rhobeg
delta = rhobeg
cpen = max(cpenmin, min(1.0E3_RP, fcratio(conmat, fval)))  ! Powell's code: CPEN = ZERO
prerec = -REALMAX
preref = -REALMAX
prerem = -REALMAX
actrem = -REALMAX
shortd = .false.
trfail = .false.
ratio = -ONE
jdrop_tr = 0
jdrop_geo = 0

! If DELTA <= GAMMA3*RHO after an update, we set DELTA to RHO. GAMMA3 must be less than GAMMA2. The
! reason is as follows. Imagine a very successful step with DENORM = the un-updated DELTA = RHO.
! Then TRRAD will update DELTA to GAMMA2*RHO. If GAMMA3 >= GAMMA2, then DELTA will be reset to RHO,
! which is not reasonable as D is very successful. See paragraph two of Sec. 5.2.5 in
! T. M. Ragonneau's thesis: "Model-Based Derivative-Free Optimization Methods and Software".
! According to test on 20230613, for COBYLA, this Powellful updating scheme of DELTA works slightly
! better than setting directly DELTA = MAX(NEW_DELTA, RHO).
gamma3 = max(ONE, min(0.75_RP * gamma2, 1.5_RP))

! MAXTR is the maximal number of trust-region iterations. Each trust-region iteration takes 1 or 2
! function evaluations unless the trust-region step is short or the trust-region subproblem solver
! fails but the geometry step is not invoked. Thus the following MAXTR is unlikely to be reached.
maxtr = max(maxfun, 2_IK * maxfun)  ! MAX: precaution against overflow, which will make 2*MAXFUN < 0.
info = MAXTR_REACHED

! Begin the iterative procedure.
! After solving a trust-region subproblem, we use three boolean variables to control the workflow.
! SHORTD - Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust-region iteration? If yes, a geometry step
! will be taken, corresponding to the "Branch (Delta)" in the COBYLA paper.
! REDUCE_RHO - Will we reduce rho after the trust-region iteration?
! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Increase the penalty parameter CPEN, if needed, so that PREREM = PREREF + CPEN * PREREC > 0.
    ! This is the first (out of two) update of CPEN, where CPEN increases or remains the same.
    ! N.B.: CPEN and the merit function PHI = FVAL + CPEN*CVAL are used at three places only.
    ! 1. In FINDPOLE/UPDATEPOLE, deciding the optimal vertex of the current simplex.
    ! 2. After the trust-region trial step, calculating the reduction radio.
    ! 3. In GEOSTEP, deciding the direction of the geometry step.
    ! They do not appear explicitly in the trust-region subproblem, though the trust-region center
    ! (i.e., the current optimal vertex) is defined by them.
    cpen = getcpen(amat, bvec, conmat, cpen, cval, delta, fval, rho, sim, simi)

    ! Switch the best vertex of the current simplex to SIM(:, N + 1).
    call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
    ! Check whether to exit due to damaging rounding in UPDATEPOLE.
    if (subinfo == DAMAGING_ROUNDING) then
        info = subinfo
        exit  ! Better action to take? Geometry step, or simply continue?
    end if

    ! Does the interpolation set have acceptable geometry? It affects IMPROVE_GEO and REDUCE_RHO.
    adequate_geo = assess_geo(delta, factor_alpha, factor_beta, sim, simi)

    ! Calculate the linear approximations to the objective and constraint functions.
    ! N.B.: TRSTLP accesses A mostly by columns, so it is more reasonable to save A instead of A^T.
    ! Zaikun 2023108: According to a test on 2023108, calculating G and A(:, M_LCON+1:M) by solving
    ! the linear systems SIM^T*G = FVAL(1:N)-FVAL(N+1) and SIM^T*A = CONMAT(:, 1:N)-CONMAT(:, N+1)
    ! does not seem to improve or worsen the performance of COBYLA in terms of the number of function
    ! evaluations. The system was solved by SOLVE in LINALG_MOD based on a QR factorization of SIM
    ! (not necessarily a good algorithm). No preconditioning was used.
    g = matprod(fval(1:n) - fval(n + 1), simi)
    A(:, 1:m_lcon) = amat
    A(:, m_lcon + 1:m) = transpose(matprod(conmat(m_lcon + 1:m, 1:n) - spread(conmat(m_lcon + 1:m, n + 1), dim=2, ncopies=n), simi))
    !!MATLAB: A(:, m_lcon+1:m) = simi'*(conmat(m_lcon+1:m, 1:n) - conmat(m_lcon+1:m, n+1))' % Implicit expansion for subtraction

    ! Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
    d = trstlp(A, -conmat(:, n + 1), delta, g)
    dnorm = min(delta, norm(d))

    ! Is the trust-region trial step short? Note that we compare DNORM with RHO, not DELTA.
    ! Powell's code essentially defines SHORTD by SHORTD = (DNORM < HALF * RHO). In our tests,
    ! TENTH seems to work better than HALF or QUART, especially for linearly constrained problems.
    ! Note that LINCOA has a slightly more sophisticated way of defining SHORTD, taking into account
    ! whether D causes a change to the active set. Should we try the same here?
    shortd = (dnorm < TENTH * rho)

    ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
    ! We have the following in precise arithmetic. They may fail to hold due to rounding errors.
    ! 1. PREREC is the reduction of the L-infinity violation of the linearized constraints achieved
    ! by D. It is nonnegative in theory; it is 0 iff CONMAT(1:M, N+1) <= 0, namely the trust-region
    ! center satisfies the constraints.
    ! 2. PREREF may be negative or 0, but it should be positive when PREREC = 0 and SHORTD is FALSE.
    ! 3. Due to 2, in theory, MAXIMUM([PREREC, PREREF]) > 0 if SHORTD is FALSE.
    preref = -inprod(d, g)  ! Can be negative.
    prerec = cval(n + 1) - maxval([conmat(:, n + 1) + matprod(d, A), ZERO])

    ! Evaluate PREREM, which is the predicted reduction in the merit function.
    ! In theory, PREREM >= 0 and it is 0 iff CPEN = 0 = PREREF. This may not be true numerically.
    prerem = preref + cpen * prerec
    trfail = (.not. prerem > 1.0E-5_RP * min(cpen, ONE) * rho**2)  ! PREREM is tiny/negative or NaN.

    if (shortd .or. trfail) then
        ! Reduce DELTA if D is short or D fails to render PREREM > 0. The latter can happen due to
        ! rounding errors. This seems important for performance.
        delta = TENTH * delta
        if (delta <= gamma3 * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if
    else
        x = sim(:, n + 1) + d

        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evaluate(calcfc_internal, x, f, constr)
        cstrv = maxval([ZERO, constr])
        nf = nf + 1_IK

        ! Print a message about the function/constraint evaluation according to IPRINT.
        call fmsg(solver, 'Trust region', iprint, nf, delta, f, x, cstrv, constr)
        ! Save X, F, CONSTR, CSTRV into the history.
        call savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
        ! Save X, F, CONSTR, CSTRV into the filter.
        call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, constr, confilt)

        ! Evaluate ACTREM, which is the actual reduction in the merit function.
        actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv)

        ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
        ratio = redrat(actrem, prerem, eta1)

        ! Update DELTA. After this, DELTA < DNORM may hold.
        ! N.B.: 1. Powell's code uses RHO as the trust-region radius and updates it as follows.
        ! Reduce RHO to GAMMA1*RHO if ADEQUATE_GEO is TRUE and either SHORTD is TRUE or RATIO < ETA1,
        ! and then revise RHO to RHOEND if its new value is not more than GAMMA3*RHOEND; RHO remains
        ! unchanged in all other cases; in particular, RHO is never increased.
        ! 2. Our implementation uses DELTA as the trust-region radius, while using RHO as a lower
        ! bound for DELTA. DELTA is updated in a way that is typical for trust-region methods, and
        ! it is revised to RHO if its new value is not more than GAMMA3*RHO. RHO reflects the current
        ! resolution of the algorithm; its update is essentially the same as the update of RHO in
        ! Powell's code (see the definition of REDUCE_RHO below). Our implementation aligns with
        ! UOBYQA/NEWUOA/BOBYQA/LINCOA and improves the performance of COBYLA.
        ! 3. The same as Powell's code, we do not reduce RHO unless ADEQUATE_GEO is TRUE. This is
        ! also how Powell updated RHO in UOBYQA/NEWUOA/BOBYQA/LINCOA. What about we also use
        ! ADEQUATE_GEO == TRUE as a prerequisite for reducing DELTA? The argument would be that the
        ! bad (small) value of RATIO may be because of a bad geometry (and hence a bad model) rather
        ! than an improperly large DELTA, and it might be good to try improving the geometry first
        ! without reducing DELTA. However, according to a test on 20230206, it does not improve the
        ! performance if we skip the update of DELTA when ADEQUATE_GEO is FALSE and RATIO < 0.1.
        ! Therefore, we choose to update DELTA without checking ADEQUATE_GEO.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= gamma3 * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if

        ! Is the newly generated X better than current best point?
        ximproved = (actrem > 0)  ! If ACTREM is NaN, then XIMPROVED should & will be FALSE.

        ! Set JDROP_TR to the index of the vertex to be replaced with X. JDROP_TR = 0 means there
        ! is no good point to replace, and X will not be included into the simplex; in this case,
        ! the geometry of the simplex likely needs improvement, which will be handled below.
        jdrop_tr = setdrop_tr(ximproved, d, delta, rho, sim, simi)

        ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_TR) is replaced with D.
        ! UPDATEXFC does nothing if JDROP_TR == 0, as the algorithm decides to discard X.
        call updatexfc(jdrop_tr, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi, subinfo)
        ! Check whether to exit due to damaging rounding in UPDATEXFC.
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit  ! Better action to take? Geometry step, or a RESCUE as in BOBYQA?
        end if

        ! Check whether to exit due to MAXFUN, FTARGET, etc.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end if  ! End of IF (SHORTD .OR. TRFAIL). The normal trust-region calculation ends.


    !----------------------------------------------------------------------------------------------!
    ! Before the next trust-region iteration, we possibly improve the geometry of simplex or
    ! reduces RHO according to IMPROVE_GEO and REDUCE_RHO. Now we decide these indicators.
    ! N.B.: We must ensure that the algorithm does not set IMPROVE_GEO = TRUE at infinitely many
    ! consecutive iterations without moving SIM(:, N+1) or reducing RHO. Otherwise, the algorithm
    ! will get stuck in repetitive invocations of GEOSTEP. This is ensured by the following facts.
    ! 1. If an iteration sets IMPROVE_GEO = TRUE, it must also reduce DELTA or set DELTA to RHO.
    ! 2. If SIM(:, N+1) and RHO remains unchanged, then ADEQUATE_GEO will become TRUE after at
    ! most N invocations of GEOSTEP.

    ! BAD_TRSTEP: Is the last trust-region step bad?
    bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. jdrop_tr == 0)
    ! IMPROVE_GEO: Should we take a geometry step to improve the geometry of the interpolation set?
    improve_geo = (bad_trstep .and. .not. adequate_geo)
    ! REDUCE_RHO: Should we enhance the resolution by reducing RHO?
    reduce_rho = (bad_trstep .and. adequate_geo .and. max(delta, dnorm) <= rho)

    ! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    ! !call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO or REDUCE_RHO are not both TRUE', srname)

    ! If SHORTD or TRFAIL is TRUE, then either IMPROVE_GEO or REDUCE_RHO is TRUE unless ADEQUATE_GEO
    ! is TRUE and MAX(DELTA, DNORM) > RHO.
    ! !call assert((.not. (shortd .or. trfail)) .or. (improve_geo .or. reduce_rho .or. &
    ! !    & (adequate_geo .and. max(delta, dnorm) > rho)), 'If SHORTD or TRFAIL is TRUE, then &
    ! !    & either IMPROVE_GEO or REDUCE_RHO is TRUE unless ADEQUATE_GEO is TRUE and MAX(DELTA, DNORM) > RHO', srname)
    !----------------------------------------------------------------------------------------------!

    ! Comments on BAD_TRSTEP:
    ! 1. Powell's definition of BAD_TRSTEP is as follows. The one used above seems to work better,
    ! especially for linearly constrained problems due to the factor TENTH (= ETA1).
    ! !bad_trstep = (shortd .or. actrem <= 0 .or. actrem < TENTH * prerem .or. jdrop_tr == 0)
    ! Besides, Powell did not check PREREM > 0 in BAD_TRSTEP, which is reasonable to do but has
    ! little impact upon the performance.
    ! 2. NEWUOA/BOBYQA/LINCOA would define BAD_TRSTEP, IMPROVE_GEO, and REDUCE_RHO as follows. Two
    ! different thresholds are used in BAD_TRSTEP. It outperforms Powell's version.
    ! !bad_trstep = (shortd .or. trfail .or. ratio <= eta1 .or. jdrop_tr == 0)
    ! !improve_geo = bad_trstep .and. .not. adequate_geo
    ! !bad_trstep = (shortd .or. trfail .or. ratio <= 0 .or. jdrop_tr == 0)
    ! !reduce_rho = bad_trstep .and. adequate_geo .and. max(delta, dnorm) <= rho
    ! 3. Theoretically, JDROP_TR > 0 when ACTREM > 0 (guaranteed by RATIO > 0). However, in Powell's
    ! implementation, JDROP_TR may be 0 even RATIO > 0 due to NaN. The modernized code has rectified
    ! this in the function SETDROP_TR. After this rectification, we can indeed simplify the
    ! definition of BAD_TRSTEP by removing the condition JDROP_TR == 0. We retain it for robustness.

    ! Comments on REDUCE_RHO:
    ! When SHORTD is TRUE, UOBYQA/NEWUOA/BOBYQA/LINCOA all set REDUCE_RHO to TRUE if the recent
    ! models are sufficiently accurate according to certain criteria. See the paragraph around (37)
    ! in the UOBYQA paper and the discussions about Box 14 in the NEWUOA paper. This strategy is
    ! crucial for the performance of the solvers. However, as of 20221111, we have not managed to
    ! make it work in COBYLA. As in NEWUOA, we recorded the errors of the recent models, and set
    ! REDUCE_RHO to true if they are small (e.g., ALL(ABS(MODERR_REC) <= 0.1 * MAXVAL(ABS(A))*RHO) or
    ! ALL(ABS(MODERR_REC) <= RHO**2)) when SHORTD is TRUE. It made little impact on the performance.


    ! Since COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously, the following
    ! two blocks are exchangeable: IF (IMPROVE_GEO) ... END IF and IF (REDUCE_RHO) ... END IF.

    ! Improve the geometry of the simplex by removing a point and adding a new one.
    ! If the current interpolation set has acceptable geometry, then we skip the geometry step.
    ! The code has a small difference from Powell's original code here: If the current geometry
    ! is acceptable, then we will continue with a new trust-region iteration; however, at the
    ! beginning of the iteration, CPEN may be updated, which may alter the pole point SIM(:, N+1)
    ! by UPDATEPOLE; the quality of the interpolation point depends on SIM(:, N + 1), meaning
    ! that the same interpolation set may have good or bad geometry with respect to different
    ! "poles"; if the geometry turns out bad with the new pole, the original COBYLA code will
    ! take a geometry step, but our code here will NOT do it but continue to take a trust-region
    ! step. The argument is this: even if the geometry step is not skipped in the first place, the
    ! geometry may turn out bad again after the pole is altered due to an update to CPEN; should
    ! we take another geometry step in that case? If no, why should we do it here? Indeed, this
    ! distinction makes no practical difference for CUTEst problems with at most 100 variables
    ! and 5000 constraints, while the algorithm framework is simplified.
    if (improve_geo .and. .not. assess_geo(delta, factor_alpha, factor_beta, sim, simi)) then
        ! Before the geometry step, UPDATEPOLE has been called either implicitly by UPDATEXFC or
        ! explicitly after CPEN is updated, so that SIM(:, N + 1) is the optimal vertex.

        ! Decide a vertex to drop from the simplex. It will be replaced with SIM(:, N + 1) + D to
        ! improve acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
        ! N.B.: COBYLA never sets JDROP_GEO = N + 1.
        jdrop_geo = setdrop_geo(delta, factor_alpha, factor_beta, sim, simi)

        ! The following JDROP_GEO comes from UOBYQA/NEWUOA/BOBYQA/LINCOA. It performs poorly!
        !!jdrop_geo = maxloc(sum(sim(:, 1:n)**2, dim=1), dim=1)

        ! JDROP_GEO is between 1 and N unless SIM and SIMI contain NaN, which should not happen
        ! at this point unless there is a bug. Nevertheless, for robustness, we include the
        ! following instruction to exit when JDROP_GEO == 0 (if JDROP_GEO does become 0, then
        ! memory error will occur if we continue, as JDROP_GEO will be used as an index of arrays.)
        if (jdrop_geo == 0) then
            info = DAMAGING_ROUNDING
            exit
        end if

        ! Calculate the geometry step D.
        ! In NEWUOA, GEOSTEP takes DELBAR = MAX(MIN(TENTH * SQRT(MAXVAL(DISTSQ)), HALF * DELTA), RHO)
        ! rather than DELTA. This should not be done here, because D should improve the geometry of
        ! the simplex when SIM(:, JDROP) is replaced with D; the quality of the geometry is defined
        ! by DELTA instead of DELBAR as in (14) of the COBYLA paper. See GEOSTEP for more detail.
        d = geostep(jdrop_geo, amat, bvec, conmat, cpen, cval, delta, fval, factor_gamma, simi)
        x = sim(:, n + 1) + d

        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evaluate(calcfc_internal, x, f, constr)
        cstrv = maxval([ZERO, constr])
        nf = nf + 1_IK

        ! Print a message about the function/constraint evaluation according to IPRINT.
        call fmsg(solver, 'Geometry', iprint, nf, delta, f, x, cstrv, constr)
        ! Save X, F, CONSTR, CSTRV into the history.
        call savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
        ! Save X, F, CONSTR, CSTRV into the filter.
        call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, constr, confilt)

        ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_GEO) is replaced with D.
        call updatexfc(jdrop_geo, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi, subinfo)
        ! Check whether to exit due to damaging rounding in UPDATEXFC.
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit  ! Better action to take? Geometry step, or simply continue?
        end if

        ! Check whether to exit due to MAXFUN, FTARGET, etc.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end if  ! End of IF (IMPROVE_GEO). The procedure of improving geometry ends.

    ! The calculations with the current RHO are complete. Enhance the resolution of the algorithm
    ! by reducing RHO; update DELTA and CPEN at the same time.
    if (reduce_rho) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        delta = max(HALF * rho, redrho(rho, rhoend))
        rho = redrho(rho, rhoend)
        ! The second (out of two) update of CPEN, where CPEN decreases or remains the same.
        ! Powell's code: CPEN = MIN(CPEN, FCRATIO(FVAL, CONMAT)), which may set CPEN to 0.
        cpen = max(cpenmin, min(cpen, fcratio(conmat, fval)))
        ! Print a message about the reduction of RHO according to IPRINT.
        call rhomsg(solver, iprint, nf, delta, fval(n + 1), rho, sim(:, n + 1), cval(n + 1), conmat(:, n + 1), cpen)
        ! Switch the best vertex of the current simplex to SIM(:, N + 1).
        call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
        ! Check whether to exit due to damaging rounding in UPDATEPOLE.
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit  ! Better action to take? Geometry step, or simply continue?
        end if
    end if  ! End of IF (REDUCE_RHO). The procedure of reducing RHO ends.

end do  ! End of DO TR = 1, MAXTR. The iterative procedure ends.

! Return from the calculation, after trying the last trust-region step if it has not been tried yet.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    ! Zaikun 20230615: UPDATEXFC or UPDATEPOLE is not called since the last trust-region step. Hence
    ! SIM(:, N + 1) remains unchanged. Otherwise, SIM(:, N + 1) + D would not make sense.
    x = sim(:, n + 1) + d
    call evaluate(calcfc_internal, x, f, constr)
    cstrv = maxval([ZERO, constr])
    nf = nf + 1_IK
    ! Print a message about the function evaluation according to IPRINT.
    ! Zaikun 20230512: DELTA has been updated. RHO is only indicative here. TO BE IMPROVED.
    ! Print a message about the function/constraint evaluation according to IPRINT.
    call fmsg(solver, 'Trust region', iprint, nf, rho, f, x, cstrv, constr)
    ! Save X, F, CONSTR, CSTRV into the history.
    call savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
    ! Save X, F, CONSTR, CSTRV into the filter.
    call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, constr, confilt)
end if

! Return the best calculated values of the variables.
! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), max(cpen, cweight), ctol)
x = xfilt(:, kopt)
f = ffilt(kopt)
constr = confilt(:, kopt)
cstrv = cfilt(kopt)

! Arrange CHIST, CONHIST, FHIST, and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist, chist, conhist)

! Print a return message according to IPRINT.
call retmsg(solver, info, iprint, nf, f, x, cstrv, constr)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nf <= maxfun, 'NF <= MAXFUN', srname)
    call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
        & 'FHIST does not contain NaN/+Inf', srname)
    call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxconhist, &
        & 'SIZE(CONHIST) == [M, MAXCONHIST]', srname)
    call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
        & is_posinf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/+Inf', srname)
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) &
        & .or. is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain negative values or NaN/+Inf', srname)
    nhist = minval([nf, maxfhist, maxchist])
    call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f, cstrv, ctol)), &
        & 'No point in the history is better than X', srname)
end if


contains


subroutine calcfc_internal(x_internal, f_internal, constr_internal)
!--------------------------------------------------------------------------------------------------!
! This internal subroutine evaluates the objective function and ALL the constraints.
! In MATLAB/Python/R/Julia, this can be implemented as a lambda function / anonymous function.
!--------------------------------------------------------------------------------------------------!b
implicit none
! Inputs
real(RP), intent(in) :: x_internal(:)
! Outputs
real(RP), intent(out) :: f_internal
real(RP), intent(out) :: constr_internal(:)
constr_internal(1:m_lcon) = matprod(x_internal, amat) - bvec
call calcfc(x_internal, f_internal, constr_internal(m_lcon + 1:m))
end subroutine calcfc_internal

end subroutine cobylb


function getcpen(amat, bvec, conmat_in, cpen_in, cval_in, delta, fval_in, rho, sim_in, simi_in) result(cpen)
!--------------------------------------------------------------------------------------------------!
! This function gets the penalty parameter CPEN so that PREREM = PREREF + CPEN * PREREC > 0.
! See the discussions around equation (9) of the COBYLA paper.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_posinf, is_nan
use, non_intrinsic :: infos_mod, only : INFO_DFT, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : matprod, inprod, isinv

! Solver-specific modules
use, non_intrinsic :: trustregion_cobyla_mod, only : trstlp
use, non_intrinsic :: update_cobyla_mod, only : findpole, updatepole

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)
real(RP), intent(in) :: bvec(:)
real(RP), intent(in) :: conmat_in(:, :)
real(RP), intent(in) :: cpen_in
real(RP), intent(in) :: cval_in(:)
real(RP), intent(in) :: delta
real(RP), intent(in) :: fval_in(:)
real(RP), intent(in) :: rho
real(RP), intent(in) :: sim_in(:, :)
real(RP), intent(in) :: simi_in(:, :)

! Outputs
real(RP) :: cpen

! Local variables
character(len=*), parameter :: srname = 'getcpen'
integer(IK) :: info
integer(IK) :: iter
integer(IK) :: m
integer(IK) :: m_lcon
integer(IK) :: n
real(RP) :: A(size(sim_in, 1), size(conmat_in, 1))
real(RP) :: conmat(size(conmat_in, 1), size(conmat_in, 2))
real(RP) :: cval(size(cval_in))
real(RP) :: d(size(sim_in, 1))
real(RP) :: fval(size(fval_in))
real(RP) :: g(size(sim_in, 1))
real(RP) :: prerec
real(RP) :: preref
real(RP) :: sim(size(sim_in, 1), size(sim_in, 2))
real(RP) :: simi(size(simi_in, 1), size(simi_in, 2))
real(RP), parameter :: itol = ONE

! Sizes
m_lcon = int(size(bvec), kind(m_lcon))
m = int(size(conmat, 1), kind(m))
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == size(bvec), 'SIZE(AMAT) == [N, SIZE(BVEC)]', srname)
    call assert(cpen_in > 0, 'CPEN > 0', srname)
    call assert(size(conmat_in, 1) == m .and. size(conmat_in, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat_in) .or. is_posinf(conmat_in)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval_in) == n + 1 .and. .not. any(cval_in < 0 .or. is_nan(cval_in) .or. is_posinf(cval_in)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fval_in) == n + 1 .and. .not. any(is_nan(fval_in) .or. is_posinf(fval_in)), &
        & 'SIZE(FVAL) == N+1 and FVAL does not contain NaN/+Inf', srname)
    call assert(size(sim_in, 1) == n .and. size(sim_in, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim_in)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim_in(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi_in, 1) == n .and. size(simi_in, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi_in)), 'SIMI is finite', srname)
    call assert(isinv(sim_in(:, 1:n), simi_in, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
    call assert(delta >= rho .and. rho > 0, 'DELTA >= RHO > 0', srname)
end if

!====================!
! Calculation starts !
!====================!

! Copy the inputs.
conmat = conmat_in
cpen = cpen_in
cval = cval_in
fval = fval_in
sim = sim_in
simi = simi_in

! Initialize INFO, PREREF, and PREREC, which are needed in the postconditions.
info = INFO_DFT
preref = ZERO
prerec = ZERO

! Increase CPEN if necessary to ensure PREREM > 0. Branch back for the next loop if this change
! alters the optimal vertex of the current simplex. Note the following.
! 1. In each loop, CPEN is changed only if PREREC > 0 > PREREF, in which case PREREM is guaranteed
! positive after the update. Note that PREREC >= 0 and MAX(PREREC, PREREF) > 0 in theory. If this
! holds numerically as well, then CPEN is not changed only if PREREC = 0 or PREREF >= 0, in which
! case PREREM is currently positive, explaining why CPEN needs no update.
! 2. Even without an upper bound for the loop counter, the loop can occur at most N+1 times. This is
! because the update of CPEN does not decrease CPEN, and hence it can make vertex J (J <= N) become
! the new optimal vertex only if CVAL(J) is less than CVAL(N+1), which can happen at most N times.
! See the paragraph below (9) in the COBYLA paper. After the "correct" optimal vertex is found,
! one more loop is needed to calculate CPEN, and hence the loop can occur at most N+1 times.
do iter = 1, n + 1_IK
    ! Switch the best vertex of the current simplex to SIM(:, N + 1).
    call updatepole(cpen, conmat, cval, fval, sim, simi, info)
    ! Check whether to exit due to damaging rounding in UPDATEPOLE.
    if (info == DAMAGING_ROUNDING) then
        exit
    end if

    ! Calculate the linear approximations to the objective and constraint functions.
    g = matprod(fval(1:n) - fval(n + 1), simi)
    A(:, 1:m_lcon) = amat
    A(:, m_lcon + 1:m) = transpose(matprod(conmat(m_lcon + 1:m, 1:n) - spread(conmat(m_lcon + 1:m, n + 1), dim=2, ncopies=n), simi))
    !!MATLAB: A(:, m_lcon+1:m) = simi'*(conmat(m_lcon+1:m, 1:n) - conmat(m_lcon+1:m, n+1))' % Implicit expansion for subtraction

    ! Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
    d = trstlp(A, -conmat(:, n + 1), delta, g)

    ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
    preref = -inprod(d, g)  ! Can be negative.
    prerec = cval(n + 1) - maxval([conmat(:, n + 1) + matprod(d, A), ZERO])

    if (.not. (prerec > 0 .and. preref < 0)) then  ! PREREC <= 0 or PREREF >= 0 or either is NaN.
        exit
    end if

    ! Powell's code defines BARMU = -PREREF / PREREC, and CPEN is increased to 2*BARMU if and
    ! only if it is currently less than 1.5*BARMU, a very "Powellful" scheme. In our implementation,
    ! however, we set CPEN directly to the maximum between its current value and 2*BARMU while
    ! handling possible overflow. This simplifies the scheme without worsening the performance.
    cpen = max(cpen, min(-TWO * (preref / prerec), REALMAX))

    if (findpole(cpen, cval, fval) == n + 1) then
        exit
    end if
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(cpen >= cpen_in .and. cpen > 0, 'CPEN >= CPEN_IN and CPEN > 0', srname)
    call assert(preref + cpen * prerec > 0 .or. info == DAMAGING_ROUNDING .or. &
        & .not. (prerec >= 0 .and. max(prerec, preref) > 0) .or. .not. is_finite(preref), &
        & 'PREREF + CPEN*PREREC > 0 unless D is short or the rounding is damaging', srname)
end if
end function getcpen


function fcratio(conmat, fval) result(r)
!--------------------------------------------------------------------------------------------------!
! This function calculates the ratio between the "typical change" of F and that of CONSTR.
! See equations (12)--(13) in Section 3 of the COBYLA paper for the definition of the ratio.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
implicit none

! Inputs
real(RP), intent(in) :: conmat(:, :)    ! CONMAT(M, N+1)
real(RP), intent(in) :: fval(:)     ! FVAL(N+1)

! Outputs
real(RP) :: r

! Local variables
real(RP) :: cmax(size(conmat, 1))
real(RP) :: cmin(size(conmat, 1))
real(RP) :: denom
real(RP) :: fmax
real(RP) :: fmin
character(len=*), parameter :: srname = 'FCRATIO'

! Preconditions
if (DEBUGGING) then
    call assert(size(fval) >= 1, 'SIZE(FVAL) >= 1', srname)
    call assert(size(conmat, 2) == size(fval), 'SIZE(CONMAT, 2) == SIZE(FVAL)', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(fval) .or. is_posinf(fval)), 'FVAL does not contain NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! N.B.: In the original version of COBYLA, Powell proposed the ratio for constraints in the form of
! CONSTR(X) >= 0, but the constraints we consider here are CONSTR(X) <= 0. Hence we need to change
! the sign of the constraints before defining CMIN and CMAX.
cmin = minval(-conmat, dim=2)
cmax = maxval(-conmat, dim=2)
fmin = minval(fval)
fmax = maxval(fval)
r = ZERO
if (any(cmin < HALF * cmax) .and. fmin < fmax) then
    denom = minval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
    ! Powell mentioned the following alternative in Section 4 of his COBYLA paper. According to a
    ! test on 20230610, it does not make much difference to the performance.
    ! !denom = maxval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
    r = (fmax - fmin) / denom
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(r >= 0, 'R >= 0', srname)
end if
end function fcratio


end module cobylb_mod
