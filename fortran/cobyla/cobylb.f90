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
! Last Modified: Friday, March 10, 2023 AM11:44:01
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobylb


contains


subroutine cobylb(calcfc, iprint, maxfilt, maxfun, ctol, cweight, eta1, eta2, ftarget, &
    & gamma1, gamma2, rhobeg, rhoend, constr, f, x, nf, chist, conhist, cstrv, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the actual calculations of COBYLA.
!
! IPRINT, MAXFILT, MAXFUN, MAXHIST, CTOL, CWEIGHT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG,
! RHOEND, X, NF, F, XHIST, FHIST, CHIST, CONHIST, CSTRV and INFO are identical to the corresponding
! arguments in subroutine COBYLA.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH, REALMAX, &
    & DEBUGGING, MIN_MAXFILT
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: infos_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: output_mod, only : retmsg, rhomsg, fmsg
use, non_intrinsic :: pintrf_mod, only : OBJCON
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: selectx_mod, only : savefilt, selectx, isbetter

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : assess_geo, setdrop_geo, setdrop_tr, geostep
use, non_intrinsic :: initialize_mod, only : initxfc, initfilt
use, non_intrinsic :: trustregion_mod, only : trstlp, trrad
use, non_intrinsic :: update_mod, only : updatexfc, updatepole, findpole

implicit none

! Inputs
procedure(OBJCON) :: calcfc ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfilt
integer(IK), intent(in) :: maxfun
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
logical :: ximproved
real(RP) :: A(size(x), size(constr) + 1)
! A(:, 1:M) contains the approximate gradient for the constraints, and A(:, M+1) is minus the
! approximate gradient for the objective function.
real(RP) :: actrem
real(RP) :: b(size(constr) + 1)
real(RP) :: barmu
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
real(RP) :: prerec  ! Predicted reduction in constraint violation
real(RP) :: preref  ! Predicted reduction in objective Function
real(RP) :: prerem  ! Predicted reduction in merit function
real(RP) :: rho
real(RP) :: sim(size(x), size(x) + 1)
real(RP) :: simi(size(x), size(x))
real(RP) :: xfilt(size(x), size(cfilt))
! N.B.: FACTOR_ALPHA, FACTOR_BETA, FACTOR_GAMMA, and FACTOR_DELTA are four factors that COBYLB uses
! when managing the simplex. Note the following.
! 1. FACTOR_ALPHA < FACTOR_GAMMA < 1 < FACTOR_DELTA <= FACTOR_BETA.
! 2. FACTOR_DELTA has nothing to do with DELTA, which is the trust-region radius.
! 3. FACTOR_GAMMA has nothing to do with GAMMA1 and GAMMA2, which are the contracting/expanding
!    factors for updating the trust-region radius DELTA.
real(RP), parameter :: factor_alpha = QUART  ! The factor alpha in the COBYLA paper
real(RP), parameter :: factor_beta = 2.1_RP  ! The factor beta in the COBYLA paper
real(RP), parameter :: factor_delta = 1.1_RP  ! The factor delta in the COBYLA paper
real(RP), parameter :: factor_gamma = HALF  ! The factor gamma in the COBYLA paper
real(RP) :: ratio  ! Reduction ratio: ACTREM/PREREM

! Sizes
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
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
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

! Initialize SIM, FVAL, CONMAT, and CVAL, together with the history.
! After the initialization, SIM(:, N+1) holds the vertex of the initial simplex with the smallest
! function value (regardless of the constraint violation), and SIM(:, 1:N) holds the displacements
! from the other vertices to SIM(:, N+1). FVAL, CONMAT, and CVAL hold the function values,
! constraint values, and constraint violations on the vertices in the order corresponding to SIM.
call initxfc(calcfc, iprint, maxfun, constr, ctol, f, ftarget, rhobeg, x, nf, chist, conhist, &
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
            & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
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
! run-time error that they are undefined. The values will not be used: when SHORTD = FALSE, they
! will be overwritten; when SHORTD = TRUE, the values are used only in BAD_TRSTEP, which is TRUE
! regardless of ACTREM or PREREM. Similar for PREREC, PREREF, PREREM, RATIO, and JDROP_TR.
rho = rhobeg
delta = rhobeg
cpen = ZERO
prerec = -REALMAX
preref = -REALMAX
prerem = -REALMAX
actrem = -REALMAX
ratio = -ONE
jdrop_tr = 0
jdrop_geo = 0

! MAXTR is the maximal number of trust-region iterations. Normally, each trust-region iteration
! takes 1 or 2 function evaluations except for the following cases:
! 1. the update of CPEN alters the optimal vertex;
! 2. the trust-region step is short or fails to reduce either the linearized objective or the
! linearized constraint violation but the geometry step is not invoked.
! The following MAXTR is unlikely to be reached.
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
    ! Before the trust-region step, UPDATEPOLE has been called either implicitly by INITXFC/UPDATEXFC
    ! or explicitly after CPEN is updated, so that SIM(:, N + 1) is the optimal vertex.

    ! Does the interpolation set have acceptable geometry? It affects IMPROVE_GEO and REDUCE_RHO.
    adequate_geo = assess_geo(delta, factor_alpha, factor_beta, sim, simi)

    ! Calculate the linear approximations to the objective and constraint functions, placing minus
    ! the objective function gradient after the constraint gradients in the array A.
    ! N.B.: TRSTLP accesses A mostly by columns, so it is more reasonable to save A instead of A^T.
    A(:, 1:m) = transpose(matprod(conmat(:, 1:n) - spread(conmat(:, n + 1), dim=2, ncopies=n), simi))
    !!MATLAB: A(:, 1:m) = simi'*(conmat(:, 1:n) - conmat(:, n+1))' % Implicit expansion for subtraction
    A(:, m + 1) = matprod(fval(n + 1) - fval(1:n), simi)

    ! Theoretically (but not numerically), the last entry of B does not affect the result of TRSTLP.
    ! We set it to -FVAL(N + 1) following Powell's code.
    b = [-conmat(:, n + 1), -fval(n + 1)]
    ! Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
    d = trstlp(A, b, delta)
    dnorm = min(delta, norm(d))

    ! Is the trust-region trial step short? N.B.: we compare DNORM with RHO, not DELTA.
    ! Powell's code especially defines SHORTD by SHORTD = (DNORM < HALF * RHO). In our tests,
    ! TENTH seems to work better than HALF or QUART, especially for linearly constrained problems.
    ! Note that LINCOA has a slightly more sophisticated way of defining SHORTD, taking into account
    ! whether D causes a change to the active set. Should we try the same here?
    shortd = (dnorm < TENTH * rho)

    ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
    ! We have the following in precise arithmetic. They may fail to hold due to rounding errors.
    ! 1. B(1:M) = -CONMAT(:, N + 1), and hence MAXVAL([B(1:M) - MATPROD(D, A(:, 1:M)), ZERO]) is the
    ! L-infinity violation of the linearized constraints corresponding to D. When D = 0, the
    ! violation is MAXVAL([B(1:M), ZERO]) = CVAL(N+1). PREREC is the reduction of this violation
    ! achieved by D, which is nonnegative in theory; PREREC = 0 iff B(1:M) <= 0, i.e., the
    ! trust-region center satisfies the linearized constraints.
    ! 2. PREREF may be negative or zero, but it is positive when PREREC = 0 and SHORTD is FALSE.
    ! 3. Due to 2, in theory, MAX(PREREC, PREREF) > 0 if SHORTD is FALSE.
    ! 4. In the code, MAX(PREREC, PREREF) is the counterpart of QRED in UOBYQA/NEWUOA/BOBYQA/LINCOA.
    prerec = cval(n + 1) - maxval([b(1:m) - matprod(d, A(:, 1:m)), ZERO])
    preref = inprod(d, A(:, m + 1))  ! Can be negative.

    if (shortd .or. .not. max(prerec, preref) > 0) then
        ! Reduce DELTA if D is short or D fails to render MAX(PREREC, PREREF) > 0, the latter can
        ! only happen due to rounding errors. This seems quite important for performance.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if
    else
        ! Increase CPEN if necessary to ensure PREREM > 0. Branch back if this change alters the
        ! optimal vertex. See the discussions around equation (9) of the COBYLA paper.
        ! This is the first (out of two) place where CPEN is updated. It can change CPEN only when
        ! PREREC > 0 > PREREF, in which case PREREM is guaranteed positive after the update.
        ! If PREREC = 0 or PREREF > 0, then PREREM is currently positive, so CPEN needs no update.
        ! However, as in Powell's implementation, if PREREC > 0 = PREREF = CPEN, then CPEN will
        ! remain zero, leaving PREREM = 0. If CPEN = 0 and PREREC > 0 > PREREF, then CPEN will
        ! become positive; if CPEN = 0, PREREC > 0, and PREREF > 0, then CPEN will remain zero.
        if (prerec > 0 .and. preref < 0) then
            ! In Powell's code, CPEN is increased to 2*BARMU if and only if it is currently less
            ! than 1.5*BARMU, a very "Powellful" scheme. In our implementation, however, we set CPEN
            ! directly to the maximum between its current value and 2*BARMU while handling possible
            ! overflow. This simplifies the scheme without worsening the performance of COBYLA.
            barmu = -preref / prerec  ! PREREF + BARMU * PREREC = 0
            cpen = max(cpen, min(TWO * barmu, REALMAX))  ! The 1st (out of 2) update of CPEN.
            if (findpole(cpen, cval, fval) <= n) then
                call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
                ! Check whether to exit due to damaging rounding in UPDATEPOLE.
                if (subinfo == DAMAGING_ROUNDING) then
                    info = subinfo
                    exit  ! Better action to take? Geometry step, or simply continue?
                end if
                cycle
                ! N.B.: The CYCLE can occur at most N times before a new function evaluation takes
                ! place. This is because the update of CPEN does not decrease CPEN, and hence it can
                ! make vertex J (J <= N) become the new optimal vertex only if CVAL(J) < CVAL(N+1),
                ! which can happen at most N times. See the paragraph below (9) in the COBYLA paper.
            end if
        end if

        x = sim(:, n + 1) + d
        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evaluate(calcfc, x, f, constr, cstrv)
        nf = nf + 1_IK

        ! Print a message about the function/constraint evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x, cstrv, constr)
        ! Save X, F, CONSTR, CSTRV into the history.
        call savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
        ! Save X, F, CONSTR, CSTRV into the filter.
        call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt, constr, confilt)

        ! Evaluate PREREM and ACTREM, which are the predicted and actual reductions in the merit
        ! function respectively.
        prerem = preref + cpen * prerec  ! Theoretically nonnegative; is 0 if CPEN = 0 = PREREF.
        actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv)
        if (cpen <= 0 .and. abs(f - fval(n + 1)) <= 0) then
            ! CPEN <= 0 indeed means CPEN == 0, while ABS(A - B) <= 0 indeed means A == B.
            prerem = prerec
            actrem = cval(n + 1) - cstrv
        end if

        ! In theory, PREREM >= 0, but this can fail due to rounding errors.
        !call assert(prerem >= 0, 'PREREM >= 0', 'COBYLA')

        !dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
        !moderrsav = [moderrsav(2:size(moderrsav)), maxval(abs([f - fval(n + 1) + preref, cstrv - cval(n + 1) + prerec]))]

        ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
        ratio = redrat(actrem, prerem, eta1)

        ! Update DELTA. After this, DELTA < DNORM may hold.
        ! N.B.: 1. Powell's code uses RHO as the trust-region radius and updates it as follows.
        ! Reduce RHO to GAMMA1*RHO if ADEQUATE_GEO is TRUE and either SHORTD is TRUE or RATIO < ETA1,
        ! and then revise RHO to RHOEND if its new value is not more than 1.5*RHOEND; RHO remains
        ! unchanged in all other cases; in particular, RHO is never increased.
        ! 2. Our implementation uses DELTA as the trust-region radius, while using RHO as a lower
        ! bound for DELTA. DELTA is updated in a way that is typical for trust-region methods, and
        ! it is revised to RHO if its new value is not more thant 1.5*RHO. RHO reflects the current
        ! resolution of the algorithm; its update is essentially the same as the update of RHO in
        ! Powell's code (see the definition of REDUCE_RHO below). Our implementation aligns with
        ! UOBYQA/NEWUOA/BOBYQA/LINCOA and improves the performance of COBYLA.
        ! 3. The same as Powell's code, we do not reduce RHO unless ADEQUATE_GEO is TRUE. This is
        ! also Powell updated RHO in UOBYQA/NEWUOA/BOBYQA/LINCOA. What about we also use
        ! ADEQUATE_GEO == TRUE as a prerequisite for reducing DELTA? The argument would be that the
        ! bad (small) value of RATIO may be because of a bad geometry (and hence a bad model) rather
        ! than an improperly large DELTA, and it might be good to try improving the geometry first
        ! without reducing DELTA. However, according to a test on 230206, it does not improve the
        ! performance if we skip the update of DELTA when ADEQUATE_GEO is FALSE and RATIO < 0.1.
        ! Therefore, we choose to update DELTA without checking ADEQUATE_GEO.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if

        ! Is the newly generated X better than current best point?
        ximproved = (actrem > 0)  ! If ACTREM is NaN, then XIMPROVED should & will be FALSE.

        ! Set JDROP_TR to the index of the vertex to be replaced with X. JDROP_TR = 0 means there
        ! is no good point to replace, and X will not be included into the simplex; in this case,
        ! the geometry of the simplex likely needs improvement, which will be handled below.
        jdrop_tr = setdrop_tr(ximproved, d, sim, simi)

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
    end if  ! End of IF (SHORTD .OR. .NOT. MAX(PREREC, PREREF) > 0). The normal trust-region calculation ends.


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
    bad_trstep = (shortd .or. (.not. max(prerec, preref) > 0) .or. ratio <= 0 .or. jdrop_tr == 0)
    ! IMPROVE_GEO: Should we take a geometry step to improve the geometry of the interpolation set?
    improve_geo = (bad_trstep .and. .not. adequate_geo)
    ! REDUCE_RHO: Should we enhance the resolution by reducing RHO?
    reduce_rho = (bad_trstep .and. adequate_geo .and. max(delta, dnorm) <= rho)

    ! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    !call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO or REDUCE_RHO are not both TRUE', srname)

    ! If SHORTD is TRUE or MAX(PREREC, PREREF) > 0 is FALSE, then either IMPROVE_GEO or REDUCE_RHO
    ! is TRUE unless ADEQUATE_GEO is TRUE and MAX(DELTA, DNORM) > RHO.
    !call assert((.not. shortd .and. max(prerec, preref) > 0) .or. (improve_geo .or. reduce_rho .or. &
    !    & (adequate_geo .and. max(delta, dnorm) > rho)), 'If SHORTD is TRUE or MAX(PREREC, PREREF) > 0 is FALSE, then&
    !    & either IMPROVE_GEO or REDUCE_RHO is TRUE unless ADEQUATE_GEO is TRUE and MAX(DELTA, DNORM) > RHO', srname)
    !----------------------------------------------------------------------------------------------!

    ! Comments on BAD_TRSTEP:
    ! 1. Powell's definition of BAD_TRSTEP is as follows. The one used above seems to work better,
    ! especially for linearly constrained problems due to the factor TENTH (= ETA1).
    ! !bad_trstep = (shortd .or. actrem <= 0 .or. actrem < TENTH * prerem .or. jdrop_tr == 0)
    ! Besides, Powell did not check MAX(PREREC, PREREF) > 0 in BAD_TRSTEP, which is reasonable to do
    ! but has little impact upon the performance.
    ! 2. NEWUOA/BOBYQA/LINCOA would define BAD_TRSTEP, IMPROVE_GEO, and REDUCE_RHO as follows. Two
    ! different thresholds are used in BAD_TRSTEP. It outperforms Powell's version.
    ! !bad_trstep = (shortd .or. (.not. max(prerec, preref) > 0) .or. ratio <= eta1 .or. jdrop_tr == 0)
    ! !improve_geo = bad_trstep .and. .not. adequate_geo
    ! !bad_trstep = (shortd .or. (.not. max(prerec, preref) > 0) .or. ratio <= 0 .or. jdrop_tr == 0)
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
    ! REDUCE_RHO to true if they are small (e.g., ALL(ABS(MODERRSAV) <= 0.1 * MAXVAL(ABS(A))*RHO) or
    ! ALL(ABS(MODERRSAV) <= RHO**2)) when SHORTD is TRUE. It made little impact on the performance.


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
        d = geostep(jdrop_geo, cpen, conmat, cval, delta, fval, factor_gamma, simi)

        x = sim(:, n + 1) + d
        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evaluate(calcfc, x, f, constr, cstrv)
        nf = nf + 1_IK

        ! Print a message about the function/constraint evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x, cstrv, constr)
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
        delta = HALF * rho
        rho = redrho(rho, rhoend)
        delta = max(delta, rho)
        cpen = min(cpen, fcratio(fval, conmat)) ! The 2nd (out of 2) update of CPEN. It may become 0
        ! Print a message about the reduction of RHO according to IPRINT.
        call rhomsg(solver, iprint, nf, fval(n + 1), rho, sim(:, n + 1), cval(n + 1), conmat(:, n + 1), cpen)
        call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
        ! Check whether to exit due to damaging rounding in UPDATEPOLE.
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit  ! Better action to take? Geometry step, or simply continue?
        end if
    end if  ! End of IF (REDUCE_RHO). The procedure of reducing RHO ends.

end do  ! End of DO TR = 1, MAXTR. The iterative procedure ends.

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
        & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) &
        & .or. is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain negative values or NaN/+Inf', srname)
    nhist = minval([nf, maxfhist, maxchist])
    call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f, cstrv, ctol)), &
        & 'No point in the history is better than X', srname)
end if

end subroutine cobylb


function fcratio(fval, conmat) result(r)
!--------------------------------------------------------------------------------------------------!
! This function calculates the ratio between the "typical change" of F and that of CONSTR.
! See equations (12)--(13) in Section 3 of the COBYLA paper for the definition of the ratio.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: CMAX(M), CMIN(M)
! Size of local arrays: REAL(RP)*(2*M)
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: fval(:)     ! FVAL(N+1)
real(RP), intent(in) :: conmat(:, :)    ! CONMAT(M, N+1)

! Outputs
real(RP) :: r

! Local variables
real(RP) :: cmax(size(conmat, 1))
real(RP) :: cmin(size(conmat, 1))
real(RP) :: denom
character(len=*), parameter :: srname = 'FCRATIO'

! Preconditions
if (DEBUGGING) then
    call assert(size(fval) >= 1, 'SIZE(FVAL) >= 1', srname)
    call assert(size(conmat, 2) == size(fval), 'SIZE(CONMAT, 2) == SIZE(FVAL)', srname)
end if

!====================!
! Calculation starts !
!====================!

cmin = minval(conmat, dim=2)
cmax = maxval(conmat, dim=2)
if (any(cmin < HALF * cmax)) then
    denom = minval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
    r = (maxval(fval) - minval(fval)) / denom
else
    r = ZERO
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
