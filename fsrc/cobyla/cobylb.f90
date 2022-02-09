module cobylb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of COBYLA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Started: July 2021
!
! Last Modified: Wednesday, February 09, 2022 AM12:30:21
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobylb


contains


subroutine cobylb(calcfc, iprint, maxfilt, maxfun, ctol, cweight, eta1, eta2, ftarget, &
    & gamma1, gamma2, rhobeg, rhoend, constr, f, x, nf, chist, conhist, cstrv, fhist, xhist, info)

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, QUART, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: consts_mod, only : MIN_MAXFILT
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: info_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: output_mod, only : retmsg, rhomsg, fmsg, cpenmsg
use, non_intrinsic :: pintrf_mod, only : OBJCON
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: selectx_mod, only : savefilt, selectx, isbetter

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : goodgeo, setdrop_geo, setdrop_tr, geostep
use, non_intrinsic :: initialize_mod, only : initxfc, initfilt
use, non_intrinsic :: trustregion_mod, only : trstlp, trrad
use, non_intrinsic :: update_mod, only : updatexfc, updatepole, findpole

implicit none

! Inputs
procedure(OBJCON) :: calcfc
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
real(RP), intent(inout) :: constr(:) ! M
real(RP), intent(inout) :: f
real(RP), intent(inout) :: x(:)  ! N

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)
real(RP), intent(out) :: conhist(:, :)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: xhist(:, :)

! Local variables
character(len=*), parameter :: solver = 'COBYLA'
character(len=*), parameter :: srname = 'COBYLB'
integer(IK) :: jdrop_geo
integer(IK) :: jdrop_tr
integer(IK) :: k
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
integer(IK) :: subinfo
integer(IK) :: tr
logical :: bad_trstep
logical :: evaluated(size(x) + 1)
logical :: good_geo
logical :: improve_geo
logical :: reduce_rho
logical :: shortd
logical :: tr_success
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
real(RP) :: prerec  ! Predicted reduction in Constraint violation
real(RP) :: preref  ! Predicted reduction in objective Function
real(RP) :: prerem  ! Predicted reduction in Merit function
real(RP) :: rho
real(RP) :: sim(size(x), size(x) + 1)  ! (n, )
real(RP) :: simi(size(x), size(x))  ! (n, )
real(RP) :: xfilt(size(x), size(cfilt))
! N.B.: FACTOR_ALPHA, FACTOR_BETA, FACTOR_GAMMA, and FACTOR_DELTA are four factors the COBYLB uses
! when managing the simplex. Note the following
! 1. FACTOR_ALPHA < FACTOR_GAMMA < 1 < FACTOR_DELTA <= FACTOR_BETA.
! 2. FACTOR_DELTA has nothing to do with DELTA, which is the trust-region radius.
! 3. FACTOR_GAMMA has nothing to do with GAMMA1 and GAMMA2, which are the contracting/expanding
!    factors for updating the trust-region radius DELTA.
real(RP), parameter :: factor_alpha = QUART  ! The factor alpha in the COBYLA paper
real(RP), parameter :: factor_beta = 2.1_RP  ! The factor better in the COBYLA paper
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
    call assert(n >= 1, 'N >= 1', srname)
    call assert(maxfun >= n + 2, 'MAXFUN >= N + 2', srname)
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
    !call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    !call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
    call assert(maxfilt >= min(MIN_MAXFILT, maxfun) .and. maxfilt <= maxfun, &
        & 'MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(conhist, 1) == m .and. maxconhist * (maxconhist - maxhist) == 0, &
        & 'SIZE(CONHIST, 1) == M, SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
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
call initfilt(conmat, ctol, cweight, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)

! Check whether to exit due to abnormal cases that may occur during the initialization.
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
        call assert(.not. any(is_nan(chist(1:min(nf, maxchist))) .or. is_posinf(chist(1:min(nf, maxchist)))), &
            & 'CHIST does not contain NaN/+Inf', srname)
        call assert(.not. any([(isbetter([fhist(k), chist(k)], [f, cstrv], ctol), &
            & k=1, minval([nf, maxfhist, maxchist]))]), 'No point in the history is better than X', srname)
    end if
    return
end if

! Initialize RHO, DELTA, and CPEN.
rho = rhobeg
delta = rhobeg
cpen = ZERO

! We must initialize ACTREM and PREREM. Otherwise, when SHORTD = TRUE, compilers may raise a
! run-time error that they are undefined. The values will not be used: when SHORTD = FALSE, they
! will be overwritten; when SHORTD = TRUE, the values are used only in BAD_TRSTEP, which is TRUE
! regardless of ACTREM or PREREM. Similar for JDROP_TR.
actrem = -HUGENUM
prerem = HUGENUM
jdrop_tr = 0_IK
jdrop_geo = 0_IK

! MAXTR is the maximal number of trust-region iterations. Normally, each trust-region iteration 
! takes 1 or 2 function evaluations unless the update of CPEN alters the optimal vertex or the
! trust-region step is short but the geometry step is not invoked. Thus the following MAXTR is 
! unlikely to be reached.
maxtr = max(maxfun, 2_IK * maxfun)  ! MAX: precaution against overflow, which will make 2*MAXFUN < 0.
info = MAXTR_REACHED

! Begin the iterative procedure.
! After solving a trust-region subproblem, COBYLA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust-region iteration? If yes, a geometry step
! will be taken, corresponding to the "Branch (Delta)" in the COBYLA paper.
! REDUCE_RHO - Will we reduce rho after the trust-region iteration?
! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Before the trust-region step, UPDATEPOLE has been called either implicitly by INITXFC/UPDATEXFC
    ! or explicitly after CPEN is updated, so that SIM(:, N + 1) is the optimal vertex.

    ! Does the current interpolation set have good geometry? It affects IMPROVE_GEO and REDUCE_RHO.
    good_geo = goodgeo(delta, factor_alpha, factor_beta, sim, simi)

    ! Calculate the linear approximations to the objective and constraint functions, placing minus
    ! the objective function gradient after the constraint gradients in the array A.
    ! N.B.:
    ! 1. When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result
    ! as Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
    ! finite-precision arithmetic.
    ! 2. TRSTLP accesses A mostly by columns, so it is not more reasonable to save A^T instead of A.
    A(:, 1:m) = transpose(matprod(conmat(:, 1:n) - spread(conmat(:, n + 1), dim=2, ncopies=n), simi))
    A(:, m + 1) = matprod(fval(n + 1) - fval(1:n), simi)

    ! Theoretically (but not numerically), the last entry of B does not affect the result of TRSTLP.
    ! We set it to -FVAL(N + 1) following Powell's code.
    b = [-conmat(:, n + 1), -fval(n + 1)]
    ! Calculate the trust-region trial step D. Note that D does NOT depend on CPEN.
    d = trstlp(A, b, delta)
    dnorm = min(delta, norm(d))

    ! Is the trust-region trial step short? N.B.: we compare DNORM with RHO, not DELTA.
    !!shortd = (dnorm < HALF * rho)  ! Powell's version
    ! TENTH seems to work better than HALF or QUART, especially for linearly constrained problems.
    shortd = (dnorm < TENTH * rho)

    if (shortd) then
        ! Reduce DELTA if D is short. This seems quite important for performance.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close.
        end if
    else
        ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
        ! We have the following in precise arithmetic. They may fail to hold due to rounding errors.
        ! 1. PREREC >= 0; PREREC == 0 iff B <= 0.
        ! 2. PREREM >= 0 because of the update of CPEN detailed in the following lines; PREREM == 0
        ! iff PREREF and CPEN are both 0.
        preref = inprod(d, A(:, m + 1))  ! Can be negative.
        prerec = cval(n + 1) - maxval([-matprod(d, A(:, 1:m)) - conmat(:, n + 1), ZERO])

        ! Increase CPEN if necessary and branch back if this change alters the optimal vertex.
        ! See the discussions around equation (9) of the COBYLA paper.
        if (prerec > 0) then
            barmu = -preref / prerec   ! PREREF + BARMU * PREREC = 0
        else  ! PREREC == 0 can happen if B <= 0.
            barmu = ZERO
        end if
        if (cpen < 1.5_RP * barmu) then
            ! This can happen only if BARMU > 0, and hence PREREC > 0 > PREREF.
            ! If CPEN == 0 and PREREC > 0 > PREREF, then CPEN will be updated to 2*BARMU > 0.
            cpen = min(TWO * barmu, HUGENUM)
            call cpenmsg(solver, iprint, cpen)
            if (findpole(cpen, cval, fval) <= n) then
                call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
                ! Check whether to exit due to damaging rounding detected in UPDATEPOLE.
                if (subinfo == DAMAGING_ROUNDING) then
                    info = subinfo
                    exit  ! Better action to take? Geometry step?
                end if
                cycle  ! Zaikun 20211111: Can this lead to infinite cycling?
            end if
        end if

        x = sim(:, n + 1) + d
        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evaluate(calcfc, x, f, constr, cstrv)
        nf = nf + 1_IK
        call fmsg(solver, iprint, nf, f, x, cstrv, constr)
        ! Save X, F, CONSTR, CSTRV into the history.
        call savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
        ! Save X, F, CONSTR, CSTRV into the filter.
        call savefilt(constr, cstrv, ctol, cweight, f, x, nfilt, cfilt, confilt, ffilt, xfilt)

        ! Begin the operations that decide whether X should replace one of the vertices of the
        ! current simplex, the change being mandatory if ACTREM is positive. PREREM and ACTREM are
        ! the predicted and actual reductions in the merit function respectively.
        prerem = preref + cpen * prerec  ! Theoretically nonnegative; equals 0 if CPEN = 0 = PREREF.
        actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv)
        if (cpen <= 0 .and. f <= fval(n + 1) .and. f >= fval(n + 1)) then
            ! CPEN <= 0 indeed means CPEN == 0, while A <= B .and. A >= B indeed mean A == B.
            ! We code in this way to avoid compilers complaining about equality comparison of reals.
            prerem = prerec
            actrem = cval(n + 1) - cstrv
        end if
        if (is_nan(actrem)) then
            actrem = -HUGENUM  ! Signify a bad trust-region step.
        end if

        !call assert(prerem >= 0, 'PREREM >= 0', 'COBYLA')
        !prerem = max(prerem, tiny(prerem))

        ratio = redrat(actrem, prerem, eta1)
        ! Update DELTA. After this, DELTA < DNORM may hold.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if

        ! Set JDROP_TR to the index of the vertex that is to be replaced by X.
        ! N.B.: COBYLA never sets JDROP_TR = N + 1.
        tr_success = (actrem > 0)
        jdrop_tr = setdrop_tr(tr_success, d, delta, factor_alpha, factor_delta, sim, simi)

        ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_TR) is replaced by D.
        ! When JDROP_TR == 0, the algorithm decides not to include X into the simplex.
        ! N.B.: UPDATEXFC does nothing when JDROP_TR == 0.
        call updatexfc(jdrop_tr, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi, subinfo)
        ! Check whether to exit due to damaging rounding detected in UPDATEPOLE (called by UPDATEXFC).
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit  ! Better action to take? Geometry step?
        end if

        ! Check whether to exit.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit  ! Better action to take? Geometry step?
        end if
    end if

    ! Is the trust-region step a bad one?
    ! N.B.:
    ! 1. THEORETICALLY, JDROP_TR > 0 when ACTREM > 0. Yet Powell's code may set JDROP_TR = 0 when
    ! ACTREM > 0 due to NaN. The modernized code has been rectified this in the function SETDROP_TR.
    ! After this rectification, we can  indeed simplify the definition of BAD_TRSTEP below by
    ! removing (JDROP_TR == 0), but we retain (JDROP_TR == 0) for robustness.
    ! 2. Powell's definition of BAD_TRSTEP is
    !!bad_trstep = (shortd .or. actrem <= 0 .or. actrem < TENTH * prerem .or. jdrop_tr == 0)
    ! But the following one seems to work better, especially for linearly constrained problems.
    bad_trstep = (shortd .or. actrem <= 0 .or. jdrop_tr == 0)

    ! Should we take a geometry step to improve the geometry of the interpolation set?
    improve_geo = bad_trstep .and. .not. good_geo

    ! Should we enhance the resolution by reducing RHO?
    reduce_rho = bad_trstep .and. good_geo .and. (max(delta, dnorm) <= rho)

    ! Improve the geometry of the simplex by removing a point and adding a new one.
    if (improve_geo) then
        ! Before the trust-region step, UPDATEPOLE has been called either implicitly by UPDATEXFC or
        ! explicitly after CPEN is updated, so that SIM(:, N + 1) is the optimal vertex.

        ! If the current interpolation set has good geometry, then we skip the geometry step.
        ! The code has a small difference from Powell's original code here: If the current geometry
        ! is good, then we will continue with a new trust-region iteration; however, at the
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
        if (.not. goodgeo(delta, factor_alpha, factor_beta, sim, simi)) then
            ! Decide a vertex to drop from the simplex. It will be replaced by SIM(:, N + 1) + D to
            ! improve acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
            ! N.B.: COBYLA never sets JDROP_GEO = N + 1.
            jdrop_geo = setdrop_geo(delta, factor_alpha, factor_beta, sim, simi)

            ! JDROP_GEO is between 1 and N unless SIM and SIMI contains NaN, which should not happen
            ! at this point unless there is a bug. Nevertheless, for robustness, we include the
            ! following instruction to exit when JDROP_GEO == 0 (if JDROP_GEO does become 0, then
            ! memory error will occur if we continue, as JDROP_GEO is used as an index of arrays.)
            if (jdrop_geo == 0) then
                info = DAMAGING_ROUNDING
                exit
            end if

            ! Calculate the geometry step D.
            d = geostep(jdrop_geo, cpen, conmat, cval, delta, fval, factor_gamma, simi)

            x = sim(:, n + 1) + d
            ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
            call evaluate(calcfc, x, f, constr, cstrv)
            nf = nf + 1_IK
            call fmsg(solver, iprint, nf, f, x, cstrv, constr)
            ! Save X, F, CONSTR, CSTRV into the history.
            call savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
            ! Save X, F, CONSTR, CSTRV into the filter.
            call savefilt(constr, cstrv, ctol, cweight, f, x, nfilt, cfilt, confilt, ffilt, xfilt)
            ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_GEO) is replaced by D.
            call updatexfc(jdrop_geo, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi, subinfo)
            ! Check whether to exit due to damaging rounding detected in UPDATEPOLE (called by UPDATEXFC).
            if (subinfo == DAMAGING_ROUNDING) then
                info = subinfo
                exit  ! Better action to take? Geometry step?
            end if

            ! Check whether to exit.
            subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
        end if
    end if

    ! Enhance the resolution of the algorithm by reducing RHO; update DELTA and CPEN at the same time.
    if (reduce_rho) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        delta = HALF * rho
        rho = redrho(rho, rhoend)
        delta = max(delta, rho)
        cpen = min(cpen, fcratio(fval, conmat))  ! It may set CPEN to 0.
        call rhomsg(solver, iprint, nf, fval(n + 1), rho, sim(:, n + 1), cval(n + 1), conmat(:, n + 1), cpen)
        call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
        ! Check whether to exit due to damaging rounding detected in UPDATEPOLE.
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit  ! Better action to take? Geometry step?
        end if
    end if
end do

! Return the best calculated values of the variables.
! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), max(cpen, cweight), ctol)
x = xfilt(:, kopt)
f = ffilt(kopt)
constr = confilt(:, kopt)
cstrv = cfilt(kopt)

! Arrange CHIST, CONHIST, FHIST, and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist, chist, conhist)

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
    call assert(.not. any(is_nan(chist(1:min(nf, maxchist))) .or. is_posinf(chist(1:min(nf, maxchist)))), &
        & 'CHIST does not contain NaN/+Inf', srname)
    call assert(.not. any([(isbetter([fhist(k), chist(k)], [f, cstrv], ctol), &
        & k=1, minval([nf, maxfhist, maxchist]))]), 'No point in the history is better than X', srname)
end if

end subroutine cobylb


function fcratio(fval, conmat) result(r)
!--------------------------------------------------------------------------------------------------!
! This function calculates the ratio between the "typical change" of F and that of CONSTR.
! See equations (12)--(13) in Section 3 of the COBYLA paper for the definition of the ratio.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: fval(:)
real(RP), intent(in) :: conmat(:, :)

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
