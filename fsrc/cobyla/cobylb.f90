module cobylb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of COBYLA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Started: July 2021
!
! Last Modified: Sunday, November 21, 2021 PM08:45:57
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobylb


contains


subroutine cobylb(calcfc, iprint, maxfun, ctol, ftarget, rhobeg, rhoend, constr, x, nf, chist, &
        & conhist, cstrv, f, fhist, xhist, info)

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evalfc
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: info_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS, NAN_MODEL, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : inprod, matprod, outprod, inv
use, non_intrinsic :: output_mod, only : retmssg, rhomssg, fmssg
use, non_intrinsic :: pintrf_mod, only : FUNCON
use, non_intrinsic :: resolution_mod, only : resenhance
use, non_intrinsic :: selectx_mod, only : savefilt, selectx, isbetter

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : goodgeo, setdrop_geo, setdrop_tr, geostep
use, non_intrinsic :: initialize_mod, only : initxfc, initfilt
use, non_intrinsic :: trustregion_mod, only : trstlp
use, non_intrinsic :: update_mod, only : updatexfc, updatepole, findpole

implicit none

! Inputs
procedure(FUNCON) :: calcfc
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ctol
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(out) :: constr(:) ! M
real(RP), intent(inout) :: x(:)  ! N

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)
real(RP), intent(out) :: conhist(:, :)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: xhist(:, :)

! Parameters
integer(IK), parameter :: maxfilt = 2000_IK  ! Must be positive. Recommended to be in [100, 10,000].

! Local variables
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
logical :: enhance_resolut
logical :: evaluated(size(x) + 1)
logical :: good_geo
logical :: improve_geo
logical :: shortd
real(RP) :: A(size(x), size(constr) + 1)
! A(:, 1:M) contains the approximate gradient for the constraints, and A(:, M+1) is minus the
! approximate gradient for the objective function.
real(RP) :: actrem
real(RP) :: b(size(constr) + 1)
real(RP) :: barmu
real(RP) :: cfilt(min(max(maxfilt, 0), maxfun))
real(RP) :: confilt(size(constr), size(cfilt))
real(RP) :: conmat(size(constr), size(x) + 1)
real(RP) :: cpen  ! Penalty parameter for constraint in merit function (PARMU in Powell's code)
real(RP) :: cval(size(x) + 1)
real(RP) :: d(size(x))
real(RP) :: factor_alpha
real(RP) :: factor_beta
real(RP) :: factor_delta
real(RP) :: factor_gamma
real(RP) :: ffilt(size(cfilt))
real(RP) :: fval(size(x) + 1)
real(RP) :: prerec  ! Predicted reduction in Constraint violation
real(RP) :: preref  ! Predicted reduction in objective Function
real(RP) :: prerem  ! Predicted reduction in Merit function
real(RP) :: rho
real(RP) :: sim(size(x), size(x) + 1)  ! (n, )
real(RP) :: simi(size(x), size(x))  ! (n, )
real(RP) :: xfilt(size(x), size(cfilt))

! Sizes
m = size(constr)
n = size(x)
maxxhist = size(xhist, 2)
maxfhist = size(fhist)
maxconhist = size(conhist, 2)
maxchist = size(chist)
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(maxfun >= n + 2, 'MAXFUN >= N + 2', srname)
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
    !call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    !call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
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

! Set the initial values of some parameters. The last column of SIM holds the optimal vertex of the
! current simplex, and the preceding N columns hold the displacements from the optimal vertex to the
! other vertices.  Further, SIMI holds the inverse of the matrix that is contained in the first N
! columns of SIM.
factor_alpha = QUART
factor_beta = 2.1E0_RP
factor_delta = 1.1E0_RP
factor_gamma = HALF
rho = rhobeg
cpen = ZERO

call initxfc(calcfc, iprint, maxfun, ctol, ftarget, rhobeg, x, nf, chist, conhist, conmat, cval, fhist,&
   & fval, sim, xhist, evaluated, subinfo)
call initfilt(conmat, ctol, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)

if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Return the best calculated values of the variables.
    ! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
    cpen = min(1.0E8_RP, HUGENUM)
    kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cpen, ctol)
    x = xfilt(:, kopt)
    f = ffilt(kopt)
    constr = confilt(:, kopt)
    cstrv = cfilt(kopt)
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
    !close (16)
    return
end if

! SIMI is the inverse of SIM(:, 1:N), which is lower triangular by Powell's initialization.
simi = inv(sim(:, 1:n))

! Initialize PREREM, ACTREM, JDROP_TR, and JDROP_GEO, or some compilers will complain that they are
! uninitialized when setting BAD_TRSTEP. Indeed, these values will not be used, because they will be
! overwritten when SHORTD = FALSE.
prerem = ONE
actrem = -ONE
jdrop_tr = 0_IK
jdrop_geo = 0_IK

! In most cases, each trust-region iteration takes at most two function evaluation. The following
! setting essentially imposes no constraint on the maximal number of trust-region iterations.
maxtr = 4_IK * maxfun
! MAXTR is unlikely to be reached, but we define the following default value for INFO for safety.
info = MAXTR_REACHED

! We must initialize ACTREM and PREREM. Otherwise, when SHORTD = TRUE, compilers may raise a
! run-time error that they are undefined. The values will not be used: when SHORTD = FALSE, they
! will be overwritten; when SHORTD = TRUE, the values are used only in BAD_TRSTEP, which is TRUE
! regardless of ACTREM or PREREM. Similar for JDROP_TR.
actrem = -HUGENUM
prerem = HUGENUM
jdrop_tr = 0_IK

! Begin the iterative procedure.
! After solving a trust-region subproblem, COBYLA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust-region iteration? If yes, a geometry step
! will be taken, corresponding to the Branch (Delta) in the COBYLA paper.
! REDUCE_RHO - Will we reduce rho after the trust-region iteration?
! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Before the trust-region step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
    !-----------------------------------------------------------------------!
    call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
    if (subinfo == DAMAGING_ROUNDING) then
        info = subinfo
        exit
        ! Instead of exiting, what about setting GOOD_GEO to FALSE in order to activate a
        ! geometry step of enhance_resolut?
    end if
    !-----------------------------------------------------------------------!

    ! Does the current interpolation set has good geometry? It affects IMPROVE_GEO and REDUCE_RHO.
    good_geo = goodgeo(factor_alpha, factor_beta, rho, sim, simi)

    ! Calculate the linear approximations to the objective and constraint functions, placing minus
    ! the objective function gradient after the constraint gradients in the array A.
    ! N.B.:
    ! 1. When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result
    ! as Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
    ! finite-precision arithmetic.
    ! 2. TRSTLP accesses A mostly by columns, so it is not more reasonable to save A^T instead of A.
    A(:, 1:m) = transpose(matprod(conmat(:, 1:n) - spread(conmat(:, n + 1), dim=2, ncopies=n), simi))
    A(:, m + 1) = matprod(fval(n + 1) - fval(1:n), simi)

    ! Exit if A contains NaN. Otherwise, TRSTLP may encounter memory errors or infinite loops.
    ! HOW EXACTLY?????
    !----------------------------------------------------------------------------------------------!
    ! POSSIBLE IMPROVEMENT: INSTEAD OF EXITING, SKIP A TRUST-REGION STEP AND PERFORM A GEOMETRY ONE!
    !----------------------------------------------------------------------------------------------!
    if (any(is_nan(A))) then
        info = NAN_MODEL
        exit
    end if

    ! Theoretically (but not numerically), the last entry of B does not affect the result of TRSTLP.
    ! We set it to -FVAL(N + 1) following Powell's code.
    b = [-conmat(:, n + 1), -fval(n + 1)]
    ! Calculate the trust-region trial step D.
    d = trstlp(A, b, rho)

    ! Is the trust-region trial step short?
    shortd = (inprod(d, d) < QUART * rho**2)

    if (.not. shortd) then
        ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
        preref = inprod(d, A(:, m + 1))
        prerec = cval(n + 1) - maxval([-matprod(d, A(:, 1:m)) - conmat(:, n + 1), ZERO])

        ! Increase CPEN if necessary and branch back if this change alters the optimal vertex.
        ! See the discussions around equation (9) of the COBYLA paper.
        barmu = -preref / prerec   ! PREREF + BARMU * PREREC = 0
        !!!!!!!!!!!!!!! Is it possible that PREREC <= 0????????????? It seems yes, but why?
        if (prerec > ZERO .and. cpen < 1.5E0_RP * barmu) then
            cpen = min(TWO * barmu, HUGENUM)
            if (findpole(cpen, cval, fval) <= n) then
                ! Zaikun 20211111: Can this lead to infinite cycling?
                !!-----------------------------------------------------------------------!
                !call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
                !if (subinfo == DAMAGING_ROUNDING) then
                !    info = subinfo
                !    exit
                !end if
                !!-----------------------------------------------------------------------!
                cycle
            end if
        end if

        x = sim(:, n + 1) + d
        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evalfc(calcfc, x, f, constr, cstrv)
        nf = nf + 1_IK
        ! Save X, F, CONSTR, CSTRV into the history.
        call savehist(nf, constr, cstrv, f, x, chist, conhist, fhist, xhist)
        ! Save X, F, CONSTR, CSTRV into the filter.
        call savefilt(constr, cstrv, ctol, f, x, nfilt, cfilt, confilt, ffilt, xfilt)

        ! Begin the operations that decide whether X should replace one of the vertices of the
        ! current simplex, the change being mandatory if ACTREM is positive.
        ! PREREM and ACTREM are the predicted and actual reductions in the merit function respectively.
        prerem = preref + cpen * prerec   ! Is it positive????
        actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv)
        if (cpen <= ZERO .and. f <= fval(n + 1) .and. f >= fval(n + 1)) then
            ! CPEN <= ZERO indeed means CPEN == ZERO, while A <= B .and. A >= B indeed mean A == B.
            ! We write them in this way to avoid compilers complaining about equality comparison
            ! between reals, which appears in the original code of Powell.
            prerem = prerec   ! Is it positive?????
            actrem = cval(n + 1) - cstrv
        end if
        if (is_nan(actrem)) then
            actrem = -HUGENUM  ! Signify a bad trust-region step.
        end if
        ! Set JDROP_TR to the index of the vertex that is to be replaced by X.
        ! N.B.: COBYLA never sets JDROP_TR = N + 1.
        jdrop_tr = setdrop_tr(actrem, d, factor_alpha, factor_delta, rho, sim, simi)
        ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_TR) is replaced by D.
        ! When JDROP_TR == 0, the algorithm decides not to include X into the simplex.
        ! N.B.:
        ! 1. UPDATEXFC does NOT manipulate the simplex so that SIM(:, N+1) is the best vertex;
        ! that is the job of UPDATEPOLE, which is called before each trust-region/geometry step.
        ! 2. UPDATEXFC does nothing when JDROP_TR == 0.
        call updatexfc(jdrop_tr, constr, cstrv, d, f, conmat, cval, fval, sim, simi)
        !!-----------------------------------------------------------------------!
        !call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
        !if (subinfo == DAMAGING_ROUNDING) then
        !    info = subinfo
        !    exit
        !end if
        !!-----------------------------------------------------------------------!

        ! Check whether to exit.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end if

    ! Should we take a geometry step to improve the geometry of the interpolation set?
    ! N.B.: THEORETICALLY, JDROP_TR > 0 when ACTREM > 0, and hence the definition of BAD_TRSTEP is
    ! mathematically equivalent to (SHORTD .OR. ACTREM <= ZERO .OR. ACTREM < TENTH * PREREM).
    ! However, Powell's code can set JDROP_TR = 0 when ACTREM >0 due to NaN. This has been rectified
    ! in the subroutine SETDROP_TR. Nevertheless, we still keep JDROP_TR == 0 for robustness.
    bad_trstep = (shortd .or. actrem <= ZERO .or. actrem < TENTH * prerem .or. jdrop_tr == 0)
    improve_geo = bad_trstep .and. .not. good_geo

    ! Should we enhance the resolution by reducing RHO?
    enhance_resolut = bad_trstep .and. good_geo

    if (improve_geo) then
        ! Before the geometry step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
        !-----------------------------------------------------------------------!
        call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit
        end if
        !-----------------------------------------------------------------------!

        ! If the current interpolation set has good geometry, then we skip the geometry step.
        ! The code has a small difference from Powell's original code here: If the current geometry
        ! is good, then we will continue with a new trust-region iteration; at the beginning of the
        ! iteration, CPEN may be updated, which may alter the pole point SIM(:, N + 1) by UPDATEPOLE;
        ! the quality of the interpolation point depends on SIM(:, N + 1), meaning that the same
        ! interpolation set may have good or bad geometry with respect to different "poles"; if the
        ! geometry turns out bad with the new pole, the original COBYLA code will take a geometry
        ! step, but our code here will NOT do it but continue to take a trust-region step.
        ! The argument is this: even if the geometry step is not skipped at the first place, the
        ! geometry may turn out bad again after the pole is altered due to an update to CPEN; should
        ! we take another geometry step in that case? If no, why should we do it here? Indeed, this
        ! distinction makes no practical difference for CUTEst problems with at most 100 variables
        ! and 5000 constraints, while the algorithm framework is simplified.
        if (.not. goodgeo(factor_alpha, factor_beta, rho, sim, simi)) then
            ! Decide a vertex to drop from the simplex. It will be replaced by SIM(:, N + 1) + D to
            ! improve acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
            ! N.B.: COBYLA never sets JDROP_GEO = N + 1.
            jdrop_geo = setdrop_geo(factor_alpha, factor_beta, rho, sim, simi)

            ! If JDROP_GEO == 0 (due to NaN in SIM or SIMI), then we exit. Without this, memory
            ! error will occur as JDROP_GEO will be used as an index of arrays.
            if (jdrop_geo == 0) then
                info = DAMAGING_ROUNDING
                exit
            end if

            ! Calculate the geometry step D.
            d = geostep(jdrop_geo, cpen, conmat, cval, fval, factor_gamma, rho, simi)
            x = sim(:, n + 1) + d
            ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
            call evalfc(calcfc, x, f, constr, cstrv)
            nf = nf + 1_IK
            ! Save X, F, CONSTR, CSTRV into the history.
            call savehist(nf, constr, cstrv, f, x, chist, conhist, fhist, xhist)
            ! Save X, F, CONSTR, CSTRV into the filter.
            call savefilt(constr, cstrv, ctol, f, x, nfilt, cfilt, confilt, ffilt, xfilt)
            ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP_GEO) is replaced by D.
            !--------------------------------------------------------------------------------------!
            ! N.B.: UPDATEXFC does NOT manipulate the simplex so that SIM(:, N+1) is the best vertex;
            ! that is the job of UPDATEPOLE, which is called before each trust-region/geometry step.
            !--------------------------------------------------------------------------------------!
            call updatexfc(jdrop_geo, constr, cstrv, d, f, conmat, cval, fval, sim, simi)
            !!-----------------------------------------------------------------------!
            !call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
            !if (subinfo == DAMAGING_ROUNDING) then
            !    info = subinfo
            !    exit
            !end if
            !!-----------------------------------------------------------------------!
            ! Check whether to exit.
            subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
        end if
    end if

    if (enhance_resolut) then  ! Enhance the resolution of the algorithm, updating RHO and CPEN.
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        call resenhance(conmat, fval, rhoend, cpen, rho)
        !!-----------------------------------------------------------------------!
        !call updatepole(cpen, conmat, cval, fval, sim, simi, subinfo)
        !if (subinfo == DAMAGING_ROUNDING) then
        !    info = subinfo
        !    exit
        !end if
        !!-----------------------------------------------------------------------!
    end if
end do

! Return the best calculated values of the variables.
! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
cpen = max(cpen, min(1.0E8_RP, HUGENUM))
kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cpen, ctol)
x = xfilt(:, kopt)
f = ffilt(kopt)
constr = confilt(:, kopt)
cstrv = cfilt(kopt)

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

!close (16)

end subroutine cobylb


end module cobylb_mod

! TODO:
! BTW, write a projection function to be used in TRSAPP, BIGLAG, BIGDEN in NEWUOA. See lines 305-306
! of trustregion.f90 and 391-392, 675-676 of geometry.f90 of NEWUOA.
!
! HYPOT is used  in trustregion.f90 for updating ZDOTA. A robust version should be written (avoid
! under/overflow in particular). This version may also be used in PLANEROT. This has been done
! previously, but it
! worsened a bit the performance of NEWUOA so it was discarded. Find it back, and test NEWUOA again.
!
! Try calculating ZDOTA from scratch (only the elements that change; there are 1 or 2 of them) instead of
! updating it. Is it more stable? Will it improve (at least not worsen) the performance of COBYLA?
! Needs tests. If it is accepted, then HYPOT is not needed.
!
! 0. In COBYLA, check what should we do with JDROP = 0, both TR and GEO
!    If ACTREM > 0 ===> JDROP > 0, why can COBYLA return sub-optimal points???
! 1. evalfc, extreme barrier, moderate excessively negative objective, which has not been done in
!    NEWUOA. Shouldn't we remove the extreme barrier in the MATLAB/Python interface after it is
!    implemented in FORTRAN?
! 3. merge UPDATEPOLE and UPDATEXFC
! 6. Do the same for NEWUOA
! 8. knew ===> jdrop
! 11.
! Enforcing programming contracts
! Programming can be thought of as requirements for correct execution of a procedure and assurances
! for the result of correct execution. The requirements and assurances might be constraints of three
! kinds:
! Preconditions (requirements): logical expressions that must evaluate to .true. when a procedure starts execution,
! Postconditions (assurances): expressions that must evaluate to .true. when a procedure finishes execution, and
! Invariants: universal pre- and postconditions that must always be true when all procedures in a class start or finish executing.
