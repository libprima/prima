! TODO:
! 1. Improve RESCUE so that it accepts an [XNEW, FNEW] that is not interpolated yet, or even accepts
! [XRESERVE, FRESERVE], which contains points that have been evaluated.
!
module bobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! N.B. (Zaikun 20230312): In Powell's code, the strategy concerning RESCUE is a bit complex.
!
! 1. Suppose that a trust-region step D is calculated. Powell's code sets KNEW_TR before evaluating
! F at the trial point XOPT+D, assuming that the value of F at this point is not better than the
! current FOPT. With this KNEW_TR, the denominator of the update is calculated. If this denominator
! is sufficiently large, then evaluate F at XOPT+D, recalculate KNEW_TR if the function value turns
! out better than FOPT, and perform the update to include XOPT+D in the interpolation. If the
! denominator is not sufficiently large, then RESCUE is called, and another trust-region step is
! taken immediately after, discarding the previously calculated trust-region step D.
!
! 2. Suppose that a geometry step D is calculated. Then KNEW_GEO must have been set before. Powell's
! code then calculates the denominator of the update. If the denominator is sufficiently large, then
! evaluate F at XOPT+D, and perform the update. If the denominator is not sufficiently large, then
! RESCUE is called; if RESCUE does not evaluate F at any new point (allowed by Powell's code but not
! ours), then take a new geometry step, or else take a trust-region step, discarding the previously
! calculated geometry step D in both cases.
!
! 3. If it turns out necessary to call RESCUE again, but no new function value has been evaluated
! after the last RESCUE, then Powell's code will terminate.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, March 13, 2023 AM09:36:58
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: bobyqb


contains


subroutine bobyqb(calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, &
    & xl, xu, x, nf, f, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the major calculations of BOBYQA.
!
! IPRINT, MAXFUN, MAXHIST, NPT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, XL, XU, X, NF,
! F, FHIST, XHIST, and INFO are identical to the corresponding arguments in subroutine BOBYQA.
!
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
!   of the model and Lagrange functions.
! SL and SU hold XL - XBASE and XU - XBASE, respectively.
! XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
!   the least calculated F so far). XOPT satisfies SL(I) <= XOPT(I) <= SU(I), with appropriate
!   equalities when XOPT is on a constraint boundary. FOPT = F(XOPT + XBASE). However, we do not
!   save XOPT and FOPT explicitly, because XOPT = XPT(:, KOPT) and FOPT = FVAL(KOPT), which is
!   explained below.
! [XPT, FVAL, KOPT] describes the interpolation set:
! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
!   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
! [GOPT, HQ, PQ] describes the quadratic model: GOPT will hold the gradient of the quadratic model
!   at XBASE + XOPT; HQ will hold the explicit second order derivatives of the quadratic model; PQ
!   will contain the parameters of the implicit second order derivatives of the quadratic model.
! [BMAT, ZMAT] describes the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of
!   the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
! ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
!   OMEGA = ZMAT*ZMAT^T, which provides both the correct rank and positive semi-definiteness. BMAT
!   will hold the last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row
!   and column of H are not saved as they are unnecessary for the calculation.
! D is reserved for trial steps from XOPT. It is chosen by subroutine TRSBOX or GEOSTEP. Usually
!   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert!, wassert, validate
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: infos_mod, only : INFO_DFT, SMALL_TR_RADIUS, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : norm
use, non_intrinsic :: output_mod, only : retmsg, rhomsg, fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calden, calvlag!, errquad
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: shiftbase_mod, only : shiftbase
use, non_intrinsic :: xinbd_mod, only : xinbd

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep, setdrop_tr
use, non_intrinsic :: initialize_mod, only : initxf, initq, inith
use, non_intrinsic :: rescue_mod, only : rescue
use, non_intrinsic :: trustregion_mod, only : trsbox, trrad
use, non_intrinsic :: update_mod, only : updatexf, updateq, tryqalt, updateh

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: npt
real(RP), intent(in) :: eta1
real(RP), intent(in) :: eta2
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: gamma1
real(RP), intent(in) :: gamma2
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend
real(RP), intent(in) :: xl(:)  ! XL(N)
real(RP), intent(in) :: xu(:)  ! XU(N)

! In-outputs
real(RP), intent(inout) :: x(:)  ! X(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'BOBYQA'
character(len=*), parameter :: srname = 'BOBYQB'
integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK)))
integer(IK) :: itest
integer(IK) :: k
integer(IK) :: knew_geo
integer(IK) :: knew_tr
integer(IK) :: kopt
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxtr
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: subinfo
integer(IK) :: tr
logical :: accurate_mod
logical :: adequate_geo
logical :: bad_trstep
logical :: close_itpset
logical :: improve_geo
logical :: reduce_rho
logical :: shortd
logical :: small_trrad
logical :: ximproved
real(RP) :: bmat(size(x), npt + size(x))
real(RP) :: crvmin
real(RP) :: d(size(x))
real(RP) :: delbar
real(RP) :: delta
real(RP) :: den(npt)
real(RP) :: distsq(npt)
real(RP) :: dnorm
real(RP) :: dnormsav(2)  ! Powell's implementation: DNORMSAV(3)
real(RP) :: ebound
real(RP) :: fval(npt)
real(RP) :: gopt(size(x))
real(RP) :: hq(size(x), size(x))
real(RP) :: moderr
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pq(npt)
real(RP) :: qred
real(RP) :: ratio
real(RP) :: rho
real(RP) :: sl(size(x))
real(RP) :: su(size(x))
real(RP) :: vlag(npt + size(x))
real(RP) :: xbase(size(x))
real(RP) :: xdrop(size(x))
real(RP) :: xosav(size(x))
real(RP) :: xpt(size(x), npt)
real(RP) :: zmat(npt, npt - size(x) - 1)

! Sizes.
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT+1', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(all(rhobeg <= (xu - xl) / TWO), 'RHOBEG <= MINVAL(XU-XL)/2', srname)
    call assert(all(is_finite(x)), 'X is finite', srname)
    call assert(all(x >= xl .and. (x <= xl .or. x >= xl + rhobeg)), 'X == XL or X >= XL + RHOBEG', srname)
    call assert(all(x <= xu .and. (x >= xu .or. x <= xu - rhobeg)), 'X == XU or X >= XU - RHOBEG', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize XBASE, XPT, FVAL, and KOPT.
call initxf(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x, ij, kopt, nf, fhist, fval, &
    & sl, su, xbase, xhist, xpt, subinfo)
x = xinbd(xbase, xpt(:, kopt), xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT.
f = fval(kopt)

! Check whether to return due to abnormal cases that may occur during the initialization.
if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Arrange FHIST and XHIST so that they are in the chronological order.
    call rangehist(nf, xhist, fhist)
    ! Print a return message according to IPRINT.
    call retmsg(solver, info, iprint, nf, f, x)
    return
end if

! Initialize GOPT, HQ, and PQ.
call initq(ij, fval, xpt, gopt, hq, pq)

! Initialize BMAT and ZMAT.
call inith(ij, xpt, bmat, zmat)

! After initializing GOPT, HQ, PQ, BMAT, ZMAT, one can also choose to return if these arrays contain
! NaN. We do not do it here. If such a model is harmful, then it will probably lead to other returns
! (NaN in X, NaN in F, trust-region subproblem fails, ...); otherwise, the code will continue to run
! and possibly recovers by geometry steps.

! Set some more initial values.
! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
rho = rhobeg
delta = rho
ebound = ZERO
shortd = .false.
ratio = -ONE
dnormsav = REALMAX
moderrsav = REALMAX
knew_tr = 0
knew_geo = 0
itest = 0

! MAXTR is the maximal number of trust-region iterations. Each trust-region iteration takes 1 or 2
! function evaluations unless the trust-region step is short or fails to reduce the trust-region
! model but the geometry step is not invoked. Thus the following MAXTR is unlikely to be reached.
maxtr = max(maxfun, 2_IK * maxfun)  ! MAX: precaution against overflow, which will make 2*MAXFUN < 0.
info = MAXTR_REACHED

! Begin the iterative procedure.
! After solving a trust-region subproblem, we use three boolean variables to control the workflow.
! SHORTD: Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO: Should we improve the geometry?
! REDUCE_RHO: Should we reduce rho?
! BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
info = INFO_DFT
do tr = 1, maxtr
    ! Generate the next trust region step D.
    call trsbox(delta, gopt, hq, pq, sl, su, xpt(:, kopt), xpt, crvmin, d)
    dnorm = min(delta, norm(d))

    shortd = (dnorm < HALF * rho)

    ! Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
    ! should be positive If it is nonpositive due to rounding errors, we will not take this step.
    qred = -quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)

    ! When D is short, make a choice between reducing RHO and improving the geometry depending
    ! on whether or not our work with the current RHO seems complete. RHO is reduced if the
    ! errors in the quadratic model at the last three interpolation points compare favourably
    ! with predictions of likely improvements to the model within distance HALF*RHO of XOPT.
    ! Why do we reduce RHO when SHORTD is true and the entries of MODERRSAV and DNORMSAV are all
    ! small? The reason is well explained by the BOBYQA paper in the paragraphs surrounding
    ! (6.8)--(6.11). Roughly speaking, in this case, a trust-region step is unlikely to decrease the
    ! objective function according to some estimations. This suggests that the current trust-region
    ! center may be an approximate local minimizer up to the current "resolution" of the algorithm.
    ! When this occurs, the algorithm takes the view that the work for the current RHO is complete,
    ! and hence it will reduce RHO, which will enhance the resolution of the algorithm in general.
    if (shortd .or. .not. qred > 0) then
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if
        ! Evaluate EBOUND. It will be used as a bound to test if the entries of MODERRSAV are small.
        ebound = errbd(crvmin, d, gopt, hq, moderrsav, pq, rho, sl, su, xpt(:, kopt), xpt)
    else
        ! Calculate the next value of the objective function.
        x = xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su)  ! X = XBASE + XOPT + D without rounding.
        call evaluate(calfun, x, f)
        nf = nf + 1_IK

        ! Print a message about the function evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x)
        ! Save X, F into the history.
        call savehist(nf, x, xhist, f, fhist)

        ! Check whether to exit
        subinfo = checkexit(maxfun, nf, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        ! Update DNORMSAV and MODERRSAV.
        ! DNORMSAV contains the DNORM of the latest 3 function evaluations with the current RHO.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
        ! MODERR is the error of the current model in predicting the change in F due to D.
        ! MODERRSAV is the prediction errors of the latest 3 models with the current RHO.
        moderr = f - fval(kopt) + qred
        moderrsav = [moderrsav(2:size(moderrsav)), moderr]

        ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
        ratio = redrat(fval(kopt) - f, qred, eta1)

        ! Update DELTA. After this, DELTA < DNORM may hold.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if

        ! Is the newly generated X better than current best point?
        ximproved = (f < fval(kopt))

        ! Call RESCUE if rounding errors have damaged the denominator corresponding to D.
        ! RESCUE is invoked sometimes though not often after a trust-region step, and it does
        ! improve the performance, especially when pursing high-precision solutions..
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        den = calden(kopt, bmat, d, xpt, zmat)
        if (ximproved .and. .not. (is_finite(sum(abs(vlag))) .and. any(den > maxval(vlag(1:npt)**2)))) then
            ! Below are some alternatives conditions for calling RESCUE. They perform fairly well.
            ! !if (.false.) then  ! Do not call RESCUE at all.
            ! !if (ximproved .and. .not. any(den > 0.25_RP * maxval(vlag(1:npt)**2))) then
            ! !if (ximproved .and. .not. any(den > HALF * maxval(vlag(1:npt)**2))) then
            ! !if (.not. any(den > HALF * maxval(vlag(1:npt)**2))) then  ! Powell's code.
            ! !if (.not. any(den > maxval(vlag(1:npt)**2))) then
            call rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, &
                & fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, subinfo)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
            dnormsav = REALMAX
            moderrsav = REALMAX

            ! RESCUE shifts XBASE to the best point before RESCUE. Update D, MODERR, and XIMPROVED.
            ! Do NOT calculate QRED according to this D, as it is not really a trust region step.
            ! Note that QRED will be used afterward for defining IMPROVE_GEO and REDUCE_RHO.
            d = max(sl, min(su, d)) - xpt(:, kopt)
            moderr = f - fval(kopt) - quadinc(d, xpt, gopt, pq, hq)
            ximproved = (f < fval(kopt))
        end if

        ! Set KNEW_TR to the index of the interpolation point to be replaced with XOPT + D.
        ! KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
        knew_tr = setdrop_tr(kopt, ximproved, bmat, d, delta, rho, xpt, zmat)

        ! Update [BMAT, ZMAT] (representing H in the BOBYQA paper), [GQ, HQ, PQ] (the quadratic
        ! model), and [FVAL, XPT, KOPT, FOPT, XOPT] so that XPT(:, KNEW_TR) becomes XOPT + D. If
        ! KNEW_TR = 0, the updating subroutines will do essentially nothing, as the algorithm
        ! decides not to include XOPT + D into XPT.
        if (knew_tr > 0) then
            xdrop = xpt(:, knew_tr)
            xosav = xpt(:, kopt)
            call updateh(knew_tr, kopt, d, xpt, bmat, zmat)
            call updatexf(knew_tr, ximproved, f, max(sl, min(su, xosav + d)), kopt, fval, xpt)
            call updateq(knew_tr, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
            ! Try whether to replace the new quadratic model with the alternative model, namely the
            ! least Frobenius norm interpolant.
            call tryqalt(bmat, fval - fval(kopt), ratio, sl, su, xpt(:, kopt), xpt, zmat, itest, gopt, hq, pq)
        end if
    end if


    !----------------------------------------------------------------------------------------------!
    ! Before the next trust-region iteration, we may improve the geometry of XPT or reduce RHO
    ! according to IMPROVE_GEO and REDUCE_RHO, which in turn depend on the following indicators.
    ! N.B.: We must ensure that the algorithm does not set IMPROVE_GEO = TRUE at infinitely many
    ! consecutive iterations without moving XOPT or reducing RHO. Otherwise, the algorithm will get
    ! stuck in repetitive invocations of GEOSTEP. To this end, make sure the following.
    ! 1. The threshold for CLOSE_ITPSET is at least DELBAR, the trust region radius for GEOSTEP.
    ! Normally, DELBAR <= DELTA <= the threshold (In Powell's UOBYQA, DELBAR = RHO < the threshold).
    ! 2. If an iteration sets IMPROVE_GEO = TRUE, it must also reduce DELTA or set DELTA to RHO.

    ! ACCURATE_MOD: Are the recent models sufficiently accurate? Used only if SHORTD is TRUE.
    accurate_mod = all(abs(moderrsav) <= ebound) .and. all(dnormsav <= rho)
    ! CLOSE_ITPSET: Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
    close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))  ! Powell's code.
    ! Below are some alternative definitions of CLOSE_ITPSET.
    ! !close_itpset = all(distsq <= (TEN * delta)**2)  ! Does not work as well as Powell's version.
    ! !close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Powell's NEWUOA code.
    ! ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
    adequate_geo = (shortd .and. accurate_mod) .or. close_itpset
    ! SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
    small_trrad = (max(delta, dnorm) <= rho)  ! Powell's code. See also (6.7) of the BOBYQA paper.
    !small_trrad = (delsav <= rho)  ! Behaves the same as Powell's version. DELSAV = unupdated DELTA.

    ! IMPROVE_GEO and REDUCE_RHO are defined as follows.
    ! N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.
    ! Powell's code does not have (.NOT. QRED>0) in BAD_TRSTEP; it terminates if QRED > 0 fails.

    ! BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= eta1 .or. knew_tr == 0)
    improve_geo = bad_trstep .and. .not. adequate_geo  ! See the text above (6.7) of the BOBYQA paper.
    ! BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    reduce_rho = bad_trstep .and. adequate_geo .and. small_trrad  ! See (6.7) of the BOBYQA paper.
    ! Zaikun 20221111: What if RESCUE has been called? Is it still reasonable to use RATIO?
    ! Zaikun 20221127: If RESCUE has been called, then KNEW_TR may be 0 even if RATIO > 0.

    ! Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    ! !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

    ! With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= eta1 .or. knew_tr == 0)
    ! !improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

    ! With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    ! !reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

    ! BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    !call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO and REDUCE_RHO are not both TRUE', srname)
    !
    ! If SHORTD is TRUE or QRED > 0 is FALSE, then either IMPROVE_GEO or REDUCE_RHO is TRUE unless
    ! CLOSE_ITPSET is TRUE but SMALL_TRRAD is FALSE.
    !call assert((.not. shortd .and. qred > 0) .or. (improve_geo .or. reduce_rho .or. &
    !    & (close_itpset .and. .not. small_trrad)), 'If SHORTD is TRUE or QRED > 0 is FALSE, then either&
    !    & IMPROVE_GEO or REDUCE_RHO is TRUE unless CLOSE_ITPSET is TRUE but SMALL_TRRAD is FALSE', srname)
    !----------------------------------------------------------------------------------------------!


    ! Since IMPROVE_GEO and REDUCE_RHO are never TRUE simultaneously, the following two blocks are
    ! exchangeable: IF (IMPROVE_GEO) ... END IF and IF (REDUCE_RHO) ... END IF.

    ! Improve the geometry of the interpolation set by removing a point and adding a new one.
    if (improve_geo) then
        ! XPT(:, KNEW_GEO) will become XOPT + D below. KNEW_GEO /= KOPT unless there is a bug.
        knew_geo = int(maxloc(distsq, dim=1), kind(knew_geo))
        delbar = max(min(TENTH * sqrt(maxval(distsq)), delta), rho)

        ! Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
        d = geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat)

        ! Call RESCUE if rounding errors have damaged the denominator corresponding to D.
        ! 1. This does make a difference, yet RESCUE seems not invoked often after a geometry step.
        ! 2. In Powell's implementation, it may happen that RESCUE only recalculates [BMAT, ZMAT]
        ! without introducing any new point into XPT. In that case, GEOSTEP will have to be called
        ! after RESCUE, without which the code may encounter an infinite cycling. We have modified
        ! RESCUE so that it introduces at least one new point into XPT and there is no need to call
        ! GEOSTEP afterward. This improves the performance a bit and simplifies the flow of the code.
        ! 3. It is temping to incorporate XOPT+D into the interpolation even if RESCUE is called.
        ! However, this cannot be done without recalculating KNEW_GEO, as XPT has been changed by
        ! RESCUE, so that it is invalid to replace XPT(:, KNEW_GEO) with XOPT+D anymore. With a new
        ! KNEW_GEO, the step D will become improper as it was chosen according to the old KNEW_GEO.
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        den = calden(kopt, bmat, d, xpt, zmat)
        if (.not. (is_finite(sum(abs(vlag))) .and. den(knew_geo) > HALF * vlag(knew_geo)**2)) then
            call rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, &
                & fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, subinfo)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if

            dnormsav = REALMAX
            moderrsav = REALMAX
        else
            ! Calculate the next value of the objective function.
            x = xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su)  ! X = XBASE + XOPT + D without rounding.
            call evaluate(calfun, x, f)
            nf = nf + 1_IK

            ! Print a message about the function evaluation according to IPRINT.
            call fmsg(solver, iprint, nf, f, x)
            ! Save X, F into the history.
            call savehist(nf, x, xhist, f, fhist)

            ! Check whether to exit
            subinfo = checkexit(maxfun, nf, f, ftarget, x)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if

            ! Update DNORMSAV and MODERRSAV.
            ! DNORMSAV contains the DNORM of the latest 3 function evaluations with the current RHO.
            ! Powell's code does not update DNORM. Therefore, DNORM is the length of the last
            ! trust-region trial step, inconsistent with MODERRSAV. The same problem exists in NEWUOA.
            dnorm = min(delbar, norm(d))
            dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
            ! MODERR is the error of the current model in predicting the change in F due to D.
            ! MODERRSAV is the prediction errors of the latest 3 models with the current RHO.
            moderr = f - fval(kopt) - quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)
            moderrsav = [moderrsav(2:size(moderrsav)), moderr]

            ! Is the newly generated X better than current best point?
            ximproved = (f < fval(kopt))

            ! Update [BMAT, ZMAT] (represents H in the BOBYQA paper), [FVAL, XPT, KOPT, FOPT, XOPT],
            ! and [GQ, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XOPT + D.
            xdrop = xpt(:, knew_geo)
            xosav = xpt(:, kopt)
            call updateh(knew_geo, kopt, d, xpt, bmat, zmat)
            call updatexf(knew_geo, ximproved, f, max(sl, min(su, xosav + d)), kopt, fval, xpt)
            call updateq(knew_geo, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)
        end if
    end if

    ! The calculations with the current value of RHO are complete. Update RHO and DELTA.
    if (reduce_rho) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        delta = HALF * rho
        rho = redrho(rho, rhoend)
        delta = max(delta, rho)
        ! Print a message about the reduction of RHO according to IPRINT.
        call rhomsg(solver, iprint, nf, fval(kopt), rho, xbase + xpt(:, kopt))
        ! DNORMSAV and MODERRSAV are corresponding to the latest 3 function evaluations with
        ! the current RHO. Update them after reducing RHO.
        dnormsav = REALMAX
        moderrsav = REALMAX
    end if

    ! Shift XBASE if XOPT may be too far from XBASE.
    ! Powell's original criteria for shifting XBASE is as follows.
    ! 1. After a trust region step that is not short, shift XBASE if SUM(XOPT**2) >= 1.0E3*DNORM**2.
    ! In this case, it seems quite important for the performance to recalculate QRED.
    ! 2. Before a geometry step, shift XBASE if SUM(XOPT**2) >= 1.0E3*DELBAR**2.
    if (sum(xpt(:, kopt)**2) >= 1.0E3_RP * delta**2) then
        ! Other possible criteria: SUM(XOPT**2) >= 1.0E4*DELTA**2, SUM(XOPT**2) >= 1.0E4*RHO**2.
        sl = min(sl - xpt(:, kopt), ZERO)
        su = max(su - xpt(:, kopt), ZERO)
        call shiftbase(kopt, xbase, xpt, zmat, bmat, pq, hq)
        xbase = max(xl, min(xu, xbase))
    end if
end do

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xinbd(xbase, xpt(:, kopt) + d, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT + D.
    call evaluate(calfun, x, f)
    nf = nf + 1_IK
    ! Print a message about the function evaluation according to IPRINT.
    call fmsg(solver, iprint, nf, f, x)
    ! Save X, F into the history.
    call savehist(nf, x, xhist, f, fhist)
end if

! Choose the [X, F] to return: either the current [X, F] or [XBASE + XOPT, FOPT].
if (fval(kopt) < f .or. is_nan(f)) then
    x = xinbd(xbase, xpt(:, kopt), xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT.
    f = fval(kopt)
end if

! Arrange FHIST and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist)

! Print a return message according to IPRINT.
call retmsg(solver, info, iprint, nf, f, x)

close (16)
!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nf <= maxfun, 'NF <= MAXFUN', srname)
    call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
    call assert(all(x >= xl) .and. all(x <= xu), 'XL <= X <= XU', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    do k = 1, min(nf, maxxhist)
        call assert(all(xhist(:, k) >= xl) .and. all(xhist(:, k) <= xu), 'XL <= XHIST <= XU', srname)
    end do
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
        & 'FHIST does not contain NaN/+Inf', srname)
    call assert(.not. any(fhist(1:min(nf, maxfhist)) < f), 'F is the smallest in FHIST', srname)
end if

end subroutine bobyqb


function errbd(crvmin, d, gopt, hq, moderrsav, pq, rho, sl, su, xopt, xpt) result(ebound)
!--------------------------------------------------------------------------------------------------!
! This function defines EBOUND, which will be used as a bound to test whether the errors in recent
! models are sufficiently small. See the elaboration on pages 30--31 of the BOBYQA paper, in the
! paragraphs surrounding (6.8)--(6.11).
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, diag, issymmetric, trueloc
use, non_intrinsic :: powalg_mod, only : hess_mul

implicit none

! Inputs
real(RP), intent(in) :: crvmin
real(RP), intent(in) :: d(:)
real(RP), intent(in) :: gopt(:)
real(RP), intent(in) :: hq(:, :)
real(RP), intent(in) :: moderrsav(:)
real(RP), intent(in) :: pq(:)
real(RP), intent(in) :: rho
real(RP), intent(in) :: sl(:)
real(RP), intent(in) :: su(:)
real(RP), intent(in) :: xopt(:)
real(RP), intent(in) :: xpt(:, :)

! Outputs
real(RP) :: ebound

! Local variables
character(len=*), parameter :: srname = 'ERRBD'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bfirst(size(d))
real(RP) :: bsecond(size(d))
real(RP) :: gnew(size(d))
real(RP) :: xnew(size(d))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(crvmin >= 0, 'CRVMIN >= 0', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(rho > 0, 'RHO > 0', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xopt >= sl .and. xopt <= su), 'SL <= XOPT <= SU', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt) .and. &
        & xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
end if

!====================!
! Calculation starts !
!====================!

xnew = xopt + d
gnew = gopt + hess_mul(d, xpt, pq, hq)
bfirst = maxval(abs(moderrsav))
bfirst(trueloc(xnew <= sl)) = gnew(trueloc(xnew <= sl)) * rho
bfirst(trueloc(xnew >= su)) = -gnew(trueloc(xnew >= su)) * rho
bsecond = HALF * (diag(hq) + matprod(xpt**2, pq)) * rho**2
ebound = minval(max(bfirst, bfirst + bsecond))
if (crvmin > 0) then
    ebound = min(ebound, 0.125_RP * crvmin * rho**2)
end if

!====================!
!  Calculation ends  !
!====================!

end function errbd


end module bobyqb_mod
