module bobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! TODO: verify that the iterates/steps respect bounds in the pre/postconditions.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, November 21, 2022 PM10:12:35
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
! The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to the
! corresponding arguments in SUBROUTINE BOBYQA.
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
! of the model and Lagrange functions.
! XPT is a 2D array that holds the coordinates of the interpolation points relative to XBASE.
! FVAL holds the values of F at the interpolation points.
! XOPT is set to the displacement from XBASE of the trust region centre.
! GOPT holds the gradient of the quadratic model at XBASE + XOPT.
! HQ holds the explicit second derivatives of the quadratic model.
! PQ contains the parameters of the implicit second derivatives of the quadratic model.
! BMAT holds the last N columns of H.
! ZMAT holds the factorization of the leading NPT by NPT submatrix of H, this factorization being
! ZMAT * ZMAT^T, which provides both the correct rank and positive semi-definiteness.
! SL and SU hold XL - XBASE and XU - XBASE, respectively. All the components of every
! XOPT are going to satisfy the bounds SL(I) <= XOPT(I) <= SU(I), with appropriate equalities when
! XOPT is on a constraint boundary.
! D is chosen by subroutine TRSBOX or GEOSTEP. Usually XBASE + XOPT + D is the vector of variables
! for the next call of CALFUN.
! VLAG contains the values of the Lagrange functions at a new point X.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert!, wassert, validate
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : INFO_DFT, SMALL_TR_RADIUS!, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : matprod, diag, trueloc, r1update!, r2update!, norm
use, non_intrinsic :: output_mod, only : retmsg, rhomsg, fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calden, calvlag, calbeta, hess_mul!, errquad
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: shiftbase_mod, only : shiftbase

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep, setdrop_tr
use, non_intrinsic :: initialize_mod, only : initxf, initq, inith
use, non_intrinsic :: rescue_mod, only : rescue
use, non_intrinsic :: trustregion_mod, only : trsbox, trrad
use, non_intrinsic :: update_mod, only : updateh
use, non_intrinsic :: xinbd_mod, only : xinbd

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
integer(IK) :: subinfo
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
real(RP) :: bmat(size(x), npt + size(x))
real(RP) :: d(size(x))
real(RP) :: fval(npt)
real(RP) :: gopt(size(x))
real(RP) :: hq(size(x), size(x))
real(RP) :: pq(npt)
real(RP) :: vlag(npt + size(x))
real(RP) :: sl(size(x))
real(RP) :: su(size(x))
real(RP) :: xbase(size(x))
real(RP) :: xnew(size(x))
real(RP) :: xopt(size(x))
real(RP) :: xpt(size(x), npt)
real(RP) :: zmat(npt, npt - size(x) - 1)
real(RP) :: gnew(size(x))
real(RP) :: delbar, bdtest(size(x)), beta, &
&        crvmin, curv(size(x)), delta, &
&        den(npt), moderr, &
&        distsq(npt), dnorm, errbd, fopt,        &
&        gisq, gqsq,       &
&        ratio, rho, qred, pqinc(npt)
real(RP) :: dnormsav(3)
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pqalt(npt), galt(size(x)), fshift(npt), pgalt(size(x)), pgopt(size(x))
integer(IK) :: itest, knew_tr, knew_geo, kopt, nfresc
integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK)))
logical :: shortd, improve_geo, tr_success, reduce_rho, small_trrad, close_itpset, accurate_mod, adequate_geo, bad_trstep, rescued


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

info = INFO_DFT

! Initialize XBASE, XPT, FVAL, and KOPT.
call initxf(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x, ij, kopt, nf, fhist, fval, &
    & sl, su, xbase, xhist, xpt, subinfo)
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = xinbd(xbase, xopt, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT.
f = fopt

! Check whether to return due to abnormal cases that may occur during the initialization.
if (subinfo /= INFO_DFT) then
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

! Set some more initial values and parameters.
rho = rhobeg
delta = rho
errbd = ZERO
dnormsav = HUGENUM
moderrsav = HUGENUM
itest = 0

! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
ratio = -ONE
shortd = .false.
improve_geo = .false.

! Begin the iterative procedure.
! After solving a trust-region subproblem, we use three boolean variables to control the workflow.
! SHORTD: Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO: Should we improve the geometry?
! REDUCE_RHO: Should we reduce rho?
! BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do while (.true.)
    ! Generate the next trust region step D.
    call trsbox(delta, gopt, hq, pq, sl, su, xopt, xpt, crvmin, d)
    dnorm = min(delta, sqrt(sum(d**2)))

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
    ! (6.9)--(6.11). Roughly speaking, in this case, a trust-region step is unlikely to decrease
    ! the objective function according to some estimations. This suggests that the current
    ! trust-region center may be an approximate local minimizer. When this occurs, the algorithm
    ! takes the view that the work for the current RHO is complete, and hence it will reduce
    ! RHO, which will enhance the resolution of the algorithm in general.
    if (shortd .or. .not. qred > 0) then
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close.
        end if

        ! Define ERRBD, which will be used in the definition of REDUCE_RHO and IMPROVE_GEO.
        gnew = gopt + hess_mul(d, xpt, pq, hq)
        bdtest = maxval(abs(moderrsav))
        xnew = max(sl, min(su, xopt + d))  ! In precise arithmetic, XNEW = XOPT + D.
        bdtest(trueloc(xnew <= sl)) = gnew(trueloc(xnew <= sl)) * rho
        bdtest(trueloc(xnew >= su)) = -gnew(trueloc(xnew >= su)) * rho
        curv = diag(hq) + matprod(xpt**2, pq)
        errbd = minval(max(bdtest, bdtest + HALF * curv * rho**2))
        if (crvmin > 0) then
            errbd = min(errbd, 0.125_RP * crvmin * rho**2)
        end if
    else
        ! Zaikun 20220528: TODO: check the shifting strategy of NEWUOA and LINCOA.
        if (sum(xopt**2) >= 1.0E3_RP * dnorm**2) then
            sl = min(sl - xopt, ZERO)
            su = max(su - xopt, ZERO)
            call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq)
            ! SHIFTBASE shifts XBASE to XBASE + XOPT and XOPT to 0.
            xbase = max(xl, min(xu, xbase))
            ! It seems important for the performance to recalculate QRED.
            qred = -quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)
        end if

        ! Calculate the next value of the objective function.
        x = xinbd(xbase, xopt + d, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT + D.
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

        fopt = fval(kopt)

        ! Update DNORMSAV and MODERRSAV.
        ! DNORMSAV contains the DNORM of the latest 3 function evaluations with the current RHO.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
        ! MODERR is the error of the current model in predicting the change in F due to D.
        ! MODERRSAV is the prediction errors of the latest 3 models with the current RHO.
        moderr = f - fopt + qred
        moderrsav = [moderrsav(2:size(moderrsav)), moderr]

        ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
        ratio = redrat(fopt - f, qred, eta1)

        ! Update DELTA. After this, DELTA < DNORM may hold.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if

        tr_success = (f < fopt)

        ! Call RESCUE if if rounding errors have damaged the denominator corresponding to D.
        ! It provides a useful safeguard, but is not invoked in most applications of BOBYQA.
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        den = calden(kopt, bmat, d, xpt, zmat)
        if (tr_success .and. .not. (is_finite(sum(abs(vlag))) .and. any(den > maxval(vlag(1:npt)**2)))) then
            ! Below are some alternatives conditions for calling RESCUE. They perform fairly well.
            ! !if (.false.) then  ! Do not call RESCUE at all.
            ! !if (tr_success .and. .not. any(den > 0.25_RP * maxval(vlag(1:npt)**2))) then
            ! !if (tr_success .and. .not. any(den > HALF * maxval(vlag(1:npt)**2))) then
            ! !if (.not. any(den > HALF * maxval(vlag(1:npt)**2))) then  ! Powell's code.
            ! !if (.not. any(den > maxval(vlag(1:npt)**2))) then
            call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fopt, &
                & fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat, subinfo)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
            dnormsav = HUGENUM
            moderrsav = HUGENUM

            ! RESCUE shifts XBASE to XBASE + XOPT. Update D, QRED, MODERR, and TR_SUCCESS.
            d = max(sl, min(su, d)) - xopt
            qred = -quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)
            moderr = f - fopt + qred
            tr_success = (f < fopt)
        end if

        ! Set KNEW_TR to the index of the interpolation point to be replaced by XOPT + D.
        ! KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
        knew_tr = setdrop_tr(kopt, tr_success, bmat, d, delta, rho, xpt, zmat)

        if (knew_tr > 0) then
            ! Update [BMAT, ZMAT] (representing H in the BOBYQA paper), [GQ, HQ, PQ] (the quadratic
            ! model), and [FVAL, XPT, KOPT, FOPT, XOPT] so that XPT(:, KNEW_TR) becomes XOPT + D. If
            ! KNEW_TR = 0, the updating subroutines will do essentially nothing, as the algorithm
            ! decides not to include XOPT + D into XPT.

            call updateh(knew_tr, kopt, d, xpt, bmat, zmat)

            call r1update(hq, pq(knew_tr), xpt(:, knew_tr))
            pq(knew_tr) = ZERO
            pqinc = matprod(zmat, moderr * zmat(knew_tr, :))
            pq = pq + pqinc
            ! Alternatives:
            ! !PQ = PQ + MATPROD(ZMAT, MODERR * ZMAT(KNEW, :))
            ! !PQ = PQ + MODERR * MATPROD(ZMAT, ZMAT(KNEW, :))

            ! Include the new interpolation point, and make the changes to GOPT at the old XOPT that
            ! are caused by the updating of the quadratic model.
            fval(knew_tr) = f
            xpt(:, knew_tr) = max(sl, min(su, xopt + d))
            gopt = gopt + moderr * bmat(:, knew_tr) + hess_mul(xopt, xpt, pqinc)

            ! Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
            if (f < fopt) then
                kopt = knew_tr
                xopt = xpt(:, kopt)
                gopt = gopt + hess_mul(d, xpt, pq, hq)
            end if

            ! Calculate the parameters of the least Frobenius norm interpolant to the current data,
            ! the gradient of this interpolant at XOPT being put into VLAG(NPT+I), I=1,2,...,N.
            fshift = fval - fval(kopt)
            pqalt = matprod(zmat, matprod(fshift, zmat))
            galt = matprod(bmat(:, 1:npt), fshift) + hess_mul(xopt, xpt, pqalt)

            pgopt = gopt
            pgopt(trueloc(xopt >= su)) = max(ZERO, gopt(trueloc(xopt >= su)))
            pgopt(trueloc(xopt <= sl)) = min(ZERO, gopt(trueloc(xopt <= sl)))
            gqsq = sum(pgopt**2)

            pgalt = galt
            pgalt(trueloc(xopt >= su)) = max(ZERO, galt(trueloc(xopt >= su)))
            pgalt(trueloc(xopt <= sl)) = min(ZERO, galt(trueloc(xopt <= sl)))
            gisq = sum(pgalt**2)

            ! Test whether to replace the new quadratic model by the least Frobenius norm interpolant,
            ! making the replacement if the test is satisfied.
            ! N.B.:
            ! 1. The replacement is done only after a trust-region step, which differs from LINCOA.
            ! 2. The replacement is done regardless of DELTA <= RHO or not, which differs from NEWUOA.
            itest = itest + 1
            if (gqsq < TEN * gisq) itest = 0
            if (itest >= 3) then
                gopt = galt
                pq = pqalt
                hq = ZERO
                itest = 0
            end if
        end if

    end if


    !----------------------------------------------------------------------------------------------!
    ! Before the next trust-region iteration, we may improve the geometry of XPT or reduce RHO
    ! according to IMPROVE_GEO and REDUCE_RHO, which in turn depend on the following indicators.
    ! ACCURATE_MOD: Are the recent models sufficiently accurate? Used only if SHORTD is TRUE.

    accurate_mod = all(abs(moderrsav) <= errbd) .and. all(dnormsav <= rho)
    ! CLOSE_ITPSET: Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xopt).^2)  % xopt should be a column! Implicit expansion
    close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))  ! Powell's code.
    ! Below are some alternative definitions of CLOSE_ITPSET.
    ! !close_itpset = all(distsq <= (TEN * rho)**2)  ! Works almost the same as Powell's version.
    ! !close_itpset = all(distsq <= (TEN * delta)**2)  ! Does not work as well as Powell's version.
    ! !close_itpset = all(distsq <= 4.0_RP * rho**2)  ! Does not work as well as Powell's version.
    ! !close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Powell's NEWUOA code.
    ! ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
    adequate_geo = (shortd .and. accurate_mod) .or. close_itpset
    ! SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
    small_trrad = (max(delta, dnorm) <= rho)  ! Powell's code.
    !small_trrad = (delsav <= rho)  ! Behaves the same as Powell's version. DELSAV = unupdated DELTA.

    ! IMPROVE_GEO and REDUCE_RHO are defined as follows.
    ! Powell's code does not have (.NOT. QRED>0) in BAD_TRSTEP; it terminates if QRED > 0 fails.
    ! BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= TENTH .or. knew_tr == 0)
    improve_geo = bad_trstep .and. .not. adequate_geo
    ! BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    reduce_rho = bad_trstep .and. adequate_geo .and. small_trrad
    ! Zaikun 20221111: What if RESCUE has been called? Is it still reasonable to use RATIO?

    ! Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    ! !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

    ! With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= TENTH .or. knew_tr == 0)
    ! !improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

    ! With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    ! !reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

    ! BOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    !call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO or REDUCE_RHO is false', srname)
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
        knew_geo = int(maxloc(distsq, dim=1), IK)
        delbar = max(min(TENTH * sqrt(maxval(distsq)), delta), rho)

        ! Zaikun 20220528: TODO: check the shifting strategy of NEWUOA and LINCOA.
        if (sum(xopt**2) >= 1.0E3_RP * delbar**2) then
            sl = min(sl - xopt, ZERO)
            su = max(su - xopt, ZERO)
            call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq)
            ! SHIFTBASE shifts XBASE to XBASE + XOPT and XOPT to 0.
            xbase = max(xl, min(xu, xbase))
        end if

        ! Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
        d = geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat)

        ! Call RESCUE if if rounding errors have damaged the denominator corresponding to D.
        ! It provides a useful safeguard, but is not invoked in most applications of BOBYQA.
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        den = calden(kopt, bmat, d, xpt, zmat)
        rescued = .false.
        if (.not. (is_finite(sum(abs(vlag))) .and. den(knew_geo) > HALF * vlag(knew_geo)**2)) then
            ! The condition below works the same as the above one, which is used in Powell's code.
            ! !if (.not. (is_finite(sum(abs(vlag))) .and. den(knew_geo) > 0.25 * vlag(knew_geo)**2)) then
            nfresc = nf
            call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fopt, &
                & fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat, subinfo)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if

            dnormsav = HUGENUM
            moderrsav = HUGENUM

            ! If NFRESC < NF, then RESCUE has updated XPT (also the model and [BMAT, ZMAT]) to
            ! improve (rescue) its geometry. Thus no geometry step is needed anymore. If NFRESC
            ! equals NF, then RESCUE did not make any change to XPT, but only recalculated the
            ! model ans [BMAT, ZMAT]; in this case, we calculate a new geometry step.
            if (nfresc == nf) then
                d = geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat)
            else
                rescued = .true.
            end if
        end if

        if (.not. rescued) then
            ! Calculate the next value of the objective function.
            x = xinbd(xbase, xopt + d, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT + D.
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

            fopt = fval(kopt)

            ! Update DNORMSAV and MODERRSAV.
            ! DNORMSAV contains the DNORM of the latest 3 function evaluations with the current RHO.
            dnorm = min(delbar, sqrt(sum(d**2)))
            dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
            ! MODERR is the error of the current model in predicting the change in F due to D.
            ! MODERRSAV is the prediction errors of the latest 3 models with the current RHO.
            moderr = f - fopt - quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)
            moderrsav = [moderrsav(2:size(moderrsav)), moderr]
            !------------------------------------------------------------------------------------------!
            ! Zaikun 20220912: Powell's code does not update DNORM. Therefore, DNORM is the length
            ! of the last  trust-region trial step, which is inconsistent with MODERRSAV. The same
            ! problem exists in NEWUOA.
            !------------------------------------------------------------------------------------------!

            ! Update [BMAT, ZMAT] (represents H in the BOBYQA paper), [FVAL, XPT, KOPT, FOPT, XOPT],
            ! and [GQ, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XOPT + D.
            !vlag = calvlag(kopt, bmat, d, xpt, zmat)
            !beta = calbeta(kopt, bmat, d, xpt, zmat)
            call updateh(knew_geo, kopt, d, xpt, bmat, zmat)

            call r1update(hq, pq(knew_geo), xpt(:, knew_geo))
            pq(knew_geo) = ZERO
            pqinc = matprod(zmat, moderr * zmat(knew_geo, :))
            pq = pq + pqinc
            ! Alternatives:
            ! !PQ = PQ + MATPROD(ZMAT, MODERR * ZMAT(KNEW, :))
            ! !PQ = PQ + MODERR * MATPROD(ZMAT, ZMAT(KNEW, :))

            ! Include the new interpolation point, and make the changes to GOPT at the old XOPT that are
            ! caused by the updating of the quadratic model.
            fval(knew_geo) = f
            xpt(:, knew_geo) = max(sl, min(su, xopt + d))

            gopt = gopt + moderr * bmat(:, knew_geo) + hess_mul(xopt, xpt, pqinc)

            ! Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
            if (f < fopt) then
                kopt = knew_geo
                xopt = xpt(:, kopt)
                gopt = gopt + hess_mul(d, xpt, pq, hq)
            end if
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
        call rhomsg(solver, iprint, nf, fopt, rho, xbase + xopt)
        ! DNORMSAV and MODERRSAV are corresponding to the latest 3 function evaluations with
        ! the current RHO. Update them after reducing RHO.
        dnormsav = HUGENUM
        moderrsav = HUGENUM
    end if
end do

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xinbd(xbase, xopt + d, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT + D.
    call evaluate(calfun, x, f)
    nf = nf + 1_IK
    ! Print a message about the function evaluation according to IPRINT.
    call fmsg(solver, iprint, nf, f, x)
    ! Save X, F into the history.
    call savehist(nf, x, xhist, f, fhist)
end if

! Choose the [X, F] to return: either the current [X, F] or [XBASE + XOPT, FOPT].
if (fval(kopt) <= f .or. is_nan(f)) then
    x = xinbd(xbase, xopt, xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XOPT.
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


end subroutine bobyqb


end module bobyqb_mod
