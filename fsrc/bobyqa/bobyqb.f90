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
! Last Modified: Friday, September 23, 2022 AM09:38:09
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
! GOPT holds the gradient of the quadratic model at XBASE+XOPT.
! HQ holds the explicit second derivatives of the quadratic model.
! PQ contains the parameters of the implicit second derivatives of the quadratic model.
! BMAT holds the last N columns of H.
! ZMAT holds the factorization of the leading NPT by NPT submatrix of H, this factorization being
! ZMAT * ZMAT^T, which provides both the correct rank and positive semi-definiteness.
! SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. All the components of every
! XOPT are going to satisfy the bounds SL(I) <= XOPT(I) <= SU(I), with appropriate equalities when
! XOPT is on a constraint boundary.
! XNEW is chosen by subroutine TRSBOX or GEOSTEP. Usually XBASE+XNEW is the vector of variables for
! the next call of CALFUN. XNEW also satisfies the SL and SU constraints in the way above-mentioned.
! XALT is an alternative to XNEW, chosen by GEOSTEP, that may replace XNEW in order to increase the
! denominator in the updating of UPDATE.
! D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
! VLAG contains the values of the Lagrange functions at a new point X.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, INFO_DFT, &
    & MAXFUN_REACHED, DAMAGING_ROUNDING, TRSUBP_FAILED, SMALL_TR_RADIUS!, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : matprod, diag, trueloc, r1update!, r2update!, norm
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calvlag, calbeta, hess_mul
!!! TODO (Zaikun 20220525): Use CALDEN instead of CALVLAG and CALBETA wherever possible!!!

use, non_intrinsic :: ieee_4dev_mod, only : ieeenan

! Solver-specific modules
use, non_intrinsic :: initialize_mod, only : initxf, initq, inith
use, non_intrinsic :: geometry_mod, only : geostep
use, non_intrinsic :: rescue_mod, only : rescue
use, non_intrinsic :: trustregion_mod, only : trsbox
use, non_intrinsic :: update_mod, only : updateh
use, non_intrinsic :: shiftbase_mod, only : shiftbase

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
real(RP) :: delbar, alpha, bdtest(size(x)), beta, &
&        biglsq, crvmin, curv(size(x)), delta,  &
&        den(npt), denom, densav, diff, &
&        dist, dsquare, distsq(npt), dnorm, dsq, errbd, fopt,        &
&        gisq, gqsq, hdiag(npt),      &
&        ratio, rho, rhosq, scaden, qred, weight(npt), pqinc(npt)
real(RP) :: dnormsav(3)
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pqalt(npt), galt(size(x)), fshift(npt), pgalt(size(x)), pgopt(size(x))
real(RP) :: score(npt), wlagsq(npt)
integer(IK) :: itest, knew, kopt, ksav, nfresc
integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK)))
logical :: shortd, improve_geo, rescue_geo


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

denom = ieeenan()  ! To be removed.

info = INFO_DFT

! Initialize XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT and ZMAT together with the corresponding values of
! of NF and KOPT, which are the number of calls of CALFUN so far and the index of the interpolation
! point at the trust region centre. Then the initial XOPT is set too. The branch to label 720 occurs
! if MAXFUN is less than NPT. GOPT will be updated if KOPT is different from KBASE.
call initxf(calfun, iprint, maxfun, ftarget, rhobeg, xl, xu, x, ij, kopt, nf, fhist, fval, &
    & sl, su, xbase, xhist, xpt, subinfo)
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = min(max(xl, xbase + xopt), xu)
x(trueloc(xopt <= sl)) = xl(trueloc(xopt <= sl))
x(trueloc(xopt >= su)) = xu(trueloc(xopt >= su))
f = fopt
if (subinfo /= INFO_DFT) then
    call rangehist(nf, xhist, fhist)
    return
end if
!write (16, *) nf, kopt, fopt

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
moderrsav = HUGENUM
dnormsav = HUGENUM
nfresc = nf
itest = 0

! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
ratio = -ONE
shortd = .false.
improve_geo = .false.
rescue_geo = .false.

do while (.true.)
    ! Generate the next point in the trust region that provides a small value of the quadratic model
    ! subject to the constraints on the variables.

    rescue_geo = .false.

    call trsbox(delta, gopt, hq, pq, sl, su, xopt, xpt, crvmin, d)

    xnew = max(min(xopt + d, su), sl)  ! In precise arithmetic, XNEW = XOPT + D.
    gnew = gopt + hess_mul(d, xpt, pq, hq)

    dsq = sum(d**2)
    dnorm = min(delta, sqrt(dsq))
    shortd = (dnorm < HALF * rho)

    ! When D is short, make a choice between reducing RHO and improving the geometry depending
    ! on whether or not our work with the current RHO seems complete. RHO is reduced if the
    ! errors in the quadratic model at the last three interpolation points compare favourably
    ! with predictions of likely improvements to the model within distance HALF*RHO of XOPT.
    ! The BOBYQA paper explains the strategy in the paragraphs between (6.7) and (6.11).
    ! Why do we reduce RHO when SHORTD is true and the entries of MODERRSAV and DNORMSAV are all
    ! small? The reason is well explained by the BOBYQA paper in the paragraphs surrounding
    ! (6.9)--(6.11). Roughly speaking, in this case, a trust-region step is unlikely to decrease
    ! the objective function according to some estimations. This suggests that the current
    ! trust-region center may be an approximate local minimizer. When this occurs, the algorithm
    ! takes the view that the work for the current RHO is complete, and hence it will reduce
    ! RHO, which will enhance the resolution of the algorithm in general.
    if (shortd) then  ! D is to short to invoke a function evaluation.
        ! Reduce DELTA.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close.
        end if

        rhosq = rho**2
        dsquare = 1.0E2_RP * rhosq
        bdtest = maxval(abs(moderrsav))
        bdtest(trueloc(xnew <= sl)) = gnew(trueloc(xnew <= sl)) * rho
        bdtest(trueloc(xnew >= su)) = -gnew(trueloc(xnew >= su)) * rho
        curv = diag(hq) + matprod(xpt**2, pq)
        errbd = minval(max(bdtest, bdtest + HALF * curv * rhosq))
        if (crvmin > 0) then
            errbd = min(errbd, 0.125_RP * crvmin * rhosq)
        end if
        improve_geo = (any(abs(moderrsav) > errbd) .or. any(dnormsav > rho))
    else  ! D is long enough to invoke a function evaluation.

        ! Zaikun 20220528: TODO: check the shifting strategy of NEWUOA and LINCOA.
        if (sum(xopt**2) >= 1.0E3_RP * dsq) then
            sl = min(sl - xopt, ZERO)
            su = max(su - xopt, ZERO)
            xnew = min(max(sl, xnew - xopt), su)  ! Needed? Will XNEW be used again later?
            call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq)
            xbase = min(max(xl, xbase), xu)
        end if

        !! Calculate VLAG and BETA for the current choice of D.
        !hdiag = sum(zmat**2, dim=2)
        !vlag = calvlag(kopt, bmat, d, xpt, zmat)
        !beta = calbeta(kopt, bmat, d, xpt, zmat)
        !den = hdiag * beta + vlag(1:npt)**2
        !distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
        !weight = max(ONE, (distsq / delta**2)**2)  ! It differs from (6.1) in the BOBYQA paper.

        !! Find the index of the interpolation point to be replaced by the trust-region trial point.
        !! The strategy here differs from that of NEWUOA: KNEW is decided before the function
        !! evaluation, and it will be recalculated after the function evaluation if the new function
        !! value is lower than the current FOPT; in NEWUOA, KNEW is calculated only after function
        !! evaluation, taking the new function value into consideration (e.g., when evaluating
        !! DISTSQ). The purpose of finding KNEW before the function evaluation is to decide whether
        !! to call RESCUE. Is there a better strategy? What about follow the scheme of NEWUOA, and
        !! call RESCUE if DENOM turns out too small (not more than SQRT(EPS)*max(ONE, VLAG(KNEW)))?
        !! Note that the condition for calling RESCUE must be very stringent, so that the newly
        !! evaluated function value is included in the next interpolation system as much as possible.
        !score = weight * den
        !score(kopt) = -ONE  ! Skip KOPT when taking the maximum of SCORE
        !knew = 0
        !scaden = ZERO
        !if (any(score > 0)) then
        !    ! SCORE(K) = NaN implies DEN(K) = NaN. We exclude such K as we want DEN to be big.
        !    knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
        !    scaden = score(knew)
        !    !!MATLAB: [scaden, knew] = max(score, [], 'omitnan');
        !    denom = den(knew)
        !end if

        !wlagsq = weight * vlag(1:npt)**2
        !wlagsq(kopt) = -ONE  ! Skip KOPT when taking the maximum of WLAGSQ
        !biglsq = ZERO
        !if (any(wlagsq > 0)) then
        !    biglsq = maxval(wlagsq, mask=(.not. is_nan(wlagsq)))
        !        !!MATLAB: biglsq = max(wlagsq, [], 'omitnan');
        !end if

        !!rescue_geo = .not. (scaden > HALF * biglsq)  ! This is the normal condition.
        !rescue_geo = .not. (scaden > biglsq)  ! This is used when verifying RESCUE.
        !improve_geo = .not. rescue_geo
        !if (.not. rescue_geo) then
        ! Put the variables for the next calculation of the objective function in XNEW, with any
        ! adjustments for the bounds. In precise arithmetic, X = XBASE + XNEW.
        x = min(max(xl, xbase + xnew), xu)
        x(trueloc(xnew <= sl)) = xl(trueloc(xnew <= sl))
        x(trueloc(xnew >= su)) = xu(trueloc(xnew >= su))
        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if
        nf = nf + 1
        if (is_nan(abs(sum(x)))) then
            f = sum(x)  ! Set F to NaN
            if (nf == 1) then
                fopt = f
                xopt = ZERO
            end if
            info = NAN_INF_X
            exit
        end if

        ! Calculate the value of the objective function at XBASE+XNEW.
        call evaluate(calfun, x, f)
        call savehist(nf, x, xhist, f, fhist)
        !write (16, *) 'tr ', nf, f, fopt, kopt

        if (is_nan(f) .or. is_posinf(f)) then
            if (nf == 1) then
                fopt = f
                xopt = ZERO
            end if
            info = NAN_INF_F
            exit
        end if
        if (f <= ftarget) then
            info = FTARGET_ACHIEVED
            exit
        end if

        ! Use the quadratic model to predict the change in F due to the step D, and set DIFF to the
        ! error of this prediction.
        fopt = fval(kopt)
        qred = -quadinc(d, xpt, gopt, pq, hq)
        diff = f - fopt + qred
        moderrsav = [moderrsav(2:size(moderrsav)), f - fopt + qred]
        ! Zaikun 20220912: If the current D is a geometry step, then DNORM is not updated. It is
        ! still the value corresponding to last trust-region step. It seems inconsistent with (6.8)
        ! of the BOBYQA paper and the elaboration below it. Is this a bug? Similar thing happened
        ! in NEWUOA, but we recognized it as a bug and fixed it.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

        ! Pick the next value of DELTA after a trust region step.
        if (.not. (qred > ZERO)) then
            !----------------------------------------------------------------------------------!
            ! Zaikun 20220405: LINCOA improves the model in this case. Try the same here?
            !----------------------------------------------------------------------------------!
            info = TRSUBP_FAILED
            exit
        end if
        ratio = (fopt - f) / qred
        if (ratio <= TENTH) then
            delta = min(HALF * delta, dnorm)
        else if (ratio <= 0.7_RP) then
            delta = max(HALF * delta, dnorm)
        else
            delta = max(HALF * delta, dnorm + dnorm)
        end if
        if (delta <= 1.5_RP * rho) delta = rho
        ! Zaikun 20220720: On the top of page 29 of the BOBYQA paper, Powell wrote: If the k-th
        ! iteration is of "alternative" type, then the (k+1)-th iteration always calculates a
        ! "trust region" step with Delta_{k+1} = Delta_k and rho_{k+1} = rho_k. This is true for
        ! rho, but isn't it a typo for Delta? NO. It simply means that the algorithm does not
        ! update Delta after taking a geometry step. In the BOBYQA paper, the iteration counter
        ! k is increase by 1 both after a trust region step and after a geometry step.

        ! Recalculate KNEW and DENOM if the new F is less than FOPT.
        if (f < fopt) then
            distsq = sum((xpt - spread(xnew, dim=2, ncopies=npt))**2, dim=1)
        else
            distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
        end if
        !ksav = knew
        !densav = denom
        hdiag = sum(zmat**2, dim=2)
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        beta = calbeta(kopt, bmat, d, xpt, zmat)
        den = hdiag * beta + vlag(1:npt)**2
        !distsq = sum((xpt - spread(xnew, dim=2, ncopies=npt))**2, dim=1)
        weight = max(ONE, (distsq / delta**2)**2)
        score = weight * den

        if (.not. (f < fopt)) then
            score(kopt) = -ONE
        end if

        knew = 0
        scaden = ZERO
        if (any(score > 0)) then
            ! SCORE(K) = NaN implies DEN(K) = NaN. We exclude such K as we want DEN to be big.
            knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
            scaden = score(knew)
                    !!MATLAB: [scaden, knew] = max(score, [], 'omitnan');
            denom = den(knew)
        end if
        !write (16, *) 'knew 1', knew

        if (knew == 0 .and. f < fopt) then
            knew = maxloc(distsq, dim=1)
        end if
        !write (16, *) 'knew 2', knew

        wlagsq = weight * vlag(1:npt)**2
        biglsq = ZERO
        if (any(wlagsq > 0)) then
            biglsq = maxval(wlagsq, mask=(.not. is_nan(wlagsq)))
                    !!MATLAB: biglsq = max(wlagsq, [], 'omitnan');
        end if

                !! KNEW > 0 is implied by SCADEN > HALF*BIGLSQ (but NOT SCADEN >= ...), yet prefer to
                !! require KNEW > 0 explicitly.
        !if (.not. (knew > 0 .and. scaden > HALF * biglsq)) then
        !    knew = ksav
        !    denom = densav
        !end if
        !end if

        if (knew > 0) then
            ! Update BMAT and ZMAT, so that the KNEW-th interpolation point can be moved. Also update
            ! the second derivative terms of the model.
            !------------------------------------------------------------------------------------------!
            call assert(.not. any(abs(vlag - calvlag(kopt, bmat, d, xpt, zmat)) > 0), 'VLAG == VLAG_TEST', srname)
            call assert(.not. abs(beta - calbeta(kopt, bmat, d, xpt, zmat)) > 0, 'BETA == BETA_TEST', srname)
            call assert(.not. abs(denom - (sum(zmat(knew, :)**2) * beta + vlag(knew)**2)) > 0, 'DENOM = DENOM_TEST', srname)
            !--------------------------------------------------------------------------------------------------!
            call updateh(knew, beta, vlag, bmat, zmat)

            call r1update(hq, pq(knew), xpt(:, knew))
            pq(knew) = ZERO
            pqinc = matprod(zmat, diff * zmat(knew, :))
            pq = pq + pqinc
            ! Alternatives:
            !!PQ = PQ + MATPROD(ZMAT, DIFF * ZMAT(KNEW, :))
            !!PQ = PQ + DIFF * MATPROD(ZMAT, ZMAT(KNEW, :))

            ! Include the new interpolation point, and make the changes to GOPT at the old XOPT that are
            ! caused by the updating of the quadratic model.
            fval(knew) = f
            xpt(:, knew) = xnew
            gopt = gopt + diff * bmat(:, knew) + hess_mul(xopt, xpt, pqinc)

            ! Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
            if (f < fopt) then
                kopt = knew
                xopt = xnew
                gopt = gopt + hess_mul(d, xpt, pq, hq)
            end if

            ! Calculate the parameters of the least Frobenius norm interpolant to the current data, the
            ! gradient of this interpolant at XOPT being put into VLAG(NPT+I), I=1,2,...,N.
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
            itest = itest + 1
            if (gqsq < TEN * gisq) itest = 0
            if (itest >= 3) then
                gopt = galt
                pq = pqalt
                hq = ZERO
                itest = 0
            end if
        end if

        ! If a trust region step has provided a sufficient decrease in F, then branch for another
        ! trust region calculation.
        !if (f <= fopt - TENTH * qred) cycle
        improve_geo = .not. (knew > 0 .and. f <= fopt - TENTH * qred)
        if (.not. improve_geo) cycle

        ! Alternatively, find out if the interpolation points are close enough to the best point so far.
        dsquare = max((TWO * delta)**2, (TEN * rho)**2)
        !end if
    end if

    if (improve_geo) then
        distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
        knew = int(maxloc([dsquare, distsq], dim=1), IK) - 1_IK ! This line cannot be exchanged with the next
        dsquare = maxval([dsquare, distsq]) ! This line cannot be exchanged with the last

        ! If KNEW is positive, then GEOSTEP finds alternative new positions for the KNEW-th
        ! interpolation point within distance DELBAR of XOPT. Otherwise, go for another trust region
        ! iteration, unless the calculations with the current RHO are complete.
        if (knew > 0) then
            dist = sqrt(dsquare)
            delbar = max(min(TENTH * dist, delta), rho)
            dsq = delbar * delbar

            ! Zaikun 20220528: TODO: check the shifting strategy of NEWUOA and LINCOA.
            if (sum(xopt**2) >= 1.0E3_RP * dsq) then
                sl = min(sl - xopt, ZERO)
                su = max(su - xopt, ZERO)
                xnew = min(max(sl, xnew - xopt), su)  ! Needed? Will XNEW be used again later?
                call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq)
                xbase = min(max(xl, xbase), xu)
            end if

            ! Calculate a geometry step.
            d = geostep(knew, kopt, bmat, delbar, sl, su, xpt, zmat)
            xnew = min(max(sl, xopt + d), su)

            ! Calculate VLAG, BETA, and DENOM for the current choice of D.
            alpha = sum(zmat(knew, :)**2)
            vlag = calvlag(kopt, bmat, d, xpt, zmat)
            beta = calbeta(kopt, bmat, d, xpt, zmat)
            denom = alpha * beta + vlag(knew)**2

            rescue_geo = .not. (denom > HALF * vlag(knew)**2) ! This is the normal condition.
            !rescue_geo = .not. (denom > vlag(knew)**2)  ! This is used when verifying RESCUE.
            if (.not. rescue_geo) then
                ! Put the variables for the next calculation of the objective function in XNEW, with any
                ! adjustments for the bounds. In precise arithmetic, X = XBASE + XNEW.
                x = min(max(xl, xbase + xnew), xu)
                x(trueloc(xnew <= sl)) = xl(trueloc(xnew <= sl))
                x(trueloc(xnew >= su)) = xu(trueloc(xnew >= su))
                if (nf >= maxfun) then
                    info = MAXFUN_REACHED
                    exit
                end if
                nf = nf + 1
                if (is_nan(abs(sum(x)))) then
                    f = sum(x)  ! Set F to NaN
                    if (nf == 1) then
                        fopt = f
                        xopt = ZERO
                    end if
                    info = NAN_INF_X
                    exit
                end if

                ! Calculate the value of the objective function at XBASE+XNEW.
                call evaluate(calfun, x, f)
                call savehist(nf, x, xhist, f, fhist)
                !write (16, *) 'geo', nf, f, fopt, kopt

                if (is_nan(f) .or. is_posinf(f)) then
                    if (nf == 1) then
                        fopt = f
                        xopt = ZERO
                    end if
                    info = NAN_INF_F
                    exit
                end if
                if (f <= ftarget) then
                    info = FTARGET_ACHIEVED
                    exit
                end if

                ! Use the quadratic model to predict the change in F due to the step D, and set DIFF to the
                ! error of this prediction.
                fopt = fval(kopt)
                qred = -quadinc(d, xpt, gopt, pq, hq)
                diff = f - fopt + qred
                moderrsav = [moderrsav(2:size(moderrsav)), f - fopt + qred]
                ! Zaikun 20220912: If the current D is a geometry step, then DNORM is not updated. It is
                ! still the value corresponding to last trust-region step. It seems inconsistent with (6.8)
                ! of the BOBYQA paper and the elaboration below it. Is this a bug? Similar thing happened
                ! in NEWUOA, but we recognized it as a bug and fixed it.
                dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

                ! Update BMAT and ZMAT, so that the KNEW-th interpolation point can be moved. Also update
                ! the second derivative terms of the model.
                !------------------------------------------------------------------------------------------!
                call assert(.not. any(abs(vlag - calvlag(kopt, bmat, d, xpt, zmat)) > 0), 'VLAG == VLAG_TEST', srname)
                call assert(.not. abs(beta - calbeta(kopt, bmat, d, xpt, zmat)) > 0, 'BETA == BETA_TEST', srname)
                call assert(.not. abs(denom - (sum(zmat(knew, :)**2) * beta + vlag(knew)**2)) > 0, 'DENOM = DENOM_TEST', srname)
                !--------------------------------------------------------------------------------------------------!
                call updateh(knew, beta, vlag, bmat, zmat)

                call r1update(hq, pq(knew), xpt(:, knew))
                pq(knew) = ZERO
                pqinc = matprod(zmat, diff * zmat(knew, :))
                pq = pq + pqinc
                ! Alternatives:
                !!PQ = PQ + MATPROD(ZMAT, DIFF * ZMAT(KNEW, :))
                !!PQ = PQ + DIFF * MATPROD(ZMAT, ZMAT(KNEW, :))

                ! Include the new interpolation point, and make the changes to GOPT at the old XOPT that are
                ! caused by the updating of the quadratic model.
                fval(knew) = f
                xpt(:, knew) = xnew
                gopt = gopt + diff * bmat(:, knew) + hess_mul(xopt, xpt, pqinc)

                ! Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
                if (f < fopt) then
                    kopt = knew
                    xopt = xnew
                    gopt = gopt + hess_mul(d, xpt, pq, hq)
                end if
                cycle
            end if
        else if ((.not. shortd) .and. (ratio > 0 .or. max(delta, dnorm) > rho)) then
            cycle
        end if
    end if
    !write (16, *) nf, kopt, fopt

! Call RESCUE if if rounding errors have damaged the denominator corresponding to D.
! 1. RESCUE is called only if rounding errors have reduced by at least a factor of TWO the
! denominator of the formula for updating the H matrix. It provides a useful safeguard, but is
! not invoked in most applications of BOBYQA.
! 2. In Powell's code, if RESCUE is called after GEO_STEP, but RESCUE invokes no function
! evaluation (so that XPT is not updated, but the model and [BMAT, ZMAT] are recalculated),
! then GEO_STEP is called again immediately. Here, however, we always call TRSTEP after RESCUE.
    if (rescue_geo) then
        if (nf <= nfresc) then
            info = DAMAGING_ROUNDING
            exit
        end if
        call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fopt, &
            & fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat, subinfo)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
        nfresc = nf
        moderrsav = HUGENUM
        dnormsav = HUGENUM
        cycle
    end if

! The calculations with the current value of RHO are complete. Update RHO and DELTA.
    if (rho <= rhoend) then
        info = SMALL_TR_RADIUS
        exit
    end if
    delta = HALF * rho
    ratio = rho / rhoend
    if (ratio <= 16.0_RP) then
        rho = rhoend
    else if (ratio <= 250.0_RP) then
        rho = sqrt(ratio) * rhoend
    else
        rho = TENTH * rho
    end if
    delta = max(delta, rho)
    moderrsav = HUGENUM
    dnormsav = HUGENUM
end do

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = min(max(xl, xbase + xnew), xu)  ! XNEW = XOPT + D??? See NEWUOA, LINCOA.
    x(trueloc(xnew <= sl)) = xl(trueloc(xnew <= sl))
    x(trueloc(xnew >= su)) = xu(trueloc(xnew >= su))
    nf = nf + 1_IK
    call evaluate(calfun, x, f)
    call savehist(nf, x, xhist, f, fhist)
end if

!--------------------------------------------------------------------------------------------------!
!  Powell: IF (FVAL(KOPT) .LE. FSAVE) THEN
!  Why update X only when FVAL(KOPT) .LE. FSAVE? This seems INCORRECT, because it may lead to
!  a return with F and X that are not the best available.
!--------------------------------------------------------------------------------------------------!
if (fval(kopt) <= f .or. is_nan(f)) then
    x = min(max(xl, xbase + xopt), xu)
    x(trueloc(xopt <= sl)) = xl(trueloc(xopt <= sl))
    x(trueloc(xopt >= su)) = xu(trueloc(xopt >= su))
    f = fval(kopt)
end if

call rangehist(nf, xhist, fhist)

!====================!
!  Calculation ends  !
!====================!

! Postconditions

close (16)

end subroutine bobyqb


end module bobyqb_mod
