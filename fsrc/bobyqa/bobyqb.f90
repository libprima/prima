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
! Last Modified: Wednesday, November 02, 2022 PM04:54:22
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
use, non_intrinsic :: debug_mod, only : assert!, wassert, validate
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, INFO_DFT, &
    & MAXFUN_REACHED, DAMAGING_ROUNDING, TRSUBP_FAILED, SMALL_TR_RADIUS!, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : matprod, diag, trueloc, r1update!, r2update!, norm
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calden, calvlag, calbeta, hess_mul!, errquad

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
real(RP) :: delbar, bdtest(size(x)), beta, &
&        crvmin, curv(size(x)), delta,  &
&        den(npt), diff, &
&        dist, dsquare, distsq(npt), dnorm, dsq, errbd, fopt,        &
&        gisq, gqsq,       &
&        ratio, rho, qred, weight(npt), pqinc(npt)
real(RP) :: dnormsav(3)
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pqalt(npt), galt(size(x)), fshift(npt), pgalt(size(x)), pgopt(size(x))
real(RP) :: score(npt)
integer(IK) :: itest, knew_tr, knew_geo, kopt, nfresc
integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK)))
logical :: shortd, improve_geo, tr_success, reduce_rho, small_trrad, close_itpset, accurate_mod, bad_trstep


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

do while (.true.)
    ! Generate the next point in the trust region that provides a small value of the quadratic model
    ! subject to the constraints on the variables.

    !rescued = .false.
    call trsbox(delta, gopt, hq, pq, sl, su, xopt, xpt, crvmin, d)

    xnew = max(min(xopt + d, su), sl)  ! In precise arithmetic, XNEW = XOPT + D.

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

        gnew = gopt + hess_mul(d, xpt, pq, hq)
        bdtest = maxval(abs(moderrsav))
        bdtest(trueloc(xnew <= sl)) = gnew(trueloc(xnew <= sl)) * rho
        bdtest(trueloc(xnew >= su)) = -gnew(trueloc(xnew >= su)) * rho
        curv = diag(hq) + matprod(xpt**2, pq)
        errbd = minval(max(bdtest, bdtest + HALF * curv * rho**2))
        if (crvmin > 0) then
            errbd = min(errbd, 0.125_RP * crvmin * rho**2)
        end if
    else  ! D is long enough to invoke a function evaluation.
        ! Zaikun 20220528: TODO: check the shifting strategy of NEWUOA and LINCOA.
        if (sum(xopt**2) >= 1.0E3_RP * dsq) then
            sl = min(sl - xopt, ZERO)
            su = max(su - xopt, ZERO)
            xnew = min(max(sl, xnew - xopt), su)
            call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq)  ! XBASE is set to XOPT, XOPT to 0.
            xbase = min(max(xl, xbase), xu)
        end if
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
        if (.not. (qred > 0)) then
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

        tr_success = (f < fopt)

        ! Call RESCUE if if rounding errors have damaged the denominator corresponding to D.
        ! It provides a useful safeguard, but is not invoked in most applications of BOBYQA.
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        den = calden(kopt, bmat, d, xpt, zmat)
        !if (tr_success .and. .not. any(den > 0.25 * maxval(vlag(1:npt)**2))) then
        !if (tr_success .and. .not. any(den > 0.5 * maxval(vlag(1:npt)**2))) then
        !if (.not. any(den > HALF * maxval(vlag(1:npt)**2))) then
        !if (.not. any(den > maxval(vlag(1:npt)**2))) then
        if (tr_success .and. .not. (is_finite(sum(abs(vlag))) .and. any(den > maxval(vlag(1:npt)**2)))) then
            if (nf == nfresc) then  ! This cannot happen.
                info = DAMAGING_ROUNDING
                exit
            end if
            nfresc = nf
            call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fopt, &
                & fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat, subinfo)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
            nfresc = nf
            moderrsav = HUGENUM
            dnormsav = HUGENUM

            ! RESCUE shifts XBASE to the pre-RESCUE value of XOPT (even if RESCUED is FALSE).
            xnew = min(max(sl, d), su)
            d = xnew - xopt
            qred = -quadinc(d, xpt, gopt, pq, hq)
            diff = f - fopt + qred
            tr_success = (f < fopt)
        end if

        ! Find the index of the interpolation point to be replaced by the trust-region trial point.

        ! Calculate the distance squares between the interpolation points and the "optimal point". When
        ! identifying the optimal point, as suggested in (7.5) of the NEWUOA paper, it is reasonable to
        ! take into account the new trust-region trial point XPT(:, KOPT) + D, which will become the optimal
        ! point in the next interpolation if TR_SUCCESS is TRUE. Strangely, considering this new point does
        ! not always lead to a better performance of BOBYQA. Here, we choose not to check TR_SUCCESS, as
        ! the performance of BOBYQA is better in this way.
        ! HOWEVER, THIS MAY WELL CHANGE IF THE OTHER PARTS OF BOBYQA ARE IMPLEMENTED DIFFERENTLY.
        !if (tr_success) then
        !    distsq = sum((xpt - spread(xopt + d, dim=2, ncopies=npt))**2, dim=1)
        !else
        !    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
        !end if
        distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)

        weight = max(ONE, distsq / delta**2)**3  ! This works better than Powell's code.
        !------------------------------------------------------------------------------------------!
        ! Other possible definitions of WEIGHT.
        !weight = max(ONE, distsq / delta**2)**2  ! Powell's original code. Works well.
        !weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3  ! NEWUOA. Better than original.
        !weight = max(ONE, distsq / rho**2)**3  ! It performs the same as the code from NEWUOA.
        !weight = distsq**3  ! This works better than Powell's code.
        !weight = distsq**4  ! This works better than Powell's code.
        !weight = max(ONE, distsq / delta**2)**4  ! This works better than Powell's code.
        !weight = max(ONE, distsq / delta**2)  ! As per (6.1) of the BOBYQA paper. It works poorly!
        !------------------------------------------------------------------------------------------!

        den = calden(kopt, bmat, d, xpt, zmat)
        score = weight * den

        ! If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
        if (.not. tr_success) then
            score(kopt) = -ONE
        end if

        ! For the first case below, NEWUOA checks ANY(SCORE>1) .OR. (TR_SUCCESS .AND. ANY(SCORE>0))
        ! instead of ANY(SCORE > 0). Such code does not seem to improve the performance of BOBYQA.
        if (any(score > 0)) then
            ! See (6.1) of the BOBYQA paper for the definition of KNEW in this case.
            ! SCORE(K) = NaN implies DEN(K) = NaN. We exclude such K as we want DEN to be big.
            knew_tr = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
            !!MATLAB: [~, knew_tr] = max(score, [], 'omitnan');
        elseif (tr_success) then
            ! Powell's code does not include the following instructions. With Powell's code, if DEN
            ! consists of only NaN, then KNEW can be 0 even when TR_SUCCESS is TRUE.
            knew_tr = int(maxloc(distsq, dim=1), IK)
        else
            knew_tr = 0_IK  ! We arrive here when TR_SUCCESS = FALSE and no entry of SCORE is positive.
        end if

        if (knew_tr > 0) then
            ! Update BMAT and ZMAT, so that the KNEW-th interpolation point can be moved. Also update
            ! the second derivative terms of the model.
            vlag = calvlag(kopt, bmat, d, xpt, zmat)
            beta = calbeta(kopt, bmat, d, xpt, zmat)
            call updateh(knew_tr, beta, vlag, bmat, zmat)

            call r1update(hq, pq(knew_tr), xpt(:, knew_tr))
            pq(knew_tr) = ZERO
            pqinc = matprod(zmat, diff * zmat(knew_tr, :))
            pq = pq + pqinc
            ! Alternatives:
            !!PQ = PQ + MATPROD(ZMAT, DIFF * ZMAT(KNEW, :))
            !!PQ = PQ + DIFF * MATPROD(ZMAT, ZMAT(KNEW, :))

            ! Include the new interpolation point, and make the changes to GOPT at the old XOPT that are
            ! caused by the updating of the quadratic model.
            fval(knew_tr) = f
            xpt(:, knew_tr) = xnew
            gopt = gopt + diff * bmat(:, knew_tr) + hess_mul(xopt, xpt, pqinc)

            ! Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
            if (f < fopt) then
                kopt = knew_tr
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
    ! ACCURATE_MOD --- Are the recent models sufficiently accurate?
    accurate_mod = all(abs(moderrsav) <= errbd) .and. all(dnormsav <= rho)
    ! SMALL_TRRAD --- Is the trust-region radius small?
    small_trrad = (max(delta, dnorm) <= rho)
    ! CLOSE_ITPSET --- Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xopt).^2)  % xopt should be a column!! Implicit expansion
    !close_itpset = all(distsq <= 4.0_RP * delta**2)
    close_itpset = all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))
    !----------------------------------------------------------------------------------------------!


    ! What if RESCUE has been called? Is it reasonable to use F and FOPT?
    improve_geo = (shortd .and. .not. (all(abs(moderrsav) <= errbd) .and. all(dnormsav <= rho))) &
        & .or. (.not. shortd .and. .not. (knew_tr > 0 .and. .not. f >= fopt - TENTH * qred))
    !if (improve_geo) then
    dsquare = max((TWO * delta)**2, (TEN * rho)**2)
    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
    knew_geo = 0_IK
    if (.not. all(distsq <= max((TWO * delta)**2, (TEN * rho)**2))) then
        knew_geo = int(maxloc(distsq, dim=1), IK)  ! This line cannot be exchanged with the next
        dsquare = distsq(knew_geo) ! This line cannot be exchanged with the last
    end if

    reduce_rho = .not. ((.not. shortd) .and. knew_tr > 0 .and. .not. f >= fopt - TENTH * qred)  &
       & .and. (.not. improve_geo .or. (knew_geo <= 0 .and. (shortd .or. (ratio <= 0 .and. max(delta, dnorm) <= rho))))
    improve_geo = improve_geo .and. (knew_geo > 0) .and. &
        & .not. ((.not. shortd) .and. knew_tr > 0 .and. .not. f >= fopt - TENTH * qred)

    bad_trstep = (shortd .or. ratio <= 0 .or. knew_tr == 0)  ! BAD_TRSTEP for REDUCE_RHO
    !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

    !bad_trstep = (shortd .or. ratio <= TENTH .or. knew_tr == 0)  ! BAD_TRSTEP for IMPROVE_GEO
    bad_trstep = (shortd .or. f >= fopt - TENTH * qred .or. knew_tr == 0)  ! BAD_TRSTEP for IMPROVE_GEO
    improve_geo = bad_trstep .and. (.not. close_itpset) .and. (.not. reduce_rho)

    ! If KNEW is positive, then GEOSTEP finds alternative new positions for the KNEW-th
    ! interpolation point within distance DELBAR of XOPT. Otherwise, go for another trust region
    ! iteration, unless the calculations with the current RHO are complete.
    if (improve_geo) then
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
        d = geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat)

        ! Call RESCUE if if rounding errors have damaged the denominator corresponding to D.
        ! It provides a useful safeguard, but is not invoked in most applications of BOBYQA.
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        den = calden(kopt, bmat, d, xpt, zmat)
        !if (.not. (is_finite(sum(abs(vlag))) .and. den(knew_geo) > HALF * vlag(knew_geo)**2)) then
        if (.not. (is_finite(sum(abs(vlag))) .and. den(knew_geo) > vlag(knew_geo)**2)) then
            if (nf == nfresc) then
                info = DAMAGING_ROUNDING
                exit
            end if
            nfresc = nf
            call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fopt, &
                & fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat, subinfo)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if

            moderrsav = HUGENUM
            dnormsav = HUGENUM

            ! If NFRESC < NF, then RESCUE has updated XPT (also the model and [BMAT, ZMAT]) to
            ! improve (rescue) its geometry. Thus no geometry step is needed anymore. If NFRESC
            ! equals NF, then RESCUE did not make any change to XPT, but only recalculated the
            ! model ans [BMAT, ZMAT]; in this case, we calculate a new geometry step.
            if (nfresc == nf) then
                d = geostep(knew_geo, kopt, bmat, delbar, sl, su, xpt, zmat)
            else
                nfresc = nf
                cycle
            end if
        end if

        ! Put the variables for the next calculation of the objective function in XNEW, with any
        ! adjustments for the bounds. In precise arithmetic, X = XBASE + XNEW.
        xnew = min(max(sl, xopt + d), su)
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
        vlag = calvlag(kopt, bmat, d, xpt, zmat)
        beta = calbeta(kopt, bmat, d, xpt, zmat)
        call updateh(knew_geo, beta, vlag, bmat, zmat)

        call r1update(hq, pq(knew_geo), xpt(:, knew_geo))
        pq(knew_geo) = ZERO
        pqinc = matprod(zmat, diff * zmat(knew_geo, :))
        pq = pq + pqinc
        ! Alternatives:
            !!PQ = PQ + MATPROD(ZMAT, DIFF * ZMAT(KNEW, :))
            !!PQ = PQ + DIFF * MATPROD(ZMAT, ZMAT(KNEW, :))

        ! Include the new interpolation point, and make the changes to GOPT at the old XOPT that are
        ! caused by the updating of the quadratic model.
        fval(knew_geo) = f
        xpt(:, knew_geo) = xnew
        gopt = gopt + diff * bmat(:, knew_geo) + hess_mul(xopt, xpt, pqinc)

        ! Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
        if (f < fopt) then
            kopt = knew_geo
            xopt = xnew
            gopt = gopt + hess_mul(d, xpt, pq, hq)
        end if
    end if

    ! The calculations with the current value of RHO are complete. Update RHO and DELTA.
    if (reduce_rho) then
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
    end if
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

close (16)
!====================!
!  Calculation ends  !
!====================!

! Postconditions


end subroutine bobyqb


end module bobyqb_mod
