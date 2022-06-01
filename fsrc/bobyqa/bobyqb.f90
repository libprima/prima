module bobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! TODO: verify that the iterates/steps respect bounds in the pre/postconditions.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, June 01, 2022 PM09:24:38
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: bobyqb


contains


subroutine bobyqb(calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, &
    & sl_in, su_in, xl, xu, x, nf, f, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the major calculations of BOBYQA.
!
! The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to the
! corresponding arguments in SUBROUTINE BOBYQA.
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
! of the model and Lagrange functions.
! XPT is a TWO-dimensional array that holds the coordinates of the interpolation points relative to XBASE.
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
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: info_mod, only : NAN_INF_X, NAN_INF_F, NAN_MODEL, FTARGET_ACHIEVED, INFO_DFT, &
    & MAXFUN_REACHED, DAMAGING_ROUNDING, TRSUBP_FAILED, SMALL_TR_RADIUS!, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : matprod, diag, trueloc, r1update!, r2update!, norm
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calvlag, calbeta, hess_mul
!!! TODO (Zaikun 20220525): Use CALDEN instead of CALVLAG and CALBETA wherever possible!!!

! Solver-specific modules
use, non_intrinsic :: initialize_mod, only : initialize
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
real(RP), intent(in) :: sl_in(:)  ! SL_IN(N)
real(RP), intent(in) :: su_in(:)  ! SU_IN(N)
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
real(RP) :: delbar, alpha, bdtest(size(x)), hqdiag(size(x)), bdtol, beta, &
&        biglsq, crvmin, curv(size(x)), delta,  &
&        den(npt), denom, densav, diff, diffa, diffb, diffc,     &
&        dist, dsquare, distsq(npt), dnorm, dsq, errbig, fopt,        &
&        frhosq, gisq, gqsq, hdiag(npt),      &
&        ratio, rho, scaden, qred, weight(npt), pqinc(npt)
real(RP) :: pqalt(npt), galt(size(x)), fshift(npt), pgalt(size(x)), pgopt(size(x))
real(RP) :: score(npt), wlagsq(npt)
integer(IK) :: itest, knew, kopt, ksav, nfsav, nresc, ntrits
logical :: shortd, improve_geo, rescued, geo_step


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
    call assert(size(sl_in) == n .and. size(su_in) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(all(sl_in <= 0 .and. (sl_in >= 0 .or. sl_in <= -rhobeg)), 'SL == 0 or SL <= -RHOBEG', srname)
    call assert(all(su_in >= 0 .and. (su_in <= 0 .or. su_in >= rhobeg)), 'SU == 0 or SU >= RHOBEG', srname)
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

sl = sl_in
su = su_in

! Initialize XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT and ZMAT together with the corresponding values of
! of NF and KOPT, which are the number of calls of CALFUN so far and the index of the interpolation
! point at the trust region centre. Then the initial XOPT is set too. The branch to label 720 occurs
! if MAXFUN is less than NPT. GOPT will be updated if KOPT is different from KBASE.
call initialize(calfun, iprint, ftarget, rhobeg, sl, su, x, xl, xu, kopt, nf, bmat, f, fhist, fval, &
    & gopt, hq, pq, xbase, xhist, xpt, zmat)
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = min(max(xl, xbase + xopt), xu)
x(trueloc(xopt <= sl)) = xl(trueloc(xopt <= sl))
x(trueloc(xopt >= su)) = xu(trueloc(xopt >= su))
f = fopt

if (is_nan(f) .or. is_posinf(f)) then
    info = NAN_INF_F
    xopt = xpt(:, kopt)
    fopt = fval(kopt)
    x = min(max(xl, xbase + xopt), xu)
    x(trueloc(xopt <= sl)) = xl(trueloc(xopt <= sl))
    x(trueloc(xopt >= su)) = xu(trueloc(xopt >= su))
    f = fopt
    return
end if
if (f <= ftarget) then
    info = FTARGET_ACHIEVED
    xopt = xpt(:, kopt)
    fopt = fval(kopt)
    x = min(max(xl, xbase + xopt), xu)
    x(trueloc(xopt <= sl)) = xl(trueloc(xopt <= sl))
    x(trueloc(xopt >= su)) = xu(trueloc(xopt >= su))
    f = fopt
    call rangehist(nf, xhist, fhist)
    return
end if

! Complete the settings that are required for the iterative procedure.
rho = rhobeg
delta = rho
nresc = nf
ntrits = 0
diffa = ZERO
diffb = ZERO
itest = 0
nfsav = nf

! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
ratio = -ONE
shortd = .false.
improve_geo = .false.
geo_step = .false.
rescued = .false.

! Generate the next point in the trust region that provides a small value of the quadratic model
! subject to the constraints on the variables. The integer NTRITS is set to the number "trust
! region" iterations that have occurred since the last "alternative" iteration. If the length of
! XNEW-XOPT is less than HALF*RHO, however, then there is a branch to label 650 or 680 with
! NTRITS=-1, instead of calculating F at XNEW.
!--------------------------------------------------------------------------------------------------!
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the models. In such a case, we
! terminate the code. Otherwise, the behavior of TRBOX, GEOSTEP, or RESCUE is not predictable, and
! Segmentation Fault or infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended  behavior of the code,
! including uninitialized indices. STILL NECESSARY???
!--------------------------------------------------------------------------------------------------!

do while (.true.)
    if (.not. geo_step) then
        if (is_nan(sum(abs(gopt)) + sum(abs(hq)) + sum(abs(pq)))) then
            info = NAN_MODEL
            exit
        end if

        call trsbox(delta, gopt, hq, pq, sl, su, xopt, xpt, crvmin, d, dsq, gnew, xnew)

        dnorm = min(delta, sqrt(dsq))
        shortd = (dnorm < HALF * rho)

        ! When D is short make a choice between improving the geometry and reducing RHO depends on
        ! whether or not our work with the current RHO seems to be complete. Either RHO is decreased
        ! or termination occurs if the errors in the quadratic model at the last three interpolation
        ! points compare favourably with predictions of likely improvements to the model within
        ! distance HALF*RHO of XOPT.
        if (shortd) then
            ntrits = -1
            dsquare = (TEN * rho)**2
            if (nf <= nfsav + 2) then
                improve_geo = .true.
            else
                errbig = max(diffa, diffb, diffc)
                frhosq = 0.125_RP * rho * rho
                if (crvmin > ZERO .and. errbig > frhosq * crvmin) then
                    improve_geo = .true.
                else
                    bdtol = errbig / rho
                    bdtest = bdtol
                    bdtest(trueloc(xnew <= sl)) = gnew(trueloc(xnew <= sl))
                    bdtest(trueloc(xnew >= su)) = -gnew(trueloc(xnew >= su))
                    hqdiag = diag(hq)
                    curv = ZERO ! Entertain Fortran compilers. No need in MATLAB/Python/Julia/R.
                    where (bdtest < bdtol) curv = hqdiag + matprod(xpt**2, pq)
                    !!MATLAB: curv(bdtest < bdtol) = hqdiag(bdtest < bdtol) + xpt(bdtest < bdtol, :).^2 * pq
                    if (any(bdtest < bdtol .and. bdtest + HALF * curv * rho < bdtol)) then
                        improve_geo = .true.
                    else
                        improve_geo = .false.
                    end if
                end if
            end if
        else
            ntrits = ntrits + 1_IK
        end if
    end if

    if (geo_step .or. .not. shortd) then
        geo_step = .false.
        ! Severe cancellation is likely to occur if XOPT is too far from XBASE. If the following
        ! test holds, then XBASE is shifted so that XOPT becomes zero. The appropriate changes are
        ! made to BMAT and to the second derivatives of the current model, beginning with the
        ! changes to BMAT that are independent of ZMAT. VLAG is used temporarily for working space.
        ! Zaikun 20220528: TODO: check the shifting strategy of NEWUOA and LINCOA.
        if (sum(xopt**2) >= 1.0E3_RP * dsq .and. .not. rescued) then
            sl = min(sl - xopt, ZERO)
            su = max(su - xopt, ZERO)
            !xnew = xnew - xopt  ! Needed? Will XNEW be used again later?
            xnew = min(max(sl, xnew - xopt), su)  ! Needed? Will XNEW be used again later?
            call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq)
            xbase = min(max(xl, xbase), xu)
        end if

        rescued = .false.

        ! Pick two alternative vectors of variables, relative to XBASE, that are suitable as new
        ! positions of the KNEW-th interpolation point. Firstly, XNEW is set to the point on a line
        ! through XOPT and another interpolation point that minimizes the predicted value of the
        ! next denominator, subject to ||XNEW - XOPT|| .LEQ. DELBAR and to the SL and SU bounds.
        ! Secondly, XALT is set to the best feasible point on a constrained version of the Cauchy
        ! step of the KNEW-th Lagrange function, the corresponding value of the square of this
        ! function being returned in CAUCHY. The choice between these alternatives is going to be
        ! made when the denominator is calculated.

        !------------------------------------------------------------------------------------------!
        !  Zaikun 23-07-2019: (STILL NECESSARY???)
        ! Although very rare, NaN can sometimes occur in BMAT or ZMAT. If it happens, we terminate
        ! the code. See the comments above line number 60. Indeed, if GEOSTEP is called with such
        ! matrices, then geostep.f90 will encounter a memory error at lines 173--174. This is
        ! because the first value of PREDSQ in GEOSTEP (see line 159 of geostep.f90) will be NaN,
        ! line 164 will not be reached, and hence no value will be assigned to IBDSAV. Such an error
        ! was observed when BOBYQA was (mistakenly) tested on CUTEst problem CONCON. CONCON is
        ! a nonlinearly constrained problem with bounds. By mistake, BOBYQA was called to solve this
        ! problem, neglecting all the constraints other than bounds. With only the bound
        ! constraints, the objective function turned to be unbounded from below, which led to
        ! abnormal values in BMAT(indeed, BETA defined in lines 366--389 took NaN/infinite values).
        !------------------------------------------------------------------------------------------!

        if (ntrits == 0) then
            if (is_nan(sum(abs(bmat)) + sum(abs(zmat)))) then
                info = NAN_MODEL
                exit
            end if

            ! Calculate a geometry step.
            d = geostep(knew, kopt, bmat, delbar, sl, su, xpt, zmat)
            !xnew = xopt + d
            xnew = min(max(sl, xopt + d), su)

            ! Calculate VLAG, BETA, and DENOM for the current choice of D.
            alpha = sum(zmat(knew, :)**2)
            vlag = calvlag(kopt, bmat, d, xpt, zmat)
            beta = calbeta(kopt, bmat, d, xpt, zmat)
            denom = alpha * beta + vlag(knew)**2

            ! Call RESCUE if if rounding errors have damaged the denominator corresponding to D.
            !if (.not. (denom > HALF * vlag(knew)**2)) then
            if (.not. (denom > vlag(knew)**2)) then  ! This is used when verifying RESCUE
                if (nf <= nresc) then
                    info = DAMAGING_ROUNDING
                    exit
                else
                    ! XBASE is also moved to XOPT by a call of RESCUE. This calculation is more
                    ! expensive than the previous shift, because new matrices BMAT and ZMAT are
                    ! generated from scratch, which may include the replacement of interpolation
                    ! points whose positions seem to be causing near linear dependence in the
                    ! interpolation conditions. Therefore RESCUE is called only if rounding errors
                    ! have reduced by at least a factor of TWO the denominator of the formula for
                    ! updating the H matrix. It provides a useful safeguard, but is not invoked in
                    ! most applications of BOBYQA.
                    !--------------------------------------------------------------------------------------!
                    ! Zaikun 2019-08-29: STILL NECESSARY?
                    if (is_nan(sum(abs(gopt)) + sum(abs(hq)) + sum(abs(pq)))) then
                        info = NAN_MODEL
                        exit
                    end if
                    if (is_nan(sum(abs(bmat)) + sum(abs(zmat)))) then
                        info = NAN_MODEL
                        exit
                    end if
                    !--------------------------------------------------------------------------------------!

                    nfsav = nf
                    call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, &
                        & fhist, fopt, fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat)

                    ! XOPT is updated now in case the loop will exit. {?}
                    ! Update GOPT?
                    if (fopt <= ftarget) then
                        info = FTARGET_ACHIEVED
                        exit
                    end if
                    if (nf >= maxfun) then
                        info = MAXFUN_REACHED
                        exit
                    end if
                    nresc = nf
                    if (nfsav < nf) then
                        nfsav = nf
                        cycle
                    end if
                    rescued = .true.
                    geo_step = .true.
                    cycle
                    !------------------------------------------------------------------------------!
                    ! After RESCUE, Powell's code takes immediately another GEOSTEP. If the geometry
                    ! becomes acceptable (i.e., DENOM > HALF * VLAG(KNEW)**2), then the algorithm
                    ! will continue.Otherwise, it will exit with INFO set to DAMAGING_ROUNDING.
                    !------------------------------------------------------------------------------!
                end if
            end if

            ! Alternatively, if NTRITS is positive, then set KNEW to the index of the next
            ! interpolation point to be deleted to make room for a trust region step. Again RESCUE
            ! may be called if rounding errors have damaged the chosen denominator, which is the
            ! reason for attempting to select KNEW before calculating the next value of the
            ! objective function.
        else
            ! Calculate VLAG and BETA for the current choice of D.
            vlag = calvlag(kopt, bmat, d, xpt, zmat)
            beta = calbeta(kopt, bmat, d, xpt, zmat)
            hdiag = sum(zmat**2, dim=2)
            den = hdiag * beta + vlag(1:npt)**2
            distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
            weight = max(ONE, (distsq / delta**2)**2)

            score = weight * den
            score(kopt) = -ONE  ! Skip KOPT when taking the maximum of SCORE
            knew = 0
            scaden = ZERO
            if (any(score > 0)) then
                ! SCORE(K) is NaN implies DENABS(K) is NaN, but we want DEN to be big. So we exclude
                ! such K.
                knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
                scaden = score(knew)
                !!MATLAB: [scaden, knew] = max(score, [], 'omitnan');
                denom = den(knew)
            end if

            wlagsq = weight * vlag(1:npt)**2
            wlagsq(kopt) = -ONE  ! Skip KOPT when taking the maximum of WLAGSQ
            biglsq = ZERO
            if (any(wlagsq > 0)) then
                biglsq = maxval(wlagsq, mask=(.not. is_nan(wlagsq)))
                !!MATLAB: biglsq = max(wlagsq, [], 'omitnan');
            end if

            ! KNEW > 0 is implied by SCADEN > HALF*BIGLSQ (but NOT SCADEN >= ...), yet we prefer to
            ! require KNEW > 0 explicitly.
            !if (.not. (knew > 0 .and. scaden > HALF * biglsq)) then
            if (.not. (knew > 0 .and. scaden > biglsq)) then  ! This is used when verifying RESCUE.
                if (nf <= nresc) then
                    info = DAMAGING_ROUNDING
                    exit
                else
                    ! XBASE is also moved to XOPT by a call of RESCUE. This calculation is more
                    ! expensive than the previous shift, because new matrices BMAT and ZMAT are
                    ! generated from scratch, which may include the replacement of interpolation
                    ! points whose positions seem to be causing near linear dependence in the
                    ! interpolation conditions. Therefore RESCUE is called only if rounding errors
                    ! have reduced by at least a factor of TWO the denominator of the formula for
                    ! updating the H matrix. It provides a useful safeguard, but is not invoked in
                    ! most applications of BOBYQA.

                    !------------------------------------------------------------------------------!
                    ! Zaikun 2019-08-29: See the comments above line number 60. STILL NECESSARY?
                    if (is_nan(sum(abs(gopt)) + sum(abs(hq)) + sum(abs(pq)))) then
                        info = NAN_MODEL
                        exit
                    end if
                    if (is_nan(sum(abs(bmat)) + sum(abs(zmat)))) then
                        info = NAN_MODEL
                        exit
                    end if
                    !------------------------------------------------------------------------------!

                    nfsav = nf
                    call rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, &
                        & fhist, fopt, fval, gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat)

                    ! XOPT is updated now in case the loop exits. (?)
                    ! Update GOPT?
                    if (fopt <= ftarget) then
                        info = FTARGET_ACHIEVED
                        exit
                    end if
                    if (nf >= maxfun) then
                        info = MAXFUN_REACHED
                        exit
                    end if
                    nresc = nf
                    nfsav = max(nfsav, nf)
                    cycle
                end if
            end if
        end if

        ! Put the variables for the next calculation of the objective function in XNEW, with any
        ! adjustments for the bounds. Calculate the value of the objective function at XBASE+XNEW.
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

        call evaluate(calfun, x, f)
        call savehist(nf, x, xhist, f, fhist)

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
        diffc = diffb
        diffb = diffa
        diffa = abs(diff)
        if (dnorm > rho) nfsav = nf

        ! Pick the next value of DELTA after a trust region step.
        if (ntrits > 0) then
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

            ! Recalculate KNEW and DENOM if the new F is less than FOPT.
            if (f < fopt) then
                ksav = knew
                densav = denom
                hdiag = sum(zmat**2, dim=2)
                den = hdiag * beta + vlag(1:npt)**2
                distsq = sum((xpt - spread(xnew, dim=2, ncopies=npt))**2, dim=1)
                weight = max(ONE, (distsq / delta**2)**2)
                score = weight * den

                knew = 0
                scaden = ZERO
                if (any(score > 0)) then
                    ! SCORE(K) is NaN implies DENABS(K) is NaN, but we want DEN to be big. So we
                    ! exclude such K.
                    knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
                    scaden = score(knew)
                    !!MATLAB: [scaden, knew] = max(score, [], 'omitnan');
                    denom = den(knew)
                end if

                wlagsq = weight * vlag(1:npt)**2
                biglsq = ZERO
                if (any(wlagsq > 0)) then
                    biglsq = maxval(wlagsq, mask=(.not. is_nan(wlagsq)))
                    !!MATLAB: biglsq = max(wlagsq, [], 'omitnan');
                end if

                ! KNEW > 0 is implied by SCADEN > HALF*BIGLSQ (but NOT SCADEN >= ...), yet prefer to
                ! require KNEW > 0 explicitly.
                if (.not. (knew > 0 .and. scaden > HALF * biglsq)) then
                    knew = ksav
                    denom = densav
                end if
            end if
        end if

        ! Update BMAT and ZMAT, so that the KNEW-th interpolation point can be moved. Also update
        ! the second derivative terms of the model.
        !------------------------------------------------------------------------------------------!
        call assert(knew >= 1, 'KNEW >= 1', srname)
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
        if (ntrits > 0) then
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
        ! trust region calculation. The case NTRITS=0 occurs when the new interpolation point was
        ! reached by an alternative step.
        if (ntrits == 0) cycle
        if (f <= fopt - TENTH * qred) cycle

        ! Alternatively, find out if the interpolation points are close enough to the best point so far.
        dsquare = max((TWO * delta)**2, (TEN * rho)**2)
        improve_geo = .true.
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
            if (ntrits == -1) then
                delta = min(TENTH * delta, HALF * dist)
                if (delta <= 1.5_RP * rho) delta = rho
            end if
            ntrits = 0
            delbar = max(min(TENTH * dist, delta), rho)
            dsq = delbar * delbar
            geo_step = .true.
            cycle
        else
            geo_step = .false.
        end if
        if (ntrits /= -1 .and. (ratio > 0 .or. max(delta, dnorm) > rho)) cycle
    end if

    ! The calculations with the current value of RHO are complete. Update RHO and DELTA.
    if (rho > rhoend) then
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
        ntrits = 0
        nfsav = nf
        cycle
    else
        info = SMALL_TR_RADIUS
        exit
    end if

end do

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
if (info == SMALL_TR_RADIUS .and. ntrits == -1 .and. nf < maxfun) then
    info = INFO_DFT  !!?? See NEWUOA
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

!close (16)

end subroutine bobyqb


end module bobyqb_mod
