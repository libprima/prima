module newuob_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of NEWUOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Saturday, October 02, 2021 PM02:22:12
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: newuob


contains


subroutine newuob(calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, &
    & rhoend, x, nf, f, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the actual calculations of NEWUOA. The arguments IPRINT, MAXFUN, MAXHIST,
! NPT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, X, NF, F, FHIST, XHIST, and INFO are
! identical to the corresponding arguments in subroutine NEWUOA.
!
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
! of the model and Lagrange functions.
! XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
! the least calculated F so far). FOPT = F(XOPT + XBASE).
! D is reserved for trial steps from XOPT.
! XNEW = XOPT + D, corresponding to the vector of variables for the next calculation of F.
! [XPT, FVAL, KOPT] describes the interpolation set:
! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
! However, there is a delay between the update of XOPT and KOPT. So they are not always consistent
! in the mid of an iteration. See the comment on the update of XOPT for details.
! values of F at the interpolation points; KOPT is the index of XOPT in XPT (XPT(:,KOPT) = XOPT).
! [GQ, HQ, PQ] describes the quadratic model: GQ will hold the gradient of the quadratic model at
! XBASE; HQ will hold the explicit second order derivatives of the quadratic model; PQ will contain
! the parameters of the implicit second order derivatives of the quadratic model.
! [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
! the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
! BMAT will hold the last N ROWs of H; ZMAT will hold the factorization of the leading NPT*NPT
! submatrix of H, this factorization being ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT)=1.
! VLAG will contain the values of the Lagrange functions at a new point X. They are part of a
! product that requires VLAG to be of length NPT + N. Both VLAG and BETA are critical for the
! updating procedure of H, which is detailed formula (4.10)--(4.12) of the NEWUOA paper.
!
! See Section 2 of the NEWUOA paper for more information about these variables.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evalf
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: info_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS
use, non_intrinsic :: linalg_mod, only : calquad, inprod, norm
use, non_intrinsic :: output_mod, only : retmssg, rhomssg, fmssg
use, non_intrinsic :: pintrf_mod, only : FUN
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: resolution_mod, only : resenhance

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : setdrop_tr, geostep
use, non_intrinsic :: initialize_mod, only : initxf, initq, inith
use, non_intrinsic :: shiftbase_mod, only : shiftbase
use, non_intrinsic :: trustregion_mod, only : trsapp, trrad
use, non_intrinsic :: update_mod, only : updateh, updateq, updatexf, tryqalt
use, non_intrinsic :: vlagbeta_mod, only : vlagbeta

implicit none

! Inputs
procedure(FUN) :: calfun
! N.B.: The INTENT attribute cannot be specified for a dummy procedure without the POINTER attribute
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

! In-outputs
real(RP), intent(inout) :: x(:)      ! X(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)   ! FHIST(MAXFHIST)
real(RP), intent(out) :: xhist(:, :)    ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'NEWUOA'
character(len=*), parameter :: srname = 'NEWUOB'
integer(IK) :: idz
integer(IK) :: ij(max(0_IK, int(npt - 2 * size(x) - 1, IK)), 2_IK)
integer(IK) :: itest
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
logical :: bad_trstep
logical :: improve_geo
logical :: reduce_rho_1
logical :: reduce_rho_2
logical :: shortd
real(RP) :: beta
real(RP) :: bmat(size(x), npt + size(x))
real(RP) :: crvmin
real(RP) :: d(size(x))
real(RP) :: delbar
real(RP) :: delta
real(RP) :: dnorm
real(RP) :: dnormsav(3)
real(RP) :: fopt
real(RP) :: fval(npt)
real(RP) :: gq(size(x))
real(RP) :: hq(size(x), size(x))
real(RP) :: moderr
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pq(npt)
real(RP) :: qred
real(RP) :: ratio
real(RP) :: rho
real(RP) :: trtol
real(RP) :: vlag(npt + size(x))
real(RP) :: xbase(size(x))
real(RP) :: xdist(npt)
real(RP) :: xnew(size(x))
real(RP) :: xopt(size(x))
real(RP) :: xpt(size(x), npt)
real(RP) :: zmat(npt, npt - size(x) - 1)

! Sizes
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', srname)
    call assert(rhobeg >= rhoend .and. rhoend > ZERO, 'RHOBEG >= RHOEND > 0', srname)
    call assert(eta1 >= ZERO .and. eta1 <= eta2 .and. eta2 < ONE, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > ZERO .and. gamma1 < ONE .and. gamma2 > ONE, &
        & '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize FVAL, XBASE, and XPT.
call initxf(calfun, iprint, maxfun, ftarget, rhobeg, x, ij, kopt, nf, fhist, fval, xbase, xhist, xpt, subinfo)
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = xbase + xopt  ! Set X.
f = fopt  ! Set F.

! Check whether to return after initialization.
if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Arrange FHIST and XHIST so that they are in the chronological order.
    call rangehist(nf, fhist, xhist)
    call retmssg(info, iprint, nf, f, x, solver)
    ! Postconditions
    if (DEBUGGING) then
        call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
        call assert(.not. any(fhist(1:min(nf, maxfhist)) < f), 'F is the smallest in FHIST', srname)
    end if
    return
end if

! Initialize GQ, HQ, and PQ.
call initq(ij, fval, xpt, gq, hq, pq, subinfo)

! Initialize BMAT and ZMAT, and IDZ.
call inith(ij, xpt, idz, bmat, zmat, subinfo)

! After initializing GQ, HQ, PQ, BMAT, ZMAT, one can also choose to return if subinfo = NAN_MODEL
! (NaN occurs in the model). We do not do it here. If such a model is harmful, then it will probably
! lead to other returns (NaN in X, NaN in F, trust-region subproblem fails, ...); otherwise, the
! code will continue to run and possibly get rid of the NaN in the model.

! Set some more initial values and parameters.
rho = rhobeg
delta = rho
moderrsav = HUGENUM
dnormsav = HUGENUM
itest = 0_IK
trtol = 1.0E-2_RP  ! Tolerance used in trsapp.
! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
ratio = -ONE

! Each trust-region iteration takes at most two function evaluation. The following setting imposes
! no constraint on the maximal number of trust-region iterations.
maxtr = 2_IK * maxfun
! MAXTR is unlikely to be reached, but we define the following default value for INFO for safety.
info = MAXTR_REACHED

! Begin the iterative procedure.
! After solving a trust-region subproblem, NEWUOA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust-region iteration?
! REDUCE_RHO - Will we reduce rho after the trust-region iteration?
! REDUCE_RHO = REDUCE_RHO_1 .OR. REDUCE_RHO_2 (see boxes 14 and 10 of Fig. 1 in the NEWUOA paper).
! NEWUOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Solve the trust-region subproblem.
    call trsapp(delta, gq, hq, pq, trtol, xopt, xpt, crvmin, d, subinfo)

    ! Calculate the length of the trial step D.
    dnorm = min(delta, norm(d))

    ! SHORTD corresponds to Box 3 of the NEWUOA paper.
    shortd = (dnorm < HALF * rho)
    ! REDUCE_RHO_1 corresponds to Box 14 of the NEWUOA paper.
    reduce_rho_1 = shortd .and. (maxval(abs(moderrsav)) <= 0.125_RP * crvmin * rho**2) .and. &
        & (maxval(dnormsav) <= rho)
    if (shortd .and. (.not. reduce_rho_1)) then
        ! Reduce DELTA. After this, DELTA < DNORM may hold.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close.
        end if
    end if

    if (.not. shortd) then  ! D is long enough.
        ! Shift XBASE if XOPT may be too far from XBASE.
        !if (inprod(d, d) <= 1.0e-3_RP*inprod(xopt, xopt)) then  ! Powell's code
        if (dnorm**2 <= 1.0E-3_RP * inprod(xopt, xopt)) then
            call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt)
        end if

        ! Calculate VLAG and BETA for D. It makes uses of XOPT, so this is done before updating XOPT.
        call vlagbeta(idz, kopt, bmat, d, xpt, zmat, beta, vlag)

        ! Use the current quadratic model to predict the change in F due to the step D.
        qred = calquad(d, gq, hq, pq, xopt, xpt)

        ! Calculate the next value of the objective function.
        xnew = xopt + d
        x = xbase + xnew
        call evalf(calfun, x, f)
        nf = nf + 1_IK
        call fmssg(iprint, nf, f, x, solver)
        ! Save X and F into the history.
        call savehist(nf, f, x, fhist, xhist)
        ! Check whether to exit
        subinfo = checkexit(maxfun, nf, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        ! DNORMSAVE constains the DNORM of the latest 3 function evaluations with the current RHO.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

        ! MODERR is the error of the current model in predicting the change in F due to D.
        moderr = f - fopt + qred
        ! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
        moderrsav = [moderrsav(2:size(moderrsav)), moderr]

        ! Calculate the reduction ratio.
        ratio = redrat(fopt - f, qred, eta1)
        ! Update DELTA. After this, DELTA < DNORM may hold.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if

        ! Set KNEW_TR to the index of the interpolation point that will be replaced by XNEW. KNEW_TR
        ! will ensure that the geometry of XPT is "good enough" after the replacement.
        ! N.B.:
        ! 1. The information of XNEW is in VLAG and BETA, which are calculated from D = XNEW - XOPT.
        ! 2. KNEW_TR = 0 means it is impossible to obtain a good interpolation set by replacing any
        ! current interpolation point with XNEW. Then XNEW and its function value will be discarded.
        ! 3. If RATIO > 0, then SETDROP_TR ensures KNEW_TR > 0 so that the XNEW is included into XPT.
        ! Otherwise, SETDROP_TR is buggy; in addition, if RATIO > 0 and KNEW_TR = 0, XOPT will
        ! differ from XPT(:, KOPT), because the former is set to XNEW but XNEW is discarded. Such
        ! a difference can lead to unexpected behaviors; for example, KNEW_GEO may equal KOPT, with
        ! which GEOSTEP will not work.
        if (f < fopt) then
            xdist = sqrt(sum((xpt - spread(xnew, dim=2, ncopies=npt))**2, dim=1))
        else
            xdist = sqrt(sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1))
        end if
        knew_tr = setdrop_tr(idz, kopt, beta, delta, ratio, rho, vlag(1:npt), xdist, zmat)

        if (knew_tr > 0) then
            ! If KNEW_TR > 0, then update BMAT, ZMAT and IDZ, so that the KNEW_TR-th interpolation
            ! point is replaced by XNEW, whose information is encoded in VLAG and BET. If KNEW_TR
            ! is 0, then probably the geometry of XPT needs improvement, which will be handled below.
            call updateh(knew_tr, beta, vlag, idz, bmat, zmat)
            ! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
            call updateq(idz, knew_tr, bmat(:, knew_tr), moderr, zmat, xpt(:, knew_tr), gq, hq, pq)
            ! Include XNEW into XPT. Then update KOPT, XOPT, and FOPT.
            call updatexf(knew_tr, f, xnew, kopt, fval, xpt, fopt, xopt)
        end if

        ! Test whether to replace the new quadratic model Q by the least-Frobenius norm interpolant
        ! Q_alt. Perform the replacement if certain criteria are satisfied.
        ! N.B.:
        ! 1. This part is OPTIONAL, but it is crucial for the performance on some problems. See
        ! Section 8 of the NEWUOA paper.
        ! 2. TRYQALT is called only after a trust-region step but not after a geometry step, maybe
        ! because the model is expected to be good after a geometry step.
        ! 3. If KNEW_TR = 0 after a trust-region step, TRYQALT is not invoked. In this case, the
        ! interpolation set is unchanged, so it seems reasonable to keep the model unchanged.
        ! 4. In theory, FVAL - FOPT in the call of TRYQALT can be replaced by FVAL + C with any
        ! constant C. This constant will not affect the result in precise arithmetic. Powell chose
        ! C = - FVAL(KOPT_OLD), where KOPT_OLD is the KOPT before the update above (Powell updated
        ! KOPT after TRYQALT). Here we use C = -FOPT, as it worked slightly better on CUTEst,
        ! although there is no difference theoretically. Note that FVAL(KOPT_OLD) may not equal
        ! FOPT_OLD --- it may happen that KNEW_TR = KOPT_OLD so that FVAL(KOPT_OLD) has been revised
        ! after the last function evaluation.
        ! 5. Question: Since TRYQALT is invoked only when DELTA equals the current RHO, why not
        ! reset ITEST to 0 when RHO is reduced?
        if (knew_tr > 0 .and. delta <= rho) then  ! DELTA = RHO.
            call tryqalt(idz, fval - fopt, ratio, bmat(:, 1:npt), zmat, itest, gq, hq, pq)
        end if
    end if  ! End of if (.not. shortd)

    ! Before next trust-region iteration, we may improve the geometry of XPT or reduce rho
    ! according to IMPROVE_GEO and REDUCE_RHO. Now we decide these two indicators.

    ! First define IMPROVE_GEO, which corresponds to Box 8 of the NEWUOA paper.
    ! The geometry of XPT likely needs improvement if the trust-region step bad --- either too short
    ! (SHORTD = TRUE) or the reduction ratio is small (RATIO < TENTH). However, if REDUCE_RHO_1 is
    ! TRUE, meaning that the step is short and the latest model errors have been small, then we do
    ! not need to improve the geometry; instead, RHO will be reduced.
    ! To improve the geometry of XPT, we will check whether the interpolation points are all close
    ! enough to the best point so far, i.e., all the points are within a ball centered at XOPT with
    ! a radius of 2*DELTA. If not, the farthest point will be replaced with a point selected by
    ! GEOSTEP, aiming to ameliorate the geometry of the interpolation set; if yes, then RHO will be
    ! reduced if MAX(DELTA, DNORM) <= RHO (if MAX(DELTA, DNORM) > RHO, then, as Powell mentioned
    ! under (2.3) of the NEWUOA paper, "RHO has not restricted the most recent choice of D", so it
    ! is not reasonable to reduce RHO).
    ! N.B.:
    ! 1. RATIO must be set even if SHORTD = TRUE. Otherwise, compilers will raise a run-time error.
    ! 2. If SHORTD = FALSE and KNEW_TR = 0, then IMPROVE_GEO = TRUE. Therefore, IMPROVE_GEO = TRUE
    ! if it is impossible to obtain a good XPT by replacing a current point with the one suggested
    ! by the trust-region step.
    ! 3. If REDUCE_RHO = FALSE and SHORTD = TRUE, then the trust-region step is not tried at all,
    ! i.e., no function evaluation is invoked at XOPT + D (when REDUCE_RHO = TRUE, the step is not
    ! tried either, but the same step will be generated again at the next trust-region iteration
    ! after RHO is reduced and DELTA is updated; see the end of Section 2 of the NEWUOA paper).
    ! 4. If SHORTD = FALSE and KNEW_TR = 0, then the trust-region step invokes a function evaluation
    ! at XOPT + D, but [XOPT + D, F(XOPT +D)] is not included into [XPT, FVAL]. In other words, this
    ! function value is discarded.
    ! 5. If SHORTD = FALSE and KNEW_TR > 0 and RATIO < TENTH, then [XPT, FVAL] is updated so that
    ! [XPT(KNEW_TR), FVAL(KNEW_TR)] = [XOPT + D, F(XOPT + D)], and the model is updated accordingly,
    ! but such a model will not be used in the next trust-region iteration, because a geometry step
    ! will be invoked to improve the geometry of the interpolation set and update the model again.
    ! 6. DELTA has been updated before arriving here: if REDUCE_RHO = FALSE and SHORTD = TRUE, then
    ! DELTA was reduced by a factor of 10; if SHORTD = FALSE, then DELTA was updated by TRRAD after
    ! the trust-region iteration.
    ! 7. If SHORTD = FALSE and KNEW_TR > 0, then XPT has been updated after the trust-region
    ! iteration; if RATIO > 0 in addition, then XOPT has been updated as well.

    xdist = sqrt(sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1))
    bad_trstep = (shortd .or. ratio < TENTH)
    improve_geo = (.not. reduce_rho_1) .and. (maxval(xdist) > TWO * delta) .and. bad_trstep

    ! If all the interpolation points are close to XOPT and the trust-region is small, but the
    ! trust-region step is "bad" (SHORTD or RATIO <= 0), then we shrink RHO (update the criterion
    ! for the "closeness" and SHORTD). REDUCE_RHO_2 corresponds to Box 10 of the NEWUOA paper.
    ! N.B.:
    ! 1. The definition of REDUCE_RHO_2 is equivalent to the following:
    ! REDUCE_RHO_2 = (.NOT. IMPROVE_GEO) .AND. (MAX(DELTA, DNORM) <= RHO) .AND. BAD_TRSTEP
    ! 2. The definition of REDUCE_RHO_2 can be moved downward below IF (IMPROVE_GEO) ... END IF.
    ! Even though DNORM gets a new value after the geometry step when IMPROVE_GEO = TRUE, this
    ! value does not affect REDUCE_RHO_2, because DNORM comes into play only if IMPROVE_GEO = FALSE.
    ! 3. DELTA < DNORM may hold due to the update of DELTA.
    bad_trstep = (shortd .or. ratio <= ZERO)
    reduce_rho_2 = (maxval(xdist) <= TWO * delta) .and. (max(delta, dnorm) <= rho) .and. bad_trstep

    ! Comments on BAD_TRSTEP:
    ! 1. Powell used different thresholds (<= 0 and < 0.1) for RATIO in the definitions of BAD_TRSTEP
    ! above. Unifying them to <= 0 makes little difference to the performance, sometimes worsening,
    ! sometimes improving, but never substantially; unifying them to 0.1 makes little difference to
    ! the performance.
    ! 2. KNEW_TR == 0 implies RATIO < 0, and hence BAD_TRSTEP = TRUE. Otherwise, SETDROP_TR is buggy.

    ! NEWUOA never sets IMPROVE_GEO and (REDUCE_RHO_1 .OR. REDUCE_RHO_2) to TRUE simultaneously. So
    ! the instructions "IF (IMPROVE_GEO) ... END IF" and "IF (REDUCE_RHO_1 .OR. REDUCE_RHO_2)" can
    ! be exchanged without changing the algorithm.
    if (improve_geo) then
        ! XPT(:, KNEW_GEO) will be dropped (replaced by XOPT + D below).
        ! KNEW_GEO should never be KOPT. Otherwise, it is a bug.
        knew_geo = int(maxloc(xdist, dim=1), kind(knew_geo))

        ! Set DELBAR, which will be used as the trust-region radius for the geometry-improving
        ! scheme GEOSTEP. We also need it to decide whether to shift XBASE or not.
        ! Note that DELTA has been updated before arriving here. See the comments above the
        ! definition of IMPROVE_GEO.
        delbar = max(min(TENTH * maxval(xdist), HALF * delta), rho)

        ! Shift XBASE if XOPT may be too far from XBASE.
        if (delbar**2 <= 1.0E-3_RP * inprod(xopt, xopt)) then
            call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt)
        end if

        ! Find a step D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) is
        ! replaced by XOPT + D. The GEOSTEP subroutine will call Powell's BIGLAG and BIGDEN. It will
        ! also calculate the VLAG and BETA for this D.
        d = geostep(idz, knew_geo, kopt, bmat, delbar, xpt, zmat)

        ! Calculate VLAG and BETA for D. It makes uses of XOPT, so this is done before updating XOPT.
        call vlagbeta(idz, kopt, bmat, d, xpt, zmat, beta, vlag)

        ! Use the current quadratic model to predict the change in F due to the step D.
        qred = calquad(d, gq, hq, pq, xopt, xpt)

        ! Calculate the next value of the objective function.
        xnew = xopt + d
        x = xbase + xnew
        call evalf(calfun, x, f)
        nf = nf + 1_IK
        call fmssg(iprint, nf, f, x, solver)
        ! Save X and F into the history.
        call savehist(nf, f, x, fhist, xhist)
        ! Check whether to exit
        subinfo = checkexit(maxfun, nf, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        ! DNORMSAVE constains the DNORM of the latest 3 function evaluations with the current RHO.
        !------------------------------------------------------------------------------------------!
        ! Powell's code does not update DNORM. Therefore, DNORM is the length of last trust-region
        ! trial step, which seems inconsistent with what is described in Section 7 (around (7.7)) of
        ! the NEWUOA paper. Seemingly we should keep DNORM = ||D|| as we do here. The value of DNORM
        ! will be used when defining REDUCE_RHO.
        dnorm = min(delbar, norm(d))  ! In theory, DNORM = DELBAR in this case.
        !------------------------------------------------------------------------------------------!
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

        ! MODERR is the error of the current model in predicting the change in F due to D.
        moderr = f - fopt + qred
        ! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
        moderrsav = [moderrsav(2:size(moderrsav)), moderr]

        ! Update BMAT, ZMAT and IDZ, so that the KNEW_GEO-th interpolation point is replaced by
        ! XNEW, whose information is encoded in VLAG and BETA.
        call updateh(knew_geo, beta, vlag, idz, bmat, zmat)
        ! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
        call updateq(idz, knew_geo, bmat(:, knew_geo), moderr, zmat, xpt(:, knew_geo), gq, hq, pq)
        ! Include XNEW into XPT. Then update KOPT, XOPT, and FOPT.
        call updatexf(knew_geo, f, xnew, kopt, fval, xpt, fopt, xopt)
    end if  ! The procedure of improving geometry ends.

    if (reduce_rho_1 .or. reduce_rho_2) then
        ! The calculations with the current RHO are complete. Pick the next values of RHO and DELTA.
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        else
            call resenhance(rhoend, delta, rho)
            call rhomssg(iprint, nf, fopt, rho, xbase + xopt, solver)
            ! DNORMSAVE and MODERRSAVE are corresponding to the latest 3 function evaluations with
            ! the current RHO. Update them after reducing RHO.
            dnormsav = HUGENUM
            moderrsav = HUGENUM
        end if
    end if  ! The procedure of reducing RHO ends.

end do  ! The iterative procedure ends.

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xbase + (xopt + d)
    call evalf(calfun, x, f)
    nf = nf + 1_IK
    call fmssg(iprint, nf, f, x, solver)
    ! Save X and F into the history.
    call savehist(nf, f, x, fhist, xhist)
end if

! Choose the [X, F] to return: either the current [X, F] or [XBASE+XOPT, FOPT].
if (is_nan(f) .or. fopt < f) then
    x = xbase + xopt
    f = fopt
end if

! Arrange FHIST and XHIST so that they are in the chronological order.
call rangehist(nf, fhist, xhist)

call retmssg(info, iprint, nf, f, x, solver)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    call assert(.not. any(fhist(1:min(nf, maxfhist)) < f), 'F is the smallest in FHIST', srname)
end if

end subroutine newuob


end module newuob_mod
