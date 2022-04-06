module newuob_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of NEWUOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Tuesday, April 05, 2022 PM05:29:29
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
! ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
! ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT-N-1)=1. BMAT will hold the last N ROWs of H
! except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th are not saved as
! they are unnecessary for the calculation.
!
! See Section 2 of the NEWUOA paper for more information about these variables.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, HALF, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: info_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS
use, non_intrinsic :: linalg_mod, only : calquad, norm
use, non_intrinsic :: output_mod, only : retmsg, rhomsg, fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : setdrop_tr, geostep
use, non_intrinsic :: initialize_mod, only : initxf, initq, inith
use, non_intrinsic :: shiftbase_mod, only : shiftbase
use, non_intrinsic :: trustregion_mod, only : trsapp, trrad
use, non_intrinsic :: update_mod, only : updateh, updateq, updatexf, tryqalt

implicit none

! Inputs
procedure(OBJ) :: calfun
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
logical :: tr_success
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
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pq(npt)
real(RP) :: qred
real(RP) :: ratio
real(RP) :: rho
real(RP) :: xbase(size(x))
real(RP) :: xdist(npt)
real(RP) :: xopt(size(x))
real(RP) :: xpt(size(x), npt)
real(RP) :: zmat(npt, npt - size(x) - 1)
real(RP), parameter :: tr_tol = 1.0E-2_RP  ! Tolerance used in TRSAPP.

! Sizes
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', srname)
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
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
    call rangehist(nf, xhist, fhist)
    call retmsg(solver, info, iprint, nf, f, x)
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
! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
ratio = -ONE
knew_tr = 0_IK
! No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
shortd = .false.

! MAXTR is the maximal number of trust-region iterations. Each trust-region iteration takes 1 or 2
! function evaluations unless the trust-region step is short but the geometry step is not invoked.
! Thus the following MAXTR is unlikely to be reached.
maxtr = max(maxfun, 2_IK * maxfun)  ! MAX: precaution against overflow, which will make 2*MAXFUN < 0.
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
    call trsapp(delta, gq, hq, pq, tr_tol, xopt, xpt, crvmin, d)
    dnorm = min(delta, norm(d))

    ! SHORTD corresponds to Box 3 of the NEWUOA paper. N.B.: we compare DNORM with RHO, not DELTA.
    shortd = (dnorm < HALF * rho)
    ! REDUCE_RHO_1 corresponds to Box 14 of the NEWUOA paper.
    reduce_rho_1 = shortd .and. (maxval(abs(moderrsav)) <= 0.125_RP * crvmin * rho**2) .and. &
        & (maxval(dnormsav) <= rho)
    if (shortd .and. (.not. reduce_rho_1)) then
        ! Reduce DELTA. After this, DELTA < DNORM may happen.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close.
        end if
    end if

    if (.not. shortd) then  ! D is long enough.
        ! DNORMSAVE contains the DNORM of the latest 3 function evaluations with the current RHO.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

        ! Shift XBASE if XOPT may be too far from XBASE.
        !if (inprod(d, d) <= 1.0e-3_RP*inprod(xopt, xopt)) then  ! Powell's code
        if (sum(xopt**2) >= 1.0E3_RP * dnorm**2) then
            call shiftbase(xbase, xopt, xpt, idz, zmat, bmat, pq, hq, gq)
        end if

        ! Calculate the next value of the objective function.
        x = xbase + (xopt + d)
        call evaluate(calfun, x, f)
        nf = nf + 1_IK
        call fmsg(solver, iprint, nf, f, x)
        ! Save X and F into the history.
        call savehist(nf, x, xhist, f, fhist)
        ! Check whether to exit
        subinfo = checkexit(maxfun, nf, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        qred = calquad(d, gq, hq, pq, xopt, xpt)  ! QRED = Q(XOPT) - Q(XOPT + D)
        ! F - FOPT + QRED is the error of the current model in predicting the change in F due to D.
        ! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
        moderrsav = [moderrsav(2:size(moderrsav)), f - fopt + qred]

        ! Calculate the reduction ratio.
        !------------------------------------------------------------------------------------------!
        ! Zaikun 20220405: REDRAT sets returns a large negative value when QRED nonnegative or NaN.
        ! This ratio will lead to a contraction of DELTA and make IMPROVE_GEO or REDUCE_RHO true.
        ! Is there a better strategy? LINCOA does something to improve the model. Applicable here?
        !------------------------------------------------------------------------------------------!
        ratio = redrat(fopt - f, qred, eta1)
        ! Update DELTA. After this, DELTA < DNORM may hold.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if

        ! Set KNEW_TR to the index of the interpolation point to be replaced by XNEW = XOPT + D.
        ! KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
        ! N.B.:
        ! 1. KNEW_TR = 0 means it is impossible to obtain a good interpolation set by replacing any
        ! current interpolation point with XNEW. Then XNEW and its function value will be discarded.
        ! In this case, the geometry of XPT likely needs improvement, which will be handled below.
        ! 2. If TR_SUCCESS = TRUE (i.e., RATIO > 0), then SETDROP_TR ensures KNEW_TR > 0 so that
        ! XNEW is included into XPT. Otherwise, SETDROP_TR is buggy; moreover, if TR_SUCCESS = TRUE
        ! and KNEW_TR = 0, XOPT will differ from XPT(:, KOPT), because the former is set to XNEW but
        ! XNEW is discarded. Such a difference can lead to unexpected behaviors; for example,
        ! KNEW_GEO may equal KOPT, with which GEOSTEP will not work.
        tr_success = (f < fopt)
        knew_tr = setdrop_tr(idz, kopt, tr_success, bmat, d, delta, rho, xpt, zmat)

        ! Update BMAT, ZMAT, IDZ (corresponding to H is the NEWUOA paper), GQ, HQ, PQ (defining the
        ! quadratic model), and FVAL, XPT, KOPT, FOPT, XOPT so that XPT(:, KNEW_TR) becomes XOPT + D.
        ! If KNEW_TR = 0, the updating subroutines will do essentially nothing, as the algorithm
        ! decides not to include XNEW into XPT.
        ! Update BMAT, ZMAT and IDZ.
        call updateh(knew_tr, kopt, idz, d, xpt, bmat, zmat)
        ! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
        call updateq(idz, knew_tr, kopt, bmat, d, f, fval, xpt, zmat, gq, hq, pq)
        ! Update XPT(:, KNEW_TR) to XOPT + D. Then update KOPT, XOPT, and FOPT.
        call updatexf(knew_tr, d, f, kopt, fval, xpt, fopt, xopt)

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
            call tryqalt(idz, fval - fopt, ratio, bmat, zmat, itest, gq, hq, pq)
        end if
    end if  ! End of IF (.NOT. SHORTD). The normal trust-region calculation ends here.


    !----------------------------------------------------------------------------------------------!
    ! Before the next trust-region iteration, we may improve the geometry of XPT or reduce RHO
    ! according to IMPROVE_GEO and REDUCE_RHO_1/2. Now we decide these indicators.
    !----------------------------------------------------------------------------------------------!

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
    !!MATLAB: xdist = sqrt(sum((xpt - xopt).^2))  % xopt should be a column!! Implicit expansion
    ! KNEW_TR == 0 implies RATIO <= 0. Therefore, we can remove KNEW_TR == 0 from the definition
    ! of BAD_TRSTEP. Nevertheless, we keep it for robustness.
    bad_trstep = (shortd .or. ratio < TENTH .or. knew_tr == 0)  ! BAD_TRSTEP for IMPROVE_GEO
    improve_geo = ((.not. reduce_rho_1) .and. maxval(xdist) > TWO * delta .and. bad_trstep)

    ! If all the interpolation points are close to XOPT and the trust-region is small, but the
    ! trust-region step is "bad" (SHORTD or RATIO <= 0), then we shrink RHO (update the criterion
    ! for the "closeness" and SHORTD). REDUCE_RHO_2 corresponds to Box 10 of the NEWUOA paper.
    ! N.B.:
    ! 0. KNEW_TR == 0 implies RATIO <= 0. Therefore, we can remove KNEW_TR == 0 from the definition
    ! of BAD_TRSTEP. Nevertheless, we keep it for robustness.
    ! 1. The definition of REDUCE_RHO_2 is equivalent to the following:
    ! REDUCE_RHO_2 = (.NOT. IMPROVE_GEO) .AND. (MAX(DELTA, DNORM) <= RHO) .AND. BAD_TRSTEP
    ! 2. The definition of REDUCE_RHO_2 can be moved downward below IF (IMPROVE_GEO) ... END IF.
    ! Even though DNORM gets a new value after the geometry step when IMPROVE_GEO = TRUE, this
    ! value does not affect REDUCE_RHO_2, because DNORM comes into play only if IMPROVE_GEO = FALSE.
    ! 3. DELTA < DNORM may hold due to the update of DELTA.
    bad_trstep = (shortd .or. ratio <= 0 .or. knew_tr == 0)  ! BAD_TRSTEP for REDUCE_RHO_2
    reduce_rho_2 = (maxval(xdist) <= TWO * delta .and. max(delta, dnorm) <= rho .and. bad_trstep)

    ! Comments on BAD_TRSTEP:
    ! 1. Powell used different thresholds (<= 0 and < 0.1) for RATIO in the definitions of BAD_TRSTEP
    ! above. Unifying them to <= 0 makes little difference to the performance, sometimes worsening,
    ! sometimes improving, but never substantially; unifying them to 0.1 makes little difference to
    ! the performance either.
    ! Update 20220204: In the current version, unifying the two thresholds to 0 seems to worsen the
    ! performance on noise-free CUTEst problems with at most 200 variables; unifying them to 0.1
    ! worsens it a bit as well.
    ! 2. KNEW_TR == 0 implies RATIO <= 0, and hence BAD_TRSTEP = TRUE. Otherwise, SETDROP_TR is buggy.
    ! Indeed, we can remove KNEW_TR == 0 from the definition of BAD_TRSTEP. It is kept for robustness.

    !----------------------------------------------------------------------------------------------!
    ! N.B.: NEWUOA never sets IMPROVE_GEO and (REDUCE_RHO_1 .OR. REDUCE_RHO_2) to TRUE
    ! simultaneously. Thus following two blocks are exchangeable:
    !!IF (IMPROVE_GEO) ... END IF
    !!IF (REDUCE_RHO_1 .OR. REDUCE_RHO_2) ... END IF
    !----------------------------------------------------------------------------------------------!

    ! Improve the geometry of the interpolation set by removing a point and adding a new one.
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
        if (sum(xopt**2) >= 1.0E3_RP * delbar**2) then
            call shiftbase(xbase, xopt, xpt, idz, zmat, bmat, pq, hq, gq)
        end if

        ! Find a step D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) is
        ! replaced by XOPT + D. The GEOSTEP subroutine will call Powell's BIGLAG and BIGDEN.
        d = geostep(idz, knew_geo, kopt, bmat, delbar, xpt, zmat)
        !------------------------------------------------------------------------------------------!
        ! Powell's code does not update DNORM. Therefore, DNORM is the length of last trust-region
        ! trial step, which seems inconsistent with what is described in Section 7 (around (7.7)) of
        ! the NEWUOA paper. Seemingly we should keep DNORM = ||D|| as we do here. The value of DNORM
        ! saved in DNORMSAVE will be used when defining REDUCE_RHO_1.
        dnorm = min(delbar, norm(d))  ! In theory, DNORM = DELBAR in this case.
        !------------------------------------------------------------------------------------------!
        ! DNORMSAVE contains the DNORM of the latest 3 function evaluations with the current RHO.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

        ! Calculate the next value of the objective function.
        x = xbase + (xopt + d)
        call evaluate(calfun, x, f)
        nf = nf + 1_IK
        call fmsg(solver, iprint, nf, f, x)
        ! Save X and F into the history.
        call savehist(nf, x, xhist, f, fhist)
        ! Check whether to exit
        subinfo = checkexit(maxfun, nf, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        qred = calquad(d, gq, hq, pq, xopt, xpt)  ! QRED = Q(XOPT) - Q(XOPT + D)
        ! F - FOPT + QRED is the error of the current model in predicting the change in F due to D.
        ! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
        moderrsav = [moderrsav(2:size(moderrsav)), f - fopt + qred]

        ! Update BMAT, ZMAT, IDZ (corresponding to H is the NEWUOA paper), GQ, HQ, PQ (defining the
        ! quadratic model), and FVAL, XPT, KOPT, FOPT, XOPT so that XPT(:, KNEW_TR) becomes XOPT + D.
        ! Update BMAT, ZMAT and IDZ.
        call updateh(knew_geo, kopt, idz, d, xpt, bmat, zmat)
        ! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
        call updateq(idz, knew_geo, kopt, bmat, d, f, fval, xpt, zmat, gq, hq, pq)
        ! Update XPT(:, KNEW_GEO) to XOPT + D. Then update KOPT, XOPT, and FOPT.
        call updatexf(knew_geo, d, f, kopt, fval, xpt, fopt, xopt)
    end if  ! End of IF (IMPROVE_GEO). The procedure of improving geometry ends.

    ! The calculations with the current RHO are complete. Enhance the resolution of the algorithm
    ! by reducing RHO; update DELTA at the same time.
    if (reduce_rho_1 .or. reduce_rho_2) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        else
            delta = HALF * rho
            rho = redrho(rho, rhoend)
            delta = max(delta, rho)
            call rhomsg(solver, iprint, nf, fopt, rho, xbase + xopt)
            ! DNORMSAVE and MODERRSAVE are corresponding to the latest 3 function evaluations with
            ! the current RHO. Update them after reducing RHO.
            dnormsav = HUGENUM
            moderrsav = HUGENUM
        end if
    end if  ! End of IF (REDUCE_RHO_1 .OR. REDUCE_RHO_2). The procedure of reducing RHO ends.

end do  ! End of Do TR = 1, MAXTR. The iterative procedure ends.

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xbase + (xopt + d)
    call evaluate(calfun, x, f)
    nf = nf + 1_IK
    call fmsg(solver, iprint, nf, f, x)
    ! Save X and F into the history.
    call savehist(nf, x, xhist, f, fhist)
end if

! Choose the [X, F] to return: either the current [X, F] or [XBASE+XOPT, FOPT].
if (is_nan(f) .or. fopt < f) then
    x = xbase + xopt
    f = fopt
end if

! Arrange FHIST and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist)

call retmsg(solver, info, iprint, nf, f, x)

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
    call assert(.not. any(fhist(1:min(nf, maxfhist)) < f), 'F is the smallest in FHIST', srname)
end if

end subroutine newuob


end module newuob_mod
