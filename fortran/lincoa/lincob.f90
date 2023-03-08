!TODO:
! 1. Check whether it is possible to change the definition of RESCON, RESNEW, RESTMP, RESACT so that
! we do not need to encode information into their signs.
!
module lincob_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, March 09, 2023 AM12:06:55
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lincob


contains


subroutine lincob(calfun, iprint, maxfilt, maxfun, npt, A_orig, amat, b_orig, bvec, ctol, cweight, &
    & eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, x, nf, chist, cstrv, f, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the actual calculations of LINCOA.
!
! The arguments IPRINT, MAXFILT, MAXFUN, MAXHIST, NPT, CTOL, CWEIGHT, ETA1, ETA2, FTARGET, GAMMA1,
! GAMMA2, RHOBEG, RHOEND, X, NF, F, XHIST, FHIST, CHIST, CSTRV and INFO are identical to the
! corresponding arguments in subroutine LINCOA.
! AMAT is a matrix whose columns are the constraint gradients, scaled so that they have unit length.
! B contains on entry the right hand sides of the constraints, scaled as above, but later B is
!   modified for variables relative to XBASE.
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
!   of the model and Lagrange functions.
! XOPT is the displacement from XBASE of the feasible vector of variables that provides the least
!   calculated F so far, this vector being the current trust region centre. FOPT = F(XOPT + XBASE).
!   However, we do not save XOPT and FOPT explicitly, because XOPT = XPT(:, KOPT) and
!   FOPT = FVAL(KOPT), which is explained below.
! [XPT, FVAL, KOPT] describes the interpolation set:
! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
!   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
! [GOPT, HQ, PQ] describes the quadratic model: GOPT will hold the gradient of the quadratic model
!   at XBASE + XOPT; HQ will hold the explicit second order derivatives of the quadratic model; PQ
!   will contain the parameters of the implicit second order derivatives of the quadratic model.
! [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
!   the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
!   ZMAT will hold a factorization of the leading NPT*NPT submatrix of H, the factorization being
!   ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT-N-1)=1. BMAT will hold the last N ROWs of H
!   except for the (NPT+1)th column. Note that the (NPT + 1)th row and column of H are not saved as
!   they are unnecessary for the calculation.
! D is reserved for trial steps from XOPT. It is chosen by subroutine TRSTEP or GEOSTEP. Usually
!   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
! IACT is an integer array for the indices of the active constraints.
! RESCON holds information about the constraint residuals at the current trust region center XOPT.
!   1. If if B(J) - AMAT(:, J)^T*XOPT <= DELTA, then RESCON(J) = B(J) - AMAT(:, J)^T*XOPT. Note that
!   RESCON >= 0 in this case, because the algorithm keeps XOPT to be feasible.
!   2. Otherwise, RESCON(J) is a negative value that B(J) - AMAT(:,J)^T*XOPT >= |RESCON(J)| >= DELTA.
!   RESCON can be updated without calculating the constraints that are far from being active, so
!   that we only need to evaluate the constraints that are nearly active.
! QFAC is the orthogonal part of the QR factorization of the matrix of active constraint gradients,
!   these gradients being ordered in accordance with IACT. When NACT is less than N, columns are
!   addedto QFAC to complete an N by N orthogonal matrix, which is important for keeping calculated
!   steps sufficiently close to the boundaries of the active constraints.
! RFAC is the upper triangular part of this QR factorization.
!--------------------------------------------------------------------------------------------------!

! Generic models
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TENTH, REALMAX, MIN_MAXFILT, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: infos_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS
use, non_intrinsic :: linalg_mod, only : matprod, maximum, eye, trueloc, linspace, norm
use, non_intrinsic :: output_mod, only : fmsg, rhomsg, retmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, omega_mul, hess_mul
use, non_intrinsic :: powalg_mod, only : updateh
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: selectx_mod, only : savefilt, selectx, isbetter
use, non_intrinsic :: shiftbase_mod, only : shiftbase

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep, setdrop_tr
use, non_intrinsic :: initialize_mod, only : initxf, inith
use, non_intrinsic :: trustregion_mod, only : trstep, trrad
use, non_intrinsic :: update_mod, only : updatexf, updateq, tryqalt, updateres

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfilt
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: npt
real(RP), intent(in) :: A_orig(:, :)  ! A_ORIG(N, M) ; Better names? necessary?
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M) ; Better names? necessary?
real(RP), intent(in) :: b_orig(:) ! B_ORIG(M) ; Better names? necessary?
real(RP), intent(in) :: bvec(:)  ! BVEC(M) ; Better names? necessary?
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
real(RP), intent(inout) :: x(:)  ! X(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)  ! CHIST(MAXCHIST)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'LINCOB'
integer(IK) :: iact(size(b_orig))
integer(IK) :: idz
integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK)))
integer(IK) :: k
integer(IK) :: knew_geo
integer(IK) :: knew_tr
integer(IK) :: kopt
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxtr
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: nact
integer(IK) :: nfilt
integer(IK) :: ngetact
integer(IK) :: nhist
integer(IK) :: subinfo
integer(IK) :: tr
logical :: accurate_mod
logical :: adequate_geo
logical :: bad_trstep
logical :: close_itpset
logical :: evaluated(npt)
logical :: feasible
logical :: improve_geo
logical :: qalt_better(3)
logical :: reduce_rho
logical :: shortd
logical :: small_trrad
logical :: ximproved
real(RP) :: b(size(b_orig))
real(RP) :: bmat(size(x), npt + size(x))
real(RP) :: cfilt(maxfilt)
real(RP) :: constr(size(b_orig))
real(RP) :: d(size(x))
real(RP) :: delbar
real(RP) :: delta
real(RP) :: distsq(npt)
real(RP) :: dnorm
real(RP) :: dnormsav(4)  ! Powell's implementation uses 5
real(RP) :: ffilt(maxfilt)
real(RP) :: fval(npt), cval(npt)
real(RP) :: galt(size(x))
real(RP) :: gopt(size(x))
real(RP) :: hq(size(x), size(x))
real(RP) :: moderr
real(RP) :: moderr_alt
real(RP) :: pq(npt)
real(RP) :: pqalt(npt)
real(RP) :: qfac(size(x), size(x))
real(RP) :: qred
real(RP) :: ratio
real(RP) :: rescon(size(b_orig))
real(RP) :: rfac(size(x), size(x))
real(RP) :: rho
real(RP) :: xbase(size(x))
real(RP) :: xdrop(size(x))
real(RP) :: xfilt(size(x), maxfilt)
real(RP) :: xosav(size(x))
real(RP) :: xpt(size(x), npt)
real(RP) :: zmat(npt, npt - size(x) - 1)

! Sizes.
m = int(size(b_orig), kind(m))
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = int(max(maxxhist, maxfhist, maxchist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT+1', srname)
    call assert(size(A_orig, 1) == n .and. size(A_orig, 2) == m, 'SIZE(A_ORIG) == [N, M]', srname)
    call assert(size(bvec) == m, 'SIZE(BVEC) == M', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    call assert(maxfilt >= min(MIN_MAXFILT, maxfun) .and. maxfilt <= maxfun, &
        & 'MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize XBASE, XPT, FVAL, and KOPT.
b = bvec
call initxf(calfun, iprint, maxfun, A_orig, amat, b_orig, ctol, ftarget, rhobeg, x, b, &
    & ij, kopt, nf, chist, cval, fhist, fval, xbase, xhist, xpt, evaluated, subinfo)
x = xbase + xpt(:, kopt)
f = fval(kopt)

! Evaluate the constraints using A_ORIG and B_ORIG. Should we do this in INITXF?
! N.B.: We must initialize CONSTR and CSTRV. Otherwise, if REDUCE_RHO is TRUE after the very first
! iteration due to SHORTD, then RHOMSG will be called with CONSTR and CSTRV uninitialized.
constr = matprod(x, A_orig) - b_orig
cstrv = maximum([ZERO, constr])

! Initialize the filter, including XFILT, FFILT, CONFILT, CFILT, and NFILT.
! N.B.: The filter is used only when selecting which iterate to return. It does not interfere with
! the iterations. LINCOA is NOT a filter method but a trust-region method. All the trust-region
! iterates are supposed to be feasible, but can be infeasible due to rounding errors; the
! geometry-improving iterates are not necessarily feasible. Powell's implementation does not use a
! filter to select the iterate, possibly returning a suboptimal iterate.
nfilt = 0
do k = 1, npt
    if (evaluated(k)) then
        call savefilt(cval(k), ctol, cweight, fval(k), xbase + xpt(:, k), nfilt, cfilt, ffilt, xfilt)
    end if
end do

! Check whether to return due to abnormal cases that may occur during the initialization.
if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Return the best calculated values of the variables. If CTOL > 0, the KOPT decided by SELECTX
    ! may not be the same as the one by INITXF.
    kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol)
    x = xfilt(:, kopt)
    f = ffilt(kopt)
    constr = matprod(x, A_orig) - b_orig
    cstrv = maximum([ZERO, constr])
    ! Arrange CHIST, FHIST, and XHIST so that they are in the chronological order.
    call rangehist(nf, xhist, fhist, chist)
    call retmsg(solver, info, iprint, nf, f, x, cstrv, constr)
    return
end if

! Initialize BMAT, ZMAT, and IDZ.
call inith(ij, xpt, idz, bmat, zmat)

! Initialize the quadratic represented by [GOPT, HQ, PQ], so that its gradient at XBASE+XOPT is
! GOPT; its Hessian is HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.
hq = ZERO
pq = omega_mul(idz, zmat, fval)
gopt = matprod(bmat(:, 1:npt), fval) + hess_mul(xpt(:, kopt), xpt, pq)
pqalt = pq
galt = gopt

! Initialize RESCON.
rescon = max(b - matprod(xpt(:, kopt), amat), ZERO)
rescon(trueloc(rescon >= rhobeg)) = -rescon(trueloc(rescon >= rhobeg))
!!MATLAB: rescon(rescon >= rhobeg) = -rescon(rescon >= rhobeg)

! Set some more initial values.
! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
! No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
rho = rhobeg
delta = rho
ratio = -ONE
dnormsav = REALMAX
shortd = .false.
qalt_better = .false.
knew_tr = 0
knew_geo = 0
qfac = eye(n)
rfac = ZERO
nact = 0
iact = linspace(1_IK, m, m)

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
! LINCOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Generate the next trust region step D by calling TRSTEP. Note that D is feasible.
    call trstep(amat, delta, gopt, hq, pq, rescon, xpt, iact, nact, qfac, rfac, d, ngetact)
    dnorm = min(delta, norm(d))

    ! A trust region step is applied whenever its length is at least 0.5*DELTA. It is also
    ! applied if its length is at least 0.1999*DELTA and if a line search of TRSTEP has caused a
    ! change to the active set, indicated by NGETACT >= 2 (note that NGETACT is at least 1).
    ! Otherwise, the trust region step is considered too short to try.
    ! N.B. The magic number 0.1999 seems to be related to the fact that a linear constraint is
    ! considered nearly active if the point under consideration is within 0.2*DELTA to the boundary
    ! of the constraint. See the subroutine GETACT and Section 3 of Powell (2015) for more details.
    shortd = ((dnorm < HALF * delta .and. ngetact < 2) .or. dnorm < 0.1999_RP * delta)
    !------------------------------------------------------------------------------------------!
    ! The SHORTD defined above needs NGETACT, which relies on Powell's trust region subproblem
    ! solver. If a different subproblem solver is used, we can take the following SHORTD adopted
    ! from UOBYQA, NEWUOA and BOBYQA.
    ! !SHORTD = (DNORM < HALF * RHO)
    !------------------------------------------------------------------------------------------!

    ! DNORMSAV saves the DNORM of last few (five) trust-region iterations. It will be used to
    ! decide whether we should improve the geometry of the interpolation set or reduce RHO when
    ! SHORTD is TRUE. Note that it does not record the geometry steps.
    dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

    ! In some cases, we reset DNORMSAV to REALMAX. This indicates a preference of improving the
    ! geometry of the interpolation set to reducing RHO in the subsequent three or more
    ! iterations. This is important for the performance of LINCOA.
    if (delta > rho .or. .not. shortd) then  ! Another possibility: IF (DELTA > RHO) THEN
        dnormsav = REALMAX
    end if

    ! Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
    ! should be positive If it is nonpositive due to rounding errors, we will not take this step.
    qred = -quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)

    if (shortd .or. .not. qred > 0) then
        ! In this case, do nothing but reducing DELTA. Afterward, DELTA < DNORM may occur.
        ! N.B.: 1. This value of DELTA will be discarded if REDUCE_RHO turns out TRUE later.
        ! 2. Powell's code does not shrink DELTA when QRED > 0 is FALSE (i.e., when VQUAD >= 0 in
        ! Powell's code, where VQUAD = -QRED). Consequently, the algorithm may be stuck in an
        ! infinite cycling, because both REDUCE_RHO and IMPROVE_GEO may end up with FALSE in this
        ! case, which did happen in tests.
        ! 3. The factor HALF works better than TENTH (used in NEWUOA/BOBYQA), 0.2, and 0.7.
        delta = HALF * delta
        ! Powell's code uses DELTA <= MIN(0.99_RP * GAMMA3, 1.5_RP) * RHO as the condition for the
        ! following IF, where GAMMA3 is a parameter in (1, GAMMA2] for updating DELTA. This aligns
        ! with the update of DELTA after a trust-region step. See the comment there for more info.
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if
    else
        ! Calculate the next value of the objective function.
        x = xbase + (xpt(:, kopt) + d)
        call evaluate(calfun, x, f)
        nf = nf + 1_IK

        ! Evaluate the constraints using A_ORIG and B_ORIG.
        constr = matprod(x, A_orig) - b_orig
        cstrv = maximum([ZERO, constr])

        ! Print a message about the function evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x, cstrv, constr)
        ! Save X, F, CSTRV into the history.
        call savehist(nf, x, xhist, f, fhist, cstrv, chist)
        ! Save X, F, CSTRV into the filter.
        call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt)

        ! Check whether to exit.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        ! QALT_BETTER is a boolean array indicating whether the recent few (three) alternative
        ! models are more accurate in predicting the function value at XOPT + D.
        ! N.B.: Do NOT change the "<" in the comparison to "<="; otherwise, the result will not be
        ! reasonable if the two values being compared are both ZERO or INF.
        moderr = f - fval(kopt) + qred
        moderr_alt = f - fval(kopt) - quadinc(d, xpt, galt, pqalt)
        qalt_better = [qalt_better(2:size(qalt_better)), abs(moderr_alt) < TENTH * abs(moderr)]

        ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
        ratio = redrat(fval(kopt) - f, qred, eta1)

        ! Update DELTA. After this, DELTA < DNORM may hold.
        ! The new DELTA lies in [GAMMA1*DNORM, GAMMA2*DNORM].
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        ! Set DELTA to RHO when it is close to or below.
        ! N.B.: Powell's code uses DELTA <= MIN(0.99_RP * GAMMA3, 1.5_RP) * RHO as the condition for
        ! the following IF, where GAMMA3 is a parameter in (1, GAMMA2] for updating DELTA. Powell's
        ! code must make sure that the multiplicative factor in the condition is less than GAMMA3.
        ! Imagine a very successful step with DENORM = the un-updated DELTA = RHO. Then, in Powell's
        ! implementation, TRRAD will update DELTA to GAMMA3*RHO. If this factor were not less than
        ! GAMMA3, then DELTA will be reset to RHO, which is not reasonable as D is very successful.
        ! See paragraph two of Sec. 5.2.5 in T. M. Ragonneau's thesis: "Model-Based Derivative-Free
        ! Optimization Methods and Software".
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if

        ! Is the newly generated X better than current best point?
        ximproved = (f < fval(kopt))

        ! Set KNEW_TR to the index of the interpolation point to be replaced with XNEW = XOPT + D.
        ! KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
        knew_tr = setdrop_tr(idz, kopt, ximproved, bmat, d, delta, rho, xpt, zmat)
        if (knew_tr > 0) then
            ! Update [BMAT, ZMAT, IDZ] (represents H in the NEWUOA paper), [XPT, FVAL, KOPT] and
            ! [GOPT, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_TR) becomes XNEW = XOPT + D.
            xdrop = xpt(:, knew_tr)
            xosav = xpt(:, kopt)
            call updateh(knew_tr, kopt, d, xpt, idz, bmat, zmat)
            call updatexf(knew_tr, ximproved, f, xosav + d, kopt, fval, xpt)
            call updateq(idz, knew_tr, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)

            ! Establish the alternative model, namely the least Frobenius norm interpolant. Replace
            ! the current model with the alternative model if the recent few (three) alternative
            ! models are more accurate in predicting the function value of XOPT + D.
            call tryqalt(idz, bmat, fval - fval(kopt), xpt(:, kopt), xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt)

            ! Update RESCON if XOPT is changed.
            ! Zaikun 20221115: Shouldn't we do it after DELTA is updated?
            call updateres(ximproved, amat, b, delta, norm(d), xpt(:, kopt), rescon)
        end if

    end if  ! End of IF (SHORTD .OR. .NOT. QRED > 0). The normal trust-region calculation ends here.

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
    ! N.B.: The ACCURATE_MOD here plays a similar role as the variable with the same name in UOBYQA,
    ! NEWUOA, and BOBYQA. However, the definition of ACCURATE_MOD here is different from that in
    ! those solvers, which do not only check whether DNORM is small in recent iterations, but also
    ! verify a curvature condition that really indicates that recent models are sufficiently
    ! accurate. Here, however, we are not really sure whether they are accurate or not. Therefore,
    ! ACCURATE_MOD is not the best name, but we keep it to align with the other solvers.
    accurate_mod = all(dnormsav <= HALF * rho) .or. all(dnormsav(3:size(dnormsav)) <= 0.2 * rho)
    ! Powell's version (note that size(dnormsav) = 5 in his implementation):
    !accurate_mod = all(dnormsav <= HALF * rho) .or. all(dnormsav(3:size(dnormsav)) <= TENTH * rho)
    ! CLOSE_ITPSET: Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
    close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Behaves the same as Powell's version.
    ! Below are some alternative definitions of CLOSE_ITPSET.
    ! !close_itpset = all(distsq <= max(delta**2, 4.0_RP * rho**2))  ! Powell's code.
    ! !close_itpset = all(distsq <= delta**2)  ! Does not work as well as Powell's version.
    ! !close_itpset = all(distsq <= 10.0_RP * delta**2)  ! Does not work as well as Powell's version.
    ! !close_itpset = all(distsq <= max((2.0_RP * delta)**2, (10.0_RP * rho)**2))  ! Powell's BOBYQA.
    ! ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
    adequate_geo = (shortd .and. accurate_mod) .or. close_itpset
    ! SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
    small_trrad = (max(delta, dnorm) <= rho)  ! Behaves the same as Powell's version.
    !small_trrad = (delsav <= rho)  ! Powell's code. DELSAV = unupdated DELTA.

    ! IMPROVE_GEO and REDUCE_RHO are defined as follows.
    ! N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.

    ! BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= eta1 .or. knew_tr == 0)
    improve_geo = bad_trstep .and. .not. adequate_geo
    ! BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    reduce_rho = bad_trstep .and. adequate_geo .and. small_trrad

    ! Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    ! !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

    ! With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= eta1 .or. knew_tr == 0)
    ! !improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

    ! With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)
    ! !reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

    ! LINCOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    !call assert(.not. (improve_geo .and. reduce_rho), 'IMPROVE_GEO and REDUCE_RHO are not both TRUE', srname)
    !
    ! If SHORTD is TRUE or QRED > 0 is FALSE, then either REDUCE_RHO or IMPROVE_GEO is TRUE unless
    ! CLOSE_ITPSET is TRUE but SMALL_TRRAD is FALSE.
    !call assert((.not. shortd .and. qred > 0) .or. (improve_geo .or. reduce_rho .or. &
    !    & (close_itpset .and. .not. small_trrad)), 'If SHORTD is TRUE or QRED > 0 is FALSE, then either&
    !    & IMPROVE_GEO or REDUCE_RHO is TRUE unless CLOSE_ITPSET is TRUE but SMALL_TRRAD is FALSE', srname)
    !----------------------------------------------------------------------------------------------!


    ! Since IMPROVE_GEO and REDUCE_RHO are never TRUE simultaneously, the following two blocks are
    ! exchangeable: IF (IMPROVE_GEO) ... END IF and IF (REDUCE_RHO) ... END IF.

    if (improve_geo) then
        ! XPT(:, KNEW_GEO) will become  XOPT + D below. KNEW_GEO /= KOPT unless there is a bug.
        knew_geo = int(maxloc(distsq, dim=1), kind(knew_geo))

        ! Set DELBAR, which will be used as the trust-region radius for the geometry-improving
        ! scheme GEOSTEP. Note that DELTA has been updated before arriving here.
        delbar = max(TENTH * delta, rho)  ! This differs from NEWUOA/BOBYQA. Possible improvement?
        ! Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
        call geostep(iact, idz, knew_geo, kopt, nact, amat, bmat, delbar, qfac, rescon, xpt, zmat, feasible, d)

        ! Calculate the next value of the objective function.
        x = xbase + (xpt(:, kopt) + d)
        call evaluate(calfun, x, f)
        nf = nf + 1_IK

        ! Evaluate the constraints using A_ORIG and B_ORIG.
        constr = matprod(x, A_orig) - b_orig
        cstrv = maximum([ZERO, constr])

        ! Print a message about the function evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x, cstrv, constr)
        ! Save X, F, CSTRV into the history.
        call savehist(nf, x, xhist, f, fhist, cstrv, chist)
        ! Save X, F, CSTRV into the filter.
        call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt)

        ! Check whether to exit.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        ! QALT_BETTER is a boolean array indicating whether the recent few (three) alternative
        ! models are more accurate in predicting the function value at XOPT + D.
        ! Powell's code takes XOPT + D into account only if it is feasible.
        ! N.B.: Do NOT change the "<" in the comparison to "<="; otherwise, the result will not be
        ! reasonable if the two values being compared are both ZERO or INF.
        moderr = f - fval(kopt) - quadinc(d, xpt, gopt, pq, hq)
        moderr_alt = f - fval(kopt) - quadinc(d, xpt, galt, pqalt)
        qalt_better = [qalt_better(2:size(qalt_better)), abs(moderr_alt) < TENTH * abs(moderr)]

        ! Is the newly generated X better than current best point?
        ximproved = (f < fval(kopt) .and. feasible)

        ! Update [BMAT, ZMAT, IDZ] (represents H in the NEWUOA paper), [XPT, FVAL, KOPT] and
        ! [GOPT, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XNEW = XOPT + D.
        xdrop = xpt(:, knew_geo)
        xosav = xpt(:, kopt)
        call updateh(knew_geo, kopt, d, xpt, idz, bmat, zmat)
        call updatexf(knew_geo, ximproved, f, xosav + d, kopt, fval, xpt)
        call updateq(idz, knew_geo, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)

        ! Establish the alternative model, namely the least Frobenius norm interpolant. Replace the
        ! current model with the alternative model if the recent few (three) alternative models are
        ! more accurate in predicting the function value of XOPT + D.
        ! N.B.: Powell's code does this only if XOPT + D is feasible.
        call tryqalt(idz, bmat, fval - fval(kopt), xpt(:, kopt), xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt)

        ! Update RESCON. Zaikun 20221115: Currently, UPDATERES does not update RESCON if XIMPROVED
        ! is FALSE. Shouldn't we do it whenever DELTA is updated? Have we MISUNDERSTOOD RESCON?
        call updateres(ximproved, amat, b, delta, norm(d), xpt(:, kopt), rescon)
    end if  ! End of IF (IMPROVE_GEO). The procedure of improving geometry ends.

    ! The calculations with the current RHO are complete. Enhance the resolution of the algorithm
    ! by reducing RHO; update DELTA at the same time.
    if (reduce_rho) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        delta = HALF * rho
        rho = redrho(rho, rhoend)
        delta = max(delta, rho)
        ! Print a message about the reduction of RHO according to IPRINT.
        call rhomsg(solver, iprint, nf, fval(kopt), rho, xbase + xpt(:, kopt), cstrv, constr)
        ! DNORMSAV is corresponding to the latest function evaluations with the current RHO.
        ! Update it after reducing RHO.
        dnormsav = REALMAX
    end if  ! End of IF (REDUCE_RHO). The procedure of reducing RHO ends.

    ! Shift XBASE if XOPT may be too far from XBASE.
    ! Powell's original criterion for shifting XBASE: before a trust region step or a geometry step,
    ! shift XBASE if SUM(XOPT**2) >= 1.0E3*DELTA**2.
    if (sum(xpt(:, kopt)**2) >= 1.0E3_RP * delta**2) then
        ! Other possible criteria: SUM(XOPT**2) >= 1.0E4*DELTA**2, SUM(XOPT**2) >= 1.0E3*RHO**2.
        b = b - matprod(xpt(:, kopt), amat)
        call shiftbase(kopt, xbase, xpt, zmat, bmat, pq, hq, idz)
        ! SHIFTBASE shifts XBASE to XBASE + XOPT and XOPT to 0.
        pqalt = omega_mul(idz, zmat, fval - fval(kopt))
        galt = matprod(bmat(:, 1:npt), fval - fval(kopt)) + hess_mul(xpt(:, kopt), xpt, pqalt)
    end if
end do  ! End of DO TR = 1, MAXTR. The iterative procedure ends.

! Return from the calculation, after trying the Newton-Raphson step if it has not been tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xbase + (xpt(:, kopt) + d)
    call evaluate(calfun, x, f)
    nf = nf + 1_IK
    constr = matprod(x, A_orig) - b_orig
    cstrv = maximum([ZERO, constr])
    ! Print a message about the function evaluation according to IPRINT.
    call fmsg(solver, iprint, nf, f, x, cstrv, constr)
    ! Save X, F, CSTRV into the history.
    call savehist(nf, x, xhist, f, fhist, cstrv, chist)
    ! Save X, F, CSTRV into the filter.
    call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt)
end if

! Return the best calculated values of the variables.
kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol)
x = xfilt(:, kopt)
f = ffilt(kopt)
constr = matprod(x, A_orig) - b_orig
cstrv = maximum([ZERO, constr])

! Arrange CHIST, FHIST, and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist, chist)

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
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) &
        & .or. is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain negative values or NaN/+Inf', srname)
    nhist = minval([nf, maxfhist, maxchist])
    call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f, cstrv, ctol)),&
        & 'No point in the history is better than X', srname)
end if

end subroutine lincob


end module lincob_mod
