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
! Last Modified: Sunday, January 08, 2023 AM12:44:55
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
! modified for variables relative to XBASE.
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
! of the model and Lagrange functions.
! XPT contains the interpolation point coordinates relative to XBASE.
! FVAL holds the values of F at the interpolation points.
! XSAV holds the best feasible vector of variables so far, without any shift of origin.
! XOPT is set to XSAV-XBASE, which is the displacement from XBASE of the feasible vector of variables
! that provides the least calculated F so far, this vector being the current trust region centre.
! GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
! HQ holds the explicit second derivatives of the quadratic model.
! PQ contains the parameters of the implicit second derivatives of the quadratic model.
! BMAT holds the last N columns of the big inverse matrix H.
! ZMAT holds the factorization of the leading NPT by NPT submatrix of H, this factorization being
! ZMAT * Diag(DZ) * ZMAT^T, where the elements of DZ are plus or minus ONE, as specified by IDZ.
! D is employed for trial steps from XOPT.
! IACT is an integer array for the indices of the active constraints.
! RESCON holds information about the constraint residuals at the current trust region center XOPT.
! 1. If if B(J) - AMAT(:, J)^T*XOPT <= DELTA, then RESCON(J) = B(J) - AMAT(:, J)^T*XOPT. Note that
! RESCON >= 0 in this case, because the algorithm keeps XOPT to be feasible.
! 2. Otherwise, RESCON(J) is a negative value that B(J) - AMAT(:, J)^T*XOPT >= |RESCON(J)| >= DELTA.
! RESCON can be updated without calculating the constraints that are far from being active, so that
! we only need to evaluate the constraints that are nearly active.
! QFAC is the orthogonal part of the QR factorization of the matrix of active constraint gradients,
! these gradients being ordered in accordance with IACT. When NACT is less than N, columns are added
! to QFAC to complete an N by N orthogonal matrix, which is important for keeping calculated steps
! sufficiently close to the boundaries of the active constraints.
! RFAC is the upper triangular part of this QR factorization, beginning with the first diagonal
! element, followed by the two elements in the upper triangular part of the second column and so on.
!--------------------------------------------------------------------------------------------------!

! Generic models
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TENTH, HUGENUM, MIN_MAXFILT, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
!use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: infos_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS
use, non_intrinsic :: linalg_mod, only : matprod, maximum, eye, trueloc
use, non_intrinsic :: output_mod, only : fmsg, rhomsg, retmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, omega_mul, hess_mul
use, non_intrinsic :: powalg_mod, only : updateh
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: selectx_mod, only : savefilt, selectx
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
integer(IK) :: iact(size(bvec))
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n, maxtr
real(RP) :: b(size(bvec))
real(RP) :: bmat(size(x), npt + size(x))
real(RP) :: fval(npt), cval(npt)
real(RP) :: gopt(size(x))
real(RP) :: hq(size(x), size(x))
real(RP) :: pq(npt)
real(RP) :: qfac(size(x), size(x))
real(RP) :: rescon(size(bvec))
real(RP) :: rfac(size(x), size(x))
real(RP) :: xfilt(size(x), maxfilt), ffilt(maxfilt), cfilt(maxfilt)
real(RP) :: d(size(x))
real(RP) :: xbase(size(x))
real(RP) :: xopt(size(x)), xosav(size(x))
real(RP) :: xpt(size(x), npt), xdrop(size(x))
real(RP) :: zmat(npt, npt - size(x) - 1)
real(RP) :: delbar, delta, &
&        distsq(npt), fopt, ratio,     &
&        rho, dnorm, &
&        qred, constr(size(bvec))
logical :: qalt_better(3)
logical :: accurate_mod, adequate_geo
logical :: bad_trstep
logical :: close_itpset
logical :: small_trrad
logical :: evaluated(npt)
logical :: feasible, shortd, improve_geo, reduce_rho, ximproved
integer(IK) :: ij(2, max(0_IK, int(npt - 2 * size(x) - 1, IK))), k
integer(IK) :: nfilt, idz, &
&           knew_tr, knew_geo, kopt, nact,      &
&           ngetact, subinfo, tr
real(RP) :: fshift(npt)
real(RP) :: pqalt(npt), galt(size(x))
real(RP) :: dnormsav(5)
real(RP) :: moderr, moderr_alt
real(RP) :: gamma3

! Sizes.
m = int(size(bvec), kind(m))
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
    call assert(size(b_orig) == m, 'SIZE(B_ORIG) == M', srname)
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
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = xbase + xopt
f = fopt
! For the output, we use A_ORIG and B_ORIG to evaluate the constraints.
cstrv = maximum([ZERO, matprod(x, A_orig) - b_orig])

nfilt = 0_IK
do k = 1, npt
    if (evaluated(k)) then
        x = xbase + xpt(:, k)
        call savefilt(cval(k), ctol, cweight, fval(k), x, nfilt, cfilt, ffilt, xfilt)
    end if
end do

if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Return the best calculated values of the variables.
    kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol)
    x = xfilt(:, kopt)
    f = ffilt(kopt)
    cstrv = cfilt(kopt)
    ! Arrange CHIST, FHIST, and XHIST so that they are in the chronological order.
    call rangehist(nf, xhist, fhist, chist)
    call retmsg(solver, info, iprint, nf, f, x, cstrv)
    !close (16)
    return
end if

! Initialize BMAT, ZMAT, and IDZ.
call inith(ij, rhobeg, xpt, idz, bmat, zmat)

! Initialize GQ, HQ, and PQ.
hq = ZERO
pq = omega_mul(idz, zmat, fval)
gopt = matprod(bmat(:, 1:npt), fval) + hess_mul(xopt, xpt, pq)
pqalt = pq
galt = gopt

! Initialize RESCON.
rescon = max(b - matprod(xopt, amat), ZERO)
rescon(trueloc(rescon >= rhobeg)) = -rescon(trueloc(rescon >= rhobeg))
!!MATLAB: rescon(rescon >= rhobeg) = -rescon(rescon >= rhobeg)

gamma3 = sqrt(gamma2)  ! Used in TRRAD; 0 < GAMMA1 < 1 < GAMMA3 <= GAMMA2.
qfac = eye(n)
rfac = ZERO
rho = rhobeg
delta = rho
qred = ZERO
ratio = -ONE
knew_tr = 0
knew_geo = 0
feasible = .false.
shortd = .false.
improve_geo = .false.
nact = 0
qalt_better = .false.
dnormsav = HUGENUM

! MAXTR is the maximal number of trust-region iterations. Each trust-region iteration takes 1 or 2
! function evaluations unless the trust-region step is short but the geometry step is not invoked.
! Thus the following MAXTR is unlikely to be reached.
maxtr = max(maxfun, 2_IK * maxfun)  ! MAX: precaution against overflow, which will make 2*MAXFUN < 0.
info = MAXTR_REACHED

! Begin the iterative procedure.
! After solving a trust-region subproblem, we use three boolean variables to control the workflow.
! SHORTD: Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO: Should we improve the geometry?
! REDUCE_RHO: Should we reduce rho?
! LINCOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Shift XBASE if XOPT may be too far from XBASE.
    ! Zaikun 20220528: The criteria is different from those in NEWUOA or BOBYQA, particularly here
    ! |XOPT| is compared with DELTA instead of DNORM. What about unifying the criteria, preferably
    ! to the one here? What about comparing with RHO? What about calling SHIFTBASE only before
    ! TRSTEP but not GEOSTEP (consider GEOSTEP as a postprocessor).
    if (sum(xopt**2) >= 1.0E4_RP * delta**2) then
        b = b - matprod(xopt, amat)
        call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq, idz)
        fshift = fval - fopt
        pqalt = omega_mul(idz, zmat, fshift)
        galt = matprod(bmat(:, 1:npt), fshift) + hess_mul(xopt, xpt, pqalt)
    end if

    ! Generate the next trust region step D by calling TRSTEP. Note that D is feasible.
    call trstep(amat, delta, gopt, hq, pq, rescon, xpt, iact, nact, qfac, rfac, ngetact, d)
    dnorm = min(delta, sqrt(sum(d**2)))

    ! A trust region step is applied whenever its length is at least 0.5*DELTA. It is also
    ! applied if its length is at least 0.1999*DELTA and if a line search of TRSTEP has caused a
    ! change to the active set, indicated by NGETACT >= 2 (note that NGETACT is at least 1).
    ! Otherwise, the trust region step is considered too short to try.
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

    ! In some cases, we reset DNORMSAV to HUGENUM. This indicates a preference of improving the
    ! geometry of the interpolation set to reducing RHO in the subsequent three or more
    ! iterations. This is important for the performance of LINCOA.
    if (delta > rho .or. .not. shortd) then  ! Another possibility: IF (DELTA > RHO) THEN
        dnormsav = HUGENUM
    end if

    ! Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
    ! should be positive If it is nonpositive due to rounding errors, we will not take this step.
    qred = -quadinc(d, xpt, gopt, pq, hq)  ! QRED = Q(XOPT) - Q(XOPT + D)

    if (shortd .or. .not. qred > 0) then
        ! In this case, do nothing but reducing DELTA. Afterward, DELTA < DNORM may occur.
        ! N.B.: 1. This value of DELTA will be discarded if REDUCE_RHO turns out TRUE later.
        ! 2. Powell's code does not shrink DELTA when QRED > 0 is FALSE (i.e., when VQUAD >= 0 in
        ! Powell's code, where VQUAD = -QRED). Consequently, the algorithm may  be stuck in an
        ! infinite cycling, because both REDUCE_RHO and IMPROVE_GEO may end up with FALSE in this
        ! case, which did happen in tests.
        ! 3. The factor HALF works better than TENTH (used in NEWUOA/BOBYQA), 0.2, and 0.7.
        ! 4. The factor 0.99*GAMMA3 aligns with the update of DELTA after a trust-region step.
        delta = HALF * delta
        if (delta <= min(0.99_RP * gamma3, 1.5_RP) * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if
    else
        ! Calculate the next value of the objective function.
        x = xbase + (xopt + d)
        call evaluate(calfun, x, f)
        nf = nf + 1_IK

        ! Evaluate the constraints using A_ORIG and B_ORIG.
        constr = matprod(x, A_orig) - b_orig
        cstrv = maximum([ZERO, constr])

        ! Print a message about the function evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x, cstrv)
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
        moderr = f - fopt + qred
        moderr_alt = f - fopt - quadinc(d, xpt, galt, pqalt)
        qalt_better = [qalt_better(2:size(qalt_better)), abs(moderr_alt) < TENTH * abs(moderr)]

        ! Calculate the reduction ratio by REDRAT, which handles Inf/NaN carefully.
        ratio = redrat(fopt - f, qred, eta1)

        ! Update DELTA. After this, DELTA < DNORM may hold.
        ! The new DELTA lies in [GAMMA1*DNORM, GAMMA3*DELTA].
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, gamma3, ratio)
        ! Set DELTA to RHO when it is close to or below. The multiplicative factor in the scheme
        ! below should be less than GAMMA3. Imagine a very successful step with DENORM = the
        ! un-updated DELTA = RHO. Then TRRAD will update DELTA to GAMMA3*RHO. If this factor were
        ! not smaller than GAMMA3, then DELTA will be reset to RHO, which is not reasonable as D is
        ! very successful. See paragraph two of Sec. 5.2.5 in T. M. Ragonneau's thesis: "Model-Based
        ! Derivative-Free Optimization Methods and Software".
        if (delta <= min(0.99_RP * gamma3, 1.5_RP) * rho) then
            delta = rho
        end if

        ! Is the newly generated X better than current best point?
        ximproved = (f < fopt)

        ! Set KNEW_TR to the index of the interpolation point to be replaced by XNEW = XOPT + D.
        ! KNEW_TR will ensure that the geometry of XPT is "good enough" after the replacement.
        knew_tr = setdrop_tr(idz, kopt, ximproved, bmat, d, xpt, zmat)
        if (knew_tr > 0) then
            ! Update [BMAT, ZMAT, IDZ] (represents H in the NEWUOA paper), [XPT, FVAL, KOPT, XOPT,
            ! FOPT] and [GQ, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_TR) becomes
            ! XNEW = XOPT + D.
            xdrop = xpt(:, knew_tr)
            xosav = xpt(:, kopt)
            call updateh(knew_tr, kopt, idz, d, xpt, bmat, zmat)
            call updatexf(knew_tr, ximproved, f, xopt + d, kopt, fval, xpt, fopt, xopt)
            call updateq(idz, knew_tr, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)

            ! Establish the alternative model, namely the least Frobenius norm interpolant. Replace
            ! the current model with the alternative model if the recent few (three) alternative
            ! models are more accurate in predicting the function value of XOPT + D.
            call tryqalt(idz, bmat, fval - fopt, xopt, xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt)

            ! Update RESCON if XOPT is changed. Zaikun 20221115: Shouldn't we do it after DELTA is updated?
            call updateres(ximproved, amat, b, delta, sqrt(sum(d**2)), xopt, rescon)
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
    accurate_mod = all(dnormsav <= HALF * rho) .or. all(dnormsav(3:size(dnormsav)) <= TENTH * rho)
    ! CLOSE_ITPSET: Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xopt).^2)  % xopt should be a column! Implicit expansion
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
        ! Shift XBASE if XOPT may be too far from XBASE.
        ! Zaikun 20220528: The criteria is different from those in NEWUOA or BOBYQA, particularly here
        ! |XOPT| is compared with DELTA instead of DNORM. What about unifying the criteria, preferably
        ! to the one here? What about comparing with RHO? What about calling SHIFTBASE only before
        ! TRSTEP but not GEOSTEP (consider GEOSTEP as a postprocessor).
        if (sum(xopt**2) >= 1.0E4_RP * delta**2) then
            b = b - matprod(xopt, amat)
            call shiftbase(xbase, xopt, xpt, zmat, bmat, pq, hq, idz)
            fshift = fval - fopt
            pqalt = omega_mul(idz, zmat, fshift)
            galt = matprod(bmat(:, 1:npt), fshift) + hess_mul(xopt, xpt, pqalt)
        end if

        ! XPT(:, KNEW_GEO) will become  XOPT + D below. KNEW_GEO /= KOPT unless there is a bug.
        knew_geo = int(maxloc(distsq, dim=1), kind(knew_geo))

        ! Set DELBAR, which will be used as the trust-region radius for the geometry-improving
        ! scheme GEOSTEP. Note that DELTA has been updated before arriving here.
        delbar = max(TENTH * delta, rho)  ! This differs from NEWUOA/BOBYQA. Possible improvement?
        ! Find D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) becomes XOPT + D.
        call geostep(iact, idz, knew_geo, kopt, nact, amat, bmat, delbar, qfac, rescon, xpt, zmat, feasible, d)

        ! Calculate the next value of the objective function.
        x = xbase + (xopt + d)
        call evaluate(calfun, x, f)
        nf = nf + 1_IK

        ! Evaluate the constraints using A_ORIG and B_ORIG.
        constr = matprod(x, A_orig) - b_orig
        cstrv = maximum([ZERO, constr])

        ! Print a message about the function evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x, cstrv)
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
        moderr = f - fopt - quadinc(d, xpt, gopt, pq, hq)
        moderr_alt = f - fopt - quadinc(d, xpt, galt, pqalt)
        qalt_better = [qalt_better(2:size(qalt_better)), abs(moderr_alt) < TENTH * abs(moderr)]

        ! Is the newly generated X better than current best point?
        ximproved = (f < fopt .and. feasible)

        ! Update [BMAT, ZMAT, IDZ] (represents H in the NEWUOA paper), [XPT, FVAL, KOPT, XOPT, FOPT]
        ! and [GQ, HQ, PQ] (the quadratic model), so that XPT(:, KNEW_GEO) becomes XNEW = XOPT + D.
        xdrop = xpt(:, knew_geo)
        xosav = xpt(:, kopt)
        call updateh(knew_geo, kopt, idz, d, xpt, bmat, zmat)
        call updatexf(knew_geo, ximproved, f, xopt + d, kopt, fval, xpt, fopt, xopt)
        call updateq(idz, knew_geo, ximproved, bmat, d, moderr, xdrop, xosav, xpt, zmat, gopt, hq, pq)

        ! Establish the alternative model, namely the least Frobenius norm interpolant. Replace the
        ! current model with the alternative model if the recent few (three) alternative models are
        ! more accurate in predicting the function value of XOPT + D.
        ! N.B.: Powell's code does this only if XOPT + D is feasible.
        call tryqalt(idz, bmat, fval - fopt, xopt, xpt, zmat, qalt_better, gopt, pq, hq, galt, pqalt)

        ! Update RESCON. Zaikun 20221115: Currently, UPDATERES does not update RESCON unless
        ! XIMPROVED is TRUE. Shouldn't we do it whenever DELTA is updated? Have we MISUNDERSTOOD
        ! RESCON?
        call updateres(ximproved, amat, b, delta, sqrt(sum(d**2)), xopt, rescon)
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
        call rhomsg(solver, iprint, nf, fopt, rho, xbase + xopt, cstrv)
        ! DNORMSAV is corresponding to the latest function evaluations with the current RHO.
        ! Update it after reducing RHO.
        dnormsav = HUGENUM
    end if  ! End of IF (REDUCE_RHO). The procedure of reducing RHO ends.
end do  ! End of DO TR = 1, MAXTR. The iterative procedure ends.

! Return from the calculation, after trying the Newton-Raphson step if it has not been tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xbase + (xopt + d)
    call evaluate(calfun, x, f)
    nf = nf + 1_IK
    constr = matprod(x, A_orig) - b_orig
    cstrv = maximum([ZERO, constr])
    ! Print a message about the function evaluation according to IPRINT.
    call fmsg(solver, iprint, nf, f, x, cstrv)
    ! Save X, F, CSTRV into the history.
    call savehist(nf, x, xhist, f, fhist, cstrv, chist)
    ! Save X, F, CSTRV into the filter.
    call savefilt(cstrv, ctol, cweight, f, x, nfilt, cfilt, ffilt, xfilt)
end if

! Return the best calculated values of the variables.
kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cweight, ctol)
x = xfilt(:, kopt)
f = ffilt(kopt)
cstrv = cfilt(kopt)

! Arrange CHIST, FHIST, and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist, chist)

! Print a return message according to IPRINT.
call retmsg(solver, info, iprint, nf, f, x, cstrv)

!====================!
!  Calculation ends  !
!====================!

! Postconditions

!close (16)

end subroutine lincob


end module lincob_mod
