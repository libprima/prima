module uobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of UOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, March 06, 2023 PM12:15:41
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: uobyqb


contains


subroutine uobyqb(calfun, iprint, maxfun, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, &
    & x, nf, f, fhist, xhist, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine performs the major calculations of UOBYQA.
!
! The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to the corresponding arguments
! in subroutine UOBYQA.
!
! XBASE will contain a shift of origin that reduces the contributions from rounding errors to values
!   of the model and Lagrange functions.
! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
!   of the model and Lagrange functions.
! XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
!   the least calculated F so far). FOPT = F(XOPT + XBASE). However, we do not save XOPT and FOPT
!   explicitly, because XOPT = XPT(:, KOPT) and FOPT = FVAL(KOPT), which is explained below.
! [XPT, FVAL, KOPT] describes the interpolation set:
! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
!   values of F at the interpolation points; KOPT is the index of XOPT in XPT.
! PQ will contain the parameters of the quadratic model.
! PL will contain the parameters of the Lagrange functions.
! D is reserved for trial steps from XOPT. It is chosen by subroutine TRSTEP or GEOSTEP. Usually
!   XBASE + XOPT + D is the vector of variables for the next call of CALFUN.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TENTH, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: infos_mod, only : INFO_DFT, SMALL_TR_RADIUS, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : vec2smat, smat_mul_vec, norm
use, non_intrinsic :: output_mod, only : fmsg, rhomsg, retmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc
use, non_intrinsic :: ratio_mod, only : redrat
use, non_intrinsic :: redrho_mod, only : redrho
use, non_intrinsic :: shiftbase_mod, only : shiftbase

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep, setdrop_tr
use, non_intrinsic :: initialize_mod, only : initxf, initq, initl
use, non_intrinsic :: trustregion_mod, only : trstep, trrad
use, non_intrinsic :: update_mod, only : update

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
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
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'UOBYQA'
character(len=*), parameter :: srname = 'UOBYQB'
integer(IK) :: knew_geo
integer(IK) :: knew_tr
integer(IK) :: kopt
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxtr
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
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
real(RP) :: crvmin
real(RP) :: d(size(x))
real(RP) :: ddmove
real(RP) :: delbar
real(RP) :: delta
real(RP) :: distsq((size(x) + 1) * (size(x) + 2) / 2)
real(RP) :: dnorm
real(RP) :: dnormsav(2)  ! Powell's implementation: DNORMSAV(3)
real(RP) :: fval(size(distsq))
real(RP) :: g(size(x))
real(RP) :: h(size(x), size(x))
real(RP) :: moderr
real(RP) :: moderrsav(size(dnormsav))
real(RP) :: pl(size(distsq) - 1, size(distsq))
real(RP) :: pq(size(distsq) - 1)
real(RP) :: qred
real(RP) :: ratio
real(RP) :: rho
real(RP) :: xbase(size(x))
real(RP) :: xpt(size(x), size(distsq))
real(RP), parameter :: trtol = 1.0E-2_RP  ! Tolerance used in TRSTEP.

! Sizes.
n = int(size(x), kind(n))
npt = (n + 1_IK) * (n + 2_IK) / 2_IK
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

! Preconditions.
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
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

call initxf(calfun, iprint, maxfun, ftarget, rhobeg, x, kopt, nf, fhist, fval, xbase, xhist, xpt, subinfo)
x = xbase + xpt(:, kopt)
f = fval(kopt)

if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Arrange FHIST and XHIST so that they are in the chronological order.
    call rangehist(nf, xhist, fhist)
    ! Print a return message according to IPRINT.
    call retmsg(solver, info, iprint, nf, f, x)
    return
end if

call initq(fval, xpt, pq)
call initl(xpt, pl)

! Set some more initial values.
! We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers may raise a run-time error that
! RATIO is undefined. The value will not be used: when SHORTD = FALSE, its value will be overwritten;
! when SHORTD = TRUE, its value is used only in BAD_TRSTEP, which is TRUE regardless of RATIO.
! Similar for KNEW_TR.
! No need to initialize SHORTD unless MAXTR < 1, but some compilers may complain if we do not do it.
rho = rhobeg
delta = rho
shortd = .false.
ratio = -ONE
ddmove = -ONE
dnormsav = REALMAX
moderrsav = REALMAX
knew_tr = 0
knew_geo = 0

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
! UOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Generate trust region step D, and also calculate a lower bound on the Hessian of Q.
    g = pq(1:n) + smat_mul_vec(pq(n + 1:npt - 1), xpt(:, kopt))
    h = vec2smat(pq(n + 1:npt - 1))
    call trstep(delta, g, h, trtol, d, crvmin)

    ! Check whether D is too short to invoke a function evaluation.
    dnorm = min(delta, norm(d))
    shortd = (dnorm < HALF * rho)

    ! Set QRED to the reduction of the quadratic model when the move D is made from XOPT. QRED
    ! should be positive If it is nonpositive due to rounding errors, we will not take this step.
    qred = -quadinc(pq, d, xpt(:, kopt))  ! QRED = Q(XOPT) - Q(XOPT + D)

    if (shortd .or. .not. qred > 0) then
        ! Powell's code does not reduce DELTA as follows. This comes from NEWUOA and works well.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close to or below.
        end if
    else
        ! Calculate the next value of the objective function.
        x = xbase + (xpt(:, kopt) + d)
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

        ! Set KNEW to the index of the next interpolation point to be deleted.
        knew_tr = setdrop_tr(kopt, ximproved, d, pl, rho, xpt)

        ! DDMOVE is norm square of DMOVE in the UOBYQA paper. See Steps 6--7 in Sec. 5 of the paper.
        ddmove = ZERO
        if (knew_tr > 0) then
            ddmove = sum((xpt(:, knew_tr) - xpt(:, kopt))**2)  ! KOPT is unupdated.
            ! Update PL, PQ, XPT, FVAL, and KOPT so that XPT(:, KNEW_TR) becomes XOPT + D.
            call update(knew_tr, d, f, moderr, kopt, fval, pl, pq, xpt)
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
    accurate_mod = all(abs(moderrsav) <= 0.125_RP * crvmin * rho**2) .and. all(dnormsav <= rho)
    ! CLOSE_ITPSET: Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
    close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Behaves the same as Powell's version.
    ! Below are some alternative definitions of CLOSE_ITPSET.
    ! !close_itpset = all(distsq <= 4.0_RP * rho**2)  ! Powell's code.
    ! !close_itpset = all(distsq <= max((2.0_RP * delta)**2, (10.0_RP * rho)**2))  ! Powell's BOBYQA.
    ! ADEQUATE_GEO: Is the geometry of the interpolation set "adequate"?
    adequate_geo = (shortd .and. accurate_mod) .or. close_itpset
    ! SMALL_TRRAD: Is the trust-region radius small? This indicator seems not impactive in practice.
    small_trrad = (max(delta, dnorm) <= rho)  ! Behaves the same as Powell's version.
    !small_trrad = (dnorm <= rho)  ! Powell's code.

    ! Comments on ACCURATE_MOD:
    ! 1. ACCURATE_MOD is needed only when SHORTD is TRUE.
    ! 2. In Powell's UOBYQA code, ACCURATE_MOD is defined according to (28), (37), and (38) in the
    ! UOBYQA paper (see also (32) of Powell 2001: "On the Lagrange functions of quadratic models
    ! that are defined by interpolation"). As elaborated in Sec. 3 of the paper (also Sec. 4 of
    ! Powell 2001), the idea is to test whether the current model is sufficiently accurate by
    ! checking whether the interpolation error bound in (28) is (sufficiently) small. If the bound
    ! is small, then set ACCURATE_MOD to TRUE. Otherwise, it identifies a "bad" interpolation point
    ! that makes a significant contribution to the bound, with a preference to the interpolation
    ! points that are a far away from the current trust-region center. Such a point will be replaced
    ! with a new point obtained by the geometry step. If all the interpolation points are close
    ! enough to the trust-region center, then they are all considered to be good.
    ! 3. Our implementation defines ACCURATE_MOD by a method from NEWUOA and BOBYQA, which is also
    ! reflected in LINCOA. It sets ACCURATE_MOD to TRUE if recent model errors and step lengths are
    ! all small. In addition, it identifies a "bad" interpolation point by simply taking the
    ! farthest point from the current trust region center, unless they are all close enough to the
    ! center. This implementation is much simpler and less costly in terms of flops yet it performs
    ! almost the same as Powell's original implementation.

    ! Powell's original definition of IMPROVE_GEO and REDUCE_RHO:
    ! !bad_trstep = (shortd .or. knew_tr == 0 .or. (ratio <= 0 .and. dnorm <= 2.0_RP*rho .and. ddmove <= 4.0_RP * rho**2))
    ! !improve_geo = bad_trstep .and. .not. (shortd .and. accurate_mod) .and. .not. close_itpset
    ! !reduce_rho = bad_trstep .and. dnorm <= rho .and. .not. improve_geo

    ! IMPROVE_GEO and REDUCE_RHO are defined as follows.
    ! N.B.: If SHORTD is TRUE at the very first iteration, then REDUCE_RHO will be set to TRUE.
    ! Powell's code does not have (.NOT. QRED>0) in BAD_TRSTEP; it terminates if QRED > 0 fails.

    ! BAD_TRSTEP (for IMPROVE_GEO): Is the last trust-region step bad? For UOBYQA, it is CRITICAL to
    ! include DMOVE <= 4.0_RP*RHO**2 in the definition of BAD_TRSTEP for IMPROVE_GEO.
    bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= eta1 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    !bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= eta1 .or. knew_tr == 0)  ! Works poorly!
    improve_geo = bad_trstep .and. .not. adequate_geo
    ! BAD_TRSTEP (for REDUCE_RHO): Is the last trust-region step bad?
    bad_trstep = (shortd .or. (.not. qred > 0) .or. ratio <= 0 .or. knew_tr == 0)  ! Performs better than the below from Powell.
    !bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    reduce_rho = bad_trstep .and. adequate_geo .and. small_trrad

    ! Equivalently, REDUCE_RHO can be set as follows. It shows that REDUCE_RHO is TRUE in two cases.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    ! !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

    ! With REDUCE_RHO properly defined, we can also set IMPROVE_GEO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= eta1 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    ! !improve_geo = bad_trstep .and. (.not. reduce_rho) .and. (.not. close_itpset)

    ! With IMPROVE_GEO properly defined, we can also set REDUCE_RHO as follows.
    ! !bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    ! !reduce_rho = bad_trstep .and. (.not. improve_geo) .and. small_trrad

    ! UOBYQA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
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

        ! DELBAR is the trust-region radius for the geometry improvement subproblem.
        ! Powell's UOBYQA code sets DELBAR = RHO, but NEWUOA/BOBYQA/LINCOA all take DELTA and/or
        ! DISTSQ into consideration. The following DELBAR is copied from NEWUOA, and it seems to
        ! improve the performance slightly according to a test on 20220720.
        delbar = max(min(TENTH * sqrt(maxval(distsq)), HALF * delta), rho)

        d = geostep(knew_geo, kopt, delbar, pl, xpt)

        ! Calculate the next value of the objective function.
        x = xbase + (xpt(:, kopt) + d)
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
        dnorm = min(delbar, norm(d))   ! In theory, DNORM = DELBAR in this case.
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
        ! MODERR is the error of the current model in predicting the change in F due to D.
        ! MODERRSAV is the prediction errors of the latest 3 models with the current RHO.
        moderr = f - fval(kopt) - quadinc(pq, d, xpt(:, kopt))  ! QUADINC = Q(XOPT + D) - Q(XOPT)
        moderrsav = [moderrsav(2:size(moderrsav)), moderr]

        ! Update PL, PQ, XPT, FVAL, and KOPT so that XPT(:, KNEW_GEO) becomes XOPT + D.
        call update(knew_geo, d, f, moderr, kopt, fval, pl, pq, xpt)
    end if

    if (reduce_rho) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if

        !! Shifting XBASE to the best point so far, and make the corresponding changes to the
        !! gradients of the Lagrange functions and the quadratic model.
        !call shiftbase(kopt, pl, pq, xbase, xpt)  ! 0918

        ! Pick the next values of RHO and DELTA.
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

    ! Shifting XBASE to the best point so far, and make the corresponding changes to the gradients
    ! of the Lagrange functions and the quadratic model. Powell's implementation does this each time
    ! after RHO is reduced. Our implementation aligns with NEWUOA/BOBYQA/LINCOA.
    if (sum(xpt(:, kopt)**2) >= 1.0E3_RP * delta**2) then
        call shiftbase(kopt, pl, pq, xbase, xpt)
    end if
end do

! Return, possibly after another Newton-Raphson step, if it is too short to have been tried before.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun) then
    x = xbase + (xpt(:, kopt) + d)
    call evaluate(calfun, x, f)
    nf = nf + 1_IK
    ! Print a message about the function evaluation according to IPRINT.
    call fmsg(solver, iprint, nf, f, x)
    ! Save X, F into the history.
    call savehist(nf, x, xhist, f, fhist)
end if

! Choose the [X, F] to return: either the current [X, F] or [XBASE + XOPT, FOPT].
if (fval(kopt) < f .or. is_nan(f)) then
    x = xbase + xpt(:, kopt)
    f = fval(kopt)
end if

! Arrange FHIST and XHIST so that they are in the chronological order.
call rangehist(nf, xhist, fhist)

! Print a return message according to IPRINT.
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

end subroutine uobyqb


end module uobyqb_mod
