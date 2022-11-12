!TODO:
!3. THE checks in rangehist cannot pass(xhist does not contain NaN)

module uobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of UOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Saturday, November 12, 2022 PM05:18:37
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
! XBASE will contain a shift of origin that reduces the contributions from rounding errors to values
! of the model and Lagrange functions.
! XOPT will be set to the displacement from XBASE of the vector of variables that provides the least
! calculated F so far.
! XNEW will be set to the displacement from XBASE of the vector of variables for the current
! calculation of F.
! XPT will contain the interpolation point coordinates relative to XBASE.
! PQ will contain the parameters of the quadratic model.
! PL will contain the parameters of the Lagrange functions.
! H will provide the second derivatives that TRSTEP and LAGMAX require.
! G will provide the first derivatives that TRSTEP and LAGMAX require.
! D is reserved for trial steps from XOPT, except that it will contain diagonal second derivatives
! during the initialization procedure.
! VLAG will contain the values of the Lagrange functions at a new point X.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: initialize_mod, only : initxf, initq, initl
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : INFO_DFT, NAN_INF_X, NAN_INF_F, FTARGET_ACHIEVED, &
    & MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS!, NAN_INF_MODEL!, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : inprod, outprod!, norm
use, non_intrinsic :: symmat_mod, only : vec2smat, smat_mul_vec
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calvlag

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep
use, non_intrinsic :: trustregion_mod, only : trstep


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
character(len=*), parameter :: srname = 'UOBYQB'
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: maxhist
integer(IK) :: maxfhist
integer(IK) :: maxxhist
real(RP) :: d(size(x))
real(RP) :: g(size(x))
real(RP) :: h(size(x), size(x))
real(RP) :: pl((size(x) + 1) * (size(x) + 2) / 2 - 1, (size(x) + 1) * (size(x) + 2) / 2)
real(RP) :: pq(size(pl, 1))
real(RP) :: vlag(size(pl, 2))
real(RP) :: xbase(size(x))
real(RP) :: xnew(size(x))
real(RP) :: xopt(size(x))
real(RP) :: xpt(size(x), size(pl, 2))
real(RP) :: ddmove, delta, diff, distsq(size(pl, 2)), delbar, &
& weight(size(pl, 2)), score(size(pl, 2)),    &
&        dnorm, errtol, crvmin, fopt,&
&        fsave, ratio, rho, &
&        trtol, &
&        qred, plknew(size(pl, 1)), fval(size(pl, 2))
integer(IK) :: k, knew_tr, knew_geo, kopt, subinfo
logical :: tr_success, shortd, improve_geo, reduce_rho, accurate_mod, adequate_geo, close_itpset, small_trrad, bad_trstep
real(RP) :: dnormsav(3), moderrsav(size(dnormsav))

! Sizes.
n = int(size(x), kind(n))
npt = (n + 1_IK) * (n + 2_IK) / 2_IK
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

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


!--------------------------------------------------------------------------------------------------!
! Temporary fix for the G95 warning that these variables are used uninitialized.
!knew_tr = 1; kopt = 1
!--------------------------------------------------------------------------------------------------!

!====================!
! Calculation starts !
!====================!

call initxf(calfun, iprint, maxfun, ftarget, rhobeg, x, kopt, nf, fhist, fval, xbase, xhist, xpt, subinfo)
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = xbase + xopt
f = fopt

if (subinfo /= INFO_DFT) then
    info = subinfo
    call rangehist(nf, xhist, fhist)
    return
end if

call initq(fval, xpt, pq)
call initl(xpt, pl)

! Set parameters to begin the iterations for the current RHO.
rho = rhobeg
delta = rho
moderrsav = HUGENUM
dnormsav = HUGENUM
ratio = -ONE
knew_tr = 0_IK
shortd = .false.
reduce_rho = .false.
trtol = 0.01_RP

! Form the gradient of the quadratic model at the trust region centre.
do while (.true.)
    xopt = xpt(:, kopt)
    g = pq(1:n) + smat_mul_vec(pq(n + 1:npt - 1), xopt)
    h = vec2smat(pq(n + 1:npt - 1))

    ! Generate the next trust region step and test its length. Set KNEW to -1 if the purpose of
    ! the next F will be to improve conditioning, and also calculate a lower bound on the
    ! Hessian term of the model Q.
    call trstep(delta, g, h, trtol, d, crvmin)
    dnorm = min(delta, sqrt(sum(d**2)))
    errtol = ZERO
    shortd = (dnorm < HALF * rho)
    ! Use the quadratic model to predict the change in F due to the step D
    qred = -quadinc(pq, d, xopt)
    if (shortd .or. .not. qred > 0) then
        ! Powell's code does not reduce DELTA as follows. This comes from NEWUOA and works well.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if
        errtol = HALF * crvmin * rho * rho
        if (nf <= npt + 9 .or. is_nan(errtol)) errtol = ZERO
    else
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]
        ! Calculate the next value of the objective function.
        xnew = xopt + d
        x = xbase + xnew

        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if

        if (is_nan(sum(abs(x)))) then
            f = sum(x) ! Set F to NaN
            if (nf == 1) then
                fopt = f
                xopt = ZERO
            end if
            info = NAN_INF_X
            exit
        end if
        call evaluate(calfun, x, f)
        nf = nf + 1
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
            xopt = xnew
            fopt = f
            exit
        end if

        moderrsav = [moderrsav(2:size(moderrsav)), f - fopt + qred]
        vlag = calvlag(pl, d, xopt, kopt)

        diff = f - fopt + qred

        ! Update FOPT and XOPT if the new F is the least value of the objective function so far.
        ! Then branch if D is not a trust region step.
        fsave = fopt
        if (f < fopt) then
            fopt = f
            xopt = xnew
        end if

        ! Pick the next value of DELTA after a trust region step.
        if (.not. (qred > 0)) then
            info = TRSUBP_FAILED
            exit
        end if
        ratio = (fsave - f) / qred
        if (ratio <= TENTH) then
            delta = HALF * dnorm
        else if (ratio <= 0.7_RP) then
            delta = max(HALF * delta, dnorm)
        else
            delta = max(delta, 1.25_RP * dnorm, dnorm + rho)
        end if
        if (delta <= 1.5_RP * rho) delta = rho

        tr_success = (f < fsave)

        ! Set KNEW to the index of the next interpolation point to be deleted.

        ! Calculate the distance squares between the interpolation points and the "optimal point".
        ! When identifying the optimal point, as suggested in (56) of the UOBYQA paper and (7.5) of
        ! the NEWUOA paper, it is reasonable to take into account the new trust-region trial point
        ! XPT(:, KOPT) + D, which will become the optimal point in the next interpolation if
        ! TR_SUCCESS is TRUE.
        distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)  ! XOPT has been updated.
        !distsq = sum((xpt - spread(xsave, dim=2, ncopies=npt))**2, dim=1)  ! XSAVE is the unupdated XOPT
        weight = max(ONE, distsq / rho**2)**4

        !------------------------------------------------------------------------------------------!
        ! Other possible definitions of WEIGHT.
        !weight = max(ONE, distsq / rho**2)**3.5_RP ! ! No better than power 4.
        !weight = max(ONE, distsq / delta**2)**3.5_RP  ! Not better than DISTSQ/RHO**2.
        !weight = max(ONE, distsq / rho**2)**1.5_RP  ! Powell's origin code: power 1.5.
        !weight = max(ONE, distsq / rho**2)**2  ! Better than power 1.5.
        !weight = max(ONE, distsq / delta**2)**2  ! Not better than DISTSQ/RHO**2.
        !weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**2  ! The same as DISTSQ/RHO**2.
        !weight = distsq**2  ! Not better than MAX(ONE, DISTSQ/RHO**2)**2
        !weight = max(ONE, distsq / rho**2)**3  ! Better than power 2.
        !weight = max(ONE, distsq / delta**2)**3  ! Similar to DISTSQ/RHO**2; not better than it.
        !weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3  ! The same as DISTSQ/RHO**2.
        !weight = distsq**3  ! Not better than MAX(ONE, DISTSQ/RHO**2)**3
        !weight = max(ONE, distsq / delta**2)**4  ! Not better than DISTSQ/RHO**2.
        !weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**4  ! The same as DISTSQ/RHO**2.
        !weight = distsq**4  ! Not better than MAX(ONE, DISTSQ/RHO**2)**4
        !------------------------------------------------------------------------------------------!

        score = weight * abs(vlag)

        ! If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
        if (.not. tr_success) then
            score(kopt) = -ONE
        end if

        knew_tr = 0_IK
        ! Changing the IF below to `IF (ANY(SCORE>0)) THEN` does not render a better performance.
        if (any(score > 1) .or. (tr_success .and. any(score > 0))) then
            ! SCORE(K) is NaN implies VLAG(K) is NaN, but we want ABS(VLAG) to be big. So we
            ! exclude such K.
            knew_tr = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
            !!MATLAB: [~, knew_tr] = max(score, [], 'omitnan');
        elseif (tr_success) then
            ! Powell's code does not include the following instructions. With Powell's code,
            ! if DENABS consists of only NaN, then KNEW can be 0 even when TR_SUCCESS is TRUE.
            knew_tr = int(maxloc(distsq, dim=1), IK)
        end if

        ! DDMOVE is norm square of DMOVE in the UOBYQA paper. See Steps 6--7 in Sec. 5 of the paper.
        ddmove = ZERO
        if (knew_tr > 0) then
            ddmove = sum((xpt(:, knew_tr) - xpt(:, kopt))**2)  ! KOPT is unupdated.
            ! Replace the interpolation point that has index KNEW by the point XNEW, and also update
            ! the Lagrange functions and the quadratic model.
            xpt(:, knew_tr) = xnew
            ! It can happen that VLAG(KNEW) = 0 due to rounding.
            pl(:, knew_tr) = pl(:, knew_tr) / vlag(knew_tr)
            plknew = pl(:, knew_tr)
            pq = pq + diff * plknew
            pl = pl - outprod(plknew, vlag)
            pl(:, knew_tr) = plknew

            ! Update KOPT if F is the least calculated value of the objective function.
            if (f < fsave) then
                kopt = knew_tr
            end if
        end if
    end if


    !----------------------------------------------------------------------------------------------!
    ! Before the next trust-region iteration, we may improve the geometry of XPT or reduce RHO
    ! according to IMPROVE_GEO and REDUCE_RHO, which in turn depend on the following indicators.
    ! ACCURATE_MOD --- Are the recent models sufficiently accurate? Used only if SHORTD is TRUE.
    accurate_mod = all(abs(moderrsav) <= 0.125_RP * crvmin * rho**2) .and. all(dnormsav <= rho)
    ! SMALL_TRRAD --- Is the trust-region radius small?  This indicator seems not impactive.
    small_trrad = (max(delta, dnorm) <= rho)  ! Behaves the same as Powell's version.
    !small_trrad = (dnorm <= rho)  ! Powell's code.
    ! CLOSE_ITPSET --- Are the interpolation points close to XOPT?
    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xopt).^2)  % xopt should be a column!! Implicit expansion
    close_itpset = all(distsq <= 4.0_RP * delta**2)  ! Behaves the same as Powell's version.
    !close_itpset = all(distsq <= 4.0_RP * rho**2)  ! Powell's code.
    !close_itpset = all(distsq <= max((2.0_RP * delta)**2, (10.0_RP * rho)**2))  ! Powell's BOBYQA code.
    !----------------------------------------------------------------------------------------------!
    adequate_geo = (shortd .and. accurate_mod) .or. close_itpset

    ! Comments on ACCURATE_MOD:
    ! 1. ACCURATE_MOD is needed only when SHORTD is TRUE.
    ! 2. In Powell's UOBYQA code, ACCURATE_MOD is defined according to (28), (37), and (38) in the
    ! UOBYQA paper. As elaborated in Sec. 3 of the paper, the idea is to test whether the current
    ! model is sufficiently accurate by checking whether the interpolation error bound in (28) is
    ! (sufficiently) small. If the bound is small, then set ACCURATE_MOD to TRUE. Otherwise, it
    ! identifies a "bad" interpolation point that makes a significant contribution to the bound,
    ! with a preference to the interpolation points that are a far away from the current
    ! trust-region center. Such a point will be replaced with a new point obtained by the geometry
    ! step. If all the interpolation points are close enough to the trust-region center, then they
    ! are all considered to be good.
    ! 3. Our implementation defines ACCURATE_MOD by a method from NEWUOA and BOBYQA, which is also
    ! reflected in LINCOA. It sets ACCURATE_MOD to TRUE if recent model errors and step lengths are
    ! all small. In addition, it identifies a "bad" interpolation point by simply taking the
    ! farthest point from the current trust region center, unless they are all close enough to the
    ! center. This implementation is much simpler and it performs almost the same as Powell's
    ! original implementation.

    ! Powell's original definition of IMPROVE_GEO and REDUCE_RHO:
    !bad_trstep = (shortd .or. knew_tr == 0 .or. (ratio <= 0 .and. dnorm <= 2.0_RP*rho .and. ddmove <= 4.0_RP * rho**2))
    !improve_geo = bad_trstep .and. .not. (shortd .and. accurate_mod) .and. .not. close_itpset
    !reduce_rho = bad_trstep .and. (dnorm <= rho) .and. (.not. improve_geo)

    ! BAD_TRSTEP for IMPROVE_GEO.
    bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= TENTH .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    improve_geo = bad_trstep .and. .not. adequate_geo
    ! BAD_TRSTEP for REDUCE_RHO.
    bad_trstep = (shortd .or. (.not. qred > 0) .or. (ratio <= 0 .and. ddmove <= 4.0_RP * delta**2) .or. knew_tr == 0)
    reduce_rho = bad_trstep .and. adequate_geo .and. small_trrad


    !bad_trstep = (shortd .or. (ratio <= TENTH .and. ddmove <= 4.0_RP * rho**2) .or. knew_tr == 0)  ! For IMPROVE_GEO
    !!bad_trstep = (shortd .or. (ratio <= TENTH .and. dnorm <= 2.0 * rho .and. ddmove <= 4.0_RP * rho**2) .or. knew_tr == 0)  ! For IMPROVE_GEO
    !improve_geo = bad_trstep .and. .not. (shortd .and. accurate_mod) .and. .not. close_itpset
    !bad_trstep = (shortd .or. (ratio <= 0 .and. ddmove <= 4.0_RP * rho**2) .or. knew_tr == 0)  ! For REDUCE_RHO
    !!bad_trstep = (shortd .or. (ratio <= 0 .and. dnorm <= 2.0 * rho .and. ddmove <= 4.0_RP * rho**2) .or. knew_tr == 0)  ! For REDUCE_RHO
    !reduce_rho = bad_trstep .and. small_trrad .and. (.not. improve_geo)



    !! The following REDUCE_RHO and IMPROVE_GEO are adopted from NEWUOA/BOBYQA/LINCOA.
    !bad_trstep = (shortd .or. (ratio <= 0 .and. ddmove <= 4.0_RP * rho**2) .or. knew_tr == 0)  ! For REDUCE_RHO
    !!bad_trstep = (shortd .or. ratio <= 0 .or. knew_tr == 0)  ! OK, but not as good as the above one.
    !reduce_rho = (shortd .and. accurate_mod) .or. (bad_trstep .and. close_itpset .and. small_trrad)

    !! It is critical to include DMOVE <= 4.0_RP*RHO**2 in the following definition of BAD_TRSTEP.
    !bad_trstep = (shortd .or. (ratio <= TENTH .and. ddmove <= 4.0_RP * rho**2) .or. knew_tr == 0)  ! For IMPROVE_GEO
    !!bad_trstep = (shortd .or. ratio <= TENTH .or. knew_tr == 0)  ! This works poorly!
    !improve_geo = bad_trstep .and. (.not. close_itpset) .and. (.not. reduce_rho)

    if (improve_geo) then
        knew_geo = int(maxloc(distsq, dim=1), kind(knew_geo))
        g = pl(1:n, knew_geo) + smat_mul_vec(pl(n + 1:npt - 1, knew_geo), xopt)
        h = vec2smat(pl(n + 1:npt - 1, knew_geo))

        ! DELBAR is the trust-region radius for the geometry improvement subproblem.
        ! Powell's UOBYQA code sets DELBAR = RHO, but NEWUOA/BOBYQA/LINCOA all take DELTA and/or
        ! DISTSQ into consideration. The following DELBAR is copied from NEWUOA, and it seems to
        ! improve the performance slightly according to a test on 20220720.
        delbar = max(min(TENTH * sqrt(maxval(distsq)), HALF * delta), rho)
        d = geostep(g, h, delbar)
        dnorm = min(delbar, sqrt(sum(d**2)))
        dnormsav = [dnormsav(2:size(dnormsav)), dnorm]

        ! Calculate the next value of the objective function.
        xnew = xopt + d
        x = xbase + xnew

        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if

        if (is_nan(sum(abs(x)))) then
            f = sum(x) ! Set F to NaN
            if (nf == 1) then
                fopt = f
                xopt = ZERO
            end if
            info = NAN_INF_X
            exit
        end if
        call evaluate(calfun, x, f)
        nf = nf + 1
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
            xopt = xnew
            fopt = f
            exit
        end if

        ! Use the quadratic model to predict the change in F due to the step D, and find the values
        ! of the Lagrange functions at the new point.
        qred = -quadinc(pq, d, xopt)
        moderrsav = [moderrsav(2:size(moderrsav)), f - fopt + qred]
        vlag = calvlag(pl, d, xopt, kopt)

        diff = f - fopt + qred

        ! Update FOPT and XOPT if the new F is the least value of the objective function so far.
        ! Then branch if D is not a trust region step.
        fsave = fopt
        if (f < fopt) then
            fopt = f
            xopt = xnew
        end if

        ! Replace the interpolation point that has index KNEW by the point XNEW, and also update
        ! the Lagrange functions and the quadratic model.
        xpt(:, knew_geo) = xnew
        ! It can happen that VLAG(KNEW) = 0 due to rounding.
        pl(:, knew_geo) = pl(:, knew_geo) / vlag(knew_geo)
        plknew = pl(:, knew_geo)
        pq = pq + diff * plknew
        pl = pl - outprod(plknew, vlag)
        pl(:, knew_geo) = plknew

        ! Update KOPT if F is the least calculated value of the objective function.
        if (f < fsave) then
            kopt = knew_geo
        end if
    end if

    if (reduce_rho) then
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        ! Prepare to reduce RHO by shifting XBASE to the best point so far, and make the
        ! corresponding changes to the gradients of the Lagrange functions and the quadratic model.
        xbase = xbase + xopt
        xpt = xpt - spread(xopt, dim=2, ncopies=npt)
        pq(1:n) = pq(1:n) + smat_mul_vec(pq(n + 1:npt - 1), xopt)  ! Model gradient
        do k = 1, npt
            pl(1:n, k) = pl(1:n, k) + smat_mul_vec(pl(n + 1:npt - 1, k), xopt)  ! Lagrange fun. gradient
        end do
        ! Pick the next values of RHO and DELTA.
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
        dnormsav = HUGENUM
        moderrsav = HUGENUM
    end if
end do

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before.
! Zaikun 20220531: For the moment, D may contain NaN. Should be avoided later.
if (info == SMALL_TR_RADIUS .and. shortd .and. nf < maxfun .and. is_finite(sum(abs(d)))) then
    x = xbase + (xopt + d)
    call evaluate(calfun, x, f)
    nf = nf + 1
    call savehist(nf, x, xhist, f, fhist)
end if

if (fopt <= f .or. is_nan(f)) then
    x = xbase + xopt
    f = fopt
end if

call rangehist(nf, xhist, fhist)

!====================!
!  Calculation ends  !
!====================!

! Postconditions

close (16)


end subroutine uobyqb


end module uobyqb_mod
