! NEWUOB_MOD is a module that performs the major calculations of NEWUOA.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Last Modified: Sunday, June 06, 2021 PM05:39:05

module newuob_mod

implicit none
private
public :: newuob


contains


subroutine newuob(calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, x, nf, f, fhist, xhist, info)
! NEWUOB performs the actual calculations of NEWUOA. The arguments IPRINT, MAXFUN, MAXHIST, NPT,
! ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, X, NF, F, FHIST, XHIST, and INFO are
! identical to the corresponding arguments in subroutine NEWUOA.

! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
! of the model and Lagrange functions.
! XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
! the least calculated F so far).
! D is reserved for trial steps from XOPT.
! XNEW = XOPT+D, corresponding to the vector of variables for the next calculation of F.
! [XPT, FVAL, KOPT] describes the interpolation set:
! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
! values of F at the interpolation points; KOPT is the index of XOPT in XPT (XPT(:,KOPT)=XOPT).
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

! Generic modules
use consts_mod, only : RP, IK, ZERO, HALF, TENTH, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : calquad, inprod

! Solver-specific modules
use pintrf_mod, only : FUNEVAL
use initialize_mod, only : initxf, initq, inith
use trustregion_mod, only : trsapp, trrad
use geometry_mod, only : setremove, geostep
use shiftbase_mod, only : shiftbase
use vlagbeta_mod, only : vlagbeta
use update_mod, only : updateh, updateq, tryqalt

implicit none

! Inputs
procedure(FUNEVAL) :: calfun
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

! In-output
real(RP), intent(inout) :: x(:)  ! SIZE(X) = N

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: xhist(:, :)

! Local variables
integer(IK) :: idz
integer(IK) :: ij(2, npt)
integer(IK) :: itest
integer(IK) :: khist
integer(IK) :: knew_geo
integer(IK) :: knew_tr
integer(IK) :: kopt
integer(IK) :: maxfhist
integer(IK) :: maxtr
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: subinfo
integer(IK) :: tr
real(RP) :: beta
real(RP) :: bmat(size(x), npt + size(x))
real(RP) :: crvmin
real(RP) :: d(size(x))
real(RP) :: delbar
real(RP) :: delta
real(RP) :: distsq
real(RP) :: dnorm
real(RP) :: dnormsave(3)
real(RP) :: fopt
real(RP) :: moderr
real(RP) :: fsave
real(RP) :: fval(npt)
real(RP) :: gq(size(x))
real(RP) :: hq(size(x), size(x))
real(RP) :: moderrsave(size(dnormsave))
real(RP) :: pq(npt)
real(RP) :: ratio
real(RP) :: rho
real(RP) :: rho_ratio
real(RP) :: trtol
real(RP) :: vlag(npt + size(x))
real(RP) :: vquad
real(RP) :: xbase(size(x))
real(RP) :: xdsq(npt)
real(RP) :: xnew(size(x))
real(RP) :: xopt(size(x))
real(RP) :: xpt(size(x), npt)
real(RP) :: zmat(npt, npt - size(x) - 1)
logical :: improve_geo
logical :: reduce_rho_1
logical :: reduce_rho_2
logical :: shortd
logical :: terminate
character(len=6), parameter :: solver = 'NEWUOA'
character(len=SRNLEN), parameter :: srname = 'NEWUOB'


! Get size.
n = int(size(x), kind(n))
maxfhist = int(size(fhist), kind(maxfhist))
maxxhist = int(size(xhist, 2), kind(maxxhist))

if (DEBUGGING) then
    if (n == 0) then
        call errstop(srname, 'SIZE(X) is invalid')
    end if
    if (size(xhist, 1) /= n .and. maxxhist > 0) then
        call errstop(srname, 'XHIST is nonempty but SIZE(XHIST, 1) /= SIZE(X)')
    end if
    if (maxfhist * maxxhist > 0 .and. maxfhist /= maxxhist) then
        call errstop(srname, 'FHIST and XHIST are both nonempty but SIZE(FHIST) /= SIZE(XHIST, 2)')
    end if
end if

maxtr = maxfun  ! Maximal number of trust region iterations.

! Initialize FVAL, XBASE, and XPT.
call initxf(calfun, iprint, x, rhobeg, ftarget, ij, kopt, nf, fhist, fval, xbase, xhist, xpt, subinfo)
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = xbase + xopt  ! Set X.
f = fopt  ! Set F.

! Check whether to return after initialization.
terminate = (subinfo == FTARGET_ACHIEVED) .or. (subinfo == NAN_X) .or. (subinfo == NAN_INF_F)

if (terminate) then
    info = subinfo
    if (abs(iprint) >= 1) then
        call retmssg(info, iprint, nf, f, x, solver)
    end if
    ! Rearrange FHIST and XHIST so that they are in the chronological order.
    if (maxfhist >= 1 .and. maxfhist < nf) then
        khist = mod(nf - 1_IK, maxfhist) + 1_IK
        fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
    end if
    if (maxxhist >= 1 .and. maxxhist < nf) then
        khist = mod(nf - 1_IK, maxxhist) + 1_IK
        xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
        ! The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
        ! order of Fortran arrays.
    end if
    return
end if

! Initialize GQ, HQ, and PQ.
call initq(ij, fval, xpt, gq, hq, pq, subinfo)

! Initialize BMAT and ZMAT, and IDZ.
call inith(ij, xpt, idz, bmat, zmat, subinfo)

! After initializing GQ, HQ, PQ, BMAT, ZMAT, one can also choose to return if subinfo = NAN_MODEL
! (NaN occurs in the model). We do not do it here. If such a model is harmful, then it will probably
! lead to other returns (NaN in X, NaN in F, trust region subproblem fails, ...); otherwise, the
! code will continue to run and possibly get rid of the NaN in the model.

! Set some more initial values and parameters.
rho = rhobeg
delta = rho
moderrsave = HUGENUM
dnormsave = HUGENUM
itest = 0
trtol = 1.0E-2_RP  ! Tolerance used in trsapp.

! Begin the iterative procedure.
! After solving a trust-region subproblem, NEWUOA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust region iteration?
! REDUCE_RHO - Will we reduce rho after the trust region iteration?
! REDUCE_RHO = REDUCE_RHO_1 .OR. REDUCE_RHO_2 (see boxes 14 and 10 of Fig. 1 in the NEWUOA paper).
! NEWUOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Solve the trust region subproblem. In Powell's NEWUOA code, VQUAD is not an output of TRSAPP.
    ! Here we output it but will NOT use it (for the moment); it will still be calculated later
    ! by CALQUAD in order to produce the same results as Powell's code.
    call trsapp(delta, gq, hq, pq, trtol, xopt, xpt, crvmin, vquad, d, subinfo)

    ! Calculate the length of the trial step D.
    dnorm = min(delta, sqrt(inprod(d, d)))

    ! SHORTD corresponds to box 3 of the NEWUOA paper.
    shortd = (dnorm < HALF * rho)
    ! REDUCE_RHO_1 corresponds to box 14 of the NEWUOA paper.
    reduce_rho_1 = shortd .and. (maxval(abs(moderrsave)) <= 0.125_RP * crvmin * rho * rho) .and. (maxval(dnormsave) <= rho)
    if (shortd .and. (.not. reduce_rho_1)) then
        ! Reduce DELTA. After this, DELTA < DNORM may hold.
        delta = TENTH * delta
        if (delta <= 1.5_RP * rho) then
            delta = rho  ! Set DELTA to RHO when it is close.
        end if
    end if

    if (.not. shortd) then  ! D is long enough.
        ! Save the current FOPT in FSAVE. It is needed later.
        fsave = fopt

        ! Shift XBASE if XOPT may be too far from XBASE.
        !if (inprod(d, d) <= 1.0e-3_RP*inprod(xopt, xopt)) then  ! Powell
        if (dnorm * dnorm <= 1.0E-3_RP * inprod(xopt, xopt)) then
            call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt)
        end if

        ! Calculate VLAG and BETA for D. It makes uses of XOPT, so this is done before updating XOPT.
        call vlagbeta(idz, kopt, bmat, d, xopt, xpt, zmat, beta, vlag)

        ! Use the current quadratic model to predict the change in F due to the step D.
        call calquad(d, gq, hq, pq, xopt, xpt, vquad)

        ! Calculate the next value of the objective function.
        xnew = xopt + d
        x = xbase + xnew
        if (any(is_nan(x))) then
            f = sum(x)  ! Set F to NaN. It is necessary.
            info = NAN_X
            exit
        end if
        call calfun(x, f)
        nf = int(nf + 1, kind(nf))
        if (abs(iprint) >= 3) then
            call fmssg(iprint, nf, f, x, solver)
        end if
        if (maxfhist >= 1) then
            khist = mod(nf - 1_IK, maxfhist) + 1_IK
            fhist(khist) = f
        end if
        if (maxxhist >= 1) then
            khist = mod(nf - 1_IK, maxxhist) + 1_IK
            xhist(:, khist) = x
        end if

        ! DNORMSAVE constains the DNORM of the latest 3 function evaluations with the current RHO.
        dnormsave = [dnorm, dnormsave(1:size(dnormsave) - 1)]

        ! MODERR is the error of the current model in predicting the change in F due to D.
        moderr = f - fsave - vquad
        ! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
        moderrsave = [moderr, moderrsave(1:size(moderrsave) - 1)]

        ! Update FOPT and XOPT
        if (f < fopt) then
            fopt = f
            xopt = xnew
        end if

        ! Check whether to exit
        if (is_nan(f) .or. is_posinf(f)) then
            info = NAN_INF_F
            exit
        end if
        if (f <= ftarget) then
            info = FTARGET_ACHIEVED
            exit
        end if
        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if

        ! Calculate the reduction ratio and update DELTA accordingly.
        if (is_nan(vquad) .or. vquad >= ZERO) then
            info = TRSUBP_FAILED
            exit
        end if
        ratio = (f - fsave) / vquad
        ! Update DELTA. After this, DELTA < DNORM may hold.
        delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)
        if (delta <= 1.5_RP * rho) then
            delta = rho
        end if

        ! Set KNEW_TR to the index of the interpolation point that will be replaced by XNEW. KNEW_TR
        ! will ensure that the geometry of XPT is "good enough" after the replacement. Note that the
        ! information of XNEW is included in VLAG and BETA, which are calculated according to
        ! D = XNEW - XOPT. KNEW_TR = 0 means it is impossible to obtain a good interpolation set
        ! by replacing any current interpolation point with XNEW.
        call setremove(idz, kopt, beta, delta, ratio, rho, vlag(1:npt), xopt, xpt, zmat, knew_tr)

        if (knew_tr > 0) then
            ! If KNEW_TR > 0, then update BMAT, ZMAT and IDZ, so that the KNEW_TR-th interpolation
            ! point is replaced by XNEW. If KNEW_TR = 0, then probably the geometry of XPT needs
            ! improvement, which will be handled below.
            call updateh(knew_tr, beta, vlag, idz, bmat, zmat)

            ! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
            call updateq(idz, knew_tr, bmat(:, knew_tr), moderr, zmat, xpt(:, knew_tr), gq, hq, pq)

            ! Include the new interpolation point. This should be done after updating the model.
            fval(knew_tr) = f
            xpt(:, knew_tr) = xnew

            ! Update KOPT to KNEW_TR if F < FSAVE (i.e., last FOPT).
            if (f < fsave) then
                kopt = knew_tr
            end if

            ! Test whether to replace the new quadratic model Q by the least-Frobenius norm
            ! interpolant Q_alt. Perform the replacement if certain criteria are satisfied. This
            ! part is OPTIONAL, but it is crucial for the performance on some problems. See
            ! Section 8 of the NEWUOA paper.
            ! In NEWUOA, TRYQALT is called only after a trust-region step but not after a geometry
            ! step. Maybe this is because the model is expected to be good after a geometry step.
            if (delta <= rho) then  ! DELTA == RHO.
                ! In theory, the FVAL - FSAVE in the following line can be replaced by FVAL + C
                ! with any constant C. This constant will not affect the result in precise
                ! arithmetic. Powell chose C = - FVAL(KOPT_ORIGINAL), where KOPT_ORIGINAL is the
                ! KOPT before the update above (i.e., Powell updated KOPT after TRYQALT). Here we
                ! use the updated KOPT, because it worked slightly better on CUTEst, although there
                ! is no difference theoretically. Note that FVAL(KOPT_ORIGINAL) may not equal FSAVE
                ! --- it may happen that KNEW_TR = KOPT_ORIGINAL so that FVAL(KOPT_ORIGINAL) has
                ! been revised after the last function evaluation.
                ! Question: Since TRYQALT is invoked only when DELTA equals the current RHO, why not
                ! reset ITEST to 0 when RHO is reduced?
                call tryqalt(idz, fval - fval(kopt), ratio, bmat(:, 1:npt), zmat, itest, gq, hq, pq)
            end if
        end if
    end if  ! End of if (.not. shortd)

    ! Before next trust region iteration, we may improve the geometry of XPT or reduce rho
    ! according to IMPROVE_GEO and REDUCE_RHO. Now we decide these two indicators.

    ! Define IMPROVE_GEO, corresponding to box 8 of the NEWUOA paper.
    improve_geo = .false.
    ! The geometry of XPT will be improved in three cases specified below. Above all,if
    ! REDUCE_RHO_1 = TRUE, meaning that the step is short and the latest model errors have been
    ! small, then we do not need to improve the geometry; instead, RHO will be reduced.
    ! 1. the trust-region step is too short (SHORTD = TRUE), or
    ! 2. it is impossible to obtain an interpolation set with good geometry by replacing a current
    ! interpolation point with the trust-region trial point (KNEW_TR = 0), or
    ! 3. the trust-region reduction ratio is small (RATION < TENTH).
    ! N.B.:
    ! 1. KNEW_TR and RATIO are both set if SHORTD = FALSE. So the expression
    ! (SHORTD .OR. KNEW_TR == 0 .OR. RATIO < TENTH) will not suffer from unset KNEW_TR or RATIO.
    ! 2. If REDUCE_RHO = FALSE and SHORTD = TRUE, then the trust-region step is not tried at all, 
    ! as no function evaluation is invoked at XOPT + D (If REDUCE_RHO = TRUE, then the trust-region
    ! step is not tried either, but the same step will be generated again at the next trust-region 
    ! iteration after RHO is reduced and DELTA is updated; see the last paragraph of Section 2 of 
    ! the NEWUOA paper).
    ! 3. If SHORTD = FALSE and KNEW_TR = 0, then the trust-region step invokes a function evaluation
    ! at XOPT + D, but [XOPT + D, F(XOPT +D)] is not included into [XPT, FVAL]. In other words, this
    ! function value is discarded. Note that KNEW_TR = 0 only if RATIO <= 0 (see SETREMOVE), so that
    ! a function value that renders a reduction is never discarded.
    ! 4. If SHORTD = FALSE and KNEW_TR > 0 and RATIO < TENTH, then [XPT, FVAL] is updated so that
    ! [XPT(KNEW_TR), FVAL(KNEW_TR)] = [XOPT + D, F(XOPT + D)], and the model is updated accordingly,
    ! but such a model will not be used in the next trust-region iteration, because a geometry step
    ! will be invoked to improve the geometry of the interpolation set and update the model again.
    if (.not. reduce_rho_1 .and. (shortd .or. knew_tr == 0 .or. ratio < TENTH)) then
        ! Find out if the interpolation points are close enough to the best point so far, i.e., all
        ! the points are within a ball centered at XOPT with a radius of 2*DELTA. If not, set
        ! KNEW_GEO to the index of the point that is the farthest away.
        distsq = 4.0_RP * delta * delta
        xdsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
        if (maxval(xdsq) > distsq) then
            knew_geo = int(maxloc(xdsq, dim=1), kind(knew_geo))
            distsq = maxval(xdsq)
        else
            knew_geo = 0
        end if

        ! If KNEW_GEO is positive (i.e., not all points are close to XOPT), then a geometry step
        ! (aka model step) will be taken to ameliorate the geometry of the interpolation set and
        ! hence improve the model.
        improve_geo = (knew_geo > 0)
    end if

    if (improve_geo) then
        ! Save the current FOPT in fsave. It is needed later.
        fsave = fopt

        ! Set DELBAR, which will be used as the trust region radius for the geometry-improving
        ! schemes GEOSTEP. We also need it to decide whether to shift XBASE or not.
        delbar = max(min(TENTH * sqrt(distsq), HALF * delta), rho)

        ! Shift XBASE if XOPT may be too far from XBASE.
        if (delbar * delbar <= 1.0E-3_RP * inprod(xopt, xopt)) then
            call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt)
        end if

        ! Find a step D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) is
        ! replaced by XOPT + D. The GEOSTEP subroutine will call Powell's BIGLAG and BIGDEN. It will
        ! also calculate the VLAG and BETA for this D.
        call geostep(idz, knew_geo, kopt, bmat, delbar, xopt, xpt, zmat, d, beta, vlag)

        ! Use the current quadratic model to predict the change in F due to the step D.
        call calquad(d, gq, hq, pq, xopt, xpt, vquad)

        ! Calculate the next value of the objective function.
        xnew = xopt + d
        x = xbase + xnew
        if (any(is_nan(x))) then
            f = sum(x)  ! Set F to NaN. It is necessary.
            info = NAN_X
            exit
        end if
        call calfun(x, f)
        nf = int(nf + 1, kind(nf))
        if (abs(iprint) >= 3) then
            call fmssg(iprint, nf, f, x, solver)
        end if
        if (maxfhist >= 1) then
            khist = mod(nf - 1_IK, maxfhist) + 1_IK
            fhist(khist) = f
        end if
        if (maxxhist >= 1) then
            khist = mod(nf - 1_IK, maxxhist) + 1_IK
            xhist(:, khist) = x
        end if

        ! DNORMSAVE constains the DNORM of the latest 3 function evaluations with the current RHO.
        !------------------------------------------------------------------------------------------!
        ! Powell's code does not update DNORM. Therefore, DNORM is the length of last trust-region
        ! trial step, which seems inconsistent with what is described in Section 7 (around (7.7)) of
        ! the NEWUOA paper. Seemingly we should keep DNORM = ||D|| as we do here. The value of DNORM
        ! will be used when defining REDUCE_RHO.
        dnorm = min(delbar, sqrt(inprod(d, d)))
        ! In theory, DNORM = DELBAR in this case.
        !------------------------------------------------------------------------------------------!
        dnormsave = [dnorm, dnormsave(1:size(dnormsave) - 1)]

        ! MODERR is the error of the current model in predicting the change in F due to D.
        moderr = f - fsave - vquad
        ! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
        moderrsave = [moderr, moderrsave(1:size(moderrsave) - 1)]

        ! Update FOPT and XOPT
        if (f < fopt) then
            fopt = f
            xopt = xnew
        end if

        ! Check whether to exit.
        if (is_nan(f) .or. is_posinf(f)) then
            info = NAN_INF_F
            exit
        end if
        if (f <= ftarget) then
            info = FTARGET_ACHIEVED
            exit
        end if
        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if

        ! Update BMAT, ZMAT and IDZ, so that the KNEW_GEO-th interpolation point can be moved.
        call updateh(knew_geo, beta, vlag, idz, bmat, zmat)

        ! Update the quadratic model.
        call updateq(idz, knew_geo, bmat(:, knew_geo), moderr, zmat, xpt(:, knew_geo), gq, hq, pq)

        ! Include the new interpolation point. This should be done after updating BMAT, ZMAT, and
        ! the model.
        fval(knew_geo) = f
        xpt(:, knew_geo) = xnew
        if (f < fsave) then
            kopt = knew_geo
        end if
    end if  ! The procedure of improving geometry ends.

    ! If all the interpolation points are close to XOPT (IMPROVE_GEO = FALSE) compared to rho, and
    ! the trust region is small, but the trust region step is "bad" (SHORTD or RATIO <= 0), then we
    ! should shrink RHO (i.e., update the standard for defining "closeness" and SHORTD).
    ! REDUCE_RHO_2 corresponds to box 10 of the NEWUOA paper. Note that DELTA < DNORM may hold due
    ! to the update of DELTA.
    reduce_rho_2 = (.not. improve_geo) .and. (max(delta, dnorm) <= rho) .and. (shortd .or. ratio <= 0)

    if (reduce_rho_1 .or. reduce_rho_2) then
        ! The calculations with the current RHO are complete. Pick the next values of RHO and DELTA.
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        else
            delta = HALF * rho
            rho_ratio = rho / rhoend
            if (rho_ratio <= 16.0_RP) then
                rho = rhoend
            else if (rho_ratio <= 250.0_RP) then
                rho = sqrt(rho_ratio) * rhoend
            else
                rho = TENTH * rho
            end if
            delta = max(delta, rho)
            ! DNORMSAVE and MODERRSAVE are corresponding to the latest 3 function evaluations with
            ! the current RHO. Update them after reducing RHO.
            dnormsave = HUGENUM
            moderrsave = HUGENUM
            if (abs(iprint) >= 2) then
                call rhomssg(iprint, nf, fopt, rho, xbase + xopt, solver)
            end if
        end if
    end if  ! The procedure of reducing RHO ends.

end do  ! The iterative procedure ends.

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before. Note that no trust region iteration has been done if MAXTR = 0, and hence we
! should not check whether SHORTD = TRUE but return immediately.
if (maxtr > 0 .and. shortd .and. nf < maxfun) then
    x = xbase + (xopt + d)
    if (any(is_nan(x))) then
        f = sum(x)  ! Set F to NaN. It is necessary.
        info = NAN_X
    else
        call calfun(x, f)
        nf = int(nf + 1, kind(nf))
        if (abs(iprint) >= 3) then
            call fmssg(iprint, nf, f, x, solver)
        end if
        if (maxfhist >= 1) then
            khist = mod(nf - 1_IK, maxfhist) + 1_IK
            fhist(khist) = f
        end if
        if (maxxhist >= 1) then
            khist = mod(nf - 1_IK, maxxhist) + 1_IK
            xhist(:, khist) = x
        end if
    end if
end if

! Note that (FOPT .LE. F) is FALSE if F is NaN; if F is NaN, it is also necessary to update X and F.
if (is_nan(f) .or. fopt <= f) then
    x = xbase + xopt
    f = fopt
end if

! Rearrange FHIST and XHIST so that they are in the chronological order.
if (maxfhist >= 1 .and. maxfhist < nf) then
    khist = mod(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if
if (maxxhist >= 1 .and. maxxhist < nf) then
    khist = mod(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
end if

if (abs(iprint) >= 1) then
    call retmssg(info, iprint, nf, f, x, solver)
end if

return
end subroutine newuob


end module newuob_mod
