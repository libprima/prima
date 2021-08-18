module cobylb_mod

implicit none
private
public :: cobylb


contains

subroutine cobylb(x, rhobeg, rhoend, iprint, maxfun, con, f, info, ftarget, cstrv, ctol, nsavmax)

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH, EPS, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, MAXTR_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F, &
   & DAMAGING_ROUNDING
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : inprod, matprod, outprod
use memory_mod, only : cstyle_sizeof
use hist_mod, only : savehist, selectx

! Solver-specific modules
use initialize_mod, only : initialize
use trustregion_mod, only : trstlp
use update_mod, only : updatepole, findpole
use geometry_mod, only : goodgeo, setdrop_geo, setdrop_tr, geostep

implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: nsavmax
real(RP), intent(in) :: ctol
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(out) :: con(:) ! M
real(RP), intent(inout) :: x(:)  ! N

! Outputs
integer(IK), intent(out) :: info
real(RP), intent(out) :: f

! Local variables
integer(IK) :: i
integer(IK) :: tr
integer(IK) :: maxtr
integer(IK) :: j
integer(IK) :: jdrop
integer(IK) :: kopt
integer(IK) :: m
integer(IK) :: n
integer(IK) :: nf
integer(IK) :: nsav
integer(IK) :: subinfo
real(RP) :: A(size(x), size(con) + 1)
! A(:, 1:M) contains the approximate gradient for the constraints, and A(:, M+1) is minus the
! approximate gradient for the objective function.
real(RP) :: b(size(con) + 1)
real(RP) :: barmu
real(RP) :: cmax(size(con))
real(RP) :: cmin(size(con))
real(RP) :: cpen  ! Penalty parameter for constraint in merit function (PARMU in Powell's code)
real(RP) :: conmat(size(con), size(x) + 1)
real(RP) :: consav(size(con), max(nsavmax, 0))
real(RP) :: csav(max(nsavmax, 0))
real(RP) :: fsav(max(nsavmax, 0))
real(RP) :: denom
real(RP) :: d(size(x))
real(RP) :: fval(size(x) + 1)
real(RP) :: cval(size(x) + 1)
real(RP) :: factor_alpha
real(RP) :: factor_beta
real(RP) :: factor_delta
real(RP) :: factor_gamma
real(RP) :: prerec  ! Predicted reduction in constraint violation
real(RP) :: preref  ! Predicted reduction in objective function
real(RP) :: prerem  ! Predicted reduction in merit function
real(RP) :: cstrv
real(RP) :: rho
real(RP) :: sim(size(x), size(x) + 1)  ! (n, )
real(RP) :: simi(size(x), size(x))  ! (n, )
real(RP) :: simid(size(x))
real(RP) :: simi_jdrop(size(x))
real(RP) :: actrem
real(RP) :: xsav(size(x), max(nsavmax, 0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY
real(RP), allocatable :: xhist(:, :)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: conhist(:, :)
real(RP), allocatable :: chist(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical :: alltrue(size(x) + 1)
logical :: bad_trstep
logical :: good_geo
logical :: improve_geo
logical :: reduce_rho
logical :: shortd
character(len=SRNLEN), parameter :: srname = 'COBYLB'

alltrue = .true.
reduce_rho = .false.

m = size(con)
n = size(x)

! Set the initial values of some parameters. The last column of SIM holds the optimal vertex of the
! current simplex, and the preceding N columns hold the displacements from the optimal vertex to the
! other vertices.  Further, SIMI holds the inverse of the matrix that is contained in the first N
! columns of SIM.
factor_alpha = QUART
factor_beta = 2.1E0_RP
factor_delta = 1.1E0_RP
factor_gamma = HALF
rho = rhobeg
cpen = ZERO

nsav = 0
fsav = HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.
consav = -HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.
csav = HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.
fval = HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.
conmat = -HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.
cval = HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.

call initialize(iprint, maxfun, ctol, ftarget, rho, x, nf, conmat, cval, fval, sim, simi, subinfo)

if (subinfo == NAN_X .or. subinfo == NAN_INF_F .or. subinfo == FTARGET_ACHIEVED .or. subinfo == MAXFUN_REACHED) then
    info = subinfo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sim(:, 1:n) = sim(:, 1:n) + spread(sim(:, n + 1), dim=2, ncopies=n)  !!! TEMPORARY
    xhist = sim
    fhist = fval
    conhist = conmat
    chist = cval
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cpen = min(1.0E8_RP, HUGENUM)
    ! Return the best calculated values of the variables.
    kopt = selectx(fhist, chist, cpen, ctol)
    x = xhist(:, kopt)
    f = fhist(kopt)
    con = conhist(:, kopt)
    cstrv = chist(kopt)
    !close (16)
    return
else
    x = sim(:, n + 1)
    f = fval(n + 1)
    con = conmat(:, n + 1)
    cstrv = cval(n + 1)
end if

! Normally, each trust-region iteration takes one function evaluation. The following setting
! essentially imposes no constraint on the maximal number of trust-region iterations.
maxtr = 10 * maxfun
! MAXTR is unlikely to be reached, but we define the following default value for INFO for safety.
info = MAXTR_REACHED

! Begin the iterative procedure.
! After solving a trust-region subproblem, COBYLA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust-region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust-region iteration? If yes, a geometry step
! will be taken, corresponding to the Branch (Delta) in the COBYLA paper.
! REDUCE_RHO - Will we reduce rho after the trust-region iteration?
! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
do tr = 1, maxtr
    ! Before the trust-region step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
    call updatepole(cpen, alltrue, conmat, cval, fval, sim, simi, subinfo)
    if (subinfo == DAMAGING_ROUNDING) then
        info = subinfo
        exit
    end if

    ! Does the current interpolation set has good geometry? It affects IMPROVE_GEO and REDUCE_RHO.
    good_geo = goodgeo(factor_alpha, factor_beta, rho, sim, simi)

    ! Calculate the linear approximations to the objective and constraint functions, placing minus
    ! the objective function gradient after the constraint gradients in the array A.
    ! N.B.:
    ! 1. When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result
    ! as Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
    ! finite-precision arithmetic.
    ! 2. TRSTLP accesses A mostly by columns, so it is not more reasonable to save A^T instead of A.
    A(:, 1:m) = transpose(matprod(conmat(:, 1:n) - spread(conmat(:, n + 1), dim=2, ncopies=n), simi))
    A(:, m + 1) = matprod(fval(n + 1) - fval(1:n), simi)

    ! Exit if A contains NaN. Otherwise, TRSTLP may encounter memory errors or infinite loops.
    !----------------------------------------------------------------------------------------------!
    ! POSSIBLE IMPROVEMENT: INSTEAD OF EXITING, SKIP A TRUST-REGION STEP AND PERFORM A GEOMETRY ONE!
    !----------------------------------------------------------------------------------------------!
    if (any(is_nan(A))) then
        info = -3
        exit
    end if

    ! In theory (but not computation), the last entry of B can be any number.
    b = [-conmat(:, n + 1), -fval(n + 1)]
    ! Calculate the trust-region trial step D.
    d = trstlp(A, b, rho)

    ! Is the trust-region trial step short?
    shortd = (inprod(d, d) < QUART * rho * rho)

    if (.not. shortd) then
        ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
        preref = inprod(d, A(:, m + 1))
        prerec = cval(n + 1) - maxval([-conmat(:, n + 1) - matprod(d, A(:, 1:m)), ZERO])

        ! Increase CPEN if necessary and branch back if this change alters the optimal vertex.
        ! Otherwise, PREREM will be set to the predicted reductions in the merit function.
        ! See the discussions around equation (9) of the COBYLA paper.
        barmu = -preref / prerec   ! PREREF + BARMU * PREREC = 0
        !!!!!!!!!!!!!!! Is it possible that PREREC <= 0????????????? It seems yes.
        if (prerec > ZERO .and. cpen < 1.5E0_RP * barmu) then
            cpen = min(TWO * barmu, HUGENUM)
            if (findpole(cpen, alltrue, cval, fval) <= n) then
                cycle
            end if
        end if

        prerem = preref + cpen * prerec   ! Is it positive????

        ! Set X.
        x = sim(:, n + 1) + d
        if (any(is_nan(x))) then
            !!!!!!!!!!!!! WHAT ABOUT NF ??? There is inconsistency. Also in NEWUOA and others. !!!!!!!!!!!
            f = sum(x)  ! Set F to NaN.
            con = f  ! Set constraint values to NaN.
            cstrv = f  ! Set constraint violation to NaN
            info = -1
            exit
        end if

        ! Evaluate the objective function and constraints at X.
        call calcfc(n, m, x, f, con)
        nf = nf + 1
        cstrv = maxval([-con, ZERO])

        ! Begin the operations that decide whether X should replace one of the vertices of the
        ! current simplex, the change being mandatory if ACTREM is positive.
        actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv)
        if (cpen <= ZERO .and. abs(f - fval(n + 1)) <= ZERO) then
            prerem = prerec   ! Is it positive?????
            actrem = cval(n + 1) - cstrv
        end if


        ! Set JDROP to the index of the vertex that is to be replaced by X.
        jdrop = setdrop_tr(actrem, d, factor_alpha, factor_delta, rho, sim, simi)

        ! When JDROP=0, the algorithm decides not to include X into the simplex.
        if (jdrop == 0) then
            call savehist(x, f, con, cstrv, xsav, fsav, consav, csav, nsav, ctol)   !?????
        else
            call savehist(sim(:, n + 1) + sim(:, jdrop), fval(jdrop), conmat(:, jdrop), cval(jdrop), &
                & xsav, fsav, consav, csav, nsav, ctol)
            ! Revise the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and CVAL.
            sim(:, jdrop) = d
            simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d)
            simi = simi - outprod(matprod(simi, d), simi_jdrop)
            simi(jdrop, :) = simi_jdrop
            fval(jdrop) = f
            conmat(:, jdrop) = con
            cval(jdrop) = cstrv
        end if

        ! TODO: The following CAN BE MOVED UPWARD to below CALCFC once XHIST etc are implemented !!!
        if (is_nan(F) .or. is_posinf(F)) then
            info = -2
            exit
        end if
        if (any(is_nan(con))) then
            cstrv = sum(con)  ! Set CSTRV to NaN
            cval(jdrop) = cstrv
            info = -2
            exit
        end if
        if (f <= ftarget .and. cstrv <= ctol) then
            info = 1
            exit
        end if
        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    end if

    ! Should we take a geometry step to improve the geometry of the interpolation set?
    ! N.B.: THEORETICALLY, JDROP > 0 when ACTREM > 0, and hence the definition of BAD_TRSTEP is
    ! mathematically equivalent to (SHORTD .OR. ACTREM <= ZERO .OR. ACTREM < TENTH * PREREM);
    ! however, JDROP may turn out 0 due to NaN even if ACTREM > 0. See SETDRTOP_TR for details.
    bad_trstep = (shortd .or. actrem <= ZERO .or. actrem < TENTH * prerem .or. jdrop == 0)
    improve_geo = bad_trstep .and. .not. good_geo

    ! Should we revise RHO (and CPEN)?
    reduce_rho = bad_trstep .and. good_geo

    if (improve_geo) then
        ! Before the geometry step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
        call updatepole(cpen, alltrue, conmat, cval, fval, sim, simi, subinfo)
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit
        end if

        ! If the current interpolation set has good geometry, then we skip the geometry step.
        ! The code has a small difference from the original COBYLA code here: If the current geometry
        ! is good, then we will continue with a new trust-region iteration; at the beginning of the
        ! iteration, CPEN may be updated, which may alter the pole point SIM(:, N + 1) by UPDATEPOLE;
        ! the quality of the interpolation point depends on SIM(:, N + 1), meaning that the same
        ! interpolation set may have good or bad geometry with respect to different "poles"; if the
        ! geometry turns out bad with the new pole, the original COBYLA code will take a geometry
        ! step, but the code here will NOT do it but continue to take a trust-region step.
        ! The argument is this: even if the geometry step is not skipped at the first place, the
        ! geometry may turn out bad again after the pole is altered due to an update to CPEN; should
        ! we take another geometry step in that case? If no, why should we do it here? Indeed, this
        ! distinction makes no practical difference for CUTEst problems with at most 100 variables
        ! and 5000 constraints, while the algorithm framework is simplified.
        if (.not. goodgeo(factor_alpha, factor_beta, rho, sim, simi)) then
            ! Decide a vertex to drop from the simplex. It will be replaced by SIM(:, N + 1) + D to
            ! improve acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
            jdrop = setdrop_geo(factor_alpha, factor_beta, rho, sim, simi)

            ! If JDROP = 0 (probably due to NaN in SIM or SIMI), then we exit. Without this, memory
            ! error may occur as JDROP will be used as an index of arrays.
            if (jdrop == 0) then
                info = DAMAGING_ROUNDING
                exit
            end if

            ! Calculate the geometry step D.
            d = geostep(jdrop, cpen, conmat, cval, fval, factor_gamma, rho, simi)

            ! Save the information of the JOPT-th vertex.
            call savehist(sim(:, n + 1) + sim(:, jdrop), fval(jdrop), conmat(:, jdrop), cval(jdrop), &
              & xsav, fsav, consav, csav, nsav, ctol)

            ! Update SIM and SIMI.
            sim(:, jdrop) = d  ! Corresponding to the new vertex SIM(:, N + 1) + D
            simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d)
            simi = simi - outprod(matprod(simi, d), simi_jdrop)
            simi(jdrop, :) = simi_jdrop

            ! Set X.
            x = sim(:, n + 1) + d
            if (any(is_nan(x))) then
                f = sum(x)  ! Set F to NaN.
                con = f  ! Set constraint values to NaN.
                cstrv = f  ! Set constraint violation to NaN.
                info = -1
                exit
            end if

            ! Evaluate the objective function and constraints at X.
            call calcfc(n, m, x, f, con)
            nf = nf + 1
            cstrv = maxval([-con, ZERO])
            fval(jdrop) = f
            conmat(:, jdrop) = con
            cval(jdrop) = cstrv

            if (is_nan(F) .or. is_posinf(F)) then
                info = -2
                exit
            end if
            if (any(is_nan(con))) then
                cstrv = sum(con)  ! Set CSTRV to NaN
                cval(jdrop) = cstrv
                info = -2
                exit
            end if
            if (f <= ftarget .and. cstrv <= ctol) then
                info = 1
                exit
            end if
            if (nf >= maxfun) then
                info = MAXFUN_REACHED
                exit
            end if
        end if
    end if

    if (reduce_rho) then  ! Update RHO and CPEN.
        if (rho <= rhoend) then
            info = 0
            exit
        end if
        ! See equation (11) in Section 3 of the COBYLA paper for the update of RHO.
        rho = HALF * rho
        if (rho <= 1.5E0_RP * rhoend) then
            rho = rhoend
        end if
        ! See equations (12)--(13) in Section 3 of the COBYLA paper for the update of CPEN.
        ! If the original CPEN = 0, then the updated CPEN is also 0.
        cmin = minval(conmat, dim=2)
        cmax = maxval(conmat, dim=2)
        if (any(cmin < HALF * cmax)) then
            denom = minval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
            cpen = min(cpen, (maxval(fval) - minval(fval)) / denom)
        else
            cpen = ZERO
        end if
    end if
end do

! Return the best calculated values of the variables.
sim(:, 1:n) = sim(:, 1:n) + spread(sim(:, n + 1), dim=2, ncopies=n)  !!! TEMPORARY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make sure that the history includes the last X.
!write (16, *) 'sim', sim
!write (16, *) 'xsav', xsav(:, 1:nsav)
!write (16, *) nf, f, x
xhist = reshape([sim, xsav(:, 1:nsav), x], [n, n + nsav + 2])
fhist = [fval, fsav(1:nsav), f]
conhist = reshape([conmat, consav(:, 1:nsav), con], [m, n + nsav + 2])
chist = [cval, csav(1:nsav), cstrv]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cpen = max(cpen, min(1.0E8_RP, HUGENUM))
kopt = selectx(fhist, chist, cpen, ctol)
x = xhist(:, kopt)
f = fhist(kopt)
con = conhist(:, kopt)
cstrv = chist(kopt)

!close (16)

end subroutine cobylb


end module cobylb_mod
