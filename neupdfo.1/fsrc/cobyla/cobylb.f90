module cobylb_mod

implicit none
private
public :: cobylb


contains


subroutine cobylb(calcfc, iprint, maxfun, ctol, ftarget, rhobeg, rhoend, constr, x, nf, chist, &
        & conhist, cstrv, f, fhist, xhist, info)

! Generic modules
use pintrf_mod, only : FUNCON
use evaluate_mod, only : evalfc
use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH, HUGENUM, DEBUGGING
use info_mod, only : INFO_DFT, MAXTR_REACHED, SMALL_TR_RADIUS, NAN_MODEL, DAMAGING_ROUNDING
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop, verisize
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : inprod, matprod, outprod, inv
use selectx_mod, only : savefilt, selectx
use history_mod, only : savehist
use checkexit_mod, only : checkexit
use resolution_mod, only : resenhance

! Solver-specific modules
use initialize_mod, only : initxfc, initfilt
use trustregion_mod, only : trstlp
use update_mod, only : updatexfc, updatepole, findpole
use geometry_mod, only : goodgeo, setdrop_geo, setdrop_tr, geostep

implicit none

! Inputs
procedure(FUNCON) :: calcfc
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ctol
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(out) :: constr(:) ! M
real(RP), intent(inout) :: x(:)  ! N

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)
real(RP), intent(out) :: conhist(:, :)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: xhist(:, :)

! Parameter
integer(IK), parameter :: maxfilt = 2000_IK  ! Must be positive. Recommended to be in [100, 10,000].

! Local variables
integer(IK) :: tr
integer(IK) :: maxtr
integer(IK) :: jdrop
integer(IK) :: kopt
integer(IK) :: m
integer(IK) :: maxxhist
integer(IK) :: maxfhist
integer(IK) :: maxconhist
integer(IK) :: maxchist
integer(IK) :: maxhist

integer(IK) :: n
integer(IK) :: nfilt
integer(IK) :: subinfo
real(RP) :: A(size(x), size(constr) + 1)
! A(:, 1:M) contains the approximate gradient for the constraints, and A(:, M+1) is minus the
! approximate gradient for the objective function.
real(RP) :: actrem
real(RP) :: b(size(constr) + 1)
real(RP) :: barmu
real(RP) :: cfilt(min(max(maxfilt, 0), maxfun))
real(RP) :: confilt(size(constr), size(cfilt))
real(RP) :: conmat(size(constr), size(x) + 1)
real(RP) :: cpen  ! Penalty parameter for constraint in merit function (PARMU in Powell's code)
real(RP) :: cval(size(x) + 1)
real(RP) :: d(size(x))
real(RP) :: factor_alpha
real(RP) :: factor_beta
real(RP) :: factor_delta
real(RP) :: factor_gamma
real(RP) :: ffilt(size(cfilt))
real(RP) :: fval(size(x) + 1)
real(RP) :: prerec  ! Predicted reduction in Constraint violation
real(RP) :: preref  ! Predicted reduction in objective Function
real(RP) :: prerem  ! Predicted reduction in Merit function
real(RP) :: rho
real(RP) :: sim(size(x), size(x) + 1)  ! (n, )
real(RP) :: simi(size(x), size(x))  ! (n, )
real(RP) :: xfilt(size(x), size(cfilt))
logical :: evaluated(size(x) + 1)
logical :: bad_trstep
logical :: good_geo
logical :: improve_geo
logical :: enhance_resolut
logical :: shortd
character(len=*), parameter :: srname = 'COBYLB'


! Get and verify the sizes.
m = size(constr)
n = size(x)
maxxhist = size(xhist, 2)
maxfhist = size(fhist)
maxconhist = size(conhist, 2)
maxchist = size(chist)
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(X) < 1')
    end if
    if (maxxhist > 0) then
        call verisize(xhist, n, maxhist)
    end if
    if (maxfhist > 0) then
        call verisize(fhist, maxhist)
    end if
    if (maxconhist > 0) then
        call verisize(conhist, m, maxhist)
    end if
    if (maxchist > 0) then
        call verisize(chist, maxhist)
    end if
end if

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

call initxfc(calcfc, iprint, maxfun, ctol, ftarget, rho, x, nf, chist, conhist, conmat, cval, fhist,&
   & fval, sim, xhist, evaluated, subinfo)
call initfilt(conmat, ctol, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)

if (subinfo /= INFO_DFT) then
    info = subinfo
    ! Return the best calculated values of the variables.
    ! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
    cpen = min(1.0E8_RP, HUGENUM)
    kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cpen, ctol)
    x = xfilt(:, kopt)
    f = ffilt(kopt)
    constr = confilt(:, kopt)
    cstrv = cfilt(kopt)
    !close (16)
    return
end if

! SIMI is the inverse of SIM(:, 1:N)
simi = inv(sim(:, 1:n), 'ltri')
! If we arrive here, the objective and constraints must have been evaluated at SIM(:, I) for all I.
evaluated = .true.

! Initialize PREREM, ACTREM, and JDROP, or some compilers will complain that they are uninitialized
! when setting BAD_TRSTEP. Indeed, these values will not be used, because they will be overwritten
! when SHORTD = FALSE.
prerem = ONE
actrem = -ONE
jdrop = 0_IK

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
    call updatepole(cpen, evaluated, conmat, cval, fval, sim, simi, subinfo)
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
    ! HOW EXACTLY?????
    !----------------------------------------------------------------------------------------------!
    ! POSSIBLE IMPROVEMENT: INSTEAD OF EXITING, SKIP A TRUST-REGION STEP AND PERFORM A GEOMETRY ONE!
    !----------------------------------------------------------------------------------------------!
    if (any(is_nan(A))) then
        info = NAN_MODEL
        exit
    end if

    ! Theoretically (but not numerically), the last entry of B does not affect the result of TRSTLP.
    ! We set it to -FVAL(N + 1) following Powell's code.
    b = [-conmat(:, n + 1), -fval(n + 1)]
    ! Calculate the trust-region trial step D.
    d = trstlp(A, b, rho)

    ! Is the trust-region trial step short?
    shortd = (inprod(d, d) < QUART * rho * rho)

    if (.not. shortd) then
        ! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
        preref = inprod(d, A(:, m + 1))
        prerec = cval(n + 1) - maxval([-matprod(d, A(:, 1:m)) - conmat(:, n + 1), ZERO])

        ! Increase CPEN if necessary and branch back if this change alters the optimal vertex.
        ! See the discussions around equation (9) of the COBYLA paper.
        barmu = -preref / prerec   ! PREREF + BARMU * PREREC = 0
        !!!!!!!!!!!!!!! Is it possible that PREREC <= 0????????????? It seems yes, but why?
        if (prerec > ZERO .and. cpen < 1.5E0_RP * barmu) then
            cpen = min(TWO * barmu, HUGENUM)
            if (findpole(cpen, evaluated, cval, fval) <= n) then
                cycle
            end if
        end if

        x = sim(:, n + 1) + d
        ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
        call evalfc(calcfc, x, f, constr, cstrv)
        nf = nf + 1_IK
        ! Save X, F, CONSTR, CSTRV into the history.
        call savehist(nf, constr, cstrv, f, x, chist, conhist, fhist, xhist)
        ! Save X, F, CONSTR, CSTRV into the filter.
        call savefilt(constr, cstrv, ctol, f, x, nfilt, cfilt, confilt, ffilt, xfilt)

        ! Begin the operations that decide whether X should replace one of the vertices of the
        ! current simplex, the change being mandatory if ACTREM is positive.
        ! PREREM and ACTREM are the predicted and actual reductions in the merit function respectively.
        prerem = preref + cpen * prerec   ! Is it positive????
        actrem = (fval(n + 1) + cpen * cval(n + 1)) - (f + cpen * cstrv)
        if (cpen <= ZERO .and. f <= fval(n + 1) .and. f >= fval(n + 1)) then
            ! CPEN <= ZERO indeed means CPEN == ZERO, while A <= B .and. A >= B indeed mean A == B.
            ! We write them in this way to avoid compilers complaining about equality comparison
            ! between reals, which appears in the original code of Powell.
            prerem = prerec   ! Is it positive?????
            actrem = cval(n + 1) - cstrv
        end if
        if (is_nan(actrem)) then
            actrem = -hugenum
        end if
        ! Set JDROP to the index of the vertex that is to be replaced by X.
        ! N.B.: COBYLA never sets JDROP = N + 1.
        jdrop = setdrop_tr(actrem, d, factor_alpha, factor_delta, rho, sim, simi)
        ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP) is replaced by D.
        ! When JDROP=0, the algorithm decides not to include X into the simplex.
        ! N.B.: UPDATEXFC does NOT manipulate the simplex so that SIM(:, N+1) is the best vertex;
        ! that is the job of UPDATEPOLE, which is called before each trust-region/geometry step.
        if (jdrop > 0) then
            call updatexfc(jdrop, constr, cstrv, d, f, conmat, cval, fval, sim, simi)
        end if

        ! Check whether to exit.
        subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end if

    ! Should we take a geometry step to improve the geometry of the interpolation set?
    ! N.B.: THEORETICALLY, JDROP > 0 when ACTREM > 0, and hence the definition of BAD_TRSTEP is
    ! mathematically equivalent to (SHORTD .OR. ACTREM <= ZERO .OR. ACTREM < TENTH * PREREM);
    ! however, JDROP may turn out 0 due to NaN even if ACTREM > 0. See SETDRTOP_TR for details.
    !bad_trstep = (shortd .or. actrem <= ZERO .or. actrem < TENTH * prerem .or. jdrop == 0)
    bad_trstep = (shortd .or. actrem <= ZERO .or. actrem < TENTH * prerem)
    improve_geo = bad_trstep .and. .not. good_geo

    ! Should we enhance the resolution by reducing RHO?
    enhance_resolut = bad_trstep .and. good_geo

    if (improve_geo) then
        ! Before the geometry step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
        call updatepole(cpen, evaluated, conmat, cval, fval, sim, simi, subinfo)
        if (subinfo == DAMAGING_ROUNDING) then
            info = subinfo
            exit
        end if

        ! If the current interpolation set has good geometry, then we skip the geometry step.
        ! The code has a small difference from Powell's original code here: If the current geometry
        ! is good, then we will continue with a new trust-region iteration; at the beginning of the
        ! iteration, CPEN may be updated, which may alter the pole point SIM(:, N + 1) by UPDATEPOLE;
        ! the quality of the interpolation point depends on SIM(:, N + 1), meaning that the same
        ! interpolation set may have good or bad geometry with respect to different "poles"; if the
        ! geometry turns out bad with the new pole, the original COBYLA code will take a geometry
        ! step, but our code here will NOT do it but continue to take a trust-region step.
        ! The argument is this: even if the geometry step is not skipped at the first place, the
        ! geometry may turn out bad again after the pole is altered due to an update to CPEN; should
        ! we take another geometry step in that case? If no, why should we do it here? Indeed, this
        ! distinction makes no practical difference for CUTEst problems with at most 100 variables
        ! and 5000 constraints, while the algorithm framework is simplified.
        if (.not. goodgeo(factor_alpha, factor_beta, rho, sim, simi)) then
            ! Decide a vertex to drop from the simplex. It will be replaced by SIM(:, N + 1) + D to
            ! improve acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
            ! N.B.: COBYLA never sets JDROP = N + 1.
            jdrop = setdrop_geo(factor_alpha, factor_beta, rho, sim, simi)

            ! If JDROP = 0 (probably due to NaN in SIM or SIMI), then we exit. Without this, memory
            ! error may occur as JDROP will be used as an index of arrays.
            if (jdrop == 0) then
                info = DAMAGING_ROUNDING
                exit
            end if

            ! Calculate the geometry step D.
            d = geostep(jdrop, cpen, conmat, cval, fval, factor_gamma, rho, simi)
            x = sim(:, n + 1) + d
            ! Evaluate the objective and constraints at X, taking care of possible Inf/NaN values.
            call evalfc(calcfc, x, f, constr, cstrv)
            nf = nf + 1_IK
            ! Save X, F, CONSTR, CSTRV into the history.
            call savehist(nf, constr, cstrv, f, x, chist, conhist, fhist, xhist)
            ! Save X, F, CONSTR, CSTRV into the filter.
            call savefilt(constr, cstrv, ctol, f, x, nfilt, cfilt, confilt, ffilt, xfilt)
            ! Update SIM, SIMI, FVAL, CONMAT, and CVAL so that SIM(:, JDROP) is replaced by D.
            ! N.B.: UPDATEXFC does NOT manipulate the simplex so that SIM(:, N+1) is the best vertex;
            ! that is the job of UPDATEPOLE, which is called before each trust-region/geometry step.
            call updatexfc(jdrop, constr, cstrv, d, f, conmat, cval, fval, sim, simi)
            ! Check whether to exit.
            subinfo = checkexit(maxfun, nf, cstrv, ctol, f, ftarget, x)
            if (subinfo /= INFO_DFT) then
                info = subinfo
                exit
            end if
        end if
    end if

    if (enhance_resolut) then  ! Enhance the resolution of the algorithm, updating RHO and CPEN.
        if (rho <= rhoend) then
            info = SMALL_TR_RADIUS
            exit
        end if
        call resenhance(conmat, fval, rhoend, cpen, rho)
    end if
end do

! Return the best calculated values of the variables.
! N.B. SELECTX and FINDPOLE choose X by different standards. One cannot replace the other.
cpen = max(cpen, min(1.0E8_RP, HUGENUM))
kopt = selectx(ffilt(1:nfilt), cfilt(1:nfilt), cpen, ctol)
x = xfilt(:, kopt)
f = ffilt(kopt)
constr = confilt(:, kopt)
cstrv = cfilt(kopt)

!close (16)

end subroutine cobylb


end module cobylb_mod

! TODO:
! 0. In COBYLA, check what should we do with JDROP = 0, both TR and GEO
!    If ACTREM > 0 ===> JDROP > 0, why can COBYLA return sub-optimal points???
! 1. evalfc, extreme barrier, moderate excessively negative objective, which has not been done in
!    NEWUOA. Shouldn't we remove the extreme barrier in the MATLAB/Python interface after it is
!    implemented in FORTRAN?
! 6. Do the same for NEWUOA
! 8. knew ===> jdrop
! 9. XPT(i, :), FVAL(:) should be indexed by j; k is for iterations (fhist and ffilt)
! 10. Check MAXVAL/MAXLOC --- are they affected by NaN?
! 11.
! Enforcing programming contracts
! Programming can be thought of as requirements for correct execution of a procedure and assurances
! for the result of correct execution. The requirements and assurances might be constraints of three
! kinds:
! Preconditions (requirements): logical expressions that must evaluate to .true. when a procedure starts execution,
! Postconditions (assurances): expressions that must evaluate to .true. when a procedure finishes execution, and
! Invariants: universal pre- and postconditions that must always be true when all procedures in a class start or finish executing.
