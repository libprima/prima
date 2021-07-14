module cobylb_mod


contains

subroutine cobylb(m, x, rhobeg, rhoend, iprint, maxfun, con, f, info, ftarget, cstrv)

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F, DAMAGING_ROUNDING
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : inprod, matprod, outprod
use memory_mod, only : cstyle_sizeof
use logging_mod, only : logging

! Solver-specific modules
!use savex_mod, only : savex
use initialize_mod, only : initialize
use trustregion_mod, only : trstlp
use update_mod, only : updatepole
use selectx_mod, only : selectx, isbetter

implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: m
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(out) :: con(:) ! m+2. Bad name; should be confr
real(RP), intent(inout) :: x(:)  ! n

! Outputs
integer(IK), intent(out) :: info
real(RP), intent(out) :: f


! Parameters
! NSAVMAX is the maximal number of "dropped X" to save
integer(IK), parameter :: nsavmax = 1000_IK
! CTOL is the tolerance for constraint violation. A point X is considered to be feasible if its
! constraint violation (CSTRV) is less than CTOL.
real(RP), parameter :: ctol = epsilon(1.0_RP)

! Local variables

integer(IK) :: i
integer(IK) :: iact(m + 1)
integer(IK) :: ibrnch
integer(IK) :: ifull
integer(IK) :: iptem
integer(IK) :: j
integer(IK) :: jdrop
integer(IK) :: jopt
integer(IK) :: n
integer(IK) :: nf
integer(IK) :: nsav
integer(IK) :: subinfo
real(RP) :: A(size(x), m + 1)  ! Better name?
! A(:, 1:m) contains the approximate gradient for the constraints, and A(:, m+1) is minus the
! approximate gradient for the objective function.
real(RP) :: barmu
real(RP) :: cmax(m)
real(RP) :: cmin(m)
real(RP) :: consav(m + 2)
real(RP) :: cpen  ! Penalty parameter for constraint in merit function (PARMU in Powell's code)
real(RP) :: cvmaxm
real(RP) :: cvmaxp
real(RP) :: datmat(m + 2, size(x) + 1)  ! CONVAL, FVAL, CVVAL
real(RP) :: datsav(m + 2, max(nsavmax, 0))
real(RP) :: denom
real(RP) :: dx(size(x))
real(RP) :: dxsign
real(RP) :: edgmax
real(RP) :: factor_alpha
real(RP) :: factor_beta
real(RP) :: factor_delta
real(RP) :: factor_gamma
real(RP) :: pareta
real(RP) :: parsig
real(RP) :: phi(size(x) + 1)  ! Merit function values
real(RP) :: phimin
real(RP) :: prerec  ! Predicted reduction in constraint violation
real(RP) :: preref  ! Predicted reduction in objective function
real(RP) :: prerem  ! Predicted reduction in merit function
real(RP) :: ratio
real(RP) :: cstrv
real(RP) :: rho
real(RP) :: sigbar(size(x))
real(RP) :: sim(size(x), size(x) + 1)  ! (n, )
real(RP) :: simi(size(x), size(x))  ! (n, )
real(RP) :: simid(size(x))
real(RP) :: tmpv(size(x))
real(RP) :: simi_jdrop(size(x))
real(RP) :: actrem
real(RP) :: veta(size(x))
real(RP) :: vmnew
real(RP) :: vmold
real(RP) :: vsig(size(x))
real(RP) :: xsav(size(x), max(nsavmax, 0))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY
real(RP), allocatable :: xhist(:, :)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: conhist(:, :)
real(RP), allocatable :: cstrvhist(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical :: terminate
logical :: improve_geo
logical :: reduce_rho
logical :: shortd


character(len=SRNLEN), parameter :: srname = 'COBYLB'

reduce_rho = .false.
terminate = .false.

n = size(x)

! Set the initial values of some parameters. The last column of SIM holds the optimal vertex of the
! current simplex, and the preceding N columns hold the displacements from the optimal vertex to the
! other vertices.  Further, SIMI holds the inverse of the matrix that is contained in the first N
! columns of SIM.
info = 2147483647
iptem = min(n, 5)
factor_alpha = QUART
factor_beta = 2.1E0_RP
factor_delta = 1.1E0_RP
factor_gamma = HALF
rho = rhobeg
cpen = ZERO
!if (iprint >= 2) then
!print 10, RHO
!10  format(/3X, 'The initial value of RHO is', 1PE13.6, 2X, 'and CPEN is set to zero.')
!end if

nsav = 0
datsav = HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.
datmat = HUGENUM  ! This is necessary; otherwise, SELECTX may return an incorrect X.

call initialize(iprint, maxfun, ctol, ftarget, rho, x, nf, datmat, sim, simi, subinfo)
x = sim(:, n + 1)
f = datmat(m + 1, n + 1)
cstrv = datmat(m + 2, n + 1)
con = datmat(:, n + 1)
consav = con

if (subinfo == NAN_X .or. subinfo == NAN_INF_F .or. subinfo == FTARGET_ACHIEVED .or. &
    & subinfo == DAMAGING_ROUNDING .or. subinfo == MAXFUN_REACHED) then
    info = subinfo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sim(:, 1:n) = sim(:, 1:n) + spread(sim(:, n + 1), dim=2, ncopies=n)  !!! TEMPORARY
    ! Make sure that the history includes the last X.
    xhist = reshape([sim, xsav(:, 1:nsav), x], [n, n + nsav + 2])
    fhist = [datmat(m + 1, :), datsav(m + 1, 1:nsav), f]
    conhist = reshape([datmat(1:m, :), datsav(1:m, 1:nsav), consav], [m, n + nsav + 2])
    cstrvhist = [datmat(m + 2, :), datsav(m + 2, 1:nsav), cstrv]
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cpen = 1.0E6_RP
    ! It is necessary to call SELECTX, because INITIALIZE chooses SIM(:, N+1) according to the
    ! function value while neglecting the constraints.
    call selectx(conhist, cstrvhist, ctol, fhist, cpen, xhist, con(1:m), cstrv, f, x)
    return
end if

ibrnch = 1
! Identify the optimal vertex of the current simplex, and switch it to SIM(:, N+1) if it is not
! there yet. Powell called SIM(:, N+1) the Pole Position of the simplex.
140 call updatepole(cpen, [(.true., i=1, n + 1)], datmat, sim, simi, subinfo)
if (subinfo == DAMAGING_ROUNDING) then
    info = subinfo
!    if (IPRINT >= 1) print 210
!210 format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
    goto 600
end if

if (info == 3) then
    !write (10, *) '420 g600'
    goto 600
end if

! Calculate the coefficients of the linear approximations to the objective and constraint functions,
! placing minus the objective function gradient after the constraint gradients in the array A.
! When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result as
! Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
! finite-precision arithmetic.
! 220
con = -datmat(:, n + 1)  ! Why put a negative sign???????????????????????????
! Is it more reasonable to save A transpose instead of A? Better name for A?
A = transpose(matprod(datmat(1:m + 1, 1:n) - spread(datmat(1:m + 1, n + 1), dim=2, ncopies=n), simi))
A(:, m + 1) = -A(:, m + 1)

if (any(is_nan(A))) then
    info = -3
    goto 600
end if

! Calculate the values of sigma and eta, and set IFLAG=0 if the current simplex is not acceptable.
parsig = factor_alpha * rho
pareta = factor_beta * rho
! For J = 1, 2, ..., n, VSIG(J) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
improve_geo = any(vsig < parsig) .or. any(veta > pareta) .or. any(is_nan([vsig, veta]))
!---------------------------------------------------------------------------------------!
!improve_geo = any(vsig < parsig) .or. any(veta > pareta)
!---------------------------------------------------------------------------------------!

if (ibrnch == 0 .and. improve_geo) then

! If a new vertex is needed to improve acceptability, then decide which vertex to drop from the simplex.
    if (maxval(veta) > pareta) then
        jdrop = int(maxloc(veta, dim=1), kind(jdrop))
    elseif (minval(vsig) < pareta) then
        jdrop = int(minloc(vsig, dim=1), kind(jdrop))
    else
        jdrop = 0_IK
    end if

! If VETA or VSIG become NaN due to rounding errors, JDROP will end up being 0. If we continue, then
! a Segmentation Fault will happen because we will access SIM(:, JDROP) and VSIG(JDROP).
    if (jdrop == 0_IK) then
!    if (IPRINT >= 1) print 286
!286 format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
        INFO = 7
        goto 600
    end if

! Save the information of the JOPT-th vertex in XSAV and DATSAV.
    call savex(sim(:, n + 1) + sim(:, jdrop), datmat(:, jdrop), xsav, datsav, nsav, ctol)

!Calculate the step to the new vertex and its sign.
    dx = factor_gamma * rho * vsig(jdrop) * simi(jdrop, :)
    cvmaxp = maxval([ZERO, -matprod(dx, A(:, 1:m)) - datmat(1:m, n + 1)])
    cvmaxm = maxval([ZERO, matprod(dx, A(:, 1:m)) - datmat(1:m, n + 1)])
    dxsign = ONE
    if (cpen * (cvmaxp - cvmaxm) > TWO * inprod(dx, a(:, m + 1))) then
        dxsign = -ONE
    end if

! Update SIM and SIMI, and set the next X.
    dx = dxsign * dx
    sim(:, jdrop) = dx
    simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), dx)
    simi = simi - outprod(matprod(simi, dx), simi_jdrop)
    simi(jdrop, :) = simi_jdrop
    x = sim(:, n + 1) + dx

!write (10, *) '624 g40'
    if (any(is_nan(x))) then
        f = sum(x)  ! Set F to NaN.
        con = sum(x)  ! Set constraint values and constraint violation to NaN.
        info = -1
        goto 600
    end if

    call calcfc(n, m, x, f, con)
    nf = nf + 1
    cstrv = maxval([ZERO, -con(1:m)])
    con(m + 1) = f
    con(m + 2) = cstrv
    consav = con

!if (nf == IPRINT - 1 .or. IPRINT == 3) then
!    print 70, nf, F, CSTRV, (X(I), I=1, IPTEM)
!70  format(/3X, 'nf =', I5, 3X, 'F =', 1PE13.6, 4X, 'MAXCV =', 1PE13.6 / 3X, 'X =', 1PE13.6, 1P4E15.6)
!    if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!80  format(1PE19.6, 1P4E15.6)
!end if

! If the objective function value or the constraints contain a NaN or an infinite value, the exit.
    if (is_nan(F) .or. is_posinf(F)) then
        info = -2
        goto 600
    end if
    if (any(is_nan(con(1:m)))) then
        cstrv = sum(con(1:m))  ! Set CSTRV to NaN
        info = -2
        goto 600
    end if
! If the objective function achieves the target value at a feasible point, then exit.
    if (f <= ftarget .and. cstrv <= ctol) then
        info = 1
        return
    end if

    if (nf >= maxfun) then
!    if (IPRINT >= 1) print 50
!50  format(/3X, 'Return from subroutine COBYLA because the ', 'MAXFUN limit has been reached.')
        info = 3
    end if

! Set the recently calculated function values in a column of DATMAT. This array has a column for
! each vertex of the current simplex, the entries of each column being the values of the constraint
! functions (if any) followed by the objective function and the greatest constraint violation at
! the vertex.
    datmat(:, jdrop) = con

    ibrnch = 1
    goto 140
end if

! Calculate DX = X(*) - X(0). Branch if the length of DX is less than 0.5*RHO.
call trstlp(n, m, A, con, rho, dx, ifull, iact)
shortd = (ifull == 0 .and. inprod(dx, dx) < QUART * rho * rho)

if (shortd) then
    ibrnch = 1
else
! Predict the change to F and to the maximum constraint violation if the variables are altered
! from X(0) to X(0)+DX.
    preref = inprod(dx, A(:, m + 1))
    prerec = datmat(m + 2, n + 1) - maxval([ZERO, con(1:m) - matprod(dx, A(:, 1:m))])

! Increase CPEN if necessary and branch back if this change alters the optimal vertex. Otherwise
! PREREM and PREREC will be set to the predicted reductions in the merit function and the maximum
! constraint violation respectively.
    barmu = zero
    if (prerec > ZERO) then
        barmu = -preref / prerec   ! What if PREREF <= 0 ??? Is it possible?
    end if
    if (cpen < 1.5E0_RP * barmu) then
        cpen = TWO * barmu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if (IPRINT >= 2) print 410, CPEN
!410 format(/3X, 'Increase in CPEN to', 1PE13.6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        phi = datmat(m + 1, :) + cpen * datmat(m + 2, :)
        phimin = minval(phi)
        if (phimin < phi(n + 1) .or. (cpen <= ZERO .and. any(datmat(m + 2, :) < datmat(m + 2, n + 1) .and. phi <= phimin))) then
            ! (CPEN <= ZERO) is indeed (CPEN == ZERO), and (PHI <= PHIMIN) is indeed (PHI == PHIMIN). We
            ! !write them in this way to avoid equality comparison of real numbers.
            goto 140
        end if
    end if

    prerem = preref + cpen * prerec


! Calculate the constraint and objective functions at X(*). Then find the actual reduction in the merit function.
    x = sim(:, n + 1) + dx
    ibrnch = 1

! Evaluate the objective function and constraints.
    if (any(is_nan(x))) then
        f = sum(x)  ! Set F to NaN.
        con = sum(x)  ! Set constraint values and constraint violation to NaN.
        info = -1
        goto 600
    end if

    call calcfc(n, m, x, f, con)
    nf = nf + 1
    cstrv = maxval([ZERO, -con(1:m)])
    con(m + 1) = f
    con(m + 2) = cstrv
    consav = con

!if (nf == IPRINT - 1 .or. IPRINT == 3) then
!    print 70, nf, F, CSTRV, (X(I), I=1, IPTEM)
!70  format(/3X, 'nf =', I5, 3X, 'F =', 1PE13.6, 4X, 'MAXCV =', 1PE13.6 / 3X, 'X =', 1PE13.6, 1P4E15.6)
!    if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!80  format(1PE19.6, 1P4E15.6)
!end if

! If the objective function value or the constraints contain a NaN or an infinite value, the exit.
    if (is_nan(F) .or. is_posinf(F)) then
        info = -2
        goto 600
    end if
    if (any(is_nan(con(1:m)))) then
        cstrv = sum(con(1:m))  ! Set CSTRV to NaN
        info = -2
        goto 600
    end if
! If the objective function achieves the target value at a feasible point, then exit.
    if (f <= ftarget .and. cstrv <= ctol) then
        info = 1
        return
    end if

    if (nf >= maxfun) then
!    if (IPRINT >= 1) print 50
!50  format(/3X, 'Return from subroutine COBYLA because the ', 'MAXFUN limit has been reached.')
        info = 3
    end if

! Set the recently calculated function values in a column of DATMAT. This array has a column for
! each vertex of the current simplex, the entries of each column being the values of the constraint
! functions (if any) followed by the objective function and the greatest constraint violation at
! the vertex.
    vmold = datmat(m + 1, n + 1) + cpen * datmat(m + 2, n + 1)
    vmnew = f + cpen * cstrv
    actrem = vmold - vmnew
    if (cpen <= ZERO .and. abs(f - datmat(m + 1, n + 1)) <= ZERO) then
        prerem = prerec
        actrem = datmat(m + 2, n + 1) - cstrv
    end if

! Begin the operations that decide whether X(*) should replace one of the vertices of the current
! simplex, the change being mandatory if ACTREM is positive. Firstly, JDROP is set to the index of
! the vertex that is to be replaced.
    if (actrem <= ZERO) then
        ratio = ONE
    else
        ratio = ZERO
    end if
    simid = matprod(simi, dx)
    sigbar = abs(simid) * vsig
    jdrop = 0
    if (maxval(abs(simid)) > ratio) then
        jdrop = int(maxloc(abs(simid), dim=1), kind(jdrop))
    end if

    edgmax = factor_delta * rho
    if (actrem > ZERO) then
        tmpv = sqrt(sum((spread(dx, dim=2, ncopies=n) - sim(:, 1:n))**2, dim=1))
    else
        tmpv = veta
    end if
    if (any(tmpv > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig))) then
        jdrop = int(maxloc(tmpv, mask=(sigbar >= parsig .or. sigbar >= vsig), dim=1), kind(jdrop))
    end if

! When jdrop=0, the algorithm decides not to include the trust-region trial point X into the
! simplex, because X is not good enough according to the merit function PHI = F + CPEN*CSTRV. In
! this case, X will simply be discarded in the original code. However, this decision depends on the
! value of CPEN. When CPEN is updated later, the discarded X might turn out better, sometimes even
! better than SIM(:, N+1), which is supposed to be the best point in the simplex. For this reason,
! we save the to-be-discarded X in XSAV and compare them with SIM(:, N+1) right before exiting. If
! a vector in XSAV turns out better than SIM(:, N+1), we replace SIM(:, N+1) by this vector. When
! jdrop > 0, SIM(:, jdrop) will be removed from the simplex according to PHI with the current CPEN.
! Similar to X, SIM(:, jdrop) may turn out better when CPEN is updated. Therefore, XSAV also takes
! SIM(:, jdrop) into account.
    if (jdrop == 0) then
        call savex(x, consav, xsav, datsav, nsav, ctol)   !?????
    else
        call savex(sim(:, n + 1) + sim(:, jdrop), datmat(:, jdrop), xsav, datsav, nsav, ctol)
        ! Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
        sim(:, jdrop) = dx
        simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), dx)
        simi = simi - outprod(matprod(simi, dx), simi_jdrop)
        simi(jdrop, :) = simi_jdrop
        datmat(:, jdrop) = con
    end if
end if

! Branch back for further iterations with the current RHO.
if ((shortd .or. actrem <= ZERO .or. actrem < TENTH * prerem) .and. improve_geo) then
    ibrnch = 0
end if

reduce_rho = (shortd .or. actrem <= ZERO .or. actrem < TENTH * prerem) .and. .not. improve_geo

if (reduce_rho) then
    if (rho <= rhoend) then
        terminate = .true.
    else
        ! Update RHO and CPEN.
        ! See equation (11) in Section 3 of the COBYLA paper for the update of RHO.
        rho = HALF * rho
        if (rho <= 1.5E0_RP * rhoend) then
            rho = rhoend
        end if
        ! See equation (12)--(13) in Section 3 of the COBYLA paper for the update of CPEN.
        if (cpen > ZERO) then
            cmin = minval(datmat(1:m, :), dim=2)
            cmax = maxval(datmat(1:m, :), dim=2)
            if (any(cmin < HALF * cmax)) then
                denom = minval(max(cmax, ZERO) - cmin, mask=(cmin < HALF * cmax))
                cpen = min(cpen, (maxval(datmat(m + 1, :)) - minval(datmat(m + 1, :))) / denom)
            else
                cpen = ZERO
            end if
        end if
        !if (IPRINT >= 2) print 580, RHO, CPEN
!580 format(/3X, 'Reduction in RHO to', 1PE13.6, '  and CPEN =', 1PE13.6)
!!    if (IPRINT == 2) then
!!        print 70, nf, DATMAT(M + 1, N + 1), DATMAT(m + 2, N + 1), (SIM(I, N + 1), I=1, IPTEM)
!!        if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!!    end if
    end if
end if

if (terminate) then
    info = 0
else
    goto 140
end if

! Return the best calculated values of the variables.
!if (iprint >= 1) print 590
!590 format(/3X, 'Normal return from subroutine COBYLA')
600 sim(:, 1:n) = sim(:, 1:n) + spread(sim(:, n + 1), dim=2, ncopies=n)  !!! TEMPORARY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make sure that the history includes the last X.
xhist = reshape([sim, xsav(:, 1:nsav), x], [n, n + nsav + 2])
fhist = [datmat(m + 1, :), datsav(m + 1, 1:nsav), f]
conhist = reshape([datmat(1:m, :), datsav(1:m, 1:nsav), consav], [m, n + nsav + 2])
cstrvhist = [datmat(m + 2, :), datsav(m + 2, 1:nsav), cstrv]
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cpen = max(cpen, 1.0E6_RP)
call selectx(conhist, cstrvhist, ctol, fhist, cpen, xhist, con(1:m), cstrv, f, x)
! We prefer SIM(:, N+1) unless the X selected above is even better.
if (.not. isbetter([f, cstrv], [datmat(m + 1, n + 1), datmat(m + 2, n + 1)], cpen, ctol)) then
    x = sim(:, n + 1)
    f = datmat(m + 1, n + 1)
    con = datmat(:, n + 1)
    cstrv = datmat(m + 2, n + 1)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!620 if (IPRINT >= 1) then
!    print 70, nf, F, CSTRV, (X(I), I=1, IPTEM)
!    if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!end if

!close (16)
return
end subroutine cobylb

subroutine savex(xdrop, datdrop, xsav, datsav, nsav, ctol)
! This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless a vector in XSAV(:, 1:NSAV) is
! better than XDROP. If XDROP is better than some vectors in XSAV(:, 1:NSAV), then these vectors
! will be removed. If XDROP is not better than any of XSAV(:, 1:NSAV) but NSAV=NSAVMAX, then we
! remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NSAV).
!
! When COBYLA calls this subroutine, XDROP is a vector to be "dropped", and  DATDROP contains its
! function/constraint information (constraint value in the first M entries, DATDROP(M+1) = F(XDROP),
! and DATDROP(M+2) = CSTRV(X)). XSAV and DATSAV save at most NSAVMAX vectors "dropped" by COBYLB
! and their function/constraint information. Only XSAV(:, 1:NSAV) and DATSAV(:, 1:NSAV) contains
! such vectors, while XSAV(:, NSAV+1:NSAVMAX) and DATSAV(:, NSAV+1:NSAVMAX) are not initialized yet.
!
! Note: We decide whether X is better than the function/constraint of Y according to the ISBETTER
! function with CPEN = -ONE. Due to the implementation of ISBETTER,
! X is better than Y with CPEN < 0
! ==> X is better than Y with any CPEN >= 0,
! ==> X is better than Y regardless of CPEN.

! Generic modules
use consts_mod, only : RP, IK, ONE
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : calquad, inprod

! Solver-specific modules
use selectx_mod, only : isbetter

implicit none

! Inputs
real(RP), intent(IN) :: ctol
real(RP), intent(IN) :: datdrop(:)  ! m+2
real(RP), intent(IN) :: xdrop(:)  ! n

! In-outputs
integer(IK), intent(INOUT) :: nsav
real(RP), intent(INOUT) :: datsav(:, :)  ! (M+2, NSAVMAX)
real(RP), intent(INOUT) :: xsav(:, :) ! (N, NSAVMAX)

! Local variables
integer(IK) :: m
integer(IK) :: n
integer(IK) :: nsavmax
integer(IK) :: i
real(RP) :: cpen
logical :: better(nsav)
logical :: keep(nsav)

m = size(datdrop) - 2
n = size(xdrop)
nsavmax = size(xsav, 2)

if (nsavmax <= 0) then
    return  ! Do nothing if NSAVMAX=0
end if

cpen = -ONE  ! See the comments above for why CPEN = -1

! Return immediately if any column of XSAV is better than XDROP.
! BETTER is defined by the array constructor with an implicit do loop.
better = [(isbetter([datsav(m + 1, i), datsav(m + 2, i)], [datdrop(m + 1), datdrop(m + 2)], cpen, ctol), i=1, nsav)]
if (any(better)) then
    return
end if

! Decide which columns of XSAV to keep. We use again the array constructor with an implicit do loop.
keep = [(.not. isbetter([datdrop(m + 1), datdrop(m + 2)], [datsav(m + 1, i), datsav(m + 2, i)], cpen, ctol), i=1, nsav)]
! If XDROP is not better than any column of XSAV, then we remove the first (oldest) column of XSAV.
if (count(keep) == nsavmax) then
    keep(1) = .false.
end if
xsav(:, 1:count(keep)) = xsav(:, pack([(i, i=1, nsav)], mask=keep))
datsav(:, 1:count(keep)) = datsav(:, pack([(i, i=1, nsav)], mask=keep))

! Update NSAV. Note that the update of XSAV and DATSAV used NSAV, so it should be updated afterward.
nsav = count(keep) + 1

! Save XDROP to XSAV(:, NSAV) and DATDROP to DATSAV(:, NSAV).
xsav(:, nsav) = xdrop(:)
datsav(:, nsav) = datdrop(:)

return

end subroutine savex


end module cobylb_mod
