!4.SIZE of d in TRSTEP is N + 1 instead of N, maybe also other arrays and other places.
!5.THE checks in rangehist cannot pass(xhist does not contain NaN)

module uobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of UOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Friday, May 06, 2022 PM12:53:04
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: uobyqb


contains


subroutine uobyqb(calfun, iprint, maxfun, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, &
    & x, nf, f, fhist, xhist, info)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: info_mod, only : NAN_INF_X, NAN_INF_F, NAN_MODEL, FTARGET_ACHIEVED, &
    & MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS!, MAXTR_REACHED
use, non_intrinsic :: linalg_mod, only : inprod, outprod!, norm
use, non_intrinsic :: symmat_mod, only : vec2smat, smat_mul_vec
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : quadinc, calvlag

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep
use, non_intrinsic :: trustregion_mod, only : trstep

! Development modules (to be removed)
use, non_intrinsic :: ieee_4dev_mod, only : ieeenan

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
!real(RP) :: d(size(x) + 1)  ! D(N+1) is accessed when N = 1
real(RP) :: d(size(x))
real(RP) :: g(size(x))
real(RP) :: h(size(x), size(x))
real(RP) :: pl((size(x) + 1) * (size(x) + 2) / 2, (size(x) + 1) * (size(x) + 2) / 2 - 1)
real(RP) :: pq(size(pl, 2))
real(RP) :: vlag(size(pl, 1))
real(RP) :: w(max(6_IK * size(x), (size(x)**2 + 3_IK * size(x) + 2_IK) / 2_IK))
real(RP) :: xbase(size(x))
real(RP) :: xnew(size(x))
real(RP) :: xopt(size(x))
real(RP) :: xpt(size(x), size(pl, 1))
!real(RP) :: xpt(size(x) + 1, size(pl, 1))  ! XPT(2, :) is accessed when N = 1
real(RP) :: ddknew, delta, diff, distsq(size(pl, 1)), weight(size(pl, 1)), score(size(pl, 1)),    &
&        distest, dnorm, errtol, estim, crvmin, fplus, fbase, fopt,&
&        fsave, ratio, rho, rhosq, sixthm, summ, &
&        sumg, temp, tempa, tol, tworsq, vmax,  &
&        qred, wmult, plknew((size(x) + 1) * (size(x) + 2) / 2 - 1)
integer(IK) :: ih, ip, iq, iw, j, k, knew, kopt, ksave
logical :: tr_success

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
knew = 1; kopt = 1
f = ieeenan()
!--------------------------------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
!       the corresponding arguments in SUBROUTINE UOBYQA.
!     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
!     XBASE will contain a shift of origin that reduces the contributions from
!       rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     PQ will contain the parameters of the quadratic model.
!     PL will contain the parameters of the Lagrange functions.
!     H will provide the second derivatives that TRSTEP and LAGMAX require.
!     G will provide the first derivatives that TRSTEP and LAGMAX require.
!     D is reserved for trial steps from XOPT, except that it will contain
!       diagonal second derivatives during the initialization procedure.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!     The array W will be used for working space. Its length must be at least
!     max [ 6*N, ( N**2 + 3*N + 2 ) / 2 ].
!
!     Set some constants.
!
tol = 0.01_RP
rho = rhobeg
rhosq = rho * rho


! Initialization. NF is the number of function calculations so far. The least function value so far,
! the corresponding X and its index are noted in FOPT, XOPT, and KOPT respectively.
xbase = x
xpt = ZERO
pl = ZERO
j = 0
ih = n
! In the following loop, FPLUS is set to F(X + RHO*e_I) when NF = 2*I, and the value of FPLUS is
! used subsequently when NF = 2*I + 1.
fplus = ZERO ! This initial value is not used but to entertain the Fortran compilers.
do nf = 1, 2_IK * n + 1_IK
    ! Pick the shift from XBASE to the next initial interpolation point that provides diagonal
    ! second derivatives.
    if (nf > 1) then
        if (modulo(nf, 2_IK) == 1_IK) then
            if (fplus < fbase) then
                w(j) = rho
                xpt(j, nf) = TWO * rho
            else
                w(j) = -rho
                xpt(j, nf) = -rho
            end if
        elseif (j < n) then
            j = j + 1
            xpt(j, nf) = rho
        end if
    end if
    x = xbase + xpt(:, nf)

    if (is_nan(sum(abs(x)))) then
        f = sum(x) ! Set F to NaN
        if (nf == 1) then
            fopt = f
            xopt = ZERO
        end if
        info = NAN_INF_X
        goto 420
    end if
    call evaluate(calfun, x, f)
    call savehist(nf, x, xhist, f, fhist)
    if (is_nan(f) .or. is_posinf(f)) then
        if (nf == 1) then
            fopt = f
            xopt = ZERO
        end if
        info = NAN_INF_F
        goto 420
    end if

    if (f <= ftarget) then
        info = FTARGET_ACHIEVED
        goto 430  ! Should not goto 420. fopt may not be defined yet
    end if

    if (nf == 1) then
        fopt = f
        kopt = nf
        xopt = xpt(:, nf)
        fbase = f
    elseif (f < fopt) then
        fopt = f
        kopt = nf
        xopt = xpt(:, nf)
    end if

    if (nf >= maxfun) then
        info = MAXFUN_REACHED
        goto 420
    end if

    if (modulo(nf, 2_IK) == 0_IK) then
        fplus = f  ! FPLUS = F(X + RHO*e_I) with I = NF/2.
        cycle
    end if

    ! Form the gradient and diagonal second derivatives of the quadratic model and Lagrange functions.
    if (j >= 1 .and. nf >= 3) then  ! NF >= 3 is implied by J >= 1. We prefer to impose it explicitly.
        ih = ih + j
        if (xpt(j, nf) > 0) then  ! XPT(J, NF) = 2*RHO
            pq(j) = (4.0_RP * fplus - 3.0_RP * fbase - f) / (TWO * rho)
            d(j) = (fbase + f - TWO * fplus) / rhosq
            pl(1, j) = -1.5_RP / rho
            pl(1, ih) = ONE / rhosq
            pl(nf - 1, j) = TWO / rho  ! Should be moved out of the loop
            pl(nf - 1, ih) = -TWO / rhosq  ! Should be moved out of the loop
        else  ! XPT(J, NF) = -RHO
            d(j) = (fplus + f - TWO * fbase) / rhosq
            pq(j) = (fplus - f) / (TWO * rho)
            pl(1, ih) = -TWO / rhosq
            pl(nf - 1, j) = HALF / rho  ! Should be moved out of the loop
            pl(nf - 1, ih) = ONE / rhosq  ! Should be moved out of the loop
        end if
        pq(ih) = d(j)
        pl(nf, j) = -HALF / rho
        pl(nf, ih) = ONE / rhosq
    end if
end do

ih = n + 1
ip = 0
iq = 2

! Form the off-diagonal second derivatives of the initial quadratic model.
do nf = 2_IK * n + 2_IK, npt
    ! Pick the shift from XBASE to the next initial interpolation point that provides off-diagonal
    ! second derivatives.
    ip = ip + 1
    if (ip == iq) then
        ih = ih + 1
        ip = 1
        iq = iq + 1
    end if
    xpt(ip, nf) = w(ip)
    ! N.B.: XPT(2, NF+1) is accessed by XPT(IQ, NF+1) even if N = 1.
    xpt(iq, nf) = w(iq)
    x = xbase + xpt(:, nf)

    if (is_nan(sum(abs(x)))) then
        f = sum(x) ! Set F to NaN
        info = NAN_INF_X
        goto 420
    end if

    call evaluate(calfun, x, f)
    call savehist(nf, x, xhist, f, fhist)
    if (is_nan(f) .or. is_posinf(f)) then
        info = NAN_INF_F
        goto 420
    end if

    if (f <= ftarget) then
        info = FTARGET_ACHIEVED
        goto 430  ! Should not goto 420. fopt may not be defined yet
    end if

    if (f < fopt) then
        fopt = f
        kopt = nf
        xopt = xpt(:, nf)
    end if
    if (nf >= maxfun) then
        info = MAXFUN_REACHED
        goto 420
    end if

    ih = ih + 1
    temp = ONE / (w(ip) * w(iq))
    tempa = f - fbase - w(ip) * pq(ip) - w(iq) * pq(iq)
    ! N.B.: D(2) is accessed by D(IQ) even if N = 1.
    pq(ih) = (tempa - HALF * rhosq * (d(ip) + d(iq))) * temp
    pl(1, ih) = temp
    iw = ip + ip
    if (w(ip) < ZERO) iw = iw + 1
    pl(iw, ih) = -temp
    iw = iq + iq
    if (w(iq) < ZERO) iw = iw + 1
    pl(iw, ih) = -temp
    pl(nf, ih) = temp
end do
!--------------------------------------------------------------------------------------------------!
! When the loop exits, the value of NF is not specified by the standard. With gfortran, it will be
! NPT+1, which is not proper for the subsequent use. !!!
nf = min(nf, npt) !!!
!--------------------------------------------------------------------------------------------------!


!
!     Set parameters to begin the iterations for the current RHO.
!
sixthm = ZERO
delta = rho
60 continue

tworsq = (TWO * rho)**2
rhosq = rho * rho
!
!     Form the gradient of the quadratic model at the trust region centre.
!
70 continue

knew = 0
xopt = xpt(:, kopt)
g = pq(1:n) + smat_mul_vec(pq(n + 1:npt - 1), xopt)
h = vec2smat(pq(n + 1:npt - 1))

if (is_nan(sum(abs(g)) + sum(abs(h)))) then
    info = NAN_MODEL
    goto 420
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
call trstep(delta, g, h, tol, d, crvmin)

dnorm = min(delta, sqrt(sum(d**2)))
errtol = -ONE
if (dnorm < HALF * rho) then
    knew = -1
    errtol = HALF * crvmin * rho * rho
    if (nf <= npt + 9) errtol = ZERO
    goto 290
end if
!
!     Calculate the next value of the objective function.


100 continue
xnew = xopt + d
x = xbase + xnew

if (nf >= maxfun) then
    info = MAXFUN_REACHED
    goto 420
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (is_nan(sum(abs(x)))) then
    f = sum(x) ! Set F to NaN
    if (nf == 1) then
        fopt = f
        xopt = ZERO
    end if
    info = NAN_INF_X
    goto 420
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------!
call evaluate(calfun, x, f)
nf = nf + 1
call savehist(nf, x, xhist, f, fhist)
!------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Zaikun (commented on 02-06-2019; implemented in 2016):
!     Exit if F has an NaN or almost infinite value.
!     If this happens at the very first function evaluation (i.e.,
!     NF=1), then it is necessary to set FOPT and XOPT before going to
!     530, because these two variables have not been set yet.
if (is_nan(f) .or. is_posinf(f)) then
    if (nf == 1) then
        fopt = f
        xopt = ZERO
    end if
    info = NAN_INF_F
    goto 420
end if

if (f <= ftarget) then
    info = FTARGET_ACHIEVED
    goto 430  ! Should not goto 420. fopt may not be defined yet
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF (KNEW .EQ. -1) GOTO 420
if (knew == -1) then
    info = SMALL_TR_RADIUS !!??
    goto 420
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use the quadratic model to predict the change in F due to the step D, and find the values of the
! Lagrange functions at the new point.
qred = -quadinc(pq, d, xopt)
vlag = calvlag(pl, d, xopt, kopt)

! Update SIXTHM, which is a lower bound on one sixth of the greatest third derivative of F. The
! lower bound is derived from (3.1) of the UOBYQA paper.
diff = f - fopt + qred
summ = ZERO
distsq = sum((xpt - spread(xnew, dim=2, ncopies=npt))**2, dim=1)
summ = inprod(distsq, sqrt(distsq) * abs(vlag))
! SUMM may become 0 due to rounding, even in double precision. Detected by ifort.
if (abs(diff) > 0 .and. summ > 0) then
    sixthm = max(sixthm, abs(diff) / summ)
end if
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. Then branch if D is not a trust region step.
!
fsave = fopt
if (f < fopt) then
    fopt = f
    xopt = xnew
end if
ksave = knew

!----------------------------------------------------------!
ddknew = ZERO ! Necessary, or DDKNEW is not always defined.
!----------------------------------------------------------!
!if (knew > 0) goto 240
if (knew <= 0) then

!
!     Pick the next value of DELTA after a trust region step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.not. (qred > 0)) then
        info = TRSUBP_FAILED
        goto 420
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
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
    distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
    weight = max(ONE, distsq / rhosq)**1.5_RP
    score = weight * abs(vlag)

    tr_success = (f < fsave)

    if (.not. tr_success) then
        ! If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
        score(kopt) = -ONE
    end if

    if (any(score > 1) .or. (tr_success .and. any(score > 0))) then
        ! SCORE(K) is NaN implies VLAG(K) is NaN, but we want ABS(VLAG) to be big. So we exclude such K.
        knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
        !!MATLAB: [~, knew] = max(score, [], 'omitnan');
        ddknew = distsq(knew)
    elseif (tr_success) then
        ! Powell's code does not include the following instructions. With Powell's code, if DENABS
        ! consists of only NaN, then KNEW can be 0 even when TR_SUCCESS is TRUE.
        knew = int(maxloc(distsq, dim=1), IK)
    else
        knew = 0_IK
    end if
    if (knew == 0) goto 290
end if
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!

xpt(:, knew) = xnew
! It can happen that VLAG(KNEW) = 0 due to rounding.
pl(knew, :) = pl(knew, :) / vlag(knew)
plknew = pl(knew, :)
pq = pq + diff * plknew
pl = pl - outprod(vlag, plknew)
pl(knew, :) = plknew

!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
if (f < fsave) then
    kopt = knew
end if

if (f < fsave .or. ksave > 0 .or. dnorm > TWO * rho .or. ddknew > tworsq) goto 70

!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
290 continue

distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
w(1:npt) = distsq

310 continue

knew = -1
distest = tworsq
if (any(w(1:npt) > distest)) then
    knew = int(maxloc(w(1:npt), mask=(.not. is_nan(w(1:npt))), dim=1), IK)
    distest = w(knew)
    !!MATLAB: [distest, knew] = max(w(1:npt), [], 'omitnan');
end if
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     HALF the summ of squares of components of the Hessian.
!
if (knew > 0) then
    !g = pl(knew, 1:n)
    g = pl(knew, 1:n) + smat_mul_vec(pl(knew, n + 1:npt - 1), xopt)
    h = vec2smat(pl(knew, n + 1:npt - 1))
    if (is_nan(sum(abs(g)) + sum(abs(h)))) then
        info = NAN_MODEL
        goto 420
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
    wmult = HUGENUM  ! Needed for Fortran, or compiler may complain that it is not defined.
    if (errtol > ZERO) then
        w(knew) = ZERO
        sumg = ZERO
        sumg = sum(g**2)
        estim = rho * (sqrt(sumg) + rho * HALF * sqrt(sum(h**2)))
        wmult = sixthm * distest**1.5_RP
        if (wmult * estim <= errtol) goto 310   ! Infinite cycling possible?
    end if
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
    call geostep(g, h, rho, d, vmax)
    !if (errtol > ZERO) then
    !    if (wmult * vmax <= errtol) goto 310
    !end if
    !goto 100
    if (errtol > 0 .and. wmult * vmax <= errtol) then
        goto 310  ! Infinite cycling possible?
        !else
    elseif (vmax > 0) then
        goto 100  ! GEOSTEP succeeds.
    else
        goto 600  ! GEOSTEP fails.
    end if
end if
if (dnorm > rho) goto 70
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
600 continue
if (rho > rhoend) then

    xbase = xbase + xopt
    xpt = xpt - spread(xopt, dim=2, ncopies=npt)
    pq(1:n) = pq(1:n) + smat_mul_vec(pq(n + 1:npt - 1), xopt)
    do k = 1, npt
        pl(k, 1:n) = pl(k, 1:n) + smat_mul_vec(pl(k, n + 1:npt - 1), xopt)
    end do

!
!     Pick the next values of RHO and DELTA.
!
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
    goto 60
end if

info = SMALL_TR_RADIUS !!??

!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.

!if (errtol >= ZERO) goto 100
if (errtol >= ZERO) then
    !----------------!
    knew = -1  ! Zaikun 20220506: Tell the algorithm to come to 420 immediately after the function evaluation.
    !----------------!
    goto 100
end if

420 continue

if (fopt <= f .or. is_nan(f)) then
    x = xbase + xopt
    f = fopt
end if

430 continue

!close (16)

call rangehist(nf, xhist, fhist)

end subroutine uobyqb


end module uobyqb_mod
