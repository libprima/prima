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
! Last Modified: Wednesday, May 04, 2022 PM07:48:20
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: uobyqb


contains


subroutine uobyqb(calfun, n, x, rhobeg, rhoend, iprint, maxfun, npt, xbase, &
     &  xopt, xnew, pq, pl, h, g, f, info, ftarget, nf, xhist, maxxhist, fhist, maxfhist)
!     &  xopt, xnew, xpt, pq, pl, h, g, d, vlag, w, f, info, ftarget, nf, xhist, maxxhist, fhist, maxfhist)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: pintrf_mod, only : OBJ

! Solver-specific modules
use, non_intrinsic :: geometry_mod, only : geostep
use, non_intrinsic :: trustregion_mod, only : trstep

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK) :: n
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfhist
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: maxxhist
integer(IK), intent(in) :: npt
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(inout) :: x(n)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
!real(RP), intent(out) :: d(n)
real(RP) :: d(n + 1)  ! D(2) may be accessed even if N = 1
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(maxfhist)
real(RP), intent(out) :: g(n)
real(RP), intent(out) :: h(n, n)
real(RP), intent(out) :: pl(npt, npt - 1_IK)
real(RP), intent(out) :: pq(npt - 1_IK)
!real(RP), intent(out) :: vlag(npt)
real(RP) :: vlag(npt)
real(RP), intent(out) :: xbase(n)
real(RP), intent(out) :: xhist(n, maxxhist)
real(RP), intent(out) :: xnew(n)
real(RP), intent(out) :: xopt(n)
!real(RP), intent(out) :: xpt(n, npt)
real(RP) :: xpt(n + 1, npt)  ! XPT(2, :) may be accessed even if N = 1
!real(RP), intent(out) :: w(max(6_IK * n, (n**2 + 3_IK * n + 2_IK) / 2_IK))
!real(RP) :: w(max(6_IK * n, (n**2 + 3_IK * n + 2_IK) / 2_IK))
! In Powell's TRSTEP code, n-dimensional vectors may be accessed at N+1, leading to memory errors.
! We assign one more place to the vectors to avoid such problems. Thus W is enlarged for this.
real(RP) :: w(max(6_IK * n, (n**2 + 3_IK * n + 2_IK) / 2_IK) + 100_IK * (n + 1_IK))

! Local variables
real(RP) :: ddknew, delta, detrat, diff,        &
&        distest, dnorm, errtol, estim, evalue, fbase, fopt,&
&        fsave, ratio, rho, rhosq, sixthm, summ, &
&        sumg, sumh, temp, tempa, tol, tworsq, vmax,  &
&        vquad, wmult, distsq(size(pl, 1)), gtmp(n), pltmp(size(pl, 1), n)
integer(IK) :: i, ih, ip, iq, iw, j, jswitch, k, knew, kopt,&
&           ksave, ktemp, nftest, nnp, nptm
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
nnp = n + n + 1
nptm = npt - 1
nftest = maxfun
!
!     Initialization. NF is the number of function calculations so far.
!
rho = rhobeg
rhosq = rho * rho
nf = 0
do i = 1, n
    xbase(i) = x(i)
    do k = 1, npt
        xpt(i, k) = ZERO
    end do
end do
do k = 1, npt
    do j = 1, nptm
        pl(k, j) = ZERO
    end do
end do
!
!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.
!
30 do i = 1, n
    x(i) = xbase(i) + xpt(i, nf + 1)
end do
goto 120
50 if (nf == 1) then
    fopt = f
    kopt = nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n
        xopt(i) = xpt(i, 1)
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fbase = f
    j = 0
    jswitch = -1
    ih = n
else
    if (f < fopt) then
        fopt = f
        kopt = nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = 1, n
            xopt(i) = xpt(i, nf)
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if
end if
!
!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.
!
if (nf <= nnp) then
    jswitch = -jswitch
    if (jswitch > 0) then
        if (j >= 1) then
            ih = ih + j
            if (w(j) < ZERO) then
                d(j) = (fsave + f - TWO * fbase) / rhosq
                pq(j) = (fsave - f) / (TWO * rho)
                pl(1, ih) = -TWO / rhosq
                pl(nf - 1, j) = HALF / rho
                pl(nf - 1, ih) = ONE / rhosq
            else
                pq(j) = (4.0_RP * fsave - 3.0_RP * fbase - f) / (TWO * rho)
                d(j) = (fbase + f - TWO * fsave) / rhosq
                pl(1, j) = -1.5_RP / rho
                pl(1, ih) = ONE / rhosq
                pl(nf - 1, j) = TWO / rho
                pl(nf - 1, ih) = -TWO / rhosq
            end if
            pq(ih) = d(j)
            pl(nf, j) = -HALF / rho
            pl(nf, ih) = ONE / rhosq
        end if
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
!
        if (j < n) then
            j = j + 1
            xpt(j, nf + 1) = rho
        end if
    else
        fsave = f
        if (f < fbase) then
            w(j) = rho
            xpt(j, nf + 1) = TWO * rho
        else
            w(j) = -rho
            xpt(j, nf + 1) = -rho
        end if
    end if
    if (nf < nnp) goto 30
!
!     Form the off-diagonal second derivatives of the initial quadratic model.
!
    ih = n
    ip = 1
    iq = 2
end if
ih = ih + 1
if (nf > nnp) then
    temp = ONE / (w(ip) * w(iq))
    tempa = f - fbase - w(ip) * pq(ip) - w(iq) * pq(iq)
    pq(ih) = (tempa - HALF * rhosq * (d(ip) + d(iq))) * temp
    pl(1, ih) = temp
    iw = ip + ip
    if (w(ip) < ZERO) iw = iw + 1
    pl(iw, ih) = -temp
    iw = iq + iq
    if (w(iq) < ZERO) iw = iw + 1
    pl(iw, ih) = -temp
    pl(nf, ih) = temp
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
!
    ip = ip + 1
end if
if (ip == iq) then
    ih = ih + 1
    ip = 1
    iq = iq + 1
end if
if (nf < npt) then
    xpt(ip, nf + 1) = w(ip)
    xpt(iq, nf + 1) = w(iq)
    goto 30
end if
!
!     Set parameters to begin the iterations for the current RHO.
!
sixthm = ZERO
delta = rho
60 tworsq = (TWO * rho)**2
rhosq = rho * rho
!
!     Form the gradient of the quadratic model at the trust region centre.
!
70 knew = 0
ih = n
do j = 1, n
    xopt(j) = xpt(j, kopt)
    !g(j) = pq(j)
    g(j) = ZERO
    do i = 1, j
        ih = ih + 1
        !g(i) = g(i) + pq(ih) * xopt(j)
        !if (i < j) g(j) = g(j) + pq(ih) * xopt(i)
        if (i < j) g(j) = g(j) + pq(ih) * xopt(i)
        g(i) = g(i) + pq(ih) * xopt(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: For ill-conditiONEd problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTEM or LAGMAX is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended
! behavior of the code, including uninitialized indices.
!   80 H(I,J)=PQ(IH)
        h(i, j) = pq(ih)
        if (h(i, j) /= h(i, j)) then
            info = -3
            goto 420
        end if
    end do
end do
!----------------!
g = g + pq(1:n)
!----------------!
do i = 1, n
    if (g(i) /= g(i)) then
        info = -3
        goto 420
    end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
!
!call trstep(n, g, h, delta, tol, d, w(1), w(n + 1), w(2 * n + 1), w(3 * n + 1), &
!& w(4 * n + 1), w(5 * n + 1), evalue)
! In Powell's TRSTEP code, n-dimensional vectors may be accessed at N+1, leading to memory errors.
! We assign one more place to the vectors to avoid such problems.
call trstep(n, g, h, delta, tol, d(1:n), w(1), w(n + 1), w(2 * n + 2), w(3 * n + 3), &
& w(4 * n + 4), w(5 * n + 5), evalue)
!trstep(n, g, h, delta, tol, d, gg, td, tn, w, piv, z, evalue)
temp = ZERO
do i = 1, n
    temp = temp + d(i)**2
end do
dnorm = min(delta, sqrt(temp))
errtol = -ONE
if (dnorm < HALF * rho) then
    knew = -1
    errtol = HALF * evalue * rho * rho
    if (nf <= npt + 9) errtol = ZERO
    goto 290
end if
!
!     Calculate the next value of the objective function.
!
100 do i = 1, n
    xnew(i) = xopt(i) + d(i)
    x(i) = xbase(i) + xnew(i)
end do
120 if (nf >= nftest) then
    info = 3
    goto 420
end if
!nf = nf + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1, n
    if (x(i) /= x(i)) then
        f = x(i) ! set f to nan
        if (nf == 1) then
            fopt = f
            do j = 1, n
                xopt(j) = ZERO
            end do
        end if
        info = -1
        goto 420
    end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------!
call evaluate(calfun, x, f)
nf = nf + 1
call savehist(nf, x, xhist, f, fhist)
!------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Zaikun (commented on 02-06-2019; implemented in 2016):
!     Exit if F has an NaN or almost infinite value.
!     If this happends at the very first function evaluation (i.e.,
!     NF=1), then it is necessary to set FOPT and XOPT before going to
!     530, because these TWO variables have not been set yet.
if (is_nan(f) .or. is_posinf(f)) then
    if (nf == 1) then
        fopt = f
        do i = 1, n
            xopt(i) = ZERO
        end do
    end if
    info = -2
    goto 420
end if
!     By Zaikun (commented on 02-06-2019; implemented in 2016):
!     Exit if F .LE. FTARGET.
if (f <= ftarget) then
    info = 1
    goto 430  ! Should not goto 420. fopt may not be defined yet
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nf <= npt) goto 50
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF (KNEW .EQ. -1) GOTO 420
if (knew == -1) then
    info = 0
    goto 420
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.
!
vquad = ZERO
ih = n
do j = 1, n
    w(j) = d(j)
    vquad = vquad + w(j) * pq(j)
    do i = 1, j
        ih = ih + 1
        w(ih) = d(i) * xnew(j) + d(j) * xopt(i)
        if (i == j) w(ih) = HALF * w(ih)
        vquad = vquad + w(ih) * pq(ih)
    end do
end do
do k = 1, npt
    temp = ZERO
    do j = 1, nptm
        temp = temp + w(j) * pl(k, j)
    end do
    vlag(k) = temp
end do
vlag(kopt) = vlag(kopt) + ONE
!
!     Update SIXTHM, which is a lower bound on ONE sixth of the greatest
!     third derivative of F.
!
diff = f - fopt - vquad
summ = ZERO
do k = 1, npt
    temp = ZERO
    do i = 1, n
        temp = temp + (xpt(i, k) - xnew(i))**2
    end do
    !temp = sqrt(temp)
    !summ = summ + abs(temp * temp * temp * vlag(k))
    summ = summ + abs(sqrt(temp) * vlag(k) * temp)
end do
!sixthm = max(sixthm, abs(diff) / summ)
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
    do i = 1, n
        xopt(i) = xnew(i)
    end do
end if
ksave = knew
if (knew > 0) goto 240
!
!     Pick the next value of DELTA after a trust region step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF (VQUAD .GE. ZERO) THEN
if (.not. (vquad < ZERO)) then
    info = 2
    goto 420
end if
ratio = (f - fsave) / vquad
if (ratio <= 0.1_RP) then
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
ktemp = 0
detrat = ZERO
if (f >= fsave) then
    ktemp = kopt
    detrat = ONE
end if
do k = 1, npt
    summ = ZERO
    do i = 1, n
        summ = summ + (xpt(i, k) - xopt(i))**2
    end do
    temp = abs(vlag(k))
    !if (summ > rhosq) temp = temp * (summ / rhosq)**1.5_RP
    temp = temp * max(ONE, summ / rhosq)**1.5_RP
    if (temp > detrat .and. k /= ktemp) then
        detrat = temp
        ddknew = summ
        knew = k
    end if
end do

!-----------------------------------------------------------------------------------------------!
distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
if (knew == 0 .and. f < fsave) then
    knew = int(maxloc(distsq, dim=1), IK)
end if
!-----------------------------------------------------------------------------------------------!

if (knew == 0) goto 290
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!
240 do i = 1, n
    xpt(i, knew) = xnew(i)
end do
temp = ONE / vlag(knew)
do j = 1, nptm
    !pl(knew, j) = temp * pl(knew, j)
    pl(knew, j) = pl(knew, j) / vlag(knew)

    pq(j) = pq(j) + diff * pl(knew, j)
end do
do k = 1, npt
    if (k /= knew) then
        temp = vlag(k)
        do j = 1, nptm
            pl(k, j) = pl(k, j) - temp * pl(knew, j)
        end do
    end if
end do
!
!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
if (f < fsave) then
    kopt = knew
    goto 70
end if
if (ksave > 0) goto 70
if (dnorm > TWO * rho) goto 70
if (ddknew > tworsq) goto 70
!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
290 do k = 1, npt
    w(k) = ZERO
    do i = 1, n
        w(k) = w(k) + (xpt(i, k) - xopt(i))**2
    end do
end do
310 knew = -1
distest = tworsq
do k = 1, npt
    if (w(k) > distest) then
        knew = k
        distest = w(k)
    end if
end do
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     HALF the summ of squares of compONEnts of the Hessian.
!
if (knew > 0) then
    ih = n
    sumh = ZERO
    do j = 1, n
        !g(j) = pl(knew, j)
        g(j) = ZERO
        do i = 1, j
            ih = ih + 1
            temp = pl(knew, ih)
            g(j) = g(j) + temp * xopt(i)
            if (i < j) then
                g(i) = g(i) + temp * xopt(j)
                !sumh = sumh + temp * temp
                sumh = sumh + TWO * temp**2
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: See the comments below line number 70
!  330     H(I,J)=TEMP
            h(i, j) = temp
            if (h(i, j) /= h(i, j)) then
                info = -3
                goto 420
            end if
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !sumh = sumh + HALF * temp * temp
        sumh = sumh + temp**2
    end do
    g = pl(knew, 1:n) + g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: See the comments below line number 70
    do i = 1, n
        if (g(i) /= g(i)) then
            info = -3
            goto 420
        end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
    if (errtol > ZERO) then
        w(knew) = ZERO
        sumg = ZERO
        do i = 1, n
            sumg = sumg + g(i)**2
        end do
        !estim = rho * (sqrt(sumg) + rho * sqrt(HALF * sumh))
        estim = rho * (sqrt(sumg) + rho * HALF * sqrt(sumh))
        wmult = sixthm * distest**1.5_RP
        if (wmult * estim <= errtol) goto 310
    end if
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
    call geostep(n, g, h, rho, d(1:n), xnew, vmax)
    if (errtol > ZERO) then
        if (wmult * vmax <= errtol) goto 310
    end if
    goto 100
end if
if (dnorm > rho) goto 70
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
if (rho > rhoend) then
    ih = n
    do j = 1, n
        xbase(j) = xbase(j) + xopt(j)
        do k = 1, npt
            xpt(j, k) = xpt(j, k) - xopt(j)
        end do
        !---------------------!
        gtmp(j) = pq(j)
        pq(j) = ZERO
        pltmp(:, j) = pl(:, j)
        pl(:, j) = ZERO
        !---------------------!
        do i = 1, j
            ih = ih + 1
            pq(i) = pq(i) + pq(ih) * xopt(j)
            if (i < j) then
                pq(j) = pq(j) + pq(ih) * xopt(i)
                do k = 1, npt
                    pl(k, j) = pl(k, j) + pl(k, ih) * xopt(i)
                end do
            end if

            !---------------------------------------------!
            do k = 1, npt
                pl(k, i) = pl(k, i) + pl(k, ih) * xopt(j)
            end do
            !---------------------------------------------!
        end do
    end do
    !---------------------!
    pq(1:n) = pq(1:n) + gtmp
    pl(:, 1:n) = pl(:, 1:n) + pltmp
    !---------------------!
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
        rho = 0.1_RP * rho
    end if
    delta = max(delta, rho)
    goto 60
end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (errtol >= ZERO) goto 100
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  420 IF (FOPT .LE. F) THEN
420 if (fopt <= f .or. f /= f) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n
        x(i) = xbase(i) + xopt(i)
    end do
    f = fopt
end if

!---------------------------------------------!
430 call rangehist(nf, xhist, fhist)
!---------------------------------------------!

end subroutine uobyqb


end module uobyqb_mod
