module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of UOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, May 09, 2022 PM01:30:39
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstep


contains


subroutine trstep(n, g, h, delta, tol, d_out, gg, td, tn, w, piv, z, evalue)

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : maximum

implicit none

! Inputs
integer(IK), intent(in) :: n
real(RP), intent(in) :: delta
real(RP), intent(in) :: g(n)
real(RP), intent(in) :: tol

! In-outputs
real(RP), intent(inout) :: d_out(n)
real(RP) :: d(n + 1)  !!! D(N+1) may be accessed
real(RP), intent(inout) :: gg(n)
real(RP), intent(inout) :: h(n, n)
!real(RP), intent(inout) :: piv(n)
!real(RP), intent(inout) :: td(n)
!real(RP), intent(inout) :: tn(n)
real(RP), intent(inout) :: piv(n + 1)  !!! PIV(N+1) may be accessed
real(RP), intent(inout) :: td(n + 1)  !!! TD(N+1) may be accessed
real(RP), intent(inout) :: tn(n + 1)  !!! TN(N+1) may be accessed
real(RP), intent(inout) :: w(n)
real(RP), intent(inout) :: z(n)

! Outputs
real(RP), intent(out) :: evalue

! Local variables
real(RP) :: delsq, dhd, dnorm, dsq, dtg, dtz, gam, gnorm,     &
&        gsq, hnorm, par, parl, parlest, paru,         &
&        paruest, phi, phil, phiu, pivksv, pivot, posdef,   &
&        scaling, shfmax, shfmin, shift, slope, summ, sumd,    &
&        tdmin, temp, tempa, tempb, wsq, wwsq, wz, zsq
real(RP) :: dsav(n)
integer(IK) :: i, iterc, j, jp, k, kp, kpp, ksav, ksave, nm
logical :: scaled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q. Only the upper triangular and diagonal
!       parts need be set. The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.
!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
!     EVALUE will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step. Then EVALUE will be positive, but otherwise it
!     will be set to ZERO.
!
!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| .LEQ. DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties
!
!             ||D|| .LEQ. DELTA  and  ACTRED .LEQ. (1-TOL)*MAXRED.
!
!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.
!
!     Initialization.
!
!

!---------------------------------------------!
! Zaikun 20220509
d_out = ZERO
evalue = ZERO

gsq = sum(g**2)
gnorm = sqrt(gsq)

if (.not. any(abs(h) > 0)) then
    if (gnorm > 0) then
        d_out = -delta * (g / gnorm)
    end if
    return
end if
!---------------------------------------------!

delsq = delta * delta
evalue = ZERO
nm = n - 1
do i = 1, n
    d(i) = ZERO
    do j = 1, i
        h(i, j) = h(j, i)
    end do
end do
!
!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H. Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.
!

!-------------------------------------------------------------------------------------------------!
! Zaikun 20220508
scaling = maxval(abs(h))
scaled = .false.
if (scaling <= 0) then
    td(1:n) = ZERO
    tn(1:n - 1) = ZERO
elseif (scaling > 1.0E6 .or. scaling < 1.0E-6) then  ! 1.0E6 and 1.0E-6 are heuristic.
    h = h / scaling
    scaled = .true.
end if
!-------------------------------------------------------------------------------------------------!

do i = 1, n
    td(i) = h(i, i)
end do

if (.not. scaling <= 0) then
    do k = 1, nm
        kp = k + 1
        summ = ZERO
        if (kp < n) then
            kpp = kp + 1
            do i = kpp, n
                summ = summ + h(i, k)**2
            end do
        end if
        if (summ == ZERO) then
            tn(k) = h(kp, k)
            h(kp, k) = ZERO
        else
            temp = h(kp, k)
            !tn(k) = sign(sqrt(summ + temp * temp), temp)
            tn(k) = sign(sqrt(summ + temp**2), temp)
            h(kp, k) = -summ / (temp + tn(k))
            temp = sqrt(TWO / (summ + h(kp, k)**2))
            do i = kp, n
                w(i) = temp * h(i, k)
                h(i, k) = w(i)
                z(i) = td(i) * w(i)
            end do
            wz = ZERO
            do j = kp, nm
                jp = j + 1
                do i = jp, n
                    z(i) = z(i) + h(i, j) * w(j)
                    z(j) = z(j) + h(i, j) * w(i)
                end do
                wz = wz + w(j) * z(j)
            end do
            wz = wz + w(n) * z(n)
            do j = kp, n
                td(j) = td(j) + w(j) * (wz * w(j) - TWO * z(j))
                if (j < n) then
                    jp = j + 1
                    do i = jp, n
                        h(i, j) = h(i, j) - w(i) * z(j) - w(j) * (z(i) - wz * w(i))
                    end do
                end if
            end do
        end if
    end do
end if


!-------------------------------------------------------------------------------------------------!
! Zaikun 20220508
if (scaled) then
    td(1:n) = td(1:n) * scaling
    tn(1:n - 1) = tn(1:n - 1) * scaling
end if
!-------------------------------------------------------------------------------------------------!



!
!     Form GG by applying the similarity transformation to G.
!
gsq = ZERO
do i = 1, n
    gg(i) = g(i)
    gsq = gsq + g(i)**2
end do
gnorm = sqrt(gsq)
do k = 1, nm
    kp = k + 1
    summ = ZERO
    do i = kp, n
        summ = summ + gg(i) * h(i, k)
    end do
    do i = kp, n
        gg(i) = gg(i) - summ * h(i, k)
    end do
end do

!---------------------------------------------------------------------------------------!
! Zaikun 20220506
tn(n) = ZERO
if (.not. is_finite(sum(abs(h)) + sum(abs(td(1:n))) + sum(abs(tn(1:n - 1))) + sum(abs(gg)))) then
    d_out = ZERO
    evalue = ZERO
    return
end if
!---------------------------------------------------------------------------------------!
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is ZERO.
!
!hnorm = abs(td(1)) + abs(tn(1))
hnorm = abs(ZERO) + abs(td(1)) + abs(tn(1))
tdmin = td(1)
tn(n) = ZERO
do i = 2, n
    temp = abs(tn(i - 1)) + abs(td(i)) + abs(tn(i))
    hnorm = max(hnorm, temp)
    tdmin = min(tdmin, td(i))
end do

!if (hnorm == ZERO) then
!    if (gnorm == ZERO) goto 400
!    scaling = delta / gnorm
!    do i = 1, n
!        !d(i) = -scaling * gg(i)
!        d(i) = -scaling * g(i)
!    end do
!    !goto 370
!    goto 400
!end if

!!--------------------------------------------------------------------------------------------------!
!! Zaikun 20220303: Exit if H, G, TD, or TN are not finite. Otherwise, the behavior of this
!! subroutine is not predictable. For example, if HNORM = GNORM = Inf, it is observed that the
!! initial value of PARL defined below will change when we add code that should not affect PARL
!! (e.g., print it, or add TD = 0, TN = 0, PIV = 0 at the beginning of this subroutine).
!! This is probably because the behavior of MAX is undefined if it receives NaN (when GNORM = HNORM
!! = Inf, GNORM/DELTA - HNORM = NaN). This also motivates us to replace the intrinsic MAX by the
!! MAXIMUM defined in LINALG_MOD. MAXIMUM will return NaN if it receives NaN, making it easier for us
!! to notice that there is a problem and hence debug.
!if (.not. (is_finite(hnorm) .and. is_finite(gnorm) .and. all(is_finite(td(1:n))) .and. all(is_finite(tn(1:n))))) then
!    goto 400
!end if
!--------------------------------------------------------------------------------------------------!

!     Set the initial values of PAR and its bounds.
!
parl = maximum([ZERO, -tdmin, gnorm / delta - hnorm])
parlest = parl
par = parl
paru = ZERO
paruest = ZERO
posdef = ZERO
iterc = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 26-06-2019: See the lines below line number 140
do i = 1, n
    dsav(i) = d(i)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
140 iterc = iterc + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 26-06-2019
! The original code can encounter infinite cycling, which did happen
! when testing the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in H (up to 10^219).
! To avoid wasting energy, we do the following
sumd = ZERO
do i = 1, n
    sumd = sumd + abs(d(i))
end do
!if (sumd >= 1.0D100 .or. sumd /= sumd) then
if (.not. is_finite(sum(abs(d(1:n))))) then
    do i = 1, n
        d(i) = dsav(i)
    end do
    goto 370
else
    do i = 1, n
        dsav(i) = d(i)
    end do
end if
!if (iterc > min(10000, 100 * n)) then
if (iterc > min(1000, 100 * n)) then
    goto 370
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ksav = 0
piv = ZERO
piv(1) = td(1) + par
k = 1
150 if (piv(k) > ZERO) then
    piv(k + 1) = td(k + 1) + par - tn(k)**2 / piv(k)
else
    if (piv(k) < ZERO .or. tn(k) /= ZERO) goto 160
    ksav = k
    piv(k + 1) = td(k + 1) + par
end if
k = k + 1
if (k < n) goto 150  ! When the loop exits, K = N.

!---------------------------------------!
! Zaikun 20220508
!if (piv(k) < ZERO) goto 160
if (.not. piv(k) >= ZERO) goto 160
!---------------------------------------!
if (piv(k) == ZERO) ksav = k
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is ZERO.
!
if (ksav == 0 .and. gsq > ZERO) goto 230
if (gsq == ZERO) then
    if (par == ZERO) goto 370
    paru = par
    paruest = par
    if (ksav == 0) goto 190
end if
k = ksav
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
160 continue

!----------------------------!
if (any(is_nan(piv))) then
    goto 370
end if
!-----------------------------!

d(k) = ONE
if (abs(tn(k)) <= abs(piv(k))) then
    dsq = ONE
    dhd = piv(k)
else
    temp = td(k + 1) + par
    if (temp <= abs(piv(k))) then
        d(k + 1) = sign(ONE, -tn(k))
        dhd = piv(k) + temp - TWO * abs(tn(k))
    else
        d(k + 1) = -tn(k) / temp
        dhd = piv(k) + tn(k) * d(k + 1)
    end if
    dsq = ONE + d(k + 1)**2
end if
170 if (k > 1) then
    k = k - 1
    if (tn(k) /= ZERO) then
        d(k) = -tn(k) * d(k + 1) / piv(k)
        dsq = dsq + d(k)**2
        goto 170
    end if
    do i = 1, k
        d(i) = ZERO
    end do
end if
parl = par
parlest = par - dhd / dsq
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
190 temp = paruest
if (gsq == ZERO) temp = temp * (ONE - tol)
if (paruest > ZERO .and. parlest >= temp) then
    dtg = ZERO
    do i = 1, n
        dtg = dtg + d(i) * gg(i)
    end do
    scaling = -sign(delta / sqrt(dsq), dtg)
    do i = 1, n
        d(i) = scaling * d(i)
    end do
    goto 370
end if
!
!     Pick the value of PAR for the next iteration.
!
220 if (paru == ZERO) then
    par = TWO * parlest + gnorm / delta
else
    par = 0.5_RP * (parl + paru)
    par = max(par, parlest)
end if
if (paruest > ZERO) par = min(par, paruest)
goto 140
!
!     Calculate D for the current PAR in the positive definite case.
!
230 w(1) = -gg(1) / piv(1)
do i = 2, n
    w(i) = (-gg(i) - tn(i - 1) * w(i - 1)) / piv(i)
end do
d(n) = w(n)
do i = nm, 1, -1
    d(i) = w(i) - tn(i) * d(i + 1) / piv(i)
end do
!
!     Branch if a Newton-Raphson step is acceptable.
!
dsq = ZERO
wsq = ZERO
do i = 1, n
    dsq = dsq + d(i)**2
    wsq = wsq + piv(i) * w(i)**2
end do
if (par == ZERO .and. dsq <= delsq) goto 320
!
!     Make the usual test for acceptability of a full trust region step.
!
dnorm = sqrt(dsq)
phi = ONE / dnorm - ONE / delta
temp = tol * (ONE + par * dsq / wsq) - dsq * phi * phi
if (temp >= ZERO) then
    !scaling = delta / dnorm
    do i = 1, n
        !d(i) = scaling * d(i)
        d(i) = delta * (d(i) / dnorm)
    end do
    goto 370
end if
if (iterc >= 2 .and. par <= parl) goto 370
if (paru > ZERO .and. par >= paru) goto 370
!
!     Complete the iteration when PHI is negative.
!
if (phi < ZERO) then
    parlest = par
    if (posdef == ONE) then
        if (phi <= phil) goto 370
        slope = (phi - phil) / (par - parl)
        parlest = par - phi / slope
    end if
    slope = ONE / gnorm
    if (paru > ZERO) slope = (phiu - phi) / (paru - par)
    temp = par - phi / slope
    if (paruest > ZERO) temp = min(temp, paruest)
    paruest = temp
    posdef = ONE
    parl = par
    phil = phi
    goto 220
end if
!
!     If required, calculate Z for the alternative test for convergence.
!
if (posdef == ZERO) then
    w(1) = ONE / piv(1)
    do i = 2, n
        temp = -tn(i - 1) * w(i - 1)
        w(i) = (sign(ONE, temp) + temp) / piv(i)
    end do
    z(n) = w(n)
    do i = nm, 1, -1
        z(i) = w(i) - tn(i) * z(i + 1) / piv(i)
    end do
    wwsq = ZERO
    zsq = ZERO
    dtz = ZERO
    do i = 1, n
        wwsq = wwsq + piv(i) * w(i)**2
        zsq = zsq + z(i)**2
        dtz = dtz + d(i) * z(i)
    end do
!
!     Apply the alternative test for convergence.
!
    tempa = abs(delsq - dsq)
    tempb = sqrt(dtz * dtz + tempa * zsq)
    gam = tempa / (sign(tempb, dtz) + dtz)
    temp = tol * (wsq + par * delsq) - gam * gam * wwsq
    if (temp >= ZERO) then
        do i = 1, n
            d(i) = d(i) + gam * z(i)
        end do
        goto 370
    end if
    parlest = max(parlest, par - wwsq / zsq)
end if
!
!     Complete the iteration when PHI is positive.
!
slope = ONE / gnorm
if (paru > ZERO) then
    if (phi >= phiu) goto 370
    slope = (phiu - phi) / (paru - par)
end if
parlest = max(parlest, par - phi / slope)
paruest = par
if (posdef == ONE) then
    slope = (phi - phil) / (par - parl)
    paruest = par - phi / slope
end if
paru = par
phiu = phi
goto 220
!
!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
!
320 shfmin = ZERO
pivot = td(1)
shfmax = pivot
do k = 2, n
    pivot = td(k) - tn(k - 1)**2 / pivot
    shfmax = min(shfmax, pivot)
end do
!
!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
ksave = 0
340 shift = 0.5_RP * (shfmin + shfmax)
k = 1
temp = td(1) - shift
350 if (temp > ZERO) then
    piv(k) = temp
    if (k < n) then
        temp = td(k + 1) - shift - tn(k)**2 / temp
        k = k + 1
        goto 350
    end if
    shfmin = shift
else
    if (k < ksave) goto 360
    if (k == ksave) then
        if (pivksv == ZERO) goto 360
        if (piv(k) - temp < temp - pivksv) then
            pivksv = temp
            shfmax = shift
        else
            pivksv = ZERO
            shfmax = (shift * piv(k) - shfmin * temp) / (piv(k) - temp)
        end if
    else
        ksave = k
        pivksv = temp
        shfmax = shift
    end if
end if
if (shfmin <= 0.99_RP * shfmax) goto 340
360 evalue = shfmin
!
!     Apply the inverse Householder transformations to D.
!
370 nm = n - 1
do k = nm, 1, -1
    kp = k + 1
    summ = ZERO
    do i = kp, n
        summ = summ + d(i) * h(i, k)
    end do
    do i = kp, n
        d(i) = d(i) - summ * h(i, k)
    end do
end do
!
!     Return from the subroutine.
!
400 d_out(1:n) = d(1:n)
return
end subroutine trstep


end module trustregion_mod
