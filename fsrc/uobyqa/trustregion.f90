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
! Last Modified: Sunday, May 08, 2022 PM04:13:59
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstep


contains


subroutine trstep(delta, g, h_in, tol, d_out, crvmin)   !!!! Possible to use D instead of D_OUT?

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : maximum, issymmetric, inprod, hessenberg

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: h_in(:, :)  ! H_IN(N, N)
real(RP), intent(in) :: tol

! In-outputs
!real(RP), intent(out) :: d(:)  ! D(N)
real(RP), intent(out) :: d_out(:)  ! D(N)  !!!! Temporary; the code below accesses D(N+1)
real(RP), intent(out) :: crvmin

! Local variables
character(len=*), parameter :: srname = 'TRSTEP'
integer(IK) :: n
real(RP) :: d(size(g) + 1) !!! D(N+1) may be accessed. !!! Possible to avoid this and spare D_OUT?
real(RP) :: gg(size(g))
real(RP) :: h(size(g), size(g))
real(RP) :: piv(size(g))
real(RP) :: td(size(g) + 1) !!! TD(N+1) may be accessed
real(RP) :: tn(size(g))
real(RP) :: w(size(g))
real(RP) :: z(size(g))
real(RP) :: dold(size(g)) !!!
real(RP) :: dnewton(size(g))  ! Newton-Raphson step; only calculated when N = 1.
real(RP) :: delsq, dhd, dnorm, dsq, dtg, dtz, gam, gnorm,     &
&        gsq, hnorm, par, parl, parlest, paru,         &
&        paruest, phi, phil, phiu, pivksv, pivot, posdef,   &
&        scaling, shfmax, shfmin, shift, slope,   &
&        tdmin, temp, tempa, tempb, wsq, wwsq, zsq
integer(IK) :: i, iter, k, ksav, ksave, maxiter

h = h_in  ! To be removed

! Sizes.
n = int(size(g), kind(n))

if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(size(h, 1) == n .and. issymmetric(h), 'H is n-by-n and symmetric', srname)
    call assert(size(d_out) == n, 'SIZE(D) == N', srname)
end if


crvmin = ZERO

! Zaikun 20220301, 20220305:
! Powell's original code requires that N >= 2.  When N = 1, the code does not work (sometimes even
! encounter memory errors). This is indeed why the original version of UOBYQA constantly terminates
! with "a trust region step has failed to reduce the quadratic model" when solving univariate problems.
if (n == 1) then
    if (h(1, 1) > 0) then
        dnewton = -g / h(1, 1)
        if (abs(dnewton(1)) <= delta) then
            d_out = dnewton
            crvmin = h(1, 1)
        else
            d_out = sign(delta, dnewton)  ! MATLAB: D_OUT = DELTA * SIGN(DNEWTON)
        end if
    else
        d_out = sign(delta, -g)  ! MATLAB: D_OUT = -DELTA * SIGN(G)
    end if
    return
end if

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
!     CRVMIN will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step. Then CRVMIN will be positive, but otherwise it
!     will be set to ZERO.
!
!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| <= DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties
!
!             ||D|| <= DELTA  and  ACTRED <= (1-TOL)*MAXRED.
!
!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.
!
!     Initialization.
!
!
delsq = delta * delta
d = ZERO

! Apply Householder transformations to obtain a tridiagonal matrix that is similar to H, and put the
! elements of the Householder vectors in the lower triangular part of H. Further, TD and TN will
! contain the diagonal and other nonzero elements of the tridiagonal matrix.
call hessenberg(h, td(1:n), tn(1:n - 1))
!!MATLAB: [P, h] = hess(h); td = diag(h); tn = diag(h, 1)
tn(n) = ZERO  ! This is necessary, as TN(N) will be accessed.


!!!!!!!!!!!! TODO: TD should have length N, and TN should have length N-1. !!!!!!!!!!!!


! Form GG by applying the similarity transformation to G.
gg = g
gsq = sum(g**2)
gnorm = sqrt(gsq)
do k = 1, n - 1
    gg(k + 1:n) = gg(k + 1:n) - inprod(gg(k + 1:n), h(k + 1:n, k)) * h(k + 1:n, k)
end do
!!MATLAB: gg = (g'*P)';  % gg = P'*g;

!--------------------------------------------------------------------------------------------------!
! Zaikun 20220303: Exit if H, G, TD, or TN are not finite. Otherwise, the behavior of this
! subroutine is not predictable. For example, if HNORM = GNORM = Inf, it is observed that the
! initial value of PARL defined below will change when we add code that should not affect PARL
! (e.g., print it, or add TD = 0, TN = 0, PIV = 0 at the beginning of this subroutine).
! This is probably because the behavior of MAX is undefined if it receives NaN (when GNORM = HNORM
! = Inf, GNORM/DELTA - HNORM = NaN). This also motivates us to replace the intrinsic MAX by the
! MAXIMUM defined in LINALG_MOD. MAXIMUM will return NaN if it receives NaN, making it easier for us
! to notice that there is a problem and hence debug.
!--------------------------------------------------------------------------------------------------!
if (.not. is_finite(sum(abs(h)) + sum(abs(td(1:n))) + sum(abs(tn)) + sum(abs(gg)))) then
    d_out = ZERO
    crvmin = ZERO
    return
end if

! Begin the trust region calculation with a tridiagonal matrix by calculating the norm of H. Then
! treat the case when H is zero.

hnorm = maxval(abs([ZERO, tn(1:n - 1)]) + abs(td(1:n)) + abs(tn))
tdmin = minval(td(1:n))  ! This leads to a difference. Why?

if (hnorm == ZERO) then
    if (gnorm == ZERO) goto 400
    scaling = delta / gnorm
    d(1:n) = -scaling * g
    goto 400
end if

! Set the initial values of PAR and its bounds.
!parl = max(ZERO, -tdmin, gnorm / delta - hnorm)
parl = maximum([ZERO, -tdmin, gnorm / delta - hnorm])
parlest = parl
par = parl
paru = ZERO
paruest = ZERO
posdef = ZERO
iter = 0
maxiter = min(1000_IK, 100_IK * int(n, IK))  ! What is the theoretical bound of iterations?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 26-06-2019: See the lines below line number 140
!do i = 1, n
!    dold(i) = d(i)
!end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
140 continue

iter = iter + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 26-06-2019
! The original code can encounter infinite cycling, which did happen
! when testing the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in H (up to 10^219).
! To avoid wasting energy, we do the following
if (.not. is_finite(sum(abs(d(1:n))))) then
    d(1:n) = dold
    goto 370
else
    dold = d(1:n)
end if
if (iter > maxiter) then
    goto 370
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ksav = 0
piv(1) = td(1) + par
do k = 1, n - 1
    if (piv(k) > ZERO) then
        piv(k + 1) = td(k + 1) + par - tn(k)**2 / piv(k)
    else
        if (piv(k) < ZERO .or. tn(k) /= ZERO) goto 160
        ksav = k
        piv(k + 1) = td(k + 1) + par
    end if
end do
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (piv(n) < ZERO) goto 160
if (piv(n) == ZERO) ksav = n
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

d(k) = ONE
if (abs(tn(k)) <= abs(piv(k))) then
    dsq = ONE
    dhd = piv(k)
else
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 20220301: The code below accesses TD(N+1), D(N+1) when K = N!!!
! Zaikun 20220507: Is K=N possible?
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
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
170 continue
if (k > 1) then
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
    par = HALF * (parl + paru)
    par = max(par, parlest)
end if
if (paruest > ZERO) par = min(par, paruest)
goto 140
!
!     Calculate D for the current PAR in the positive definite case.
!
230 continue
w(1) = -gg(1) / piv(1)
do i = 2, n
    w(i) = (-gg(i) - tn(i - 1) * w(i - 1)) / piv(i)
end do
d(n) = w(n)
do i = n - 1, 1, -1
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
    scaling = delta / dnorm
    do i = 1, n
        d(i) = scaling * d(i)
    end do
    goto 370
end if
if (iter >= 2 .and. par <= parl) goto 370
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
    do i = n - 1, 1, -1
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
!     Set CRVMIN to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on CRVMIN.
!
320 shfmin = ZERO
pivot = td(1)
shfmax = pivot
do k = 2, n
    pivot = td(k) - tn(k - 1)**2 / pivot
    shfmax = min(shfmax, pivot)
end do
!
!     Find CRVMIN by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
ksave = 0
340 shift = HALF * (shfmin + shfmax)
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

360 continue

crvmin = shfmin
!
!     Apply the inverse Householder transformations to D.
!
370 continue

do k = n - 1, 1, -1
    d(k + 1:n) = d(k + 1:n) - inprod(d(k + 1:n), h(k + 1:n, k)) * h(k + 1:n, k)
end do
!!MATLAB: d = P*d;

400 continue

d_out(1:n) = d(1:n)  !!!! Temporary

end subroutine trstep


end module trustregion_mod
