subroutine trstep(n, g, h, delta, tol, d, gg, td, tn, w, piv, z, evalue)
implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
dimension g(*), h(n, *), d(*), gg(*), td(*), tn(*), w(*), piv(*), z(*), dsav(n)
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
!     will be set to zero.
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
one = 1.0D0
two = 2.0D0
zero = 0.0D0
delsq = delta * delta
evalue = zero
nm = n - 1
do i = 1, n
    d(i) = zero
    td(i) = h(i, i)
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
do k = 1, nm
    kp = k + 1
    sum = zero
    if (kp < n) then
        kpp = kp + 1
        do i = kpp, n
            sum = sum + h(i, k)**2
        end do
    end if
    if (sum == zero) then
        tn(k) = h(kp, k)
        h(kp, k) = zero
    else
        temp = h(kp, k)
        tn(k) = dsign(sqrt(sum + temp * temp), temp)
        h(kp, k) = -sum / (temp + tn(k))
        temp = sqrt(two / (sum + h(kp, k)**2))
        do i = kp, n
            w(i) = temp * h(i, k)
            h(i, k) = w(i)
            z(i) = td(i) * w(i)
        end do
        wz = zero
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
            td(j) = td(j) + w(j) * (wz * w(j) - two * z(j))
            if (j < n) then
                jp = j + 1
                do i = jp, n
                    h(i, j) = h(i, j) - w(i) * z(j) - w(j) * (z(i) - wz * w(i))
                end do
            end if
        end do
    end if
end do
!
!     Form GG by applying the similarity transformation to G.
!
gsq = zero
do i = 1, n
    gg(i) = g(i)
    gsq = gsq + g(i)**2
end do
gnorm = sqrt(gsq)
do k = 1, nm
    kp = k + 1
    sum = zero
    do i = kp, n
        sum = sum + gg(i) * h(i, k)
    end do
    do i = kp, n
        gg(i) = gg(i) - sum * h(i, k)
    end do
end do
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.
!
hnorm = abs(td(1)) + abs(tn(1))
tdmin = td(1)
tn(n) = zero
do i = 2, n
    temp = abs(tn(i - 1)) + abs(td(i)) + abs(tn(i))
    hnorm = max(hnorm, temp)
    tdmin = min(tdmin, td(i))
end do
if (hnorm == zero) then
    if (gnorm == zero) goto 400
    scale = delta / gnorm
    do i = 1, n
        d(i) = -scale * gg(i)
    end do
    goto 370
end if
!
!     Set the initial values of PAR and its bounds.
!
parl = max(zero, -tdmin, gnorm / delta - hnorm)
parlest = parl
par = parl
paru = zero
paruest = zero
posdef = zero
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
sumd = zero
do i = 1, n
    sumd = sumd + abs(d(i))
end do
if (sumd >= 1.0D100 .or. sumd /= sumd) then
    do i = 1, n
        d(i) = dsav(i)
    end do
    goto 370
else
    do i = 1, n
        dsav(i) = d(i)
    end do
end if
if (iterc > min(10000, 100 * n)) then
    goto 370
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ksav = 0
piv(1) = td(1) + par
k = 1
150 if (piv(k) > zero) then
    piv(k + 1) = td(k + 1) + par - tn(k)**2 / piv(k)
else
    if (piv(k) < zero .or. tn(k) /= zero) goto 160
    ksav = k
    piv(k + 1) = td(k + 1) + par
end if
k = k + 1
if (k < n) goto 150
if (piv(k) < zero) goto 160
if (piv(k) == zero) ksav = k
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
!
if (ksav == 0 .and. gsq > zero) goto 230
if (gsq == zero) then
    if (par == zero) goto 370
    paru = par
    paruest = par
    if (ksav == 0) goto 190
end if
k = ksav
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
160 d(k) = one
if (abs(tn(k)) <= abs(piv(k))) then
    dsq = one
    dhd = piv(k)
else
    temp = td(k + 1) + par
    if (temp <= abs(piv(k))) then
        d(k + 1) = dsign(one, -tn(k))
        dhd = piv(k) + temp - two * abs(tn(k))
    else
        d(k + 1) = -tn(k) / temp
        dhd = piv(k) + tn(k) * d(k + 1)
    end if
    dsq = one + d(k + 1)**2
end if
170 if (k > 1) then
    k = k - 1
    if (tn(k) /= zero) then
        d(k) = -tn(k) * d(k + 1) / piv(k)
        dsq = dsq + d(k)**2
        goto 170
    end if
    do i = 1, k
        d(i) = zero
    end do
end if
parl = par
parlest = par - dhd / dsq
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
190 temp = paruest
if (gsq == zero) temp = temp * (one - tol)
if (paruest > zero .and. parlest >= temp) then
    dtg = zero
    do i = 1, n
        dtg = dtg + d(i) * gg(i)
    end do
    scale = -dsign(delta / sqrt(dsq), dtg)
    do i = 1, n
        d(i) = scale * d(i)
    end do
    goto 370
end if
!
!     Pick the value of PAR for the next iteration.
!
220 if (paru == zero) then
    par = two * parlest + gnorm / delta
else
    par = 0.5D0 * (parl + paru)
    par = max(par, parlest)
end if
if (paruest > zero) par = min(par, paruest)
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
dsq = zero
wsq = zero
do i = 1, n
    dsq = dsq + d(i)**2
    wsq = wsq + piv(i) * w(i)**2
end do
if (par == zero .and. dsq <= delsq) goto 320
!
!     Make the usual test for acceptability of a full trust region step.
!
dnorm = sqrt(dsq)
phi = one / dnorm - one / delta
temp = tol * (one + par * dsq / wsq) - dsq * phi * phi
if (temp >= zero) then
    scale = delta / dnorm
    do i = 1, n
        d(i) = scale * d(i)
    end do
    goto 370
end if
if (iterc >= 2 .and. par <= parl) goto 370
if (paru > zero .and. par >= paru) goto 370
!
!     Complete the iteration when PHI is negative.
!
if (phi < zero) then
    parlest = par
    if (posdef == one) then
        if (phi <= phil) goto 370
        slope = (phi - phil) / (par - parl)
        parlest = par - phi / slope
    end if
    slope = one / gnorm
    if (paru > zero) slope = (phiu - phi) / (paru - par)
    temp = par - phi / slope
    if (paruest > zero) temp = min(temp, paruest)
    paruest = temp
    posdef = one
    parl = par
    phil = phi
    goto 220
end if
!
!     If required, calculate Z for the alternative test for convergence.
!
if (posdef == zero) then
    w(1) = one / piv(1)
    do i = 2, n
        temp = -tn(i - 1) * w(i - 1)
        w(i) = (dsign(one, temp) + temp) / piv(i)
    end do
    z(n) = w(n)
    do i = nm, 1, -1
        z(i) = w(i) - tn(i) * z(i + 1) / piv(i)
    end do
    wwsq = zero
    zsq = zero
    dtz = zero
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
    gam = tempa / (dsign(tempb, dtz) + dtz)
    temp = tol * (wsq + par * delsq) - gam * gam * wwsq
    if (temp >= zero) then
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
slope = one / gnorm
if (paru > zero) then
    if (phi >= phiu) goto 370
    slope = (phiu - phi) / (paru - par)
end if
parlest = max(parlest, par - phi / slope)
paruest = par
if (posdef == one) then
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
320 shfmin = zero
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
340 shift = 0.5D0 * (shfmin + shfmax)
k = 1
temp = td(1) - shift
350 if (temp > zero) then
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
        if (pivksv == zero) goto 360
        if (piv(k) - temp < temp - pivksv) then
            pivksv = temp
            shfmax = shift
        else
            pivksv = zero
            shfmax = (shift * piv(k) - shfmin * temp) / (piv(k) - temp)
        end if
    else
        ksave = k
        pivksv = temp
        shfmax = shift
    end if
end if
if (shfmin <= 0.99D0 * shfmax) goto 340
360 evalue = shfmin
!
!     Apply the inverse Householder transformations to D.
!
370 nm = n - 1
do k = nm, 1, -1
    kp = k + 1
    sum = zero
    do i = kp, n
        sum = sum + d(i) * h(i, k)
    end do
    do i = kp, n
        d(i) = d(i) - sum * h(i, k)
    end do
end do
!
!     Return from the subroutine.
!
400 return
end
