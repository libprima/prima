subroutine altmov(n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, knew, adelt, xnew, xalt, alpha, cauchy, glag, hcol, w)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension xpt(npt, *), xopt(*), bmat(ndim, *), zmat(npt, *), sl(*), su(*), xnew(*), xalt(*), glag(*), hcol(*), w(*)
!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!
half = 0.5D0
one = 1.0D0
zero = 0.0D0
const = one + dsqrt(2.0D0)
do k = 1, npt
    hcol(k) = zero
end do
do j = 1, npt - n - 1
    temp = zmat(knew, j)
    do k = 1, npt
        hcol(k) = hcol(k) + temp * zmat(k, j)
    end do
end do
alpha = hcol(knew)
ha = half * alpha
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
do i = 1, n
    glag(i) = bmat(knew, i)
end do
do k = 1, npt
    temp = zero
    do j = 1, n
        temp = temp + xpt(k, j) * xopt(j)
    end do
    temp = hcol(k) * temp
    do i = 1, n
        glag(i) = glag(i) + temp * xpt(k, i)
    end do
end do
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
presav = zero
do k = 1, npt
    if (k == kopt) cycle
    dderiv = zero
    distsq = zero
    do i = 1, n
        temp = xpt(k, i) - xopt(i)
        dderiv = dderiv + glag(i) * temp
        distsq = distsq + temp * temp
    end do
    subd = adelt / dsqrt(distsq)
    slbd = -subd
    ilbd = 0
    iubd = 0
    sumin = dmin1(one, subd)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
    do i = 1, n
        temp = xpt(k, i) - xopt(i)
        if (temp > zero) then
            if (slbd * temp < sl(i) - xopt(i)) then
                slbd = (sl(i) - xopt(i)) / temp
                ilbd = -i
            end if
            if (subd * temp > su(i) - xopt(i)) then
                subd = dmax1(sumin, (su(i) - xopt(i)) / temp)
                iubd = i
            end if
        else if (temp < zero) then
            if (slbd * temp > su(i) - xopt(i)) then
                slbd = (su(i) - xopt(i)) / temp
                ilbd = i
            end if
            if (subd * temp < sl(i) - xopt(i)) then
                subd = dmax1(sumin, (sl(i) - xopt(i)) / temp)
                iubd = -i
            end if
        end if
    end do
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
    if (k == knew) then
        diff = dderiv - one
        step = slbd
        vlag = slbd * (dderiv - slbd * diff)
        isbd = ilbd
        temp = subd * (dderiv - subd * diff)
        if (dabs(temp) > dabs(vlag)) then
            step = subd
            vlag = temp
            isbd = iubd
        end if
        tempd = half * dderiv
        tempa = tempd - diff * slbd
        tempb = tempd - diff * subd
        if (tempa * tempb < zero) then
            temp = tempd * tempd / diff
            if (dabs(temp) > dabs(vlag)) then
                step = tempd / diff
                vlag = temp
                isbd = 0
            end if
        end if
!
!     Search along each of the other lines through XOPT and another point.
!
    else
        step = slbd
        vlag = slbd * (one - slbd)
        isbd = ilbd
        temp = subd * (one - subd)
        if (dabs(temp) > dabs(vlag)) then
            step = subd
            vlag = temp
            isbd = iubd
        end if
        if (subd > half) then
            if (dabs(vlag) < 0.25D0) then
                step = half
                vlag = 0.25D0
                isbd = 0
            end if
        end if
        vlag = vlag * dderiv
    end if
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
    temp = step * (one - step) * distsq
    predsq = vlag * vlag * (vlag * vlag + ha * temp * temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: With the original code, if either PREDSQ or PRESAV
! is NaN, KSAV/STPSAV/IBDSAV will not get a value. This may cause
! Segmentation Fault.
!      IF (PREDSQ .GT. PRESAV) THEN
    if (.not. (predsq <= presav)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        presav = predsq
        ksav = k
        stpsav = step
        ibdsav = isbd
    end if
end do
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
do i = 1, n
    temp = xopt(i) + stpsav * (xpt(ksav, i) - xopt(i))
    xnew(i) = dmax1(sl(i), dmin1(su(i), temp))
end do
if (ibdsav < 0) xnew(-ibdsav) = sl(-ibdsav)
if (ibdsav > 0) xnew(ibdsav) = su(ibdsav)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
bigstp = adelt + adelt
iflag = 0
100 wfixsq = zero
ggfree = zero
do i = 1, n
    w(i) = zero
    tempa = dmin1(xopt(i) - sl(i), glag(i))
    tempb = dmax1(xopt(i) - su(i), glag(i))
    if (tempa > zero .or. tempb < zero) then
        w(i) = bigstp
        ggfree = ggfree + glag(i)**2
    end if
end do
if (ggfree == zero) then
    cauchy = zero
    goto 200
end if
!
!     Investigate whether more components of W can be fixed.
!
120 temp = adelt * adelt - wfixsq
if (temp > zero) then
    wsqsav = wfixsq
    step = dsqrt(temp / ggfree)
    ggfree = zero
    do i = 1, n
        if (w(i) == bigstp) then
            temp = xopt(i) - step * glag(i)
            if (temp <= sl(i)) then
                w(i) = sl(i) - xopt(i)
                wfixsq = wfixsq + w(i)**2
            else if (temp >= su(i)) then
                w(i) = su(i) - xopt(i)
                wfixsq = wfixsq + w(i)**2
            else
                ggfree = ggfree + glag(i)**2
            end if
        end if
    end do
    if (wfixsq > wsqsav .and. ggfree > zero) goto 120
end if
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
gw = zero
do i = 1, n
    if (w(i) == bigstp) then
        w(i) = -step * glag(i)
        xalt(i) = dmax1(sl(i), dmin1(su(i), xopt(i) + w(i)))
    else if (w(i) == zero) then
        xalt(i) = xopt(i)
    else if (glag(i) > zero) then
        xalt(i) = sl(i)
    else
        xalt(i) = su(i)
    end if
    gw = gw + glag(i) * w(i)
end do
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
curv = zero
do k = 1, npt
    temp = zero
    do j = 1, n
        temp = temp + xpt(k, j) * w(j)
    end do
    curv = curv + hcol(k) * temp * temp
end do
if (iflag == 1) curv = -curv
if (curv > -gw .and. curv < -const * gw) then
    scale = -gw / curv
    do i = 1, n
        temp = xopt(i) + scale * w(i)
        xalt(i) = dmax1(sl(i), dmin1(su(i), temp))
    end do
    cauchy = (half * gw * scale)**2
else
    cauchy = (gw + half * curv)**2
end if
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
if (iflag == 0) then
    do i = 1, n
        glag(i) = -glag(i)
        w(n + i) = xalt(i)
    end do
    csave = cauchy
    iflag = 1
    goto 100
end if
if (csave > cauchy) then
    do i = 1, n
        xalt(i) = w(n + i)
    end do
    cauchy = csave
end if
200 return
end
