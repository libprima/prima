module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Sunday, February 27, 2022 AM12:49:05
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, knew, adelt, xnew, xalt, alpha, cauchy, glag, hcol, w)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF

implicit none

! Inputs
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
integer(IK), intent(in) :: n
integer(IK), intent(in) :: ndim
integer(IK), intent(in) :: npt
real(RP), intent(in) :: adelt
real(RP), intent(in) :: bmat(n, npt + n)
real(RP), intent(in) :: sl(n)
real(RP), intent(in) :: su(n)
real(RP), intent(in) :: xopt(n)
real(RP), intent(in) :: xpt(n, npt)
real(RP), intent(in) :: zmat(npt, npt - n - 1_IK)

! In-outputs
real(RP), intent(inout) :: alpha
real(RP), intent(inout) :: cauchy
real(RP), intent(inout) :: glag(n)
real(RP), intent(inout) :: hcol(npt)
real(RP), intent(inout) :: w(2_IK * n)
real(RP), intent(inout) :: xalt(n)

! Outputs
real(RP), intent(out) :: xnew(n)

! Local variables
real(RP) :: bigstp, const, csave, curv, dderiv, diff, distsq,  &
&        ggfree, gw, ha, predsq, presav, scaling, &
&        slbd, step, stpsav, subd, sumin, temp, tempa,      &
&        tempb, tempd, vlag, wfixsq, wsqsav
integer(IK) :: i, ibdsav, iflag, ilbd, isbd, iubd, j, k, ksav

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
!       function subject to the constraints that have been mentiONEd, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all compONEnts of GLAG (see below) are ZERO.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to ZERO if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT compONEnts of W to the leading elements of the
!     KNEW-th column of the H matrix.
!

const = ONE + sqrt(TWO)
do k = 1, npt
    hcol(k) = ZERO
end do
do j = 1, npt - n - 1
    temp = zmat(knew, j)
    do k = 1, npt
        hcol(k) = hcol(k) + temp * zmat(k, j)
    end do
end do
alpha = hcol(knew)
ha = HALF * alpha
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
do i = 1, n
    glag(i) = bmat(i, knew)
end do
do k = 1, npt
    temp = ZERO
    do j = 1, n
        temp = temp + xpt(j, k) * xopt(j)
    end do
    temp = hcol(k) * temp
    do i = 1, n
        glag(i) = glag(i) + temp * xpt(i, k)
    end do
end do
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
presav = ZERO
do k = 1, npt
    if (k == kopt) cycle
    dderiv = ZERO
    distsq = ZERO
    do i = 1, n
        temp = xpt(i, k) - xopt(i)
        dderiv = dderiv + glag(i) * temp
        distsq = distsq + temp * temp
    end do
    subd = adelt / sqrt(distsq)
    slbd = -subd
    ilbd = 0
    iubd = 0
    sumin = min(ONE, subd)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
    do i = 1, n
        temp = xpt(i, k) - xopt(i)
        if (temp > ZERO) then
            if (slbd * temp < sl(i) - xopt(i)) then
                slbd = (sl(i) - xopt(i)) / temp
                ilbd = -i
            end if
            if (subd * temp > su(i) - xopt(i)) then
                subd = max(sumin, (su(i) - xopt(i)) / temp)
                iubd = i
            end if
        else if (temp < ZERO) then
            if (slbd * temp > su(i) - xopt(i)) then
                slbd = (su(i) - xopt(i)) / temp
                ilbd = i
            end if
            if (subd * temp < sl(i) - xopt(i)) then
                subd = max(sumin, (sl(i) - xopt(i)) / temp)
                iubd = -i
            end if
        end if
    end do
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
    if (k == knew) then
        diff = dderiv - ONE
        step = slbd
        vlag = slbd * (dderiv - slbd * diff)
        isbd = ilbd
        temp = subd * (dderiv - subd * diff)
        if (abs(temp) > abs(vlag)) then
            step = subd
            vlag = temp
            isbd = iubd
        end if
        tempd = HALF * dderiv
        tempa = tempd - diff * slbd
        tempb = tempd - diff * subd
        if (tempa * tempb < ZERO) then
            temp = tempd * tempd / diff
            if (abs(temp) > abs(vlag)) then
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
        vlag = slbd * (ONE - slbd)
        isbd = ilbd
        temp = subd * (ONE - subd)
        if (abs(temp) > abs(vlag)) then
            step = subd
            vlag = temp
            isbd = iubd
        end if
        if (subd > HALF) then
            if (abs(vlag) < 0.25_RP) then
                step = HALF
                vlag = 0.25_RP
                isbd = 0
            end if
        end if
        vlag = vlag * dderiv
    end if
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
    temp = step * (ONE - step) * distsq
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
    temp = xopt(i) + stpsav * (xpt(i, ksav) - xopt(i))
    xnew(i) = max(sl(i), min(su(i), temp))
end do
if (ibdsav < 0) xnew(-ibdsav) = sl(-ibdsav)
if (ibdsav > 0) xnew(ibdsav) = su(ibdsav)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed compONEnts of W is formed in
!     WFIXSQ, and the free compONEnts of W are set to BIGSTP.
!
bigstp = adelt + adelt
iflag = 0
100 wfixsq = ZERO
ggfree = ZERO
do i = 1, n
    w(i) = ZERO
    tempa = min(xopt(i) - sl(i), glag(i))
    tempb = max(xopt(i) - su(i), glag(i))
    if (tempa > ZERO .or. tempb < ZERO) then
        w(i) = bigstp
        ggfree = ggfree + glag(i)**2
    end if
end do
if (ggfree == ZERO) then
    cauchy = ZERO
    goto 200
end if
!
!     Investigate whether more compONEnts of W can be fixed.
!
120 temp = adelt * adelt - wfixsq
if (temp > ZERO) then
    wsqsav = wfixsq
    step = sqrt(temp / ggfree)
    ggfree = ZERO
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
    if (wfixsq > wsqsav .and. ggfree > ZERO) goto 120
end if
!
!     Set the remaining free compONEnts of W and all compONEnts of XALT,
!     except that W may be scaled later.
!
gw = ZERO
do i = 1, n
    if (w(i) == bigstp) then
        w(i) = -step * glag(i)
        xalt(i) = max(sl(i), min(su(i), xopt(i) + w(i)))
    else if (w(i) == ZERO) then
        xalt(i) = xopt(i)
    else if (glag(i) > ZERO) then
        xalt(i) = sl(i)
    else
        xalt(i) = su(i)
    end if
    gw = gw + glag(i) * w(i)
end do
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than ONE if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
curv = ZERO
do k = 1, npt
    temp = ZERO
    do j = 1, n
        temp = temp + xpt(j, k) * w(j)
    end do
    curv = curv + hcol(k) * temp * temp
end do
if (iflag == 1) curv = -curv
if (curv > -gw .and. curv < -const * gw) then
    scaling = -gw / curv
    do i = 1, n
        temp = xopt(i) + scaling * w(i)
        xalt(i) = max(sl(i), min(su(i), temp))
    end do
    cauchy = (HALF * gw * scaling)**2
else
    cauchy = (gw + HALF * curv)**2
end if
!
!     If IFLAG is ZERO, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The ONE that
!     is chosen is the ONE that gives the larger value of CAUCHY.
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
end subroutine geostep


end module geometry_mod
