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
! Last Modified: Wednesday, April 27, 2022 PM09:30:16
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(knew, kopt, adelt, bmat, sl, su, xopt, xpt, zmat, alpha, cauchy, xalt, xnew)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, inprod, trueloc

implicit none

! Inputs
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: adelt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: sl(:)  ! SL(N)
real(RP), intent(in) :: su(:)  ! SU(N)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
real(RP), intent(out) :: alpha
real(RP), intent(out) :: cauchy
real(RP), intent(out) :: xalt(:)  ! XALT(N)
real(RP), intent(out) :: xnew(:)  ! XNEW(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: glag(size(xpt, 1))
real(RP) :: hcol(size(xpt, 2))
real(RP) :: w(2 * size(xpt, 1))
real(RP) :: bigstp, csave, curv, dderiv, diff, distsq,  &
&        ggfree, gw, ha, predsq, presav, scaling, &
&        slbd, step, stpsav, subd, sumin, temp, tempa,      &
&        tempb, tempd, vlag, wfixsq, wsqsav, xtest(size(xpt, 1))
logical :: mask_fixl(size(xpt, 1)), mask_fixu(size(xpt, 1)), mask_free(size(xpt, 1))
integer(IK) :: i, ibdsav, iflag, ilbd, isbd, iubd, j, k, ksav


! Sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt, 'KNEW /= KOPT', srname)
    call assert(adelt > 0, 'ADELT > 0', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(xopt) == n, 'SIZE(XOPT) == N', srname)
    call assert(size(xalt) == n, 'SIZE(XALT) == N', srname)
    call assert(size(xnew) == n, 'SIZE(XNEW) == N', srname)
end if

!--------------------------------------------------------------------------------------------------!
! Zaikun 20220305: Temporary fix for G95 warning: ‘ksav’/'ibdsav' may be used uninitialized in this function
ksav = 1_IK; ibdsav = 1_IK
!--------------------------------------------------------------------------------------------------!

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
!       not calculated if all components of GLAG (see below) are ZERO.
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
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!

hcol = matprod(zmat, zmat(knew, :))
alpha = hcol(knew)
ha = HALF * alpha

!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
!!glag = bmat(:, knew) + hess_mul(xopt, xpt, hcol)
glag = bmat(:, knew)
do k = 1, npt
    glag = glag + hcol(k) * inprod(xopt, xpt(:, k)) * xpt(:, k)
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
    dderiv = inprod(glag, xpt(:, k) - xopt)
    distsq = sum((xpt(:, k) - xopt)**2)
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
            if (abs(vlag) < QUART) then
                step = HALF
                vlag = QUART
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
xnew = max(sl, min(su, xopt + stpsav * (xpt(:, ksav) - xopt)))
if (ibdsav < 0) xnew(-ibdsav) = sl(-ibdsav)
if (ibdsav > 0) xnew(ibdsav) = su(ibdsav)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
bigstp = adelt + adelt
iflag = 0

100 continue

w(1:n) = ZERO
mask_free = (min(xopt - sl, glag) > 0 .or. max(xopt - su, glag) < 0)
w(trueloc(mask_free)) = bigstp
ggfree = sum(glag(trueloc(mask_free))**2)

if (ggfree <= ZERO) then
    cauchy = ZERO
    return
end if
!
!     Investigate whether more components of W can be fixed.
!
wfixsq = ZERO
do k = 1, n
    temp = adelt**2 - wfixsq
    if (.not. temp > 0) exit
    wsqsav = wfixsq
    step = sqrt(temp / ggfree)
    xtest = xopt - step * glag
    mask_fixl = (w(1:n) == bigstp .and. xtest <= sl)
    mask_fixu = (w(1:n) == bigstp .and. xtest >= su)
    mask_free = (w(1:n) == bigstp .and. .not. (mask_fixl .or. mask_fixu))
    w(trueloc(mask_fixl)) = sl(trueloc(mask_fixl)) - xopt(trueloc(mask_fixl))
    w(trueloc(mask_fixu)) = su(trueloc(mask_fixu)) - xopt(trueloc(mask_fixu))
    wfixsq = wfixsq + sum(w(trueloc(mask_fixl .or. mask_fixu))**2)
    ggfree = sum(glag(trueloc(mask_free))**2)
    if (.not. (wfixsq > wsqsav .and. ggfree > ZERO)) exit
end do
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
where (w == bigstp)
    w = -step * glag
    xalt = max(sl, min(su, xopt + w))
elsewhere(w == ZERO)
    xalt = xopt
elsewhere(glag > ZERO)
    xalt = sl
elsewhere
    xalt = su
end where

gw = inprod(glag, w(1:n))



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
if (curv > -gw .and. curv < -(ONE + sqrt(TWO)) * gw) then
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
!     If IFLAG is 0, then XALT is calculated as before after reversing
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
return
end subroutine geostep


end module geometry_mod
