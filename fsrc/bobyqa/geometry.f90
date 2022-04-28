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
! Last Modified: Thursday, April 28, 2022 PM12:18:41
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(knew, kopt, adelt, bmat, sl, su, xopt, xpt, zmat, alpha, cauchy, xalt, xnew)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
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
real(RP) :: s(size(xpt, 1)), xsav(size(xpt, 1))
real(RP) :: bigstp, csave, curv, dderiv(size(xpt, 2)), diff, distsq(size(xpt, 2)),  &
&        ggfree, gs, predsq, presav, scaling, &
&        resis, slbd, stplen, stpsav, subd, sumin, temp, tempa,      &
&        tempb, tempd, vlag, sfixsq, ssqsav, xtemp(size(xpt, 1)), sxpt(size(xpt, 2)),  &
&        subd_test(size(xpt, 1)), slbd_test(size(xpt, 1)), &
&        ubd_test(size(xpt, 1)), lbd_test(size(xpt, 1)), xdiff(size(xpt, 1))
logical :: mask_fixl(size(xpt, 1)), mask_fixu(size(xpt, 1)), mask_free(size(xpt, 1))
integer(IK) :: i, ibdsav, iflag, ilbd, isbd, iubd, k, ksav


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

hcol = matprod(zmat, zmat(knew, :))
alpha = hcol(knew)

!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
!!glag = bmat(:, knew) + hess_mul(xopt, xpt, hcol)
glag = bmat(:, knew)
do k = 1, npt
    glag = glag + hcol(k) * inprod(xopt, xpt(:, k)) * xpt(:, k)
end do

if (any(is_nan(glag))) then
    alpha = ZERO
    cauchy = ZERO
    xnew = xopt
    xalt = xopt
    return
end if

!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
presav = ZERO
dderiv = matprod(glag, xpt - spread(xopt, dim=2, ncopies=npt))
distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
do k = 1, npt
    if (k == kopt) cycle
    subd = adelt / sqrt(distsq(k))
    slbd = -subd
    ilbd = 0
    iubd = 0
    sumin = min(ONE, subd)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
    !do i = 1, n
    !    temp = xpt(i, k) - xopt(i)
    !    if (temp > ZERO) then
    !        if (slbd < (sl(i) - xopt(i)) / temp) then
    !            slbd = (sl(i) - xopt(i)) / temp
    !            ilbd = -i
    !        end if
    !        if (subd > (su(i) - xopt(i)) / temp) then
    !            subd = max(sumin, (su(i) - xopt(i)) / temp)
    !            iubd = i
    !        end if
    !    else if (temp < ZERO) then
    !        if (slbd < (su(i) - xopt(i)) / temp) then
    !            slbd = (su(i) - xopt(i)) / temp
    !            ilbd = i
    !        end if
    !        if (subd > (sl(i) - xopt(i)) / temp) then
    !            subd = max(sumin, (sl(i) - xopt(i)) / temp)
    !            iubd = -i
    !        end if
    !    end if
    !end do

    xdiff = xpt(:, k) - xopt
    where (xdiff /= 0)
        lbd_test = (sl - xopt) / xdiff
        ubd_test = (su - xopt) / xdiff
    end where

    slbd_test = slbd
    slbd_test(trueloc(xdiff > 0)) = lbd_test(trueloc(xdiff > 0))
    slbd_test(trueloc(xdiff < 0)) = ubd_test(trueloc(xdiff < 0))
    if (any(slbd_test > slbd)) then
        ilbd = maxloc(slbd_test, mask=(.not. is_nan(slbd_test)), dim=1)
        slbd = slbd_test(ilbd)
        ilbd = -ilbd * int(sign(ONE, xdiff(ilbd)), IK)
    end if

    subd_test = subd
    subd_test(trueloc(xdiff > 0)) = ubd_test(trueloc(xdiff > 0))
    subd_test(trueloc(xdiff < 0)) = lbd_test(trueloc(xdiff < 0))
    if (any(subd_test < subd)) then
        iubd = minloc(subd_test, mask=(.not. is_nan(subd_test)), dim=1)
        subd = max(sumin, subd_test(iubd))
        iubd = iubd * int(sign(ONE, xdiff(iubd)), IK)
    end if
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
    if (k == knew) then
        diff = dderiv(k) - ONE
        stplen = slbd
        vlag = slbd * (dderiv(k) - slbd * diff)
        isbd = ilbd
        temp = subd * (dderiv(k) - subd * diff)
        if (abs(temp) > abs(vlag)) then
            stplen = subd
            vlag = temp
            isbd = iubd
        end if
        tempd = HALF * dderiv(k)
        tempa = tempd - diff * slbd
        tempb = tempd - diff * subd
        if (tempa * tempb < ZERO) then
            temp = tempd * tempd / diff
            if (abs(temp) > abs(vlag)) then
                stplen = tempd / diff
                vlag = temp
                isbd = 0
            end if
        end if
!
!     Search along each of the other lines through XOPT and another point.
!
    else
        stplen = slbd
        vlag = slbd * (ONE - slbd)
        isbd = ilbd
        temp = subd * (ONE - subd)
        if (abs(temp) > abs(vlag)) then
            stplen = subd
            vlag = temp
            isbd = iubd
        end if
        if (subd > HALF) then
            if (abs(vlag) < QUART) then
                stplen = HALF
                vlag = QUART
                isbd = 0
            end if
        end if
        vlag = vlag * dderiv(k)
    end if
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
    temp = stplen * (ONE - stplen) * distsq(k)
    predsq = vlag * vlag * (vlag * vlag + HALF * alpha * temp * temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: With the original code, if either PREDSQ or PRESAV
! is NaN, KSAV/STPSAV/IBDSAV will not get a value. This may cause
! Segmentation Fault.
!      IF (PREDSQ .GT. PRESAV) THEN
    if (.not. (predsq <= presav)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        presav = predsq
        ksav = k
        stpsav = stplen
        ibdsav = isbd
    end if
end do


! Construct XNEW in a way that satisfies the bound constraints exactly.
xnew = max(sl, min(su, xopt + stpsav * (xpt(:, ksav) - xopt)))
if (ibdsav < 0) then
    xnew(-ibdsav) = sl(-ibdsav)
end if
if (ibdsav > 0) then
    xnew(ibdsav) = su(ibdsav)
end if

! Prepare for the method that assembles the constrained Cauchy step in S. The sum of squares of the
! fixed components of S is formed in SFIXSQ, and the free components of S are set to BIGSTP.
! When IFLAG = 0, the method calculates the downhill version of XALT, which intends to minimize the
! KNEW-th Lagrange function; when IFLAG = 1, it calculates the uphill version that intends to
! maximize the Lagrange function.
bigstp = adelt + adelt
do iflag = 0, 1
    s = ZERO
    mask_free = (min(xopt - sl, glag) > 0 .or. max(xopt - su, glag) < 0)
    s(trueloc(mask_free)) = bigstp
    ggfree = sum(glag(trueloc(mask_free))**2)
    if (ggfree <= ZERO) then
        cauchy = ZERO
        return
    end if

    ! Investigate whether more components of S can be fixed. Note that the loop counter K does not
    ! appear in the loop body. The purpose of K is only to impose an explicit bound on the number of
    ! loops. Powell's code does not have such a bound. The bound is not a true restriction, because
    ! we can check that (SFIXSQ > SSQSAV .AND. GGFREE > ZERO) must fail within N loops.
    sfixsq = ZERO
    do k = 1, n
        resis = adelt**2 - sfixsq
        if (resis <= 0) exit
        ssqsav = sfixsq
        stplen = sqrt(resis / ggfree)
        xtemp = xopt - stplen * glag
        mask_fixl = (s == bigstp .and. xtemp <= sl)
        mask_fixu = (s == bigstp .and. xtemp >= su)
        mask_free = (s == bigstp .and. .not. (mask_fixl .or. mask_fixu))
        s(trueloc(mask_fixl)) = sl(trueloc(mask_fixl)) - xopt(trueloc(mask_fixl))
        s(trueloc(mask_fixu)) = su(trueloc(mask_fixu)) - xopt(trueloc(mask_fixu))
        sfixsq = sfixsq + sum(s(trueloc(mask_fixl .or. mask_fixu))**2)
        ggfree = sum(glag(trueloc(mask_free))**2)
        if (.not. (sfixsq > ssqsav .and. ggfree > ZERO)) exit
    end do

    ! Set the remaining free components of S and all components of XALT, except that S may be
    ! scaled later.
    xalt(trueloc(glag > 0)) = sl(trueloc(glag > 0))
    xalt(trueloc(glag <= 0)) = su(trueloc(glag <= 0))
    xalt(trueloc(s == 0)) = xopt(trueloc(s == 0))
    xtemp = max(sl, min(su, xopt - stplen * glag))
    xalt(trueloc(s == bigstp)) = xtemp(trueloc(s == bigstp))
    s(trueloc(s == bigstp)) = -stplen * glag(trueloc(s == bigstp))
    gs = inprod(glag, s)

    ! Set CURV to the curvature of the KNEW-th Lagrange function along S. Scale S by a factor less
    ! than ONE if that can reduce the modulus of the Lagrange function at XOPT+S. Set CAUCHY to the
    ! final value of the square of this function.
    sxpt = matprod(s, xpt)
    curv = inprod(sxpt, hcol * sxpt)
    if (iflag == 1) then
        curv = -curv
    end if
    if (curv > -gs .and. curv < -(ONE + sqrt(TWO)) * gs) then
        scaling = -gs / curv
        xalt = max(sl, min(su, xopt + scaling * s))
        cauchy = (HALF * gs * scaling)**2
    else
        cauchy = (gs + HALF * curv)**2
    end if

    ! If IFLAG is 0, then XALT is calculated as before after reversing the sign of GLAG. Thus two
    ! XALT vectors become available. The one that is chosen is the one that gives the larger value
    ! of CAUCHY.
    if (iflag == 0) then
        glag = -glag
        xsav = xalt
        csave = cauchy
    end if
end do

if (csave > cauchy) then
    xalt = xsav
    cauchy = csave
end if
end subroutine geostep


end module geometry_mod
