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
! Last Modified: Friday, April 29, 2022 AM08:09:25
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(knew, kopt, adelt, bmat, sl, su, xopt, xpt, zmat, cauchy, xalt, xnew)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
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
real(RP) :: bigstp, csave, curv, dderiv(size(xpt, 2)), distsq(size(xpt, 2)),  &
&        ggfree, gs, predsq(3, size(xpt, 2)), scaling, &
&        resis, slbd, stplen(3, size(xpt, 2)), grdstp, subd, sumin, betabd(3, size(xpt, 2)), &
&         vlag(3, size(xpt, 2)), sfixsq, ssqsav, xtemp(size(xpt, 1)), sxpt(size(xpt, 2)),  &
&        subd_test(size(xpt, 1)), slbd_test(size(xpt, 1)), &
&        ufrac(size(xpt, 1)), lfrac(size(xpt, 1)), xdiff(size(xpt, 1)), stpm
real(RP) :: alpha, stpsiz
logical :: mask_fixl(size(xpt, 1)), mask_fixu(size(xpt, 1)), mask_free(size(xpt, 1))
integer(IK) :: ibd, uphill, ilbd, isbd(3, size(xpt, 2)), iubd, k, ksqs(3), ksq, isq


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

cauchy = ZERO
xnew = xopt
xalt = xopt
if (any(is_nan(glag))) then
    return
end if


! Search for a large denominator along the straight lines through XOPT and another interpolation
! point. SLBD and SUBD will be lower and upper bounds on the step along each of these lines in turn.
! PREDSQ will be set to the square of the predicted denominator for each line.
!dderiv = matprod(glag, xpt - spread(xopt, dim=2, ncopies=npt))
!distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)

!do k = 1, npt
!    if (k == kopt .or. is_nan(dderiv(k))) then
!        predsq(k) = ZERO
!        stplen(k) = ZERO
!        isbd(k) = 0_IK
!        cycle
!    end if

!    subd = adelt / sqrt(distsq(k))
!    slbd = -subd
!    ilbd = 0
!    iubd = 0
!    sumin = min(ONE, subd)

!    ! Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!    xdiff = xpt(:, k) - xopt
!    where (xdiff /= 0)
!        lfrac = (sl - xopt) / xdiff
!        ufrac = (su - xopt) / xdiff
!    end where

!    slbd_test = slbd
!    slbd_test(trueloc(xdiff > 0)) = lfrac(trueloc(xdiff > 0))
!    slbd_test(trueloc(xdiff < 0)) = ufrac(trueloc(xdiff < 0))
!    if (any(slbd_test > slbd)) then
!        ilbd = maxloc(slbd_test, mask=(.not. is_nan(slbd_test)), dim=1)
!        slbd = slbd_test(ilbd)
!        ilbd = -ilbd * int(sign(ONE, xdiff(ilbd)), IK)
!        !!MATLAB:
!        !![slbd, ilbd] = max(slbd_test, [], 'omitnan');
!        !!ilbd = -ilbd * sign(xdiff(ilbd));
!    end if

!    subd_test = subd
!    subd_test(trueloc(xdiff > 0)) = ufrac(trueloc(xdiff > 0))
!    subd_test(trueloc(xdiff < 0)) = lfrac(trueloc(xdiff < 0))
!    if (any(subd_test < subd)) then
!        iubd = minloc(subd_test, mask=(.not. is_nan(subd_test)), dim=1)
!        subd = max(sumin, subd_test(iubd))
!        iubd = iubd * int(sign(ONE, xdiff(iubd)), IK)
!        !!MATLAB:
!        !![subd, iubd] = min(subd_test, [], 'omitnan');
!        !!subd = max(sumin, subd);
!        !!iubd = iubd * sign(xdiff(iubd));
!    end if

!    !slbds(k) = slbd
!    !ilbds(k) = ilbd
!    !subds(k) = subd
!    !iubds(k) = iubd
!!end do
!!
!!do k = 1, npt
!    !if (k == kopt .or. .not. abs(dderiv(k)) > 0) cycle
!    if (k == knew) then
!        ! Seek a large modulus of the KNEW-th Lagrange function when the index of the other
!        ! interpolation point on the line through XOPT is KNEW.
!        diff = dderiv(k) - ONE
!        vlagl = slbd * (dderiv(k) - slbd * diff)
!        vlagu = subd * (dderiv(k) - subd * diff)
!        !vlagm = stpm * (dderiv(k) - stpm * diff)
!        vlagm = ZERO
!        stpm = slbd
!        if (diff /= 0) then
!            stpm = HALF * dderiv(k) / diff
!        end if
!        if (stpm > slbd .and. stpm < subd) then
!            vlagm = (HALF * dderiv(k)) * (HALF * dderiv(k)) / diff
!        end if
!    else
!        ! Search along each of the other lines through XOPT and another point.
!        vlagl = slbd * (ONE - slbd) * dderiv(k)
!        vlagu = subd * (ONE - subd) * dderiv(k)
!        vlagm = ZERO
!        stpm = HALF
!        if (stpm > slbd .and. stpm < subd) then
!            vlagm = QUART * dderiv(k)
!        end if
!    end if

!    vlags = [vlagl, vlagu, vlagm]
!    stplens = [slbd, subd, stpm]
!    isbds = [ilbd, iubd, 0_IK]
!    !ilag = int(maxloc(abs(vlags), dim=1), IK)
!    !vlag = vlags(ilag)
!    !stplen(k) = stplens(ilag)
!    !isbd(k) = isbds(ilag)

!    ! Calculate PREDSQ for the current line search.
!    !distprod = stplen(k) * (ONE - stplen(k)) * distsq(k)
!    !predsq(k) = vlag * vlag * (vlag * vlag + HALF * alpha * distprod * distprod)

!    distprodl = stplens(1) * (ONE - stplens(1)) * distsq(k)
!    distprodu = stplens(2) * (ONE - stplens(2)) * distsq(k)
!    distprodm = stplens(3) * (ONE - stplens(3)) * distsq(k)
!    predsql = vlagl * vlagl * (vlagl * vlagl + HALF * alpha * distprodl * distprodl)
!    predsqu = vlagu * vlagu * (vlagu * vlagu + HALF * alpha * distprodu * distprodu)
!    predsqm = vlagm * vlagm * (vlagm * vlagm + HALF * alpha * distprodm * distprodm)

!    predsqs = [predsql, predsqu, predsqm]
!    ilag = int(maxloc(predsqs, dim=1), IK)
!    predsq(k) = predsqs(ilag)
!    stplen(k) = stplens(ilag)
!    isbd(k) = isbds(ilag)
!end do
!kden = kopt
!if (any(predsq > 0)) then
!    kden = maxloc(predsq, mask=(.not. is_nan(predsq)), dim=1)
!end if

!! Construct XNEW in a way that satisfies the bound constraints exactly.
!xnew = max(sl, min(su, xopt + stplen(kden) * (xpt(:, kden) - xopt)))
!ibd = isbd(kden)
!if (ibd < 0) then
!    xnew(-ibd) = sl(-ibd)
!end if
!if (ibd > 0) then
!    xnew(ibd) = su(ibd)
!end if

! Search for a large denominator along the straight lines through XOPT and another interpolation
! point. SLBD and SUBD will be lower and upper bounds on the step along each of these lines in turn.
! On each line, we will evaluate the value of the KNEW-th Lagrange function at 3 trial points, and
! estimate the denominator accordingly. The three points take the form (1-t)*XOPT + t*XPT(:, K)
! with step lengths t = SLBD, SUBD, and STPM. In total, 3*(NPT-1) trial points will be considered.
! On the K-th line, we intend to maximize the modulus of PHI_K(t) = LFUN((1-t)*XOPT + t*XPT(:,K)),
! where LFUN is the KNEW-th Lagrange function; overall, we attempt to find a trial point that
! renders a large value of the quantity PREDSQ defined in (3.11) of the BOBYQA paper.
!
! We start with the following DO loop, the purpose of which is to define two 3-by-NPT arrays STPLEN
! and ISBD. For each K, STPLEN(1:3, K) and ISBD(1:3, K) corresponds to the straight line through
! XOPT and XPT(:, K). STPLEN(1:3, K) contains SLBD. SUBD, and STPM in this order, which are the step
! lengths for the three trial points on the K-th line. The three entries of SBDI(1:3, K) indicate
! whether the corresponding trial points lie on bounds; SBDI(I, K) = J > 0 means that the I-th trail
! point on the K-th line attainins the J-th upper bound, SBDI(I, K) = -J indicates reaching the J-th
! lower bound, and SBDI(I, K) = 0 means not touching any bound.
dderiv = matprod(glag, xpt - spread(xopt, dim=2, ncopies=npt))  ! It contains derivatives PHI_K'(0).
distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
do k = 1, npt

    if (k == kopt .or. is_nan(dderiv(k))) then
        stplen(:, k) = ZERO
        isbd(:, k) = 0_IK
        cycle
    end if

    subd = adelt / sqrt(distsq(k))
    slbd = -subd
    ilbd = 0
    iubd = 0
    sumin = min(ONE, subd)

    ! Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
    xdiff = xpt(:, k) - xopt
    where (xdiff /= 0)
        lfrac = (sl - xopt) / xdiff
        ufrac = (su - xopt) / xdiff
    end where

    ! First, revise SLBD.
    slbd_test = slbd
    slbd_test(trueloc(xdiff > 0)) = lfrac(trueloc(xdiff > 0))
    slbd_test(trueloc(xdiff < 0)) = ufrac(trueloc(xdiff < 0))
    if (any(slbd_test > slbd)) then
        ilbd = maxloc(slbd_test, mask=(.not. is_nan(slbd_test)), dim=1)
        slbd = slbd_test(ilbd)
        ilbd = -ilbd * int(sign(ONE, xdiff(ilbd)), IK)
        !!MATLAB:
        !![slbd, ilbd] = max(slbd_test, [], 'omitnan');
        !!ilbd = -ilbd * sign(xdiff(ilbd));
    end if

    ! Second, revise SUBD.
    subd_test = subd
    subd_test(trueloc(xdiff > 0)) = ufrac(trueloc(xdiff > 0))
    subd_test(trueloc(xdiff < 0)) = lfrac(trueloc(xdiff < 0))
    if (any(subd_test < subd)) then
        iubd = minloc(subd_test, mask=(.not. is_nan(subd_test)), dim=1)
        subd = max(sumin, subd_test(iubd))
        iubd = iubd * int(sign(ONE, xdiff(iubd)), IK)
        !!MATLAB:
        !![subd, iubd] = min(subd_test, [], 'omitnan');
        !!subd = max(sumin, subd);
        !!iubd = iubd * sign(xdiff(iubd));
    end if

    ! Now, define the step length STPM between SLBD and SUBD by finding the critical point of the
    ! function PHI_K(t) = LFUN((1-t)*XOPT + t*XPT(:,K)) mentioned above. It is a quadratic since
    ! LFUN is the KNEW-th Lagrange function. For K /= KNEW, the critical point is 0.5, as
    ! PHI_K(0) = 1 = PHI_K(1); when K = KNEW, it is -0.5*PHI_K'(0) / (1 - PHI_K'(0)), because
    ! PHI_K(0) = 0 and PHI_K(1) = 1.
    stpm = HALF
    if (k == knew) then
        stpm = slbd
        if (ONE - dderiv(k) /= 0) then
            stpm = -HALF * dderiv(k) / (ONE - dderiv(k))
        end if
    end if
    stpm = min(max(slbd, stpm), subd)

    stplen(:, k) = [slbd, subd, stpm]
    isbd(:, k) = [ilbd, iubd, 0_IK]
end do

! The following lines calculate PREDSQ for all the 3*(NPT-1) trial points.
! First, compute VLAG = PHI(STPLEN). Using the fact that PHI_K(0) = 0, PHI_K(1) = delta_{K, KNEW}
! (Kroneker delta), and recalling the PHI_K is quadratic, we can find that
! PHI_K(t) = t*(1-t)*PHI_K'(0) for K /= KNEW, and PHI_KNEW = t*[t*(1-PHI_K'(0)) + PHI_K'(0)].
vlag = stplen * (ONE - stplen) * spread(dderiv, dim=1, ncopies=3)
vlag(:, knew) = stplen(:, knew) * (stplen(:, knew) * (ONE - dderiv(knew)) + dderiv(knew))
! Second, BETABD is the upper bound of BETA given in (3.10) of the BOBYQA paper.
betabd = HALF * (stplen * (ONE - stplen) * spread(distsq, dim=1, ncopies=3))**2
! Finally, PREDSQ is the quantity defined in (3.11) of the BOBYQA paper.
predsq = vlag * vlag * (vlag * vlag + alpha * betabd)
predsq(:, kopt) = ZERO

! Locate the trial point the renders the maximum of PREDSQ. It is the ISQ-th trial point on the
! straight line through XOPT and XPT(:, KSQ).
! N.B.: The strategy is a bit different from Powell's original code. In Powell's code and the BOBYQA
! paper, we first select the trial point that gives the largest value of ABS(VLAG) on each straight,
! line, and then maximize PREDSQ among the (NPT-1) selected points. Here we maximize PREDSQ among
! all the trial points. It works slightly better than Powell's version in a test on 20220428.
ksqs = int(maxloc(predsq, mask=(.not. is_nan(predsq)), dim=2), IK)
isq = int(maxloc([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))], dim=1), IK)
ksq = ksqs(isq)

! Construct XNEW in a way that satisfies the bound constraints exactly.
stpsiz = stplen(isq, ksq)
ibd = isbd(isq, ksq)
xnew = max(sl, min(su, xopt + stpsiz * (xpt(:, ksq) - xopt)))
if (ibd < 0) then
    xnew(-ibd) = sl(-ibd)
end if
if (ibd > 0) then
    xnew(ibd) = su(ibd)
end if

! Prepare for the method that assembles the constrained Cauchy step in S. The sum of squares of the
! fixed components of S is formed in SFIXSQ, and the free components of S are set to BIGSTP. When
! UPHILL = 0, the method calculates the downhill version of XALT, which intends to minimize the
! KNEW-th Lagrange function; when UPHILL = 1, it calculates the uphill version that intends to
! maximize the Lagrange function.
bigstp = adelt + adelt
do uphill = 0, 1
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
        grdstp = sqrt(resis / ggfree)
        xtemp = xopt - grdstp * glag
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
    xtemp = max(sl, min(su, xopt - grdstp * glag))
    xalt(trueloc(s == bigstp)) = xtemp(trueloc(s == bigstp))
    s(trueloc(s == bigstp)) = -grdstp * glag(trueloc(s == bigstp))
    gs = inprod(glag, s)

    ! Set CURV to the curvature of the KNEW-th Lagrange function along S. Scale S by a factor less
    ! than ONE if that can reduce the modulus of the Lagrange function at XOPT+S. Set CAUCHY to the
    ! final value of the square of this function.
    sxpt = matprod(s, xpt)
    curv = inprod(sxpt, hcol * sxpt)
    if (uphill == 1) then
        curv = -curv
    end if
    if (curv > -gs .and. curv < -(ONE + sqrt(TWO)) * gs) then
        scaling = -gs / curv
        xalt = max(sl, min(su, xopt + scaling * s))
        cauchy = (HALF * gs * scaling)**2
    else
        cauchy = (gs + HALF * curv)**2
    end if

    ! If UPHILL is 0, then XALT is calculated as before after reversing the sign of GLAG. Thus two
    ! XALT vectors become available. The one that is chosen is the one that gives the larger value
    ! of CAUCHY.
    if (uphill == 0) then
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
