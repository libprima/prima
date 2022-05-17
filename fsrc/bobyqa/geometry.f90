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
! Last Modified: Tuesday, May 17, 2022 PM03:52:03
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


function geostep(knew, kopt, adelt, bmat, sl, su, xopt, xpt, zmat) result(d)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : matprod, inprod, trueloc
use, non_intrinsic :: powalg_mod, only : hess_mul, calvlag, calbeta

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
real(RP) :: d(size(xpt, 1))  ! D(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: vlagsq, denom_cauchy, beta_cauchy, vlag_cauchy(size(xpt, 1) + size(xpt, 2))
real(RP) :: vlag_line(size(xpt, 1) + size(xpt, 2))
real(RP) :: beta_line
real(RP) :: denom_line
real(RP) :: xcauchy(size(xpt, 1))
real(RP) :: x(size(xpt, 1))
real(RP) :: xline(size(xpt, 1))
real(RP) :: glag(size(xpt, 1))
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: s(size(xpt, 1))
real(RP) :: bigstp, vlagsq_cauchy, curv, dderiv(size(xpt, 2)), distsq(size(xpt, 2)),  &
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
end if

!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XLINE will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XLINE-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XCAUCHY also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XLINE being that XCAUCHY-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XCAUCHY is
!       not calculated if all components of GLAG (see below) are ZERO.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XCAUCHY-XOPT from XOPT for the vector XCAUCHY that is returned,
!       except that CAUCHY is set to ZERO if XCAUCHY is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.


! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC, which is the KNEW-th Lagrange function.
pqlag = matprod(zmat, zmat(knew, :))
alpha = pqlag(knew)

! Calculate the gradient of the KNEW-th Lagrange function at XOPT.
glag = bmat(:, knew) + hess_mul(xopt, xpt, pqlag)
if (any(is_nan(glag)) .or. is_nan(adelt)) then  ! ADELT is not NaN if the input is correct.
    d = ZERO
    return
end if

! Search for a large denominator along the straight lines through XOPT and another interpolation
! point, subject to the bound constraints and the trust region. According to these constraints, SLBD
! and SUBD will be lower and upper bounds on the step along each of these lines in turn. On each
! line, we will evaluate the value of the KNEW-th Lagrange function at 3 trial points, and estimate
! the denominator accordingly. The three points take the form (1-t)*XOPT + t*XPT(:, K) with step
! lengths t = SLBD, SUBD, and STPM, corresponding the upper (U) bound of t, the lower (L) bound of
! t, and a medium (M) step. In total, 3*(NPT-1) trial points will be considered. On the K-th line,
! we intend to maximize the modulus of PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:,K)); overall, we intend
! to find a trial point rendering a large value of the PREDSQ defined in (3.11) of the BOBYQA paper.
!
! We start with the following DO loop, the purpose of which is to define two 3-by-NPT arrays STPLEN
! and ISBD. For each K, STPLEN(1:3, K) and ISBD(1:3, K) corresponds to the straight line through
! XOPT and XPT(:, K). STPLEN(1:3, K) contains SLBD, SUBD, and STPM in this order, which are the step
! lengths for the three trial points on this line. The three entries of SBDI(1:3, K) indicate
! whether the corresponding trial points lie on bounds; SBDI(I, K) = J > 0 means that the I-th trail
! point on the K-th line attains the J-th upper bound, SBDI(I, K) = -J < 0 indicates reaching the
! J-th lower bound, and SBDI(I, K) = 0 means not touching any bound.
dderiv = matprod(glag, xpt) - inprod(glag, xopt) ! The derivatives PHI_K'(0).
distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
do k = 1, npt
    ! It does not make sense to consider "straight line through XOPT and XPT(:, KOPT)". Hence set
    ! STPLEN(:, KOPT) = 0 and ISBD(:, KOPT) = 0 so that VLAG(:, K) and PREDSQ(:, K) obtained after
    ! this loop will be both zero and the search will skip K = KOPT. To avoid undesired/unpredictable
    ! behavior due to possible NaN, set DDERIV(K) = 0 if K = KOPT or if DDERIV(K) is originally NaN.
    if (k == kopt .or. is_nan(dderiv(k))) then
        dderiv(k) = ZERO
        stplen(:, k) = ZERO
        isbd(:, k) = 0_IK
        cycle
    end if

    subd = adelt / sqrt(distsq(k))  ! DISTSQ(K) > 0 unless K == KOPT as long as the input is correct.
    slbd = -subd
    ilbd = 0
    iubd = 0
    sumin = min(ONE, subd)

    ! Revise SLBD and SUBD if necessary because of the bounds in SL and SU according to LFRAC, UFRAC.
    ! N.B.: We calculate LFRAC only at the positions where SL - XOPT > -ABS(XDIFF) * SUBD, because
    ! the values of LFRAC are relevant only at the positions where ABS(LFRAC) < SUBD. Powell's code
    ! does not check this inequality before evaluating LFRAC, and overflow may occur due to large
    ! entries of SL. Note that SL - XOPT > -ABS(XDIFF) * SUBD implies that DIFF /= 0, as long as
    ! SL <= XOPT is ensured. Similar things can be said about UFRAC.
    xdiff = xpt(:, k) - xopt
    lfrac = sign(subd, -xdiff)
    where (sl - xopt > -abs(xdiff) * subd) lfrac = (sl - xopt) / xdiff
    ufrac = sign(subd, xdiff)
    where (su - xopt < abs(xdiff) * subd) ufrac = (su - xopt) / xdiff
    !!MATLAB code for LFRAC and UFRAC (the code is simpler as we are not concerned about overflow):
    !!xdiff = xpt(:, k) - xopt;
    !!lfrac = (sl - xopt) / xdiff;
    !!ufrac = (su - xopt) / xdiff;

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
    ! function PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:,K)) mentioned above. It is a quadratic since
    ! LFUNC is the KNEW-th Lagrange function. For K /= KNEW, the critical point is 0.5, as
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
!!MATLAB: vlag = stplen .* (1 - stplen) .* dderiv; % Implicit expansion; dderiv is a row!!
vlag(:, knew) = stplen(:, knew) * (stplen(:, knew) * (ONE - dderiv(knew)) + dderiv(knew))
! Set NaNs in VLAG to 0 so that the behavior of MAXVAL(ABS(VLAG)) is predictable. VLAG does not have
! NaN unless XPT does, which would be a bug. MAXVAL(ABS(VLAG)) appears in Powell's code, not here.
where (is_nan(vlag)) vlag = ZERO  !!MATLAB: vlag(isnan(vlag)) = 0;
!
! Second, BETABD is the upper bound of BETA given in (3.10) of the BOBYQA paper.
betabd = HALF * (stplen * (ONE - stplen) * spread(distsq, dim=1, ncopies=3))**2
!!MATLAB: betabd = 0.5 * (stplen .* (1-stplen) .* distsq).^2 % Implicit expansion; distsq is a row!!
!
! Finally, PREDSQ is the quantity defined in (3.11) of the BOBYQA paper.
predsq = vlag * vlag * (vlag * vlag + alpha * betabd)
! Set NaNs in PREDSQ to 0 so that the behavior of MAXLOC(PREDSQ) is predictable. PREDSQ does not
! have NaN unless XPT does, which would be a bug.
where (is_nan(predsq)) predsq = ZERO  !!MATLAB: predsq(isnan(predsq)) = 0

! Locate the trial point the renders the maximum of PREDSQ. It is the ISQ-th trial point on the
! straight line through XOPT and XPT(:, KSQ).
! N.B.: 1. The strategy is a bit different from Powell's original code. In Powell's code and the
! BOBYQA paper, we first select the trial point that gives the largest value of ABS(VLAG) on each
! straight line, and then maximize PREDSQ among the (NPT-1) selected points. Here we maximize PREDSQ
! among all the trial points. It works slightly better than Powell's version in a test on 20220428.
! Powell's version is as follows.
!--------------------------------------------------------------!
!isqs = int(maxloc(abs(vlag), dim=1), IK)  ! SIZE(ISQS) = NPT
!ksq = int(maxloc([(predsq(isqs(k), k), k=1, npt)], dim=1), IK)
!isq = isqs(ksq)
!--------------------------------------------------------------!
! 2. Recall that we have set the NaN entries of PREDSQ to zero, if there is any. Thus the KSQS below
! is a well defined integer array, all the three entries lying between 1 and NPT.
ksqs = int(maxloc(predsq, dim=2), IK)
isq = int(maxloc([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))], dim=1), IK)
ksq = ksqs(isq)
!!MATLAB:
!![~, ksqs] = max(predsq, [], 'omitnan');
!![~, isq] = max([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))]);
!!ksq = ksqs(isq);

! Construct XLINE in a way that satisfies the bound constraints exactly.
stpsiz = stplen(isq, ksq)
ibd = isbd(isq, ksq)

xline = max(sl, min(su, xopt + stpsiz * (xpt(:, ksq) - xopt)))
if (ibd < 0) then
    xline(-ibd) = sl(-ibd)
end if
if (ibd > 0) then
    xline(ibd) = su(ibd)
end if

! Calculate DENOM for the current choice of D.
! N.B.: Indeed, only VLAG_LINE(KNEW) is needed. Moreover, VLAG_LINE(KNEW) = VLAG(ISQ, KSQ), which
! can be obtained directly from the array VLAG without any calculation. Therefore, evaluating the
! whole VLAG_LINE from scratch is not really needed here. However, we still choose to do so, because
! VLAG_LINE and BETA_LINE can be calculated in one single function if the language supports allow a
! functions to have multiple output (e.g., MATLAB/Python/Julia/R, but NOT Fortran and C), and the
! calculation of BETA_LINE needs the full VLAG_LINE. For instance, in MATLAB, the calculation will
! look like [vlag_line, beta_line] = vlagbeta(kopt, bmat, d, xpt, zmat). See the code of CALVLAG and
! CALBETA for more information. In other words, our code here is a suboptimal implementation but it
! intends to illustrate what should be done in MATLAB/Python/Julia/R.
d = xline - xopt
vlag_line = calvlag(kopt, bmat, d, xpt, zmat)
beta_line = calbeta(kopt, bmat, d, xpt, zmat)
denom_line = alpha * beta_line + vlag_line(knew)**2

! It works surprisingly well to try the Cauchy step only in the late stage of the algorithm,
! e.g., when ADELT is small. Why? How to find an adaptive rule?
! 1.0E-2 works well if we use DENOM_CAUCHY to decide whether to take the Cauchy step;
! 1.0E-3 works well if we use VLAGSQ instead.
!if (adelt > 1.0E-3) then
if (adelt > 1.0E-2) then
    return
end if

! Prepare for the method that assembles the constrained Cauchy step in S. The sum of squares of the
! fixed components of S is formed in SFIXSQ, and the free components of S are set to BIGSTP. When
! UPHILL = 0, the method calculates the downhill version of XCAUCHY, which intends to minimize the
! KNEW-th Lagrange function; when UPHILL = 1, it calculates the uphill version that intends to
! maximize the Lagrange function.
bigstp = adelt + adelt
xcauchy = xopt
vlagsq_cauchy = ZERO
do uphill = 0, 1
    if (uphill == 1) then
        glag = -glag
    end if
    s = ZERO
    mask_free = (min(xopt - sl, glag) > 0 .or. max(xopt - su, glag) < 0)
    s(trueloc(mask_free)) = bigstp
    ggfree = sum(glag(trueloc(mask_free))**2)
    if (ggfree <= ZERO) then
        cycle
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

    ! Set the remaining free components of S and all components of XCAUCHY. S may be scaled later.
    x(trueloc(glag > 0)) = sl(trueloc(glag > 0))
    x(trueloc(glag <= 0)) = su(trueloc(glag <= 0))
    x(trueloc(s == 0)) = xopt(trueloc(s == 0))
    xtemp = max(sl, min(su, xopt - grdstp * glag))
    x(trueloc(s == bigstp)) = xtemp(trueloc(s == bigstp))
    s(trueloc(s == bigstp)) = -grdstp * glag(trueloc(s == bigstp))
    gs = inprod(glag, s)

    ! Set CURV to the curvature of the KNEW-th Lagrange function along S. Scale S by a factor less
    ! than ONE if that can reduce the modulus of the Lagrange function at XOPT+S. Set CAUCHY to the
    ! final value of the square of this function.
    sxpt = matprod(s, xpt)
    curv = inprod(sxpt, pqlag * sxpt)  ! CURV = INPROD(S, HESS_MUL(S, XPT, PQLAG))
    if (uphill == 1) then
        curv = -curv
    end if
    if (curv > -gs .and. curv < -(ONE + sqrt(TWO)) * gs) then
        scaling = -gs / curv
        x = max(sl, min(su, xopt + scaling * s))
        vlagsq = (HALF * gs * scaling)**2
    else
        vlagsq = (gs + HALF * curv)**2
    end if

    if (vlagsq > vlagsq_cauchy) then
        xcauchy = x
        vlagsq_cauchy = vlagsq
    end if
end do

s = xcauchy - xopt
vlag_cauchy = calvlag(kopt, bmat, s, xpt, zmat)
beta_cauchy = calbeta(kopt, bmat, s, xpt, zmat)
denom_cauchy = alpha * beta_cauchy + vlag_cauchy(knew)**2

! Take the Cauchy step if it is likely to render a larger denominator.
!if (vlagsq_cauchy > max(denom_line, ZERO) .or. is_nan(denom_line)) then  ! Powell's version
if (denom_cauchy > max(denom_line, ZERO) .or. is_nan(denom_line)) then  ! This works slightly better
    d = s
end if

end function geostep


end module geometry_mod
