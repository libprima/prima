module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, March 06, 2023 PM05:20:07
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: setdrop_tr, geostep


contains


function setdrop_tr(kopt, ximproved, bmat, d, delta, rho, xpt, zmat) result(knew)
!--------------------------------------------------------------------------------------------------!
! This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
! REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
! XPT(:, KNEW) is replaced with XNEW = XOPT + D, where D is the trust-region step. See discussions
! around (6.1) of the BOBYQA paper.
! N.B.:
! 1. If XIMPROVED = TRUE, then KNEW > 0 so that XNEW is included into XPT. Otherwise, it is a bug.
! 2. If XIMPROVED = FALSE, then KNEW /= KOPT so that XPT(:, KOPT) stays. Otherwise, it is a bug.
! 3. It is tempting to take the function value into consideration when defining KNEW, for example,
! set KNEW so that FVAL(KNEW) = MAX(FVAL) as long as F(XNEW) < MAX(FVAL), unless there is a better
! choice. However, this is not a good idea, because the definition of KNEW should benefit the
! quality of the model that interpolates f at XPT. A set of points with low function values is not
! necessarily a good interpolation set. In contrast, a good interpolation set needs to include
! points with relatively high function values; otherwise, the interpolant will unlikely reflect the
! landscape of the function sufficiently.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: HDIAG(NPT), DENABS(NPT), SCORE(NPT) VLAG(N+NPT), XDIST(NPT)
! Size of local arrays: REAL(RP)*(4*NPT+N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : issymmetric
use, non_intrinsic :: powalg_mod, only : calden

implicit none

! Inputs
integer(IK), intent(in) :: kopt
logical, intent(in) :: ximproved
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)        ! D(N)
real(RP), intent(in) :: delta
real(RP), intent(in) :: rho
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
integer(IK) :: knew

! Local variables
character(len=*), parameter :: srname = 'SETDROP_TR'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: den(size(xpt, 2))
real(RP) :: distsq(size(xpt, 2))
real(RP) :: score(size(xpt, 2))
real(RP) :: weight(size(xpt, 2))

! Sizes
n = int(size(xpt, 1), kind(npt))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(delta >= rho .and. rho > 0, 'DELTA >= RHO > 0', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the distance squares between the interpolation points and the "optimal point". When
! identifying the optimal point, it is reasonable to take into account the new trust-region trial
! point XPT(:, KOPT) + D, which will become the optimal point in the next iteration if XIMPROVED
! is TRUE. Powell suggested this in
! - (56) of the UOBYQA paper, lines 276--297 of uobyqb.f,
! - (7.5) and Box 5 of the NEWUOA paper, lines 383--409 of newuob.f,
! - the last paragraph of page 26 of the BOBYQA paper, lines 435--465 of bobyqb.f.
! However, Powell's LINCOA code is different. In his code, the KNEW after a trust-region step is
! picked in lines 72--96 of the update.f for LINCOA, where DISTSQ is calculated as the square of the
! distance to XPT(KOPT, :) (Powell recorded the interpolation points in rows). However, note that
! the trust-region trial point has not been included in to XPT yet --- it can not be included
! without knowing KNEW (see lines 332-344 and 404--431 of lincob.f). Hence Powell's LINCOA code
! picks KNEW based on the distance to the un-updated "optimal point", which is unreasonable.
! This has been corrected in our implementation of LINCOA, yet it does not boost the performance.
if (ximproved) then
    distsq = sum((xpt - spread(xpt(:, kopt) + d, dim=2, ncopies=npt))**2, dim=1)
else
    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
end if

weight = max(ONE, distsq / rho**2)**4
! Other possible definitions of WEIGHT.
! !weight = max(ONE, distsq / rho**2)**3.5  ! Quite similar to power 4
! !weight = max(ONE, distsq / rho**2)**3  ! Not bad
! !weight = max(ONE, distsq / delta**2)**2  ! Powell's code. Does not works as well as the above.
! !weight = max(ONE, distsq / rho**2)**2  ! Similar to Powell's code, not better.
! !weight = max(ONE, distsq / delta**2)  ! Defined in (6.1) of the BOBYQA paper. It works poorly!
! !weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3.5  ! The same as DISTSQ/RHO**2.
! The following WEIGHT all perform a bit worse than the above one.
! !weight = max(ONE, distsq / delta**2)**3.5
! !weight = max(ONE, distsq / delta**2)**2.5
! !weight = max(ONE, distsq / delta**2)**3
! !weight = max(ONE, distsq / delta**2)**4
! !weight = max(ONE, distsq / delta**2)**4.5
! !weight = max(ONE, distsq / rho**2)**2.5
! !weight = max(ONE, distsq / rho**2)**3
! !weight = max(ONE, distsq / rho**2)**4.5

! Different from NEWUOA/LINCOA, the possibility that entries in DEN become negative is handled by
! RESCUE. Hence the SCORE here uses DEN in contrast to ABS(DEN) in NEWUOA/LINCOA.
den = calden(kopt, bmat, d, xpt, zmat)
score = weight * den

! If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
if (.not. ximproved) then
    score(kopt) = -ONE
end if

! For the first case below, NEWUOA checks ANY(SCORE>1) .OR. (XIMPROVED .AND. ANY(SCORE>0))
! instead of ANY(SCORE > 0). This seems to improve the performance of BOBYQA very slightly.
if (any(score > 1) .or. (ximproved .and. any(score > 0))) then  ! Condition in NEWUOA.
    ! !if (any(score > 0)) then  ! Powell's original condition in BOBYQA.
    ! See (6.1) of the BOBYQA paper for the definition of KNEW in this case.
    ! SCORE(K) = NaN implies DEN(K) = NaN. We exclude such K as we want DEN to be big.
    knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), kind(knew))
    !!MATLAB: [~, knew] = max(score, [], 'omitnan');
elseif (ximproved) then
    ! Powell's code does not include the following instructions. With Powell's code, if DEN
    ! consists of only NaN, then KNEW can be 0 even when XIMPROVED is TRUE.
    knew = int(maxloc(distsq, dim=1), kind(knew))
else
    knew = 0  ! We arrive here when XIMPROVED = FALSE and no entry of SCORE is positive.
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(knew /= kopt .or. ximproved, 'KNEW /= KOPT unless XIMPROVED = TRUE', srname)
    call assert(knew >= 1 .or. .not. ximproved, 'KNEW >= 1 unless XIMPROVED = FALSE', srname)
    ! KNEW >= 1 when XIMPROVED = TRUE unless NaN occurs in DISTSQ, which should not happen if the
    ! starting point does not contain NaN and the trust-region/geometry steps never contain NaN.
end if

end function setdrop_tr


function geostep(knew, kopt, bmat, delbar, sl, su, xpt, zmat) result(d)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds a step D that intends to improve the geometry of the interpolation set
! when XPT(:, KNEW) is changed to XOPT + D, where XOPT = XPT(:, KOPT). See Section 3 of the BOBYQA
! paper, particularly the discussions starting from (3.7).
!
! The arguments XPT, BMAT, ZMAT, SL and SU all have the same meanings as in BOBYQB.
! KOPT is the index of the optimal interpolation point.
! KNEW is the index of the interpolation point that is going to be moved.
! DELBAR is the trust region bound for the geometry step.
! XLINE will be a suitable new position for the interpolation point XPT(:, KNEW). Specifically, it
!   satisfies the SL, SU and trust region bounds and it should provide a large denominator in the
!   next call of UPDATE. The step XLINE-XOPT from XOPT is restricted to moves along the straight
!   lines through XOPT and another interpolation point.
! XCAUCHY also provides a large value of the modulus of the KNEW-th Lagrange function subject to
!   the constraints that have been mentioned, its main difference from XLINE being that XCAUCHY-XOPT
!   is a bound-constrained version of the Cauchy step within the trust region.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, trueloc, norm, issymmetric
use, non_intrinsic :: powalg_mod, only : hess_mul, calden

implicit none

! Inputs
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: delbar
real(RP), intent(in) :: sl(:)  ! SL(N)
real(RP), intent(in) :: su(:)  ! SU(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
real(RP) :: d(size(xpt, 1))  ! D(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: ibd
integer(IK) :: ilbd
integer(IK) :: isbd(3, size(xpt, 2))
integer(IK) :: isq
integer(IK) :: iubd
integer(IK) :: k
integer(IK) :: ksq
integer(IK) :: ksqs(3)
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: uphill
logical :: mask_fixl(size(xpt, 1))
logical :: mask_fixu(size(xpt, 1))
logical :: mask_free(size(xpt, 1))
real(RP) :: alpha, stpsiz
real(RP) :: betabd(3, size(xpt, 2))
real(RP) :: bigstp
real(RP) :: curv
real(RP) :: dderiv(size(xpt, 2))
real(RP) :: den_cauchy(size(xpt, 2))
real(RP) :: den_line(size(xpt, 2))
real(RP) :: distsq(size(xpt, 2))
real(RP) :: ggfree
real(RP) :: glag(size(xpt, 1))
real(RP) :: grdstp
real(RP) :: gs
real(RP) :: lfrac(size(xpt, 1))
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: predsq(3, size(xpt, 2))
real(RP) :: resis
real(RP) :: s(size(xpt, 1))
real(RP) :: scaling
real(RP) :: sfixsq
real(RP) :: slbd
real(RP) :: slbd_test(size(xpt, 1))
real(RP) :: ssqsav
real(RP) :: stplen(3, size(xpt, 2))
real(RP) :: stpm
real(RP) :: subd
real(RP) :: subd_test(size(xpt, 1))
real(RP) :: sumin
real(RP) :: sxpt(size(xpt, 2))
real(RP) :: ufrac(size(xpt, 1))
real(RP) :: vlag(3, size(xpt, 2))
real(RP) :: vlagsq
real(RP) :: vlagsq_cauchy
real(RP) :: x(size(xpt, 1))
real(RP) :: xcauchy(size(xpt, 1))
real(RP) :: xdiff(size(xpt, 1))
real(RP) :: xline(size(xpt, 1))
real(RP) :: xopt(size(xpt, 1))
real(RP) :: xtemp(size(xpt, 1))

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
    call assert(delbar > 0, 'DELBAR > 0', srname)
    call assert(size(sl) == n .and. all(sl <= 0), 'SIZE(SL) == N, SL <= 0', srname)
    call assert(size(su) == n .and. all(su >= 0), 'SIZE(SU) == N, SU >= 0', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt)) .and. &
        & all(xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC, which is the KNEW-th Lagrange function. ALPHA will is the KNEW-th
! diagonal element of the H matrix.
pqlag = matprod(zmat, zmat(knew, :))
alpha = pqlag(knew)

! Read XOPT.
xopt = xpt(:, kopt)

! Calculate the gradient GLAG of the KNEW-th Lagrange function at XOPT.
glag = bmat(:, knew) + hess_mul(xopt, xpt, pqlag)

! In case GLAG contains NaN, set D to a displacement from XOPT to XPT(:, KNEW) and return. Powell's
! code does not have this, and D may be NaN in the end. Note that it is crucial to ensure that a
! geometry step is nonzero.
if (is_nan(sum(abs(glag)))) then
    d = xpt(:, knew) - xopt
    d = min(HALF, delbar / norm(d)) * d  ! Since XPT respects the bounds, so does XOPT + D.
    return
end if

! Search for a large denominator along the straight lines through XOPT and another interpolation
! point, subject to the bound constraints and the trust region. According to these constraints, SLBD
! and SUBD will be lower and upper bounds on the step along each of these lines in turn. On each
! line, we will evaluate the value of the KNEW-th Lagrange function at 3 trial points, and estimate
! the denominator accordingly. The three points take the form (1-t)*XOPT + t*XPT(:, K) with step
! lengths t = SLBD, SUBD, and STPM, corresponding to the upper (U) bound of t, the lower (L) bound
! of t, and a medium (M) step. In total, 3*(NPT-1) trial points will be considered. On the K-th line,
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
        isbd(:, k) = 0
        cycle
    end if

    subd = delbar / sqrt(distsq(k))  ! DISTSQ(K) > 0 unless K == KOPT or the input is incorrect.
    slbd = -subd
    ilbd = 0
    iubd = 0
    sumin = min(ONE, subd)

    ! Revise SLBD and SUBD if necessary because of the bounds in SL and SU according to LFRAC, UFRAC.
    ! N.B.: We calculate LFRAC only at the positions where SL - XOPT > -ABS(XDIFF) * SUBD, because
    ! the values of LFRAC are relevant only at the positions where ABS(LFRAC) < SUBD. Powell's code
    ! does not check this inequality before evaluating LFRAC, and overflow may occur due to large
    ! entries of SL. Note that SL - XOPT > -ABS(XDIFF) * SUBD implies that XDIFF /= 0, as long as
    ! SL <= XOPT is ensured. In addition, when initializing LFRAC to SIGN(SUBD, -XDIFF), we do not
    ! need to worry about the case where XDIFF = 0, because we only use LFRAC when XDIFF /= 0.
    ! Similar things can be said about UFRAC.
    xdiff = xpt(:, k) - xopt
    lfrac = sign(subd, -xdiff)
    where (sl - xopt > -abs(xdiff) * subd) lfrac = (sl - xopt) / xdiff
    ufrac = sign(subd, xdiff)
    where (su - xopt < abs(xdiff) * subd) ufrac = (su - xopt) / xdiff
    !!MATLAB code for LFRAC and UFRAC (the code is simpler as we are not concerned about overflow):
    !!xdiff = xpt(:, k) - xopt;
    !!lfrac = (sl - xopt) / xdiff;
    !!ufrac = (su - xopt) / xdiff;

    ! First, revise SLBD. Note that SLBD_TEST <= 0 unless the input violates XOPT >= SL.
    slbd_test = slbd
    slbd_test(trueloc(xdiff > 0)) = lfrac(trueloc(xdiff > 0))
    slbd_test(trueloc(xdiff < 0)) = ufrac(trueloc(xdiff < 0))
    if (any(slbd_test > slbd)) then
        ilbd = int(maxloc(slbd_test, mask=(.not. is_nan(slbd_test)), dim=1), kind(ilbd))
        slbd = slbd_test(ilbd)
        ilbd = -ilbd * nint(sign(ONE, xdiff(ilbd)), kind(ilbd))
        !!MATLAB:
        !![slbd, ilbd] = max(slbd_test, [], 'omitnan');
        !!ilbd = -ilbd * sign(xdiff(ilbd));
    end if

    ! Second, revise SUBD. Note that SUBD_TEST >= 0 unless the input violates XOPT <= SU.
    subd_test = subd
    subd_test(trueloc(xdiff > 0)) = ufrac(trueloc(xdiff > 0))
    subd_test(trueloc(xdiff < 0)) = lfrac(trueloc(xdiff < 0))
    if (any(subd_test < subd)) then
        iubd = int(minloc(subd_test, mask=(.not. is_nan(subd_test)), dim=1), kind(iubd))
        subd = max(sumin, subd_test(iubd))
        iubd = iubd * nint(sign(ONE, xdiff(iubd)), kind(iubd))
        !!MATLAB:
        !![subd, iubd] = min(subd_test, [], 'omitnan');
        !!subd = max(sumin, subd);
        !!iubd = iubd * sign(xdiff(iubd));
    end if

    if (DEBUGGING) then
        call assert(slbd <= 0 .and. subd >= 0, 'SLBD <= 0 <= SUBD', srname)
    end if

    ! Now, define the step length STPM between SLBD and SUBD by finding the critical point of the
    ! function PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:,K)) mentioned above. It is a quadratic since
    ! LFUNC is the KNEW-th Lagrange function. For K /= KNEW, the critical point is 0.5, as
    ! PHI_K(0) = 1 = PHI_K(1); when K = KNEW, it is -0.5*PHI_K'(0) / (1 - PHI_K'(0)), because
    ! PHI_K(0) = 0 and PHI_K(1) = 1.
    stpm = HALF
    if (k == knew) then
        stpm = slbd
        if (abs(ONE - dderiv(k)) > 0) then
            stpm = -HALF * dderiv(k) / (ONE - dderiv(k))
        end if
    end if
    stpm = max(slbd, min(subd, stpm))

    stplen(:, k) = [slbd, subd, stpm]
    isbd(:, k) = [ilbd, iubd, 0_IK]
end do

! The following lines calculate PREDSQ for all the 3*(NPT-1) trial points.
! First, compute VLAG = PHI(STPLEN). Using the fact that PHI_K(0) = 0, PHI_K(1) = delta_{K, KNEW}
! (Kroneker delta), and recalling the PHI_K is quadratic, we can find that
! PHI_K(t) = t*(1-t)*PHI_K'(0) for K /= KNEW, and PHI_KNEW = t*[t*(1-PHI_K'(0)) + PHI_K'(0)].
vlag = stplen * (ONE - stplen) * spread(dderiv, dim=1, ncopies=3)
!!MATLAB: vlag = stplen .* (1 - stplen) .* dderiv; % Implicit expansion; dderiv is a row!
vlag(:, knew) = stplen(:, knew) * (stplen(:, knew) * (ONE - dderiv(knew)) + dderiv(knew))
! Set NaNs in VLAG to 0 so that the behavior of MAXVAL(ABS(VLAG)) is predictable. VLAG does not have
! NaN unless XPT does, which would be a bug. MAXVAL(ABS(VLAG)) appears in Powell's code, not here.
where (is_nan(vlag)) vlag = ZERO  !!MATLAB: vlag(isnan(vlag)) = 0;
!
! Second, BETABD is the upper bound of BETA given in (3.10) of the BOBYQA paper.
betabd = HALF * (stplen * (ONE - stplen) * spread(distsq, dim=1, ncopies=3))**2
!!MATLAB: betabd = 0.5 * (stplen .* (1-stplen) .* distsq).^2 % Implicit expansion; distsq is a row!
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
!---------------------------------------------------------------------!
!isqs = int(maxloc(abs(vlag), dim=1), kind(isqs))  ! SIZE(ISQS) = NPT
!ksq = int(maxloc([(predsq(isqs(k), k), k=1, npt)], dim=1), kind(ksq))
!isq = isqs(ksq)
!---------------------------------------------------------------------!
! 2. Recall that we have set the NaN entries of PREDSQ to zero, if there is any. Thus the KSQS below
! is a well defined integer array, all the three entries lying between 1 and NPT.
ksqs = int(maxloc(predsq, dim=2), kind(ksqs))
isq = int(maxloc([predsq(1, ksqs(1)), predsq(2, ksqs(2)), predsq(3, ksqs(3))], dim=1), kind(isq))
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

! Calculate DENOM for the current choice of D. Indeed, only DEN_LINE(KNEW) is needed.
d = xline - xopt
den_line = calden(kopt, bmat, d, xpt, zmat)

!--------------------------------------------------------------------------------------------------!
! The following IF ... END IF does not exist in Powell's code. SURPRISINGLY, the performance of
! BOBYQA on bound constrained problems (but NOT unconstrained ones) is evidently improved by this IF
! ... END IF, which means to try the Cauchy step only in the late stage of the algorithm, e.g., when
! DELBAR is relatively small. WHY? In the following condition, 1.0E-2 works well if we use
! DEN_CAUCHY to decide whether to take the Cauchy step; 1.0E-3 works well if we use VLAGSQ instead.
! How to make this condition adaptive? A naive idea is to replace the thresholds to,
! e.g.,1.0E-2*RHOBEG. However, in a test on 20220517, this adaptation worsened the performance. In
! such a test, RHOBEG must take a value that is quite different from one. We tried RHOBEG = 0.9E-2.
!if (delbar > 1.0E-3) then
!if (delbar > 1.0E-1) then
if (delbar > 1.0E-2) then
    return
end if
!--------------------------------------------------------------------------------------------------!

! Prepare for the method that assembles the constrained Cauchy step in S. The sum of squares of the
! fixed components of S is formed in SFIXSQ, and the free components of S are set to BIGSTP. When
! UPHILL = 0, the method calculates the downhill version of XCAUCHY, which intends to minimize the
! KNEW-th Lagrange function; when UPHILL = 1, it calculates the uphill version that intends to
! maximize the Lagrange function.
bigstp = delbar + delbar  ! N.B.: In the sequel, S <= BIGSTP.
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
    ! In Powell's code, the subroutine returns immediately if GGFREE is 0. However, GGFREE depends
    ! on GLAG, which in turn depends on UPHILL. It can happen that GGFREE is 0 when UPHILL = 0 but
    ! not so when UPHILL= 1. Thus we skip the iteration for the current UPHILL but do not return.
    if (ggfree <= 0) then
        cycle
    end if

    ! Investigate whether more components of S can be fixed. Note that the loop counter K does not
    ! appear in the loop body. The purpose of K is only to impose an explicit bound on the number of
    ! loops. Powell's code does not have such a bound. The bound is not a true restriction, because
    ! we can check that (SFIXSQ > SSQSAV .AND. GGFREE > 0) must fail within N loops.
    sfixsq = ZERO
    grdstp = ZERO
    do k = 1, n
        resis = delbar**2 - sfixsq
        if (resis <= 0) then
            exit
        end if
        ssqsav = sfixsq
        grdstp = sqrt(resis / ggfree)
        xtemp = xopt - grdstp * glag
        mask_fixl = (s >= bigstp .and. xtemp <= sl)  ! S == BIGSTP & XTEMP == SL
        mask_fixu = (s >= bigstp .and. xtemp >= su)  ! S == BIGSTP & XTEMP == SU
        mask_free = (s >= bigstp .and. .not. (mask_fixl .or. mask_fixu))
        s(trueloc(mask_fixl)) = sl(trueloc(mask_fixl)) - xopt(trueloc(mask_fixl))
        s(trueloc(mask_fixu)) = su(trueloc(mask_fixu)) - xopt(trueloc(mask_fixu))
        sfixsq = sfixsq + sum(s(trueloc(mask_fixl .or. mask_fixu))**2)
        ggfree = sum(glag(trueloc(mask_free))**2)
        if (.not. (sfixsq > ssqsav .and. ggfree > 0)) then
            exit
        end if
    end do

    ! Set the remaining free components of S and all components of XCAUCHY. S may be scaled later.
    x(trueloc(glag > 0)) = sl(trueloc(glag > 0))
    x(trueloc(glag <= 0)) = su(trueloc(glag <= 0))
    x(trueloc(abs(s) <= 0)) = xopt(trueloc(abs(s) <= 0))
    xtemp = max(sl, min(su, xopt - grdstp * glag))
    x(trueloc(s >= bigstp)) = xtemp(trueloc(s >= bigstp))  ! S == BIGSTP
    s(trueloc(s >= bigstp)) = -grdstp * glag(trueloc(s >= bigstp))  ! S == BIGSTP
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

! Calculate the denominator rendered by the Cauchy step. Indeed, only DEN_CAUCHY(KNEW) is needed.
s = xcauchy - xopt
den_cauchy = calden(kopt, bmat, s, xpt, zmat)

! Take the Cauchy step if it is likely to render a larger denominator.
!IF (VLAGSQ_CAUCHY > MAX(DEN_LINE(KNEW), ZERO) .OR. IS_NAN(DEN_LINE(KNEW))) THEN  ! Powell's version
if (den_cauchy(knew) > max(den_line(knew), ZERO) .or. is_nan(den_line(knew))) then  ! Works better
    d = s
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n, 'SIZE(D) == N', srname)
    call assert(all(is_finite(d)), 'D is finite', srname)
    ! In theory, ||D|| <= DELBAR, which may be false due to rounding, but ||D|| >= 2*DELBAR is unlikely.
    ! It is crucial to ensure that the geometry step is nonzero, which holds in theory. However, due
    ! to the bound constraints, ||D|| may be much smaller than DELBAR.
    call assert(norm(d) > 0 .and. norm(d) < TWO * delbar, '0 < ||D|| < 2*DELBAR', srname)
    ! D is supposed to satisfy the bound constraints SL <= XOPT + D <= SU.
    call assert(all(xopt + d >= sl - TEN * EPS * max(ONE, abs(sl)) .and. &
        & xopt + d <= su + TEN * EPS * max(ONE, abs(su))), 'SL <= XOPT + D <= SU', srname)
end if

end function geostep


end module geometry_mod
