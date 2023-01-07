module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, December 28, 2022 AM01:39:06
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: setdrop_tr, geostep


contains


function setdrop_tr(kopt, ximproved, d, pl, rho, xpt) result(knew)
!--------------------------------------------------------------------------------------------------!
! This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
! REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
! XPT(:, KNEW) is replaced by XNEW = XOPT + D, where D is the trust-region step.
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

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: powalg_mod, only : calvlag

implicit none

! Inputs
integer(IK), intent(in) :: kopt
logical, intent(in) :: ximproved
real(RP), intent(in) :: d(:)  ! D(N)
real(RP), intent(in) :: pl(:, :)  ! PL(NPT-1, NPT)
real(RP), intent(in) :: rho
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK) :: knew

! Local variables
character(len=*), parameter :: srname = 'SETDROP_TR'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: distsq(size(xpt, 2))
real(RP) :: score(size(xpt, 2))
real(RP) :: vlag(size(xpt, 2))
real(RP) :: weight(size(xpt, 2))

! Sizes
n = int(size(xpt, 1), kind(npt))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(pl, 1) == npt - 1 .and. size(pl, 2) == npt, 'SIZE(PL)==[NPT - 1, NPT]', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the distance squares between the interpolation points and the "optimal point". When
! identifying the optimal point, as suggested in (56) of the UOBYQA paper and (7.5) of the NEWUOA
! paper, it is reasonable to take into account the new trust-region trial point XPT(:, KOPT) + D,
! which will become the optimal point in the next interpolation if XIMPROVED is TRUE.
if (ximproved) then
    distsq = sum((xpt - spread(xpt(:, kopt) + d, dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - (xpt(:, kopt) + d)).^2)  % d should be a column! Implicit expansion
else
    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
    !!MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
end if

weight = max(ONE, distsq / rho**2)**4

!------------------------------------------------------------------------------------------!
! Other possible definitions of WEIGHT.
!weight = max(ONE, distsq / rho**2)**3.5_RP ! ! No better than power 4.
!weight = max(ONE, distsq / delta**2)**3.5_RP  ! Not better than DISTSQ/RHO**2.
!weight = max(ONE, distsq / rho**2)**1.5_RP  ! Powell's origin code: power 1.5.
!weight = max(ONE, distsq / rho**2)**2  ! Better than power 1.5.
!weight = max(ONE, distsq / delta**2)**2  ! Not better than DISTSQ/RHO**2.
!weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**2  ! The same as DISTSQ/RHO**2.
!weight = distsq**2  ! Not better than MAX(ONE, DISTSQ/RHO**2)**2
!weight = max(ONE, distsq / rho**2)**3  ! Better than power 2.
!weight = max(ONE, distsq / delta**2)**3  ! Similar to DISTSQ/RHO**2; not better than it.
!weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3  ! The same as DISTSQ/RHO**2.
!weight = distsq**3  ! Not better than MAX(ONE, DISTSQ/RHO**2)**3
!weight = max(ONE, distsq / delta**2)**4  ! Not better than DISTSQ/RHO**2.
!weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**4  ! The same as DISTSQ/RHO**2.
!weight = distsq**4  ! Not better than MAX(ONE, DISTSQ/RHO**2)**4
!------------------------------------------------------------------------------------------!

vlag = calvlag(pl, d, xpt(:, kopt), kopt)
score = weight * abs(vlag)

! If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
if (.not. ximproved) then
    score(kopt) = -ONE
end if

knew = 0_IK
! Changing the IF below to `IF (ANY(SCORE>0)) THEN` does not render a better performance.
if (any(score > 1) .or. (ximproved .and. any(score > 0))) then
    ! SCORE(K) is NaN implies VLAG(K) is NaN, but we want ABS(VLAG) to be big. So we
    ! exclude such K.
    knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
            !!MATLAB: [~, knew] = max(score, [], 'omitnan');
elseif (ximproved) then
    ! Powell's code does not include the following instructions. With Powell's code,
    ! if DENABS consists of only NaN, then KNEW can be 0 even when XIMPROVED is TRUE.
    knew = int(maxloc(distsq, dim=1), IK)
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


function geostep(g, h, delbar, xopt, xpt, knew) result(d)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : issymmetric, matprod, inprod, norm, trueloc

implicit none

! Inputs
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: h(:, :)  ! H(N, N)
real(RP), intent(in) :: delbar
real(RP) :: xpt(:, :), xopt(:)
integer(IK) :: knew

! Outputs
real(RP) :: d(size(g))  ! D(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: n
real(RP) :: v(size(g))
real(RP) :: dcauchy(size(g)), dd, dhd, dlin, gd, gg, ghg, gnorm, &
&        ratio, scaling, temp, &
&        tempa, tempb, tempc, tempd, tempv, vhg, vhv, vhd, &
&        vlin, vmu, vnorm, vv, wcos, wsin, hv(size(g))
integer(IK) :: k


! Sizes.
n = int(size(g), kind(n))

if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(delbar > 0, 'DELBAR > 0', srname)
    call assert(size(h, 1) == n .and. issymmetric(h), 'H is n-by-n and symmetric', srname)
    call assert(size(d) == n, 'SIZE(D) == N', srname)
end if

!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     DELBAR is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.
!
!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. DELBAR
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.
!

! Calculate the Cauchy step as a backup.
gg = sum(g**2)
ghg = inprod(g, matprod(h, g))
dcauchy = ZERO
if (gg > 0) then
    dcauchy = (delbar / sqrt(gg)) * g
    if (ghg < 0) then
        dcauchy = -dcauchy
    end if
    !dcauchy(trueloc(is_nan(dcauchy))) = ZERO  ! DCAUCHY may contain NaN if the problem is ill-conditioned.
end if

!if (.not. any(abs(dcauchy) > 0)) then
!if (sum(dcauchy**2) <= 0) then
if (any(is_nan(dcauchy)) .or. .not. any(abs(dcauchy) > 0)) then
    dcauchy = xpt(:, knew) - xopt
    dd = sum(dcauchy**2)
    !dcauchy = min(HALF, delbar / sqrt(dd)) * dcauchy
    dcauchy = max(0.6_RP * (delbar / sqrt(dd)), min(HALF, delbar / sqrt(dd))) * dcauchy
    if (inprod(g, dcauchy) * inprod(dcauchy, matprod(h, dcauchy)) < 0) then
        dcauchy = -dcauchy
    end if
end if

if (is_nan(sum(abs(h)) + sum(abs(g))) .or. sum(abs(h)) <= 0) then
    d = dcauchy
    return
end if

if (n == 1) then
    if (g(1) * h(1, 1) > 0) then
        d = delbar
    else
        d = -delbar
    end if
    return
end if


! Pick V such that ||HV|| / ||V|| is large.
k = int(maxloc(sum(h**2, dim=1), dim=1), IK)
v = h(:, k)
v = v / sqrt(sum(v**2))

! Set D to a vector in the subspace span{V, HV} that maximizes |(D, HD)|/(D, D), except that we set
! D = HV if V and HV are nearly parallel.
vv = sum(v**2)
d = matprod(h, v)
vhv = inprod(v, d)
if (vhv * vhv <= 0.9999_RP * sum(d**2) * vv) then
    d = d - (vhv / vv) * v
    dd = sum(d**2)
    ratio = sqrt(dd / vv)
    dhd = inprod(d, matprod(h, d))
    v = ratio * v
    vhv = ratio * ratio * vhv
    vhd = ratio * dd
    temp = HALF * (dhd - vhv)
    if (dhd + vhv < 0) then
        d = vhd * v + (temp - sqrt(temp**2 + vhd**2)) * d
    else
        d = vhd * v + (temp + sqrt(temp**2 + vhd**2)) * d
    end if
end if

! We now turn our attention to the subspace span{G, D}. A multiple of the current D is returned if
! that choice seems to be adequate.
gd = inprod(g, d)
dd = sum(d**2)
dhd = inprod(d, matprod(h, d))

! Zaikun 20220504: GG and DD can become 0 at this point due to rounding. Detected by IFORT.
if (.not. (gg > 0 .and. dd > 0)) then
    d = dcauchy
    return
end if

v = d - (gd / gg) * g
vv = sum(v**2)
if (gd * dhd < 0) then
    scaling = -delbar / sqrt(dd)
else
    scaling = delbar / sqrt(dd)
end if
d = scaling * d
!d(trueloc(is_nan(d))) = ZERO
gnorm = sqrt(gg)

if (.not. (gnorm * dd > 0.5E-2_RP * delbar * abs(dhd) .and. vv > 1.0E-4_RP * dd)) then
    !if (sum(d**2) <= 0) then
    !if (any(is_nan(d))) then
    if (.not. sum(abs(d)) > 0) then
        d = dcauchy
    end if
    return
end if

! G and V are now orthogonal in the subspace span{G, D}. Hence we generate an orthonormal basis of
! this subspace such that (D, HV) is negligible or 0, where D and V will be the basis vectors.
hv = matprod(h, v)
vhg = inprod(g, hv)
vhv = inprod(v, hv)
vnorm = sqrt(vv)
ghg = ghg / gg
vhg = vhg / (vnorm * gnorm)
vhv = vhv / vv
if (abs(vhg) <= 0.01_RP * max(abs(ghg), abs(vhv))) then
    vmu = ghg - vhv
    wcos = ONE
    wsin = ZERO
else
    temp = HALF * (ghg - vhv)
    if (temp < 0) then
        vmu = temp - sqrt(temp**2 + vhg**2)
    else
        vmu = temp + sqrt(temp**2 + vhg**2)
    end if
    temp = sqrt(vmu**2 + vhg**2)
    wcos = vmu / temp
    wsin = vhg / temp
end if
tempa = wcos / gnorm
tempb = wsin / vnorm
tempc = wcos / vnorm
tempd = wsin / gnorm
d = tempa * g + tempb * v
v = tempc * v - tempd * g

! The final D is a multiple of the current D, V, D + V or D - V. We make the choice from these
! possibilities that is optimal.
dlin = wcos * gnorm / delbar
vlin = -wsin * gnorm / delbar
tempa = abs(dlin) + HALF * abs(vmu + vhv)
tempb = abs(vlin) + HALF * abs(ghg - vmu)
tempc = sqrt(HALF) * (abs(dlin) + abs(vlin)) + QUART * abs(ghg + vhv)
if (tempa >= tempb .and. tempa >= tempc) then
    if (dlin * (vmu + vhv) < 0) then
        tempd = -delbar
    else
        tempd = delbar
    end if
    tempv = ZERO
else if (tempb >= tempc) then
    tempd = ZERO
    if (vlin * (ghg - vmu) < 0) then
        tempv = -delbar
    else
        tempv = delbar
    end if
else
    if (dlin * (ghg + vhv) < 0) then
        tempd = -sqrt(HALF) * delbar
    else
        tempd = sqrt(HALF) * delbar
    end if
    if (vlin * (ghg + vhv) < 0) then
        tempv = -sqrt(HALF) * delbar
    else
        tempv = sqrt(HALF) * delbar
    end if
end if
d = tempd * d + tempv * v
!d(trueloc(is_nan(d))) = ZERO  ! D may contain NaN if the problem is ill-conditioned.
!if (sum(d**2) <= 0) then
if (any(is_nan(d))) then
    d = dcauchy
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    ! Due to rounding, it may happen that |D| > DELBAR, but |D| > 2*DELBAR is highly improbable.
    call assert(norm(d) <= TWO * delbar, '|D| <= 2*DELBAR', srname)
end if
end function geostep


end module geometry_mod
