module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Friday, May 06, 2022 PM09:41:35
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(g, h, rho, d, vmax)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, QUART, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : issymmetric, matprod, inprod

implicit none

! Inputs
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: h(:, :)  ! H(N, N)
real(RP), intent(in) :: rho  ! NEWUOA etc uses DELBAR, which is NOT RHO; possible improvement?

! Outputs
real(RP), intent(out) :: d(:)  ! D(N)
real(RP), intent(out) :: vmax

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: n
real(RP) :: v(size(g))
real(RP) :: dd, dhd, dlin, gd, gg, ghg, gnorm, &
&        ratio, scaling, temp, &
&        tempa, tempb, tempc, tempd, tempv, vhg, vhv, vhd, &
&        vlin, vmu, vnorm, vv, wcos, wsin, hv(size(g))
integer(IK) :: k


! Sizes.
n = int(size(g), kind(n))

if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(rho > 0, 'RHO > 0', srname)
    call assert(size(h, 1) == n .and. issymmetric(h), 'H is n-by-n and symmetric', srname)
    call assert(size(d) == n, 'SIZE(D) == N', srname)
end if

!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     RHO is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.
!
!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. RHO
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.
!
!     Preliminary calculations.
!


if (is_nan(sum(abs(h)) + sum(abs(g)))) then
    d = ZERO
    vmax = ZERO
    return
end if

! Pick V such that ||HV|| / ||V|| is large.
k = int(maxloc(sum(h**2, dim=1), dim=1), IK)
v = h(:, k)

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
    temp = temp + sign(sqrt(temp**2 + vhd**2), dhd + vhv)
    d = vhd * v + temp * d
end if

! We now turn our attention to the subspace span{G, D}. A multiple of the current D is returned if
! that choice seems to be adequate.
gg = sum(g**2)
gd = inprod(g, d)
dd = sum(d**2)
dhd = inprod(d, matprod(h, d))

! Zaikun 20220504: GG and DD can become 0 at this point due to rounding. Detected by IFORT.
if (.not. (gg > 0 .and. dd > 0)) then
    d = ZERO
    vmax = ZERO
    return
end if

scaling = sign(rho / sqrt(dd), gd * dhd)
v = d - (gd / gg) * g
vv = sum(v**2)
d = scaling * d
gnorm = sqrt(gg)

if (.not. (gnorm * dd > 0.5E-2_RP * rho * abs(dhd) .and. vv > 1.0E-4_RP * dd)) then
    vmax = abs(scaling * (gd + HALF * scaling * dhd))
    return
end if

! G and V are now orthogonal in the subspace span{G, D}. Hence we generate an orthonormal basis of
! this subspace such that (D, HV) is negligible or 0, where D and V will be the basis vectors.
ghg = inprod(g, matprod(h, g))
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
    vmu = temp + sign(sqrt(temp**2 + vhg**2), temp)
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
dlin = wcos * gnorm / rho
vlin = -wsin * gnorm / rho
tempa = abs(dlin) + HALF * abs(vmu + vhv)
tempb = abs(vlin) + HALF * abs(ghg - vmu)
tempc = sqrt(HALF) * (abs(dlin) + abs(vlin)) + QUART * abs(ghg + vhv)
if (tempa >= tempb .and. tempa >= tempc) then
    tempd = sign(rho, dlin * (vmu + vhv))
    tempv = ZERO
else if (tempb >= tempc) then
    tempd = ZERO
    tempv = sign(rho, vlin * (ghg - vmu))
else
    tempd = sign(sqrt(HALF) * rho, dlin * (ghg + vhv))
    tempv = sign(sqrt(HALF) * rho, vlin * (ghg + vhv))
end if
d = tempd * d + tempv * v
vmax = rho * rho * max(tempa, tempb, tempc)

end subroutine geostep


end module geometry_mod
