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
! Last Modified: Thursday, March 03, 2022 PM08:45:26
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(n, g, h, rho, d, v, vmax)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF

implicit none

! Inputs
integer(IK), intent(in) :: n
real(RP), intent(in) :: g(n)
real(RP), intent(in) :: rho

! In-outputs
real(RP), intent(out) :: d(n)
real(RP), intent(inout) :: h(n, n)
real(RP), intent(out) :: v(n)

! Outputs
real(RP), intent(out) :: vmax

! Local variables
real(RP) :: dd, dhd, dlin, dsq, gd, gg, ghg, gnorm, &
&        halfrt, hmax, ratio, scaling, summ, sumv, temp, &
&        tempa, tempb, tempc, tempd, tempv, vhg, vhv, vhw, &
&        vlin, vmu, vnorm, vsq, vv, wcos, whw, wsin, wsq
integer(IK) :: i, j, k
!
!     N is the number of variables of a quadratic objective function, q say.
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

halfrt = sqrt(HALF)

!
!     Pick V such that ||HV|| / ||V|| is large.
!
!--------------------------------------------------------------------------------------------------!
! Zaikun 20220303: The procedure below may leave K uninitialized due to NaN, leading to SEGFAULT.
k = 1_IK
!--------------------------------------------------------------------------------------------------!
hmax = ZERO
do i = 1, n
    summ = ZERO
    do j = 1, n
        h(j, i) = h(i, j)
        summ = summ + h(i, j)**2
    end do
    if (summ > hmax) then
        hmax = summ
        k = i
    end if
end do
do j = 1, n
    v(j) = h(k, j)
end do
!
!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.
!
vsq = ZERO
vhv = ZERO
dsq = ZERO
do i = 1, n
    vsq = vsq + v(i)**2
    d(i) = ZERO
    do j = 1, n
        d(i) = d(i) + h(i, j) * v(j)
    end do
    vhv = vhv + v(i) * d(i)
    dsq = dsq + d(i)**2
end do
if (vhv * vhv <= 0.9999D0 * dsq * vsq) then
    temp = vhv / vsq
    wsq = ZERO
    do i = 1, n
        d(i) = d(i) - temp * v(i)
        wsq = wsq + d(i)**2
    end do
    whw = ZERO
    ratio = sqrt(wsq / vsq)
    do i = 1, n
        temp = ZERO
        do j = 1, n
            temp = temp + h(i, j) * d(j)
        end do
        whw = whw + temp * d(i)
        v(i) = ratio * v(i)
    end do
    vhv = ratio * ratio * vhv
    vhw = ratio * wsq
    temp = HALF * (whw - vhv)
    temp = temp + sign(sqrt(temp**2 + vhw**2), whw + vhv)
    do i = 1, n
        d(i) = vhw * v(i) + temp * d(i)
    end do
end if
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
gg = ZERO
gd = ZERO
dd = ZERO
dhd = ZERO
do i = 1, n
    gg = gg + g(i)**2
    gd = gd + g(i) * d(i)
    dd = dd + d(i)**2
    summ = ZERO
    do j = 1, n
        summ = summ + h(i, j) * d(j)
    end do
    dhd = dhd + summ * d(i)
end do
temp = gd / gg
vv = ZERO
scaling = sign(rho / sqrt(dd), gd * dhd)
do i = 1, n
    v(i) = d(i) - temp * g(i)
    vv = vv + v(i)**2
    d(i) = scaling * d(i)
end do
gnorm = sqrt(gg)
if (gnorm * dd <= 0.5D-2 * rho * abs(dhd) .or. vv / dd <= 1.0D-4) then
    vmax = abs(scaling * (gd + HALF * scaling * dhd))
    goto 170
end if
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or ZERO, where D and V will be the basis vectors.
!
ghg = ZERO
vhg = ZERO
vhv = ZERO
do i = 1, n
    summ = ZERO
    sumv = ZERO
    do j = 1, n
        summ = summ + h(i, j) * g(j)
        sumv = sumv + h(i, j) * v(j)
    end do
    ghg = ghg + summ * g(i)
    vhg = vhg + sumv * g(i)
    vhv = vhv + sumv * v(i)
end do
vnorm = sqrt(vv)
ghg = ghg / gg
vhg = vhg / (vnorm * gnorm)
vhv = vhv / vv
if (abs(vhg) <= 0.01D0 * max(abs(ghg), abs(vhv))) then
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
do i = 1, n
    d(i) = tempa * g(i) + tempb * v(i)
    v(i) = tempc * v(i) - tempd * g(i)
end do
!
!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.
!
dlin = wcos * gnorm / rho
vlin = -wsin * gnorm / rho
tempa = abs(dlin) + HALF * abs(vmu + vhv)
tempb = abs(vlin) + HALF * abs(ghg - vmu)
tempc = halfrt * (abs(dlin) + abs(vlin)) + 0.25D0 * abs(ghg + vhv)
if (tempa >= tempb .and. tempa >= tempc) then
    tempd = sign(rho, dlin * (vmu + vhv))
    tempv = ZERO
else if (tempb >= tempc) then
    tempd = ZERO
    tempv = sign(rho, vlin * (ghg - vmu))
else
    tempd = sign(halfrt * rho, dlin * (ghg + vhv))
    tempv = sign(halfrt * rho, vlin * (ghg + vhv))
end if
do i = 1, n
    d(i) = tempd * d(i) + tempv * v(i)
end do
vmax = rho * rho * max(tempa, tempb, tempc)
170 return
end subroutine geostep


end module geometry_mod
