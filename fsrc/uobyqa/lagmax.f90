subroutine lagmax(n, g, h, rho, d, v, vmax)

use, non_intrinsic :: consts_mod, only : RP, IK
implicit real(RP) (a - h, o - z)
implicit integer(IK) (i - n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension g(*), h(n, *), d(*), v(*)
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
half = 0.5D0
halfrt = sqrt(half)
one = 1.0D0
zero = 0.0D0
!
!     Pick V such that ||HV|| / ||V|| is large.
!
hmax = zero
do i = 1, n
    sum = zero
    do j = 1, n
        h(j, i) = h(i, j)
        sum = sum + h(i, j)**2
    end do
    if (sum > hmax) then
        hmax = sum
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
vsq = zero
vhv = zero
dsq = zero
do i = 1, n
    vsq = vsq + v(i)**2
    d(i) = zero
    do j = 1, n
        d(i) = d(i) + h(i, j) * v(j)
    end do
    vhv = vhv + v(i) * d(i)
    dsq = dsq + d(i)**2
end do
if (vhv * vhv <= 0.9999D0 * dsq * vsq) then
    temp = vhv / vsq
    wsq = zero
    do i = 1, n
        d(i) = d(i) - temp * v(i)
        wsq = wsq + d(i)**2
    end do
    whw = zero
    ratio = sqrt(wsq / vsq)
    do i = 1, n
        temp = zero
        do j = 1, n
            temp = temp + h(i, j) * d(j)
        end do
        whw = whw + temp * d(i)
        v(i) = ratio * v(i)
    end do
    vhv = ratio * ratio * vhv
    vhw = ratio * wsq
    temp = half * (whw - vhv)
    temp = temp + sign(sqrt(temp**2 + vhw**2), whw + vhv)
    do i = 1, n
        d(i) = vhw * v(i) + temp * d(i)
    end do
end if
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
gg = zero
gd = zero
dd = zero
dhd = zero
do i = 1, n
    gg = gg + g(i)**2
    gd = gd + g(i) * d(i)
    dd = dd + d(i)**2
    sum = zero
    do j = 1, n
        sum = sum + h(i, j) * d(j)
    end do
    dhd = dhd + sum * d(i)
end do
temp = gd / gg
vv = zero
scale = sign(rho / sqrt(dd), gd * dhd)
do i = 1, n
    v(i) = d(i) - temp * g(i)
    vv = vv + v(i)**2
    d(i) = scale * d(i)
end do
gnorm = sqrt(gg)
if (gnorm * dd <= 0.5D-2 * rho * abs(dhd) .or. vv / dd <= 1.0D-4) then
    vmax = abs(scale * (gd + half * scale * dhd))
    goto 170
end if
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.
!
ghg = zero
vhg = zero
vhv = zero
do i = 1, n
    sum = zero
    sumv = zero
    do j = 1, n
        sum = sum + h(i, j) * g(j)
        sumv = sumv + h(i, j) * v(j)
    end do
    ghg = ghg + sum * g(i)
    vhg = vhg + sumv * g(i)
    vhv = vhv + sumv * v(i)
end do
vnorm = sqrt(vv)
ghg = ghg / gg
vhg = vhg / (vnorm * gnorm)
vhv = vhv / vv
if (abs(vhg) <= 0.01D0 * max(abs(ghg), abs(vhv))) then
    vmu = ghg - vhv
    wcos = one
    wsin = zero
else
    temp = half * (ghg - vhv)
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
tempa = abs(dlin) + half * abs(vmu + vhv)
tempb = abs(vlin) + half * abs(ghg - vmu)
tempc = halfrt * (abs(dlin) + abs(vlin)) + 0.25D0 * abs(ghg + vhv)
if (tempa >= tempb .and. tempa >= tempc) then
    tempd = sign(rho, dlin * (vmu + vhv))
    tempv = zero
else if (tempb >= tempc) then
    tempd = zero
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
end
