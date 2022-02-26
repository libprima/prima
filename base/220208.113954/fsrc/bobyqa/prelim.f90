subroutine prelim(n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, &
& xpt, fval, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, kopt, f, ftarget)

implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
dimension x(*), xl(*), xu(*), xbase(*), xpt(npt, *), fval(*), gopt(*), &
& hq(*), pq(*), bmat(ndim, *), zmat(npt, *), sl(*), su(*)
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
half = 0.5D0
one = 1.0D0
two = 2.0D0
zero = 0.0D0
rhosq = rhobeg * rhobeg
recip = one / rhosq
np = n + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
almost_infinity = huge(0.0D0) / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
do j = 1, n
    xbase(j) = x(j)
    do k = 1, npt
        xpt(k, j) = zero
    end do
    do i = 1, ndim
        bmat(i, j) = zero
    end do
end do
do ih = 1, (n * np) / 2
    hq(ih) = zero
end do
do k = 1, npt
    pq(k) = zero
    do j = 1, npt - np
        zmat(k, j) = zero
    end do
end do
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
nf = 0
50 nfm = nf
nfx = nf - n
nf = nf + 1
if (nfm <= 2 * n) then
    if (nfm >= 1 .and. nfm <= n) then
        stepa = rhobeg
        if (su(nfm) == zero) stepa = -stepa
        xpt(nf, nfm) = stepa
    else if (nfm > n) then
        stepa = xpt(nf - n, nfx)
        stepb = -rhobeg
        if (sl(nfx) == zero) stepb = dmin1(two * rhobeg, su(nfx))
        if (su(nfx) == zero) stepb = dmax1(-two * rhobeg, sl(nfx))
        xpt(nf, nfx) = stepb
    end if
else
    itemp = (nfm - np) / n
    jpt = nfm - itemp * n - n
    ipt = jpt + itemp
    if (ipt > n) then
        itemp = jpt
        jpt = ipt - n
        ipt = itemp
    end if
    xpt(nf, ipt) = xpt(ipt + 1, ipt)
    xpt(nf, jpt) = xpt(jpt + 1, jpt)
end if
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
do j = 1, n
    x(j) = dmin1(dmax1(xl(j), xbase(j) + xpt(nf, j)), xu(j))
    if (xpt(nf, j) == sl(j)) x(j) = xl(j)
    if (xpt(nf, j) == su(j)) x(j) = xu(j)
end do
call calfun(n, x, f)
fval(nf) = f
if (nf == 1) then
    fbeg = f
    kopt = 1
else if (f < fval(kopt)) then
    kopt = nf
end if
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
if (nf <= 2 * n + 1) then
    if (nf >= 2 .and. nf <= n + 1) then
        gopt(nfm) = (f - fbeg) / stepa
        if (npt < nf + n) then
            bmat(1, nfm) = -one / stepa
            bmat(nf, nfm) = one / stepa
            bmat(npt + nfm, nfm) = -half * rhosq
        end if
    else if (nf >= n + 2) then
        ih = (nfx * (nfx + 1)) / 2
        temp = (f - fbeg) / stepb
        diff = stepb - stepa
        hq(ih) = two * (temp - gopt(nfx)) / diff
        gopt(nfx) = (gopt(nfx) * stepb - temp * stepa) / diff
        if (stepa * stepb < zero) then
            if (f < fval(nf - n)) then
                fval(nf) = fval(nf - n)
                fval(nf - n) = f
                if (kopt == nf) kopt = nf - n
                xpt(nf - n, nfx) = stepb
                xpt(nf, nfx) = stepa
            end if
        end if
        bmat(1, nfx) = -(stepa + stepb) / (stepa * stepb)
        bmat(nf, nfx) = -half / xpt(nf - n, nfx)
        bmat(nf - n, nfx) = -bmat(1, nfx) - bmat(nf, nfx)
        zmat(1, nfx) = dsqrt(two) / (stepa * stepb)
        zmat(nf, nfx) = dsqrt(half) / rhosq
        zmat(nf - n, nfx) = -zmat(1, nfx) - zmat(nf, nfx)
    end if
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
else
    ih = (ipt * (ipt - 1)) / 2 + jpt
    zmat(1, nfx) = recip
    zmat(nf, nfx) = recip
    zmat(ipt + 1, nfx) = -recip
    zmat(jpt + 1, nfx) = -recip
    temp = xpt(nf, ipt) * xpt(nf, jpt)
    hq(ih) = (fbeg - fval(ipt + 1) - fval(jpt + 1) + f) / temp
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 04-06-2019):
!     If the evaluation returns an NaN or an infinity value, this
!     subroutine is stopped.
if (f /= f .or. f > almost_infinity) goto 80
!     By Tom (on 04-06-2019):
!     If the target value is reached, stop the algorithm.
if (f <= ftarget) goto 80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nf < npt .and. nf < maxfun) goto 50
80 return
end
