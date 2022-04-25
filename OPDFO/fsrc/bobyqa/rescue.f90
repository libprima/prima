module rescue_mod
!--------------------------------------------------------------------------------------------------!
! This module contains the `rescue` subroutine.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, April 25, 2022 AM09:42:34
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: rescue


contains


subroutine rescue(calfun, n, npt, xl, xu, iprint, maxfun, xbase, xpt, &
     &  fval, xopt, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, delta, &
     &  kopt, vlag, ptsaux, ptsid, w, f, ftarget, &
     & xhist, maxxhist, fhist, maxfhist)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: pintrf_mod, only : OBJ

! Solver-specif modules
use, non_intrinsic :: update_mod, only : update

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK) :: n
integer(IK) :: ndim
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfhist
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: maxxhist
integer(IK), intent(in) :: npt
real(RP), intent(in) :: delta
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: xl(n)
real(RP), intent(in) :: xu(n)

! In-outputs
integer(IK), intent(inout) :: kopt
integer(IK), intent(inout) :: nf
real(RP), intent(inout) :: f
real(RP), intent(inout) :: bmat(n, npt + n)
real(RP), intent(inout) :: fval(npt)
real(RP), intent(inout) :: gopt(n)
real(RP), intent(inout) :: hq(n * (n + 1_IK) / 2_IK)
real(RP), intent(inout) :: pq(npt)
real(RP), intent(inout) :: ptsaux(2, n)
real(RP), intent(inout) :: ptsid(npt)
real(RP), intent(inout) :: sl(n)
real(RP), intent(inout) :: su(n)
real(RP), intent(inout) :: vlag(npt + n)
real(RP), intent(inout) :: w(2_IK * npt + n)
real(RP), intent(inout) :: xbase(n)
real(RP), intent(inout) :: xopt(n)
real(RP), intent(inout) :: xpt(n, npt)
real(RP), intent(inout) :: zmat(npt, npt - n - 1_IK)

! Outputs
real(RP), intent(out) :: fhist(maxfhist)
real(RP), intent(out) :: xhist(n, maxxhist)

! Local variables
real(RP) :: x(n)
real(RP) :: beta, bsum, den, denom, diff,      &
&        distsq, dsqmin, fbase, hdiag, sfrac,    &
&        summ, sumpq, temp, vlmxsq, vquad, vtemp, winc, xp, xq, v1(npt), v2(npt)
integer(IK) :: i, ih, ihp, ihq, ip, iq, iw, j, jp, jpn, k, &
&           knew, kold, kpt, np, nptm, nrem

!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonZERO step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonZERO. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!
np = n + 1
sfrac = HALF / real(np, RP)
nptm = npt - np
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to ZERO. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
sumpq = ZERO
winc = ZERO
do k = 1, npt
    distsq = ZERO
    do j = 1, n
        xpt(j, k) = xpt(j, k) - xopt(j)
        distsq = distsq + xpt(j, k)**2
    end do
    sumpq = sumpq + pq(k)
    w(ndim + k) = distsq
    winc = max(winc, distsq)
    do j = 1, nptm
        zmat(k, j) = ZERO
    end do
end do
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
ih = 0
do j = 1, n
    w(j) = HALF * sumpq * xopt(j)
!-------------------------------------------------------!
!   Zaikun 20220424
    !do k = 1, npt
    !    w(j) = w(j) + pq(k) * xpt(j, k)
    !end do
    w(j) = w(j) + inprod(xpt(j, :), pq(1:npt))
!-------------------------------------------------------!

    do i = 1, j
        ih = ih + 1
        hq(ih) = hq(ih) + w(i) * xopt(j) + w(j) * xopt(i)
    end do
end do
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to ZERO, and
!     also set the elements of PTSAUX.
!
do j = 1, n
    xbase(j) = xbase(j) + xopt(j)
    sl(j) = sl(j) - xopt(j)
    su(j) = su(j) - xopt(j)
    xopt(j) = ZERO
    ptsaux(1, j) = min(delta, su(j))
    ptsaux(2, j) = max(-delta, sl(j))
    if (ptsaux(1, j) + ptsaux(2, j) < ZERO) then
        temp = ptsaux(1, j)
        ptsaux(1, j) = ptsaux(2, j)
        ptsaux(2, j) = temp
    end if
    if (abs(ptsaux(2, j)) < HALF * abs(ptsaux(1, j))) then
        ptsaux(2, j) = HALF * ptsaux(1, j)
    end if
    do i = 1, ndim
        bmat(j, i) = ZERO
    end do
end do
fbase = fval(kopt)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonZERO elements of BMAT and ZMAT.
!
ptsid(1) = sfrac
do j = 1, n
    jp = j + 1
    jpn = jp + n
    ptsid(jp) = real(j, RP) + sfrac
    if (jpn <= npt) then
        ptsid(jpn) = real(j, RP) / real(np, RP) + sfrac
        temp = ONE / (ptsaux(1, j) - ptsaux(2, j))
        bmat(j, jp) = -temp + ONE / ptsaux(1, j)
        bmat(j, jpn) = temp + ONE / ptsaux(2, j)
        bmat(j, 1) = -bmat(j, jp) - bmat(j, jpn)
        zmat(1, j) = sqrt(2.0_RP) / abs(ptsaux(1, j) * ptsaux(2, j))
        zmat(jp, j) = zmat(1, j) * ptsaux(2, j) * temp
        zmat(jpn, j) = -zmat(1, j) * ptsaux(1, j) * temp
    else
        bmat(j, 1) = -ONE / ptsaux(1, j)
        bmat(j, jp) = ONE / ptsaux(1, j)
        bmat(j, j + npt) = -HALF * ptsaux(1, j)**2
    end if
end do
!
!     Set any remaining identifiers with their nonZERO elements of ZMAT.
!
if (npt >= n + np) then
    do k = 2 * np, npt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
        iw = int((real(k - np, RP) - HALF) / real(n, RP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ip = k - np - iw * n
        iq = ip + iw
        if (iq > n) iq = iq - n
        ptsid(k) = real(ip, RP) + real(iq, RP) / real(np, RP) + sfrac
        temp = ONE / (ptsaux(1, ip) * ptsaux(1, iq))
        zmat(1, k - np) = temp
        zmat(ip + 1, k - np) = -temp
        zmat(iq + 1, k - np) = -temp
        zmat(k, k - np) = temp
    end do
end if
nrem = npt
kold = 1
knew = kopt
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
80 do j = 1, n
    temp = bmat(j, kold)
    bmat(j, kold) = bmat(j, knew)
    bmat(j, knew) = temp
end do
do j = 1, nptm
    temp = zmat(kold, j)
    zmat(kold, j) = zmat(knew, j)
    zmat(knew, j) = temp
end do
ptsid(kold) = ptsid(knew)
ptsid(knew) = ZERO
w(ndim + knew) = ZERO
nrem = nrem - 1
if (knew /= kopt) then
    temp = vlag(kold)
    vlag(kold) = vlag(knew)
    vlag(knew) = temp
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
    call update(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w)
    if (nrem == 0) goto 350
    do k = 1, npt
        w(ndim + k) = abs(w(ndim + k))
    end do
end if
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced ONE of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
120 dsqmin = ZERO
do k = 1, npt
    if (w(ndim + k) > ZERO) then
        if (dsqmin == ZERO .or. w(ndim + k) < dsqmin) then
            knew = k
            dsqmin = w(ndim + k)
        end if
    end if
end do
if (dsqmin == ZERO) goto 260
!
!     Form the W-vector of the chosen original interpolation point.
!
do j = 1, n
    w(npt + j) = xpt(j, knew)
end do
do k = 1, npt
    summ = ZERO
    if (k == kopt) then
        continue
    else if (ptsid(k) == ZERO) then
        do j = 1, n
            summ = summ + w(npt + j) * xpt(j, k)
        end do
    else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IP=PTSID(K)
        ip = int(ptsid(k))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ip > 0) summ = w(npt + ip) * ptsaux(1, ip)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
        iq = int(real(np, RP) * ptsid(k) - real(ip * np, RP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (iq > 0) then
            iw = 1
            if (ip == 0) iw = 2
            summ = summ + w(npt + iq) * ptsaux(iw, iq)
        end if
    end if
    w(k) = HALF * summ * summ
end do
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
do k = 1, npt
    summ = ZERO
    do j = 1, n
        summ = summ + bmat(j, k) * w(npt + j)
    end do
    vlag(k) = summ
end do

!------------------------------------------!
! Zaikun 20220425
v1 = vlag(1:npt); vlag(1:npt) = ZERO
!------------------------------------------!

beta = ZERO
do j = 1, nptm
    summ = ZERO
    do k = 1, npt
        summ = summ + zmat(k, j) * w(k)
    end do
    !beta = beta - summ * summ
    beta = beta - summ**2
    do k = 1, npt
        vlag(k) = vlag(k) + summ * zmat(k, j)
    end do
end do

!-----------------------------------------------!
! Zaikun 20220425
v2 = vlag(1:npt); vlag(1:npt) = v1 + v2
!-----------------------------------------------!

bsum = ZERO
distsq = ZERO
do j = 1, n
    summ = ZERO
    do k = 1, npt
        summ = summ + bmat(j, k) * w(k)
    end do

    jp = j + npt
    bsum = bsum + summ * w(jp)
    do ip = npt + 1, ndim
        summ = summ + bmat(j, ip) * w(ip)
    end do

    bsum = bsum + summ * w(jp)
    vlag(jp) = summ
    distsq = distsq + xpt(j, knew)**2
end do

!-------------------------------------------------!
! Zaikun 20220425
bsum = inprod(w(1:n), matprod(bmat(:, 1:npt), w(1:npt)) + matprod(bmat, w(1:npt + n)))
!beta = HALF * distsq * distsq + beta - bsum
beta = HALF * distsq**2 + beta - bsum

!beta = HALF * distsq**2 - inprod(vlag, w(1:npt + n))
!-------------------------------------------------!

vlag(kopt) = vlag(kopt) + ONE
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
denom = ZERO
vlmxsq = ZERO
do k = 1, npt
    if (ptsid(k) /= ZERO) then
        hdiag = ZERO
        do j = 1, nptm
            hdiag = hdiag + zmat(k, j)**2
        end do
        !den = beta * hdiag + vlag(k)**2
        den = hdiag * beta + vlag(k)**2
        if (den > denom) then
            kold = k
            denom = den
        end if
    end if
    vlmxsq = max(vlmxsq, vlag(k)**2)
end do
if (denom <= 1.0E-2_RP * vlmxsq) then
    w(ndim + knew) = -w(ndim + knew) - winc
    goto 120
end if
goto 80
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be dONE. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to ZERO,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
260 do kpt = 1, npt
    if (ptsid(kpt) == ZERO) cycle
    if (nf >= maxfun) then
        nf = -1
        goto 350
    end if
    ih = 0
    do j = 1, n
        w(j) = xpt(j, kpt)
        xpt(j, kpt) = ZERO
        temp = pq(kpt) * w(j)
        do i = 1, j
            ih = ih + 1
            hq(ih) = hq(ih) + temp * w(i)
        end do
    end do
    pq(kpt) = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IP=PTSID(KPT)
!      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
    ip = int(ptsid(kpt))
    iq = int(real(np, RP) * ptsid(kpt) - real(ip * np, RP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ip > 0) then
        xp = ptsaux(1, ip)
        xpt(ip, kpt) = xp
    end if
    if (iq > 0) then
        xq = ptsaux(1, iq)
        if (ip == 0) xq = ptsaux(2, iq)
        xpt(iq, kpt) = xq
    end if
!
!     Set VQUAD to the value of the current model at the new point.
!
    vquad = fbase
    if (ip > 0) then
        ihp = (ip + ip * ip) / 2
        vquad = vquad + xp * (gopt(ip) + HALF * xp * hq(ihp))
    end if
    if (iq > 0) then
        ihq = (iq + iq * iq) / 2
        vquad = vquad + xq * (gopt(iq) + HALF * xq * hq(ihq))
        if (ip > 0) then
            iw = max(ihp, ihq) - abs(ip - iq)
            vquad = vquad + xp * xq * hq(iw)
        end if
    end if
    vtemp = ZERO
    do k = 1, npt
        temp = ZERO
        if (ip > 0) temp = temp + xp * xpt(ip, k)
        if (iq > 0) temp = temp + xq * xpt(iq, k)
        vtemp = vtemp + pq(k) * temp * temp
    end do
    vquad = vquad + HALF * vtemp
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
    do i = 1, n
        w(i) = min(max(xl(i), xbase(i) + xpt(i, kpt)), xu(i))
        if (xpt(i, kpt) == sl(i)) w(i) = xl(i)
        if (xpt(i, kpt) == su(i)) w(i) = xu(i)
    end do
    nf = nf + 1

!-------------------------------------------------------------------!
    !call calfun(n, w, f)
    x = w(1:n)
    call evaluate(calfun, x, f)
    call savehist(nf, x, xhist, f, fhist)
!-------------------------------------------------------------------!

    fval(kpt) = f
    if (f < fval(kopt)) kopt = kpt
    diff = f - vquad
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
    do i = 1, n
        gopt(i) = gopt(i) + diff * bmat(i, kpt)
    end do
    do k = 1, npt
        summ = ZERO
        do j = 1, nptm
            summ = summ + zmat(k, j) * zmat(kpt, j)
        end do
        temp = diff * summ
        if (ptsid(k) == ZERO) then
            pq(k) = pq(k) + temp
        else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IP=PTSID(K)
!          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
            ip = int(ptsid(k))
            iq = int(real(np, RP) * ptsid(k) - real(ip * np, RP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ihq = (iq * iq + iq) / 2
            if (ip == 0) then
                hq(ihq) = hq(ihq) + temp * ptsaux(2, iq)**2
            else
                ihp = (ip * ip + ip) / 2
                hq(ihp) = hq(ihp) + temp * ptsaux(1, ip)**2
                if (iq > 0) then
                    hq(ihq) = hq(ihq) + temp * ptsaux(1, iq)**2
                    iw = max(ihp, ihq) - abs(iq - ip)
                    hq(iw) = hq(iw) + temp * ptsaux(1, ip) * ptsaux(1, iq)
                end if
            end if
        end if
    end do
    ptsid(kpt) = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 03-06-2019):
!     If a NaN or an infinite value has been reached during the
!     evaluation of the objective function, the loop exit after setting
!     all the parameters, not to raise an exception. KOPT is set to KPT
!     to check in BOBYQB weather FVAL(KOPT) is NaN or infinite value or
!     not.
    if (is_nan(f) .or. is_posinf(f)) then
        exit
    end if
!     By Tom (on 04-06-2019):
!     If the target function value is reached, the loop exit and KOPT is
!     set to KPT to check in BOBYQB weather FVAL(KOPT) .LE. FTARGET
    if (f <= ftarget) then
        exit
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end do
350 return
end subroutine rescue


end module rescue_mod
