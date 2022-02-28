module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, February 28, 2022 PM05:44:06
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initialize


contains


subroutine initialize(calfun, iprint, ftarget, rhobeg, sl, su, xl, xu, x, &
    & kopt, nf, bmat, f, fhist, fval, gopt, hq, pq, xbase, xhist, xpt, zmat)

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: pintrf_mod, only : OBJ

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK), intent(in) :: iprint
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: sl(:)  ! SL(N)
real(RP), intent(in) :: su(:)  ! SU(N)
real(RP), intent(in) :: xl(:)  ! XL(N)
real(RP), intent(in) :: xu(:)  ! XU(N)

! In-outputs
real(RP), intent(inout) :: x(:)  ! X(N)

! Outputs
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: fval(:)  ! FVAL(NPT)
real(RP), intent(out) :: gopt(:)  ! GOPT(N)
real(RP), intent(out) :: hq(:)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)
real(RP), intent(out) :: xbase(:)  ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(out) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Local variables
character(len=*), parameter :: srname = 'INITIIALIZE'
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
real(RP) :: diff, fbeg, recip, rhosq, stepa, stepb, temp
integer(IK) :: i, ih, ipt, itemp, j, jpt, k, nfm, nfx, np

! Sizes.
n = int(size(x), kind(n))
npt = int(size(fval), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)

    !----------------------------------------------------------------------------------------------!
    !call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)
    call assert(size(hq) == n * (n + 1_IK) / 2_IK, 'HQ is n-by-n and symmetric', srname)
    !----------------------------------------------------------------------------------------------!

    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
end if

!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonZERO, BOBYQB will change it to its usual value later.
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
rhosq = rhobeg * rhobeg
recip = ONE / rhosq
np = n + 1
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to ZERO.
!
do j = 1, n
    xbase(j) = x(j)
    do k = 1, npt
        xpt(j, k) = ZERO
    end do
    do i = 1, npt + n
        bmat(j, i) = ZERO
    end do
end do
do ih = 1, (n * np) / 2
    hq(ih) = ZERO
end do
do k = 1, npt
    pq(k) = ZERO
    do j = 1, npt - np
        zmat(k, j) = ZERO
    end do
end do
!
!     Begin the initialization procedure. NF becomes ONE more than the number
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
        if (su(nfm) == ZERO) stepa = -stepa
        xpt(nfm, nf) = stepa
    else if (nfm > n) then
        stepa = xpt(nfx, nf - n)
        stepb = -rhobeg
        if (sl(nfx) == ZERO) stepb = min(TWO * rhobeg, su(nfx))
        if (su(nfx) == ZERO) stepb = max(-TWO * rhobeg, sl(nfx))
        xpt(nfx, nf) = stepb
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
    xpt(ipt, nf) = xpt(ipt, ipt + 1)
    xpt(jpt, nf) = xpt(jpt, jpt + 1)
end if
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
do j = 1, n
    x(j) = min(max(xl(j), xbase(j) + xpt(j, nf)), xu(j))
    if (xpt(j, nf) == sl(j)) x(j) = xl(j)
    if (xpt(j, nf) == su(j)) x(j) = xu(j)
end do


!-------------------------------------------------------------------!
!call calfun(n, x, f)
call evaluate(calfun, x, f)
call savehist(nf, x, xhist, f, fhist)
!-------------------------------------------------------------------!


fval(nf) = f
if (nf == 1) then
    fbeg = f
    kopt = 1
else if (f < fval(kopt)) then
    kopt = nf
end if
!
!     Set the nonZERO initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
if (nf <= 2 * n + 1) then
    if (nf >= 2 .and. nf <= n + 1) then
        gopt(nfm) = (f - fbeg) / stepa
        if (npt < nf + n) then
            bmat(nfm, 1) = -ONE / stepa
            bmat(nfm, nf) = ONE / stepa
            bmat(nfm, npt + nfm) = -HALF * rhosq
        end if
    else if (nf >= n + 2) then
        ih = (nfx * (nfx + 1)) / 2
        temp = (f - fbeg) / stepb
        diff = stepb - stepa
        hq(ih) = TWO * (temp - gopt(nfx)) / diff
        gopt(nfx) = (gopt(nfx) * stepb - temp * stepa) / diff
        if (stepa * stepb < ZERO) then
            if (f < fval(nf - n)) then
                fval(nf) = fval(nf - n)
                fval(nf - n) = f
                if (kopt == nf) kopt = nf - n
                xpt(nfx, nf - n) = stepb
                xpt(nfx, nf) = stepa
            end if
        end if
        bmat(nfx, 1) = -(stepa + stepb) / (stepa * stepb)
        bmat(nfx, nf) = -HALF / xpt(nfx, nf - n)
        bmat(nfx, nf - n) = -bmat(nfx, 1) - bmat(nfx, nf)
        zmat(1, nfx) = sqrt(TWO) / (stepa * stepb)
        zmat(nf, nfx) = sqrt(HALF) / rhosq
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
    temp = xpt(ipt, nf) * xpt(jpt, nf)
    hq(ih) = (fbeg - fval(ipt + 1) - fval(jpt + 1) + f) / temp
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 04-06-2019):
!     If the evaluation returns an NaN or an infinity value, this
!     subroutine is stopped.
if (is_nan(f) .or. is_posinf(f)) goto 80
!     By Tom (on 04-06-2019):
!     If the target value is reached, stop the algorithm.
if (f <= ftarget) goto 80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nf < npt) goto 50
80 nf = min(nf, npt)  ! nf = npt + 1 at exit of the loop

end subroutine initialize


end module initialize_mod
