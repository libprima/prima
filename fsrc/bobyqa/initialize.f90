module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, June 06, 2022 AM12:36:23
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initialize


contains


subroutine initialize(calfun, iprint, maxfun, ftarget, rhobeg, sl, su, x0, xl, xu, &
    & kopt, nf, bmat, f, fhist, fval, gopt, hq, pq, xbase, xhist, xpt, zmat, info)

use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_INF_F, FTARGET_ACHIEVED
use, non_intrinsic :: linalg_mod, only : matprod, trueloc, issymmetric!, norm
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: sl(:)  ! SL(N)
real(RP), intent(in) :: su(:)  ! SU(N)
real(RP), intent(in) :: x0(:)  ! X(N)
real(RP), intent(in) :: xl(:)  ! XL(N)
real(RP), intent(in) :: xu(:)  ! XU(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: fval(:)  ! FVAL(NPT)
real(RP), intent(out) :: gopt(:)  ! GOPT(N)
real(RP), intent(out) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)
real(RP), intent(out) :: xbase(:)  ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(out) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)


! Local variables
character(len=*), parameter :: solver = 'NEWUOA'
character(len=*), parameter :: srname = 'INITIIALIZE'
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: subinfo
real(RP) :: x(size(xpt, 1))
real(RP) :: fbase, recip, rhosq, stepa, stepb
integer(IK) :: ipt(size(xpt, 2)), itemp, jpt(size(xpt, 2)), k
logical :: evaluated(size(xpt, 2))

! Sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(size(fval) == npt, 'SIZE(FVAL) == NPT', srname)
    call assert(size(sl) == n .and. all(sl <= 0), 'SIZE(SL) == N, SL <= 0', srname)
    call assert(size(su) == n .and. all(su >= 0), 'SIZE(SU) == N, SU >= 0', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(size(x0) == n .and. all(is_finite(x0)), 'SIZE(X0) == N, X0 is finite', srname)
    call assert(all(x0 >= xl .and. x0 <= xu), 'XL <= X0 <= XU', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) == [N, N]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
end if

!--------------------------------------------------------------------------------------------------!
ipt = 1; jpt = 1  ! Temporary fix for G95 warning about these variables used uninitialized
!--------------------------------------------------------------------------------------------------!
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

info = INFO_DFT

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously.
evaluated = .false.
! Initialize FVAL to HUGENUM. Otherwise, compilers may complain that FVAL is not (completely)
! initialized if the initialization aborts due to abnormality (see CHECKEXIT).
fval = HUGENUM
rhosq = rhobeg * rhobeg
recip = ONE / rhosq

! Set XBASE to the initial vector of variables, and set the initial elements of XPT, BMAT, HQ, PQ
! and ZMAT to ZERO.
xbase = x0
xpt = ZERO
hq = ZERO
pq = ZERO
bmat = ZERO
zmat = ZERO


! Set XPT(:, 1 : 2*N + 1).
xpt(:, 1) = ZERO
do k = 1, min(npt, int(2 * n + 1, kind(npt)))
    if (k >= 2 .and. k <= n + 1) then
        xpt(k - 1, k) = rhobeg
        if (su(k - 1) <= 0) then  ! SU(NF - 1) == 0
            xpt(k - 1, k) = -rhobeg
        end if
    else if (k >= n + 2) then
        xpt(k - n - 1, k) = -rhobeg
        if (sl(k - n - 1) >= 0) then  ! SL(NF - N - 1) == 0
            xpt(k - n - 1, k) = min(TWO * rhobeg, su(k - n - 1))
        end if
        if (su(k - n - 1) <= 0) then  ! SU(NF - N - 1) == 0
            xpt(k - n - 1, k) = max(-TWO * rhobeg, sl(k - n - 1))
        end if
    end if
end do

! Set FVAL(1 : 2*N + 1) by evaluating F. Totally parallelizable except for FMSG.
do k = 1, min(npt, int(2 * n + 1, kind(npt)))
    x = min(max(xl, xbase + xpt(:, k)), xu)
    x(trueloc(xpt(:, k) <= sl)) = xl(trueloc(xpt(:, k) <= sl))
    x(trueloc(xpt(:, k) >= su)) = xu(trueloc(xpt(:, k) >= su))
    call evaluate(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call savehist(k, x, xhist, f, fhist)
    call fmsg(solver, iprint, k, f, x)
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do

! For each K between 2 and N + 1, switch XPT(:, K) and XPT(:, K+1) if XPT(K-1, K) * XPT(K-1, K+N) is
! negative and FVAL(K) <= FVAL(K+N).
do k = 2, min(npt - n, int(n + 1, kind(npt)))
    if (xpt(k - 1, k) * xpt(k - 1, k + n) < 0 .and. fval(k + n) < fval(k)) then
        fval([k, k + n]) = fval([k + n, k])
        xpt(:, [k, k + n]) = xpt(:, [k + n, k])
        ! Indeed, only XPT(K-1, [K, K+N]) needs switching, as the other entries are zero.
    end if
end do

if (info == INFO_DFT) then
    do k = int(2 * n + 2, kind(k)), npt
        itemp = (k - n - 2) / n
        jpt(k) = k - (itemp + 1) * n - 1
        ipt(k) = jpt(k) + itemp
        if (ipt(k) > n) then
            itemp = jpt(k)
            jpt(k) = ipt(k) - n
            ipt(k) = itemp
        end if

        xpt(ipt(k), k) = xpt(ipt(k), ipt(k) + 1)
        xpt(jpt(k), k) = xpt(jpt(k), jpt(k) + 1)

        x = min(max(xl, xbase + xpt(:, k)), xu)
        x(trueloc(xpt(:, k) <= sl)) = xl(trueloc(xpt(:, k) <= sl))
        x(trueloc(xpt(:, k) >= su)) = xu(trueloc(xpt(:, k) >= su))
        call evaluate(calfun, x, f)
        call savehist(k, x, xhist, f, fhist)
        evaluated(k) = .true.
        fval(k) = f
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end do
end if

nf = count(evaluated)
kopt = int(minloc(fval, mask=evaluated, dim=1), kind(kopt))

! Set the nonzero initial elements of BMAT and the quadratic model in the cases when NF is at most
! 2*N+1. If NF exceeds N+1, then the positions of the NF-th and (NF-N)-th interpolation points may
! be switched, in order that the function value at the first of them contributes to the off-diagonal
! second derivative terms of the initial quadratic model.
if (all(evaluated)) then
    fbase = fval(1)
    do k = 1, npt
        if (k <= 2 * n + 1) then
            if (k >= 2 .and. k <= n + 1) then
                gopt(k - 1) = (fval(k) - fbase) / xpt(k - 1, k)
                if (npt < k + n) then
                    bmat(k - 1, 1) = -ONE / xpt(k - 1, k)
                    bmat(k - 1, k) = ONE / xpt(k - 1, k)
                    bmat(k - 1, npt + k - 1) = -HALF * rhosq
                end if
            else if (k >= n + 2) then
                stepa = xpt(k - n - 1, k - n)
                stepb = xpt(k - n - 1, k)
                hq(k - n - 1, k - n - 1) = TWO * ((fval(k) - fbase) / stepb - (fval(k - n) - fbase) / stepa) / (stepb - stepa)
                gopt(k - n - 1) = (((fval(k - n) - fbase) / stepa) * stepb - ((fval(k) - fbase) / stepb) * stepa) / (stepb - stepa)
                bmat(k - n - 1, 1) = -(stepa + stepb) / (stepa * stepb)
                bmat(k - n - 1, k) = -HALF / xpt(k - n - 1, k - n)
                bmat(k - n - 1, k - n) = -bmat(k - n - 1, 1) - bmat(k - n - 1, k)
                zmat(1, k - n - 1) = sqrt(TWO) / (stepa * stepb)
                zmat(k, k - n - 1) = sqrt(HALF) / rhosq
                zmat(k - n, k - n - 1) = -zmat(1, k - n - 1) - zmat(k, k - n - 1)
            end if
! Set the off-diagonal second derivatives of the Lagrange functions and the initial quadratic model.
        else
            zmat(1, k - n - 1) = recip
            zmat(k, k - n - 1) = recip
            zmat(ipt(k) + 1, k - n - 1) = -recip
            zmat(jpt(k) + 1, k - n - 1) = -recip
            hq(ipt(k), jpt(k)) = (fbase - fval(ipt(k) + 1) - fval(jpt(k) + 1) + fval(k)) / (xpt(ipt(k), k) * xpt(jpt(k), k))
            hq(jpt(k), ipt(k)) = hq(ipt(k), jpt(k))
        end if
    end do

    if (kopt /= 1) then
        gopt = gopt + matprod(hq, xpt(:, kopt))
    end if
end if

!write (16, *) nf, fval
!write (16, *) xpt
!write (16, *) bmat
!write (16, *) zmat
!write (16, *) gopt, hq

! Postconditions
if (DEBUGGING) then
    call assert(nf <= npt, 'NF <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(all(xbase >= xl .and. xbase <= xu), 'XL <= XBASE <= XU', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt)) .and. &
        & all(xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    call assert(size(fval) == npt .and. .not. any(evaluated .and. (is_nan(fval) .or. is_posinf(fval))), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN or +Inf', srname)
    call assert(.not. any(fval < fval(kopt) .and. evaluated), 'FVAL(KOPT) = MINVAL(FVAL)', srname)
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
end if

end subroutine initialize


end module initialize_mod
