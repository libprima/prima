module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Friday, February 25, 2022 PM11:21:10
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initialize


contains


subroutine initialize(calfun, n, npt, m, amat, b, x, rhobeg, iprint, xbase, &
     &  xpt, fval, xsav, xopt, gopt, kopt, hq, pq, bmat, zmat, idz, ndim, &
     &  rsp, rescon, step, pqw, w, f, ftarget, &
     &  A_orig, b_orig, &
     & cstrv, nf, xhist, maxxhist, fhist, maxfhist, chist, maxchist)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm, maximum
use, non_intrinsic :: pintrf_mod, only : OBJ

! Solver-specific modules
use, non_intrinsic :: update_mod, only : update

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: m
integer(IK), intent(in) :: maxchist
integer(IK), intent(in) :: maxfhist
integer(IK), intent(in) :: maxxhist
integer(IK), intent(in) :: n
integer(IK), intent(in) :: ndim
integer(IK), intent(in) :: npt
real(RP), intent(in) :: A_orig(n, m)
real(RP), intent(in) :: amat(n, m)
real(RP), intent(in) :: b_orig(m)
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg

! In-outputs
integer(IK), intent(inout) :: idz
integer(IK), intent(inout) :: kopt
real(RP), intent(inout) :: b(m)
real(RP), intent(inout) :: f
real(RP), intent(inout) :: fval(npt)
real(RP), intent(inout) :: gopt(n)
real(RP), intent(inout) :: pq(npt)
real(RP), intent(inout) :: pqw(npt)
real(RP), intent(inout) :: rsp(2_IK * npt)
real(RP), intent(inout) :: step(n)
real(RP), intent(inout) :: w(n + npt)
real(RP), intent(inout) :: x(n)
real(RP), intent(inout) :: xbase(n)
real(RP), intent(inout) :: xopt(n)
real(RP), intent(inout) :: xpt(npt, n)
real(RP), intent(inout) :: bmat(npt + size(x), size(x))
real(RP), intent(inout) :: zmat(npt, npt - size(x) - 1)

! Outputs
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(maxchist)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: fhist(maxfhist)
real(RP), intent(out) :: hq(n * (n + 1_IK) / 2_IK)
real(RP), intent(out) :: rescon(m)
real(RP), intent(out) :: xhist(n, maxxhist)
real(RP), intent(out) :: xsav(n)

! Local variables
real(RP) :: bigv, feas, recip, reciq, resid, rhosq, temp, test
integer(IK) :: i, ipt, itemp, j, jp, jpt, jsav, k, kbase, knew, nptm


!*++
!*++ End of declarations rewritten by SPAG
!
!     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
!       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
!       same as the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
!       region centre.
!     IDZ is going to be set to ONE, so that every element of Diag(DZ) is
!       ONE in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
!       factorization of the leading NPT by NPT submatrix of H.
!     STEP, PQW and W are used for working space, the arrays STEP and PQW
!       being taken from LINCOB. The length of W must be at least N+NPT.
!
!     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
!       for the first iteration, an important feature being that, if any of
!       of the columns of XPT is an infeasible point, then the largest of
!       the constraint violations there is at least 0.2*RHOBEG. It also sets
!       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
!
!     Set some constants.
!
nptm = npt - n - 1
rhosq = rhobeg * rhobeg
recip = ONE / rhosq
reciq = sqrt(HALF) / rhosq
test = 0.2D0 * rhobeg
idz = 1
kbase = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the initial elements of XPT, BMAT, SP and ZMAT to ZERO.
!
do j = 1, n
    xbase(j) = x(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xopt(j) = ZERO
    xsav(j) = xbase(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1, npt
        xpt(k, j) = ZERO
    end do
    do i = 1, ndim
        bmat(i, j) = ZERO
    end do
end do
do k = 1, npt
    rsp(k) = ZERO
    do j = 1, npt - n - 1
        zmat(k, j) = ZERO
    end do
end do
!
!     Set the nonZERO coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
!       but they may be altered later to make a constraint violation
!       sufficiently large. The initial nonZERO elements of BMAT and of
!       the first min[N,NPT-N-1] columns of ZMAT are set also.
!
do j = 1, n
    xpt(j + 1, j) = rhobeg
    if (j < npt - n) then
        jp = n + j + 1
        xpt(jp, j) = -rhobeg
        bmat(j + 1, j) = HALF / rhobeg
        bmat(jp, j) = -HALF / rhobeg
        zmat(1, j) = -reciq - reciq
        zmat(j + 1, j) = reciq
        zmat(jp, j) = reciq
    else
        bmat(1, j) = -ONE / rhobeg
        bmat(j + 1, j) = ONE / rhobeg
        bmat(npt + j, j) = -HALF * rhosq
    end if
end do
!
!     Set the remaining initial nonZERO elements of XPT and ZMAT when the
!       number of interpolation points exceeds 2*N+1.
!
if (npt > 2 * n + 1) then
    do k = n + 1, npt - n - 1
        itemp = (k - 1) / n
        ipt = k - itemp * n
        jpt = ipt + itemp
        if (jpt > n) jpt = jpt - n
        xpt(n + k + 1, ipt) = rhobeg
        xpt(n + k + 1, jpt) = rhobeg
        zmat(1, k) = recip
        zmat(ipt + 1, k) = -recip
        zmat(jpt + 1, k) = -recip
        zmat(n + k + 1, k) = recip
    end do
end if
!
!     Update the constraint right hand sides to allow for the shift XBASE.
!
if (m > 0) then
    do j = 1, m
        temp = ZERO
        do i = 1, n
            temp = temp + amat(i, j) * xbase(i)
        end do
        b(j) = b(j) - temp
    end do
end if
!
!     Go through the initial points, shifting every infeasible point if
!       necessary so that its constraint violation is at least 0.2*RHOBEG.
!
do nf = 1, npt
    feas = ONE
    bigv = ZERO
    j = 0
80  j = j + 1
    if (j <= m .and. nf >= 2) then
        resid = -b(j)
        do i = 1, n
            resid = resid + xpt(nf, i) * amat(i, j)
        end do
        if (resid <= bigv) goto 80
        bigv = resid
        jsav = j
        if (resid <= test) then
            feas = -ONE
            goto 80
        end if
        feas = ZERO
    end if
    if (feas < ZERO) then
        do i = 1, n
            step(i) = xpt(nf, i) + (test - bigv) * amat(i, jsav)
        end do
        do k = 1, npt
            rsp(npt + k) = ZERO
            do j = 1, n
                rsp(npt + k) = rsp(npt + k) + xpt(k, j) * step(j)
            end do
        end do
        knew = nf
        call update(n, npt, xpt, bmat, zmat, idz, ndim, rsp, step, kbase, knew, pqw, w)
        do i = 1, n
            xpt(nf, i) = step(i)
        end do
    end if
!
!     Calculate the objective function at the current interpolation point,
!       and set KOPT to the index of the first trust region centre.
!
    do j = 1, n
        x(j) = xbase(j) + xpt(nf, j)
    end do
    f = feas
    !---------------------------------------------------!
    call evaluate(calfun, x, f)
    cstrv = maximum([ZERO, matprod(x, A_orig) - b_orig])
    call savehist(nf, x, xhist, f, fhist, cstrv, chist)
    !---------------------------------------------------!
    if (nf == 1) then
        kopt = 1
    else if (f < fval(kopt) .and. feas > ZERO) then
        kopt = nf
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  150 FVAL(NF)=F
    fval(nf) = f
!     By Tom/Zaikun (on 03-06-2019/07-06-2019):
!     If the objective function reached a NaN or infinite value, or if
!     the value is under the target value, the algorithm go back to
!     LINCOB with updated KOPT and XSAV.
!     Note that we should NOT compare F and FTARGET, because X may not
!     be feasible.
    if (is_nan(f) .or. is_posinf(f) .or. fval(kopt) <= ftarget) then
        exit
    end if
end do
!----------------------------------------------------------!
nf = min(nf, npt)  ! At exit of the loop, nf = npt + 1
!----------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set PQ for the first quadratic model.
!
do j = 1, nptm
    w(j) = ZERO
    do k = 1, npt
        w(j) = w(j) + zmat(k, j) * fval(k)
    end do
end do
do k = 1, npt
    pq(k) = ZERO
    do j = 1, nptm
        pq(k) = pq(k) + zmat(k, j) * w(j)
    end do
end do
!
!     Set XOPT, SP, GOPT and HQ for the first quadratic model.
!
do j = 1, n
    xopt(j) = xpt(kopt, j)
    xsav(j) = xbase(j) + xopt(j)
    gopt(j) = ZERO
end do
do k = 1, npt
    rsp(k) = ZERO
    do j = 1, n
        rsp(k) = rsp(k) + xpt(k, j) * xopt(j)
    end do
    temp = pq(k) * rsp(k)
    do j = 1, n
        gopt(j) = gopt(j) + fval(k) * bmat(k, j) + temp * xpt(k, j)
    end do
end do
do i = 1, (n * n + n) / 2
    hq(i) = ZERO
end do
!
!     Set the initial elements of RESCON.
!
do j = 1, m
    temp = b(j)
    do i = 1, n
        temp = temp - xopt(i) * amat(i, j)
    end do
    temp = max(temp, ZERO)
    if (temp >= rhobeg) temp = -temp
    rescon(j) = temp
end do


end subroutine initialize


end module initialize_mod
