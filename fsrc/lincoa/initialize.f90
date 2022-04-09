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
! Last Modified: Sunday, April 10, 2022 AM01:35:36
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initialize


contains


subroutine initialize(calfun, iprint, A_orig, amat, b_orig, ftarget, rhobeg, x0, b, &
    & idz, kopt, nf, bmat, chist, cstrv, f, fhist, fval, gopt, hq, pq, rescon, xsxpt, &
    & step, vlag, xbase, xhist, xopt, xpt, xsav, zmat)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : matprod, maximum
use, non_intrinsic :: pintrf_mod, only : OBJ

! Solver-specific modules
use, non_intrinsic :: update_mod, only : update

! Development modules (to be removed)
use, non_intrinsic :: ieee_4dev_mod, only : ieeenan

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
real(RP), intent(in) :: A_orig(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: b_orig(:)  ! B_ORIG(M)
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)  ! X0(N)

! In-outputs
real(RP), intent(inout) :: b(:)  ! B(M)

! Outputs
integer(IK), intent(out) :: idz
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(out) :: chist(:)  ! CHIST(MAXCHIST)
real(RP), intent(out) :: cstrv
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out) :: fval(:)  ! FVAL(NPT)
real(RP), intent(out) :: gopt(:)  ! GOPT(N)
real(RP), intent(out) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(out) :: pq(:)  ! PQ(NPT)
real(RP), intent(out) :: rescon(:)  ! RESCON(M)
real(RP), intent(out) :: xsxpt(:)  ! RXSXPT(2*NPT)
real(RP), intent(out) :: step(:)  ! STEP(N)
real(RP), intent(out) :: vlag(:)  ! VLAG(NPT+N) The size is NPT + N instead of NPT
real(RP), intent(out) :: xbase(:)  ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xopt(:)  ! XOPT(N)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(out) :: xsav(:)  ! XSAV(N)  ; necessary ???
real(RP), intent(out) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Local variables
character(len=*), parameter :: srname = 'INITIIALIZE'
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
real(RP) :: w(size(zmat, 2))
real(RP) :: x(size(x0))
real(RP) :: bigv, feas, recip, reciq, resid, rhosq, temp, test
integer(IK) :: i, ipt, itemp, j, jp, jpt, jsav, k, kbase, knew


! Sizes.
m = int(size(b), kind(m))
n = int(size(x), kind(n))
npt = int(size(fval), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = int(max(maxxhist, maxfhist, maxchist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(size(A_orig, 1) == n .and. size(A_orig, 2) == m, 'SIZE(A_ORIG) == [N, M]', srname)
    call assert(size(b_orig) == m, 'SIZE(B_ORIG) == M', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    call assert(size(xsxpt) == 2_IK * npt, 'SIZE(RXSXPT) == 2*NPT', srname)
    call assert(size(step) == n, 'SIZE(STEP) == N', srname)
    call assert(size(vlag) == npt + n, 'SIZE(VLAG) == NPT+N', srname)
    call assert(size(xbase) == n, 'SIZE(XBASE) == N', srname)
    call assert(size(xopt) == n, 'SIZE(XOPT) == N', srname)
    call assert(size(xsav) == n, 'SIZE(XSAV) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(hq, 1) == n .and. size(hq, 2) == n, 'SIZE(HQ) = [N, N]', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
end if


!*++
!*++ End of declarations rewritten by XSXPTAG
!
!     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
!       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, XSXPT and RESCON are the
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
!       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, XSXPT and RESCON.
!
!     Set some constants.
!
!--------------------------------------------------------------------------------------------------!
! Temporary fix for uninitialized variables when initialization terminates because of f < ftarget etc
bmat = ieeenan(); zmat = ieeenan(); fval = ieeenan()
!--------------------------------------------------------------------------------------------------!

rhosq = rhobeg * rhobeg
recip = ONE / rhosq
reciq = sqrt(HALF) / rhosq
test = 0.2_RP * rhobeg
idz = 1
kbase = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the initial elements of XPT, BMAT, XSXPT and ZMAT to ZERO.
!
do j = 1, n
    xbase(j) = x0(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xopt(j) = ZERO
    xsav(j) = xbase(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1, npt
        xpt(j, k) = ZERO
    end do
    do i = 1, npt + n
        bmat(j, i) = ZERO
    end do
end do
do k = 1, npt
    xsxpt(k) = ZERO
    do j = 1, npt - n - 1
        zmat(k, j) = ZERO
    end do
end do
!
!     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
!       but they may be altered later to make a constraint violation
!       sufficiently large. The initial nonzero elements of BMAT and of
!       the first min[N,NPT-N-1] columns of ZMAT are set also.
!
do j = 1, n
    xpt(j, j + 1) = rhobeg
    if (j < npt - n) then
        jp = n + j + 1
        xpt(j, jp) = -rhobeg
        bmat(j, j + 1) = HALF / rhobeg
        bmat(j, jp) = -HALF / rhobeg
        zmat(1, j) = -reciq - reciq
        zmat(j + 1, j) = reciq
        zmat(jp, j) = reciq
    else
        bmat(j, 1) = -ONE / rhobeg
        bmat(j, j + 1) = ONE / rhobeg
        bmat(j, npt + j) = -HALF * rhosq
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
        xpt(ipt, n + k + 1) = rhobeg
        xpt(jpt, n + k + 1) = rhobeg
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
!--------------------------------------------------------------------------------------------------!
jsav = 1_IK  ! Temporary fix for attention: «jsav» may be used uninitialized in this function from g95
!--------------------------------------------------------------------------------------------------!
do nf = 1, npt
    feas = ONE
    bigv = ZERO
    j = 0
80  j = j + 1
    if (j <= m .and. nf >= 2) then
        resid = -b(j)
        do i = 1, n
            resid = resid + xpt(i, nf) * amat(i, j)
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
            step(i) = xpt(i, nf) + (test - bigv) * amat(i, jsav)
        end do
        do k = 1, npt
            xsxpt(npt + k) = ZERO
            do j = 1, n
                xsxpt(npt + k) = xsxpt(npt + k) + xpt(j, k) * step(j)
            end do
        end do
        knew = nf

        call update(kbase, step, xpt, idz, knew, bmat, zmat, vlag)

        do i = 1, n
            xpt(i, nf) = step(i)
        end do
    end if
!
!     Calculate the objective function at the current interpolation point,
!       and set KOPT to the index of the first trust region centre.
!
    do j = 1, n
        x(j) = xbase(j) + xpt(j, nf)
    end do
    f = feas
    !---------------------------------------------------!
    call evaluate(calfun, x, f)  ! What if X contains NaN?
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
do j = 1, npt - n - 1
    w(j) = ZERO
    do k = 1, npt
        w(j) = w(j) + zmat(k, j) * fval(k)
    end do
end do
do k = 1, npt
    pq(k) = ZERO
    do j = 1, npt - n - 1
        pq(k) = pq(k) + zmat(k, j) * w(j)
    end do
end do
!
!     Set XOPT, XSXPT, GOPT and HQ for the first quadratic model.
!
do j = 1, n
    xopt(j) = xpt(j, kopt)
    xsav(j) = xbase(j) + xopt(j)
    gopt(j) = ZERO
end do
do k = 1, npt
    xsxpt(k) = ZERO
    do j = 1, n
        xsxpt(k) = xsxpt(k) + xpt(j, k) * xopt(j)
    end do
    temp = pq(k) * xsxpt(k)
    do j = 1, n
        gopt(j) = gopt(j) + fval(k) * bmat(j, k) + temp * xpt(j, k)
    end do
end do

hq = ZERO

!!
!!     Set the initial elements of RESCON.
!!
!do j = 1, m
!    temp = b(j)
!    do i = 1, n
!        temp = temp - xopt(i) * amat(i, j)
!    end do
!    temp = max(temp, ZERO)
!    if (temp >= rhobeg) temp = -temp
!    rescon(j) = temp
!end do

!
!     Set the initial elements of RESCON.
!
!     RESCON holds useful information about the constraint residuals. Every
!       nonnegative RESCON(J) is the residual of the J-th constraint at the
!       current trust region centre. Otherwise, if RESCON(J) is negative, the
!       J-th constraint holds as a strict inequality at the trust region
!       centre, its residual being at least |RESCON(J)|; further, the value
!       of |RESCON(J)| is at least the current trust region radius DELTA.
!
! 1. Normally, RESCON = B - AMAT^T*XOPT (theoretically, B - AMAT^T*XOPT >= 0 since XOPT is feasible)
! 2. If RESCON(J) >= DELTA (current trust-region radius), its sign is flipped: RESCON(J) = -RESCON(J).
rescon = max(b - matprod(xopt, amat), ZERO)  ! Calculation changed
where (rescon >= rhobeg)
    rescon = -rescon
end where
!!MATLAB: rescon(rescon >= rhobeg) = -rescon


end subroutine initialize


end module initialize_mod
