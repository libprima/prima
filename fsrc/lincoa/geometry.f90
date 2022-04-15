module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
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
! Last Modified: Friday, April 15, 2022 PM02:56:14
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(iact, knew, kopt, nact, amat, del, gl_in, pqw, qfac, rescon, xopt, xpt, ifeas, step)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TEN, TENTH, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : inprod, isorth

implicit none

! Inputs
integer(IK), intent(in) :: iact(:)  ! IACT(M)
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
integer(IK), intent(in) :: nact
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: del
real(RP), intent(in) :: gl_in(:)  ! GL_IN(N)
real(RP), intent(in) :: pqw(:)  ! PQW(NPT)  ; better name?
real(RP), intent(in) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(in) :: rescon(:)  ! RESCON(M)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! Outputs
integer(IK), intent(out) :: ifeas
real(RP), intent(out) :: step(:)  ! STEP(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: m
integer(IK) :: n
integer(IK) :: npt
real(RP) :: rstat(size(amat, 2))
real(RP) :: w(size(xopt))
real(RP) :: gl(size(gl_in))
real(RP) :: bigv, ctol, gg, ghg, resmax, sp, ss, tol, &
&        stp, stplen(size(pqw)), stpsav, summ, temp, mincv, vbig, vgrad, vlag(size(pqw)), vnew, ww
integer(IK) :: i, j, jsav, k, ksav

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(xopt), kind(n))
npt = int(size(pqw), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt, 'KNEW /= KOPT', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(del > 0, 'DEL > 0', srname)
    call assert(size(gl_in) == n, 'SIZE(GL_IN) == N', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
end if

ifeas = 0_IK !??? Added by Zaikun 20220227
gl = gl_in

!
!     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
!       same as the terms with these names in SUBROUTINE LINCOB.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DEL is the current restriction on the length of STEP, which is never
!       greater than the current trust region radius DELTA.
!     STEP will be set to the required step from XOPT to the new point.
!     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
!       is the KNEW-th Lagrange function. It is used also for some other
!       gradients of LFUNC.
!     PQW provides the second derivative parameters of LFUNC.
!     RSTAT and W are used for working space. Their lengths must be at least
!       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
!       J-th constraint is irrelevant, active, or both inactive and relevant,
!       respectively.
!     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
!
!     STEP is chosen to provide a relatively large value of the modulus of
!       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
!       calculated too, within the trust region, that does not alter the
!       residuals of the active constraints. The projected step is preferred
!       if its value of |LFUNC(XOPT+STEP)| is at least ONE fifth of the
!       original ONE, but the greatest violation of a linear constraint must
!       be at least MINCV = 0.2*DEL, in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes MINCV.
!
!     Set some constants.
!
!--------------------------------------------------------------------------------------------------!
ksav = 1_IK; jsav = 1_IK  ! Temporary fix for ksav/jsav may be uninitialized from G95
!--------------------------------------------------------------------------------------------------!

mincv = 0.2_RP * del  ! Is this really better than 0? According to an experiment of Tom on 20220225, NO

! Replace GL by the gradient of LFUNC at the trust region centre, and set the elements of RSTAT.
do k = 1, npt
    !temp = ZERO
    !do j = 1, n
    !    temp = temp + xpt(j, k) * xopt(j)
    !end do
    !temp = pqw(k) * temp
    !do i = 1, n
    !    gl(i) = gl(i) + temp * xpt(i, k)
    !end do
    gl = gl + pqw(k) * inprod(xopt, xpt(:, k)) * xpt(:, k)
end do

rstat = ONE
where (abs(rescon) >= del)
    rstat = -ONE
end where
!!MATLAB: rstat(abs(rescon) >= del) = -1;
rstat(iact(1:nact)) = ZERO

! Maximize |LFUNC| within the trust region on the lines through XOPT and other interpolation points.
vlag = ZERO
do k = 1, npt
    if (k == kopt) then
        cycle
    end if
    ss = sum((xpt(:, k) - xopt)**2)
    sp = inprod(gl, xpt(:, k) - xopt)
    stplen(k) = -del / sqrt(ss)
    if (k == knew) then
        if (sp * (sp - ONE) < ZERO) stplen(k) = -stplen(k)
        vlag(k) = abs(stplen(k) * sp) + stplen(k)**2 * abs(sp - ONE)
    else
        vlag(k) = abs(stplen(k) * (ONE - stplen(k)) * sp)
    end if
end do

! N.B.: We define KSAV slightly differently from Powell's code, which sets KSAV = MAXLOC(VLAG,DIM=1)
! by comparing the entries of VLAG one by one.
! 1. If VLAG contains only NaN, which can happen, Powell's code leaves KSAV uninitialized.
! 2. If VLAG(KNEW) = MINVAL(VLAG) = VLAG(K) with K < KNEW, Powell's code does not set KSAV = KNEW.
ksav = knew
if (maxval(vlag) > vlag(knew)) then
    ksav = maxloc(vlag, dim=1)
end if
vbig = vlag(ksav)
stpsav = stplen(ksav)

!
!     Set STEP to the move that gives the greatest modulus calculated above.
!       This move may be replaced by a steepest ascent step from XOPT.
!
gg = ZERO
do i = 1, n
    gg = gg + gl(i)**2
    step(i) = stpsav * (xpt(i, ksav) - xopt(i))
end do
vgrad = del * sqrt(gg)
if (vgrad <= TENTH * vbig) goto 220
!
!     Make the replacement if it provides a larger value of VBIG.
!
ghg = ZERO
do k = 1, npt
    temp = ZERO
    do j = 1, n
        temp = temp + xpt(j, k) * gl(j)
    end do
    ghg = ghg + pqw(k) * temp * temp
end do
vnew = vgrad + abs(HALF * del * del * ghg / gg)
if (vnew > vbig .or. (is_nan(vbig) .and. .not. is_nan(vnew))) then
    vbig = vnew
    stp = del / sqrt(gg)
    if (ghg < ZERO) stp = -stp
    do i = 1, n
        step(i) = stp * gl(i)
    end do
end if
if (nact == 0 .or. nact == n) goto 220
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
do k = nact + 1, n
    w(k) = ZERO
    do i = 1, n
        w(k) = w(k) + gl(i) * qfac(i, k)
    end do
end do
gg = ZERO
do i = 1, n
    gl(i) = ZERO
    do k = nact + 1, n
        gl(i) = gl(i) + qfac(i, k) * w(k)
    end do
    gg = gg + gl(i)**2
end do
vgrad = del * sqrt(gg)
if (vgrad <= TENTH * vbig) goto 220
ghg = ZERO
do k = 1, npt
    temp = ZERO
    do j = 1, n
        temp = temp + xpt(j, k) * gl(j)
    end do
    ghg = ghg + pqw(k) * temp * temp
end do
vnew = vgrad + abs(HALF * del * del * ghg / gg)
!
!     Set W to the possible move along the projected gradient.
!
stp = del / sqrt(gg)
if (ghg < ZERO) stp = -stp
ww = ZERO
do i = 1, n
    w(i) = stp * gl(i)
    ww = ww + w(i)**2
end do
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least MINCV. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
if (vnew / vbig >= 0.2_RP) then
    ifeas = 1
    bigv = ZERO
    !j = 0
!170 j = j + 1
    !if (j <= m) then
    !    if (rstat(j) == ONE) then
    !        temp = -rescon(j)
    !        do i = 1, n
    !            temp = temp + w(i) * amat(i, j)
    !        end do
    !        bigv = max(bigv, temp)
    !    end if
    !    if (bigv < mincv) goto 170
    !    ifeas = 0
    !end if
    do j = 1, m
        if (rstat(j) == ONE) then
            !temp = -rescon(j)
            !do i = 1, n
            !    temp = temp + w(i) * amat(i, j)
            !end do
            !bigv = max(bigv, temp)

            bigv = max(bigv, inprod(amat(:, j), w) - rescon(j))  ! Calculation changed
        end if
        if (.not. bigv < mincv) then
            ifeas = 0
            exit
        end if
    end do

    !bigv = max(ZERO, cummax(matrod(w, amat) - rescon, mask=(rstat == ONE)))

    ctol = ZERO
    temp = 0.01_RP * sqrt(ww)
    if (bigv > ZERO .and. bigv < temp) then
        do k = 1, nact
            j = iact(k)
            summ = ZERO
            do i = 1, n
                summ = summ + w(i) * amat(i, j)
            end do
            ctol = max(ctol, abs(summ))
        end do
    end if
    if (bigv <= TEN * ctol .or. bigv >= mincv) then
        do i = 1, n
            step(i) = w(i)
        end do
        goto 260
    end if
end if
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
220 ifeas = 1
bigv = ZERO
resmax = ZERO
j = 0
230 j = j + 1
if (j <= m) then
    if (rstat(j) < 0) goto 230
    temp = -rescon(j)
    do i = 1, n
        temp = temp + step(i) * amat(i, j)
    end do
    resmax = max(resmax, temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IF (TEMP .LT. TEST) THEN
    if (.not. (temp >= mincv)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (temp <= bigv) goto 230
        bigv = temp
        jsav = j
        ifeas = -1
        goto 230
    end if
    ifeas = 0
end if
if (ifeas == -1) then
    do i = 1, n
        step(i) = step(i) + (mincv - bigv) * amat(i, jsav)
    end do
    ifeas = 0
end if
!
!     Return the calculated STEP and the value of IFEAS.
!
260 return
end subroutine geostep


end module geometry_mod
