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
! Last Modified: Friday, April 15, 2022 PM10:01:07
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(n, npt, m, amat, xpt, xopt, nact, iact, &
     &  rescon, qfac, kopt, knew, del, step, gl, pqw, rstat, w, ifeas)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TENTH
use, non_intrinsic :: linalg_mod, only : maximum
use, non_intrinsic :: infnan_mod, only : is_nan

implicit none

! Inputs
integer(IK), intent(in) :: iact(:)
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n
integer(IK), intent(in) :: nact
integer(IK), intent(in) :: npt
real(RP), intent(in) :: amat(n, m)
real(RP), intent(in) :: del
real(RP), intent(in) :: pqw(npt)
real(RP), intent(in) :: qfac(n, n)
real(RP), intent(in) :: rescon(m)
real(RP), intent(in) :: xopt(n)
real(RP), intent(in) :: xpt(n, npt)

! In-outputs
integer(IK), intent(inout) :: ifeas
real(RP), intent(inout) :: gl(n)
real(RP), intent(inout) :: rstat(m)
real(RP), intent(inout) :: step(n)
real(RP), intent(inout) :: w(n)

! Local variables
real(RP) :: bigv, ctol, gg, ghg, resmax, sp, ss,  &
&        stp, stpsav, summ, temp, test, vbig, vgrad, &
&        vlag, vnew, ww
integer(IK) :: i, j, jsav, k, ksav


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
!       if its value of | LFUNC(XOPT+STEP) | is at least ONE fifth of the
!       original ONE, but the greatest violation of a linear constraint must
!       be at least 0.2*DEL, in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes 0.2*DEL.
!
!     Set some constants.
!
test = 0.2_RP * del  ! Is this really better than 0? According to an experiment of Tom on 20220225, NO
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
do k = 1, npt
    temp = ZERO
    do j = 1, n
        temp = temp + xpt(j, k) * xopt(j)
    end do
    temp = pqw(k) * temp
    do i = 1, n
        gl(i) = gl(i) + temp * xpt(i, k)
    end do
end do
if (m > 0) then
    do j = 1, m
        rstat(j) = ONE
        if (abs(rescon(j)) >= del) rstat(j) = -ONE
    end do
    do k = 1, nact
        rstat(iact(k)) = ZERO
    end do
end if
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-15: IFLAG is never used
!      IFLAG=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------!
! Zaikun 20220415
!vbig = ZERO
ss = ZERO
sp = ZERO
do i = 1, n
    temp = xpt(i, knew) - xopt(i)
    !ss = ss + temp * temp
    ss = ss + temp**2
    sp = sp + gl(i) * temp
end do
stp = -del / sqrt(ss)
if (sp * (sp - ONE) < ZERO) stp = -stp
vlag = abs(stp * sp) + stp**2 * abs(sp - ONE)
vbig = vlag
ksav = knew
stpsav = stp
!-----------------------------------------------------------!
do k = 1, npt
    if (k == kopt) cycle
    ss = ZERO
    sp = ZERO
    do i = 1, n
        temp = xpt(i, k) - xopt(i)
        !ss = ss + temp * temp
        ss = ss + temp**2
        sp = sp + gl(i) * temp
    end do
    stp = -del / sqrt(ss)
    if (k == knew) then
        if (sp * (sp - ONE) < ZERO) stp = -stp
        vlag = abs(stp * sp) + stp**2 * abs(sp - ONE)
    else
        vlag = abs(stp * (ONE - stp) * sp)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: With the original code, if either VLAG or VBIG is
! NaN, KSAV will not get a value. This may cause Segmentation Fault
! because XPT(KSAV, :) will later be accessed.
!      IF (VLAG .GT. VBIG) THEN
!    if (.not. (vlag <= vbig)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if (vlag > vbig .and. vlag > 0) then
    if (vlag > vbig) then
        ksav = k
        stpsav = stp
        vbig = vlag
    end if
end do
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
!if (vnew > vbig) then
!--------------------------------------------------------------------!
! Zaikun 20220415
if (vnew > vbig .or. (is_nan(vbig) .and. .not. is_nan(vnew))) then
!--------------------------------------------------------------------!
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
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
if (vnew / vbig >= 0.2_RP) then
    ifeas = 1
    bigv = ZERO
    j = 0
170 j = j + 1
    if (j <= m) then
        if (rstat(j) == ONE) then
            !temp = -rescon(j)
            !do i = 1, n
            !    temp = temp + w(i) * amat(i, j)
            !end do
            temp = ZERO
            do i = 1, n
                temp = temp + w(i) * amat(i, j)
            end do
            temp = temp - rescon(j)
            !-------------------------------------!
            ! Zaikun 20220415
            !bigv = max(bigv, temp)
            bigv = maximum([bigv, temp])
            !-------------------------------------!
        end if
        if (bigv < test) goto 170
        ifeas = 0
    end if
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
    if (bigv <= 10.0_RP * ctol .or. bigv >= test) then
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
    write (17, *) j, jsav, ifeas, bigv, resmax
    if (rstat(j) < 0) goto 230

    !temp = -rescon(j)
    temp = ZERO
    do i = 1, n
        temp = temp + step(i) * amat(i, j)
    end do
    temp = temp - rescon(j)

    resmax = max(resmax, temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (temp < test) then
!    if (.not. (temp >= test)) then
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
        step(i) = step(i) + (test - bigv) * amat(i, jsav)
    end do
    ifeas = 0
end if
!
!     Return the calculated STEP and the value of IFEAS.
!
260 return
end subroutine geostep


end module geometry_mod
