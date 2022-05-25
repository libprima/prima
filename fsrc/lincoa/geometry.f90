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
! Last Modified: Wednesday, May 25, 2022 PM11:16:07
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(iact, idz, knew, kopt, nact, amat, bmat, delbar, qfac, rescon, xopt, xpt, zmat, ifeas, step)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds STEP such that intends to improve the geometry of the interpolation set
! when XPT(:, KNEW) is changed to XOPT + STEP, where XOPT = XPT(:, KOPT).
!  STEP is chosen to provide a relatively large value of the modulus of
!       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DELBAR. A projected STEP is
!       calculated too, within the trust region, that does not alter the
!       residuals of the active constraints. The projected step is preferred
!       if its value of |LFUNC(XOPT+STEP)| is at least one fifth of the
!       original one but the greatest violation of a linear constraint must
!       be at least MINCV = 0.2*DELBAR in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes MINCV.
!     IFEAS will be set to TRUE or FALSE if XOPT+STEP is feasible or infeasible.
!
!     AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the same as the terms with these names
!     in subroutine LINCOB. KNEW is the index of the interpolation point that is going to be moved.
!     DELBAR is the current restriction on the length of STEP, which is never greater than the
!     current trust region radius DELTA.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TEN, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, isorth, maximum, trueloc, norm
use, non_intrinsic :: powalg_mod, only : hess_mul, omega_col, calden

implicit none

! Inputs
integer(IK), intent(in) :: iact(:)  ! IACT(M)
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
integer(IK), intent(in) :: nact
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: delbar
real(RP), intent(in) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(in) :: rescon(:)  ! RESCON(M)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
logical, intent(out) :: ifeas
real(RP), intent(out) :: step(:)  ! STEP(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: rstat(size(amat, 2))
real(RP) :: constr(size(amat, 2))
real(RP) :: cv_prjg
real(RP) :: cvtol
real(RP) :: cvtol_prjg
real(RP) :: den_grad(size(xpt, 2))
real(RP) :: den_line(size(xpt, 2))
real(RP) :: den_prjg(size(xpt, 2))
real(RP) :: denabs
real(RP) :: gg
real(RP) :: ghg
real(RP) :: gl(size(xpt, 1))
real(RP) :: gxpt(size(xpt, 2))
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: sp
real(RP) :: ss
real(RP) :: step_grad(size(xpt, 1))
real(RP) :: step_line(size(xpt, 1))
real(RP) :: step_prjg(size(xpt, 1))
real(RP) :: stp
real(RP) :: stplen(size(xpt, 2))
real(RP) :: tol
real(RP) :: vlag(size(xpt, 2))

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt, 'KNEW /= KOPT', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(delbar > 0, 'DELBAR> 0', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N and XOPT is finite', srname)
    call wassert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(abs(xopt - xpt(:, kopt)) <= 0), 'XOPT = XPT(:, KOPT)', srname) ! If this is true, XOPT should not be passed.
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, 'SIZE(ZMAT) == [NPT, NPT- N-1]', srname)
end if

gl = bmat(:, knew)

! RSTAT(J) = -1, 0, or 1 respectively means constraint J is irrelevant, active, or inactive&relevant.
! RSTAT never changes after being set below.
rstat = 1_IK
rstat(trueloc(abs(rescon) >= delbar)) = -1_IK
rstat(iact(1:nact)) = 0_IK

! Replace GL by the gradient of LFUNC at the trust region centre, and set the elements of RSTAT.
! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC.
pqlag = omega_col(idz, zmat, knew)
gl = gl + hess_mul(xopt, xpt, pqlag)

! Maximize |LFUNC| within the trust region on the lines through XOPT and other interpolation points.
! In the following loop, VLAG(K) is set to the maximum of |PHI_K(t)| subject to the trust-region
! constraint with PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:, K)). Note that PHI_K(t) is a quadratic
! function with PHI_K'(0) = SP = INPROD(GL, XPT(:, K) - XOPT) and PHI_K(0) = 0; PHI_K(1) = 0 if
! K /= KNEW and PHI_K(1) = 1 if K = KNEW. The linear constraints are not considered.
vlag = ZERO
do k = 1, npt
    if (k == kopt) then
        cycle
    end if
    sp = inprod(gl, xpt(:, k) - xopt)
    ss = sum((xpt(:, k) - xopt)**2)
    stplen(k) = -delbar / sqrt(ss)
    if (k == knew) then
        if (sp * (sp - ONE) < 0) then
            stplen(k) = -stplen(k)
        end if
        vlag(k) = abs(stplen(k) * sp) + stplen(k)**2 * abs(sp - ONE)
    else
        vlag(k) = abs(stplen(k) * (ONE - stplen(k)) * sp)
    end if
end do

! N.B.: 1. We define K in a way slightly different from Powell's code, which sets K to MAXLOC(VLAG)
! by comparing the entries of VLAG sequentially.
! 1. If VLAG contains only NaN, which can happen, Powell's code leaves K uninitialized.
! 2. If VLAG(KNEW) = MINVAL(VLAG) = VLAG(K) with K < KNEW, Powell's code does not set K = KNEW.
k = knew
if (any(vlag > vlag(knew))) then
    k = maxloc(vlag, mask=(.not. is_nan(vlag)), dim=1)
    !!MATLAB: [~, k] = max(vlag, [], 'omitnan');
end if
step_line = stplen(k) * (xpt(:, k) - xopt)

! Replace STEP by a steepest ascent step from XOPT if the latter provides a larger value of VBIG.
gg = sum(gl**2)
!vgrad = delbar * sqrt(gg)
!if (vgrad <= TENTH * vbig) goto 220
gxpt = matprod(gl, xpt)
ghg = inprod(gxpt, pqlag * gxpt)  ! GHG = INPROD(G, HESS_MUL(G, XPT, PQLAG))
stp = delbar / sqrt(gg)  ! GG can be ZERO!
if (ghg < 0) then
    stp = -stp
end if
step_grad = stp * gl

! Overwrite GL by its projection to the column space of QFAC(:, NACT+1:N). Then set VNEW to the
! greatest value of |LFUNC| on the projected gradient from XOPT subject to the trust region bound.
! If VNEW is sufficiently large, then STEP may be changed to a move along the projected gradient,
! which is the STMP below. STMP does not alter the residuals of the active constraints.
gl = matprod(qfac(:, nact + 1:n), matprod(gl, qfac(:, nact + 1:n)))
!!MATLAB: gl = qfac(:, nact+1:n) * (gl' * qfac(:, nact+1:n))';
gg = sum(gl**2)
!vgrad = delbar * sqrt(gg)
!if (vgrad <= TENTH * vbig) goto 220
gxpt = matprod(gl, xpt)
ghg = inprod(gxpt, pqlag * gxpt)  ! GHG = INPROD(G, HESS_MUL(G, XPT, PQLAG))
stp = delbar / sqrt(gg)  ! GG can be ZERO!
if (ghg < 0) then
    stp = -stp
end if
step_prjg = stp * gl

!--------------------------------------------------------------------------------------------------!
if (nact == n) then
    step_prjg = ZERO
else if (nact == 0) then
    step_prjg = step_grad
end if
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!N.B.: LINCOB tests whether xdiff > 0.1*RHO. Check whether it should be turned off.
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

den_line = calden(kopt, bmat, step_line, xpt, zmat, idz)
den_grad = calden(kopt, bmat, step_grad, xpt, zmat, idz)
den_prjg = calden(kopt, bmat, step_prjg, xpt, zmat, idz)

if (abs(den_grad(knew)) > abs(den_line(knew))) then
    denabs = abs(den_grad(knew))
    step = step_grad
else
    denabs = abs(den_line(knew))
    step = step_line
end if
cvtol = ZERO

cv_prjg = maximum(matprod(step_prjg, amat(:, trueloc(rstat == 1))) - rescon(trueloc(rstat == 1)))
cvtol_prjg = min(0.01_RP * sqrt(sum(step_prjg**2)), TEN * norm(matprod(step_prjg, amat(:, iact(1:nact))), 'inf'))
if (abs(den_prjg(knew)) > 0.1_RP * denabs .and. cv_prjg <= cvtol_prjg) then
    step = step_prjg
    cvtol = cvtol_prjg
end if

constr = ZERO
constr(trueloc(rstat >= 0)) = matprod(step, amat(:, trueloc(rstat >= 0))) - rescon(trueloc(rstat >= 0))
ifeas = all(constr <= cvtol)

end subroutine geostep


end module geometry_mod

!return
!!--------------------------------------------------------------------------------------------------!

!! Set STEP to STMP if STMP gives a sufficiently large value of the modulus of the Lagrange function,
!! and if STMP either preserves feasibility or gives a constraint violation of at least MINCV. The
!! purpose of CVTOL below is to provide a check on feasibility that includes a tolerance for
!! contributions from computer rounding errors. Note that CVTOL equals 0 in precise arithmetic.
!! As commented by Powell, "the projected step is preferred if its value of |LFUNC(XOPT+STEP)| is at
!! least one fifth of the original one but the greatest violation of a linear constraint must be at
!! least MINCV = 0.2*DELBAR in order to keep the interpolation points apart." WHY THE PREFERENCE?
!if (vnew >= 0.2_RP * vbig .or. (is_nan(vbig) .and. .not. is_nan(vnew))) then
!    bigcv = maximum(matprod(stmp, amat(:, trueloc(rstat == 1))) - rescon(trueloc(rstat == 1)))
!    cvtol = min(0.01_RP * sqrt(sstmp), TEN * norm(matprod(stmp, amat(:, iact(1:nact))), 'inf'))
!    ifeas = (bigcv <= cvtol)
!    if (ifeas .or. bigcv >= mincv) then
!        step = stmp
!        return
!    end if
!end if

!220 continue

!! Calculate the greatest constraint violation at XOPT+STEP with STEP at its original value. Modify
!! STEP if this violation is unacceptable.
!! RSTAT(J) = -1, 0, or 1 means constraint J is irrelevant, active, or inactive&relevant, respectively.
!constr = ZERO
!constr(trueloc(rstat >= 0)) = matprod(step, amat(:, trueloc(rstat >= 0))) - rescon(trueloc(rstat >= 0))
!ifeas = all(constr <= 0)
!if (all(constr < mincv) .and. any(constr > 0)) then
!    jsav = int(maxloc(constr, mask=(.not. is_nan(constr)), dim=1), IK)
!    bigcv = constr(jsav)
!    !!MATLAB: [bigcv, jsav] = max(constr, [], 'omitnan');
!    step = step + (mincv - bigcv) * amat(:, jsav)
!end if

!!--------------------------------------------------------------------------------------------------!
!! Zaikun: For the projected gradient, the schemes below work evidently worse than the one above as
!! tested on 20220417. Why?
!!------------------------------------------------------------------------!
!! VESION 1:
!!!gl = gl - matprod(qfac(:, 1:nact), matprod(gl, qfac(:, 1:nact)))
!!------------------------------------------------------------------------!
!! VERSION 2:
!!!if (2 * nact < n) then
!!!    gl = gl - matprod(qfac(:, 1:nact), matprod(gl, qfac(:, 1:nact)))
!!!else
!!!    gl = matprod(qfac(:, nact + 1:n), matprod(gl, qfac(:, nact + 1:n)))
!!!end if
!!------------------------------------------------------------------------!
!!--------------------------------------------------------------------------------------------------!
