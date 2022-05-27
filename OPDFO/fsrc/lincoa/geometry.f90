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
! Last Modified: Friday, May 27, 2022 PM10:36:20
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(iact, idz, knew, kopt, nact, amat, del, gl_in, qfac, rescon, xopt, xpt, zmat, ifeas, step, bmat)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF, TEN, TENTH, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, isorth, maximum, trueloc, norm
use, non_intrinsic :: powalg_mod, only : hess_mul, omega_col, calvlag, calbeta

implicit none

! Inputs
integer(IK), intent(in) :: iact(:)  ! IACT(M)
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
integer(IK), intent(in) :: nact
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: del
real(RP), intent(in) :: gl_in(:)  ! GL_IN(N)
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(in) :: rescon(:)  ! RESCON(M)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
logical, intent(out) :: ifeas
real(RP), intent(out) :: step(:)  ! STEP(N)

! Local variables
logical :: prjg
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: m
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: rstat(size(amat, 2))
real(RP) :: stmp(size(xopt))
real(RP) :: gl(size(gl_in))
real(RP) :: constr(size(amat, 2))
real(RP) :: bigcv, cvtol, cvtol_prjg, gg, gxpt(size(xpt, 2)), ghg, sp, ss, tol, &
&        stp, stplen(size(xpt, 2)), stpsav, mincv, vbig, vgrad, vlag(size(xpt, 2)), vnew, sstmp
real(RP) :: pqlag(size(xpt, 2))  ! PQLAG(NPT)
integer(IK) :: jsav, k, ksav
real(RP) :: step_prjg(size(xpt, 1)), step_grad(size(xpt, 1)), step_line(size(xpt, 1))
real(RP) :: denabs, beta_line, beta_grad, beta_prjg, denom_line, denom_grad, denom_prjg, denom, denmax, &
    & vlag_line(size(xpt, 1) + size(xpt, 2)), &
    & vlag_grad(size(xpt, 1) + size(xpt, 2)), vlag_prjg(size(xpt, 1) + size(xpt, 2)), alpha

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
    call assert(del > 0, 'DEL > 0', srname)
    call assert(size(gl_in) == n, 'SIZE(GL_IN) == N', srname)
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

ifeas = .false. !??? Added by Zaikun 20220227
prjg = .false.
alpha = -sum(zmat(knew, 1:idz - 1)**2) + sum(zmat(knew, idz:size(zmat, 2))**2)
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
!     IFEAS will be set to TRUE or FALSE if XOPT+STEP is feasible or infeasible.
!
!     STEP is chosen to provide a relatively large value of the modulus of
!       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
!       calculated too, within the trust region, that does not alter the
!       residuals of the active constraints. The projected step is preferred
!       if its value of |LFUNC(XOPT+STEP)| is at least one fifth of the
!       original one but the greatest violation of a linear constraint must
!       be at least MINCV = 0.2*DEL, in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes MINCV.
!
!     Set some constants.
!

mincv = 0.2_RP * del  ! Is this really better than 0? According to an experiment of Tom on 20220225, NO

! Replace GL by the gradient of LFUNC at the trust region centre, and set the elements of RSTAT.
! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC.
pqlag = omega_col(idz, zmat, knew)
gl = gl + hess_mul(xopt, xpt, pqlag)

! RSTAT(J) = -1, 0, or 1 respectively means constraint J is irrelevant, active, or inactive&relevant.
! RSTAT never changes after being set below.
rstat = 1_IK
rstat(trueloc(abs(rescon) >= del)) = -1_IK
rstat(iact(1:nact)) = 0_IK

! Maximize |LFUNC| within the trust region on the lines through XOPT and other interpolation points.
! In the following loop, VLAG(K) is set to the maximum of PHI_K(t) with subject to the trust-region
! constraint with PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:, K)). Note that PHI_K(t) is a quadratic
! function with PHI_K(0) = 0, PHI_K(1) = 1, and PHI_K'(0) = SP = INPROD(GL, XPT(:, K) - XOPT).
vlag = ZERO
do k = 1, npt
    if (k == kopt) then
        cycle
    end if
    ss = sum((xpt(:, k) - xopt)**2)
    !sp = inprod(gl, xpt(:, k) - xopt)
    sp = inprod(gl, xpt(:, k)) - inprod(gl, xopt)
    stplen(k) = -del / sqrt(ss)
    if (k == knew) then
        if (sp * (sp - ONE) < ZERO) then
            stplen(k) = -stplen(k)
        end if
        vlag(k) = abs(stplen(k) * sp) + stplen(k)**2 * abs(sp - ONE)
    else
        vlag(k) = abs(stplen(k) * (ONE - stplen(k)) * sp)
    end if
end do

! N.B.: 0. We define KSAV in a way slightly different from Powell's code, which sets KSAV to
! MAXLOC(VLAG,DIM=1) by comparing the entries of VLAG sequentially.
! 1. If VLAG contains only NaN, which can happen, Powell's code leaves KSAV uninitialized.
! 2. If VLAG(KNEW) = MINVAL(VLAG) = VLAG(K) with K < KNEW, Powell's code does not set KSAV = KNEW.
ksav = knew
if (any(vlag > vlag(knew))) then
    ksav = maxloc(vlag, mask=(.not. is_nan(vlag)), dim=1)
    !!MATLAB: [~, ksav] = max(vlag, [], 'omitnan');
end if
vbig = vlag(ksav)
stpsav = stplen(ksav)

! Set STEP to the move that gives the greatest modulus calculated above. This move may be replaced
! by a steepest ascent step from XOPT.
step = stpsav * (xpt(:, ksav) - xopt)
step_line = step
vlag_line = calvlag(kopt, bmat, step_line, xpt, zmat, idz)
beta_line = calbeta(kopt, bmat, step_line, xpt, zmat, idz)
denom_line = alpha * beta_line + vlag_line(knew)**2
denabs = abs(denom_line)
cvtol = ZERO

gg = sum(gl**2)
vgrad = del * sqrt(gg)
!if (vgrad <= TENTH * vbig) goto 220
if (gg > 0) then
! Make the replacement if it provides a larger value of VBIG.
    gxpt = matprod(gl, xpt)
    ghg = inprod(gxpt, pqlag * gxpt)  ! GHG = INPROD(G, HESS_MUL(G, XPT, PQLAG))
    vnew = vgrad + abs(HALF * del * del * ghg / gg)
!if (vnew > vbig .or. (is_nan(vbig) .and. .not. is_nan(vnew))) then
    vbig = vnew
    stp = del / sqrt(gg)
    if (ghg < ZERO) then
        stp = -stp
    end if
!end if
    step_grad = stp * gl
    vlag_grad = calvlag(kopt, bmat, step_grad, xpt, zmat, idz)
    beta_grad = calbeta(kopt, bmat, step_grad, xpt, zmat, idz)
    denom_grad = alpha * beta_grad + vlag_grad(knew)**2
    if (abs(denom_grad) > denabs) then
        denabs = abs(denom_grad)
        step = step_grad
    end if
end if

constr = ZERO
constr(trueloc(rstat >= 0)) = matprod(step, amat(:, trueloc(rstat >= 0))) - rescon(trueloc(rstat >= 0))
ifeas = all(constr <= 0)

! Overwrite GL by its projection to the column space of QFAC(:, NACT+1:N). Then set VNEW to the
! greatest value of |LFUNC| on the projected gradient from XOPT subject to the trust region bound.
! If VNEW is sufficiently large, then STEP may be changed to a move along the projected gradient.
gl = matprod(qfac(:, nact + 1:n), matprod(gl, qfac(:, nact + 1:n)))
!!MATLAB: gl = qfac(:, nact+1:n) * (gl' * qfac(:, nact+1:n))';
!--------------------------------------------------------------------------------------------------!
! Zaikun: The schemes below work evidently worse than the one above as tested on 20220417. Why?
!------------------------------------------------------------------------!
! VESION 1:
!!gl = gl - matprod(qfac(:, 1:nact), matprod(gl, qfac(:, 1:nact)))
!------------------------------------------------------------------------!
! VERSION 2:
!!if (2 * nact < n) then
!!    gl = gl - matprod(qfac(:, 1:nact), matprod(gl, qfac(:, 1:nact)))
!!else
!!    gl = matprod(qfac(:, nact + 1:n), matprod(gl, qfac(:, nact + 1:n)))
!!end if
!------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

gg = sum(gl**2)
if (nact > 0 .and. nact < n .and. gg > 0) then
    vgrad = del * sqrt(gg)
!if (vgrad <= TENTH * vbig) goto 220
    gxpt = matprod(gl, xpt)
    ghg = inprod(gxpt, pqlag * gxpt)  ! GHG = INPROD(G, HESS_MUL(G, XPT, PQLAG))
    vnew = vgrad + abs(HALF * del * del * ghg / gg)

! Set STMP to the possible move along the projected gradient.
    stp = del / sqrt(gg)
    if (ghg < ZERO) then
        stp = -stp
    end if
    stmp = stp * gl
    sstmp = sum(stmp**2)
    step_prjg = stmp
    vlag_prjg = calvlag(kopt, bmat, step_prjg, xpt, zmat, idz)
    beta_prjg = calbeta(kopt, bmat, step_prjg, xpt, zmat, idz)
    denom_prjg = alpha * beta_prjg + vlag_prjg(knew)**2

    bigcv = maximum(matprod(stmp, amat(:, trueloc(rstat == 1))) - rescon(trueloc(rstat == 1)))
    cvtol_prjg = min(0.01_RP * sqrt(sstmp), TEN * norm(matprod(stmp, amat(:, iact(1:nact))), 'inf'))
    if (abs(denom_prjg) > 0.1_RP * denabs .and. bigcv <= cvtol_prjg) then
        step = step_prjg
        cvtol = cvtol_prjg
        prjg = .true.
        ifeas = (bigcv <= cvtol_prjg)
    end if
end if
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


!! Set STEP to STMP if STMP gives a sufficiently large value of the modulus of the Lagrange function,
!! andif STMP either preserves feasibility or gives a constraint violation of at least MINCV. The
!! purpose of CVTOL below is to provide a check on feasibility that includes a tolerance for
!! contributions from computer rounding errors.
!! As commented by Powell, "the projected step is preferred if its value of |LFUNC(XOPT+STEP)| is at
!! least one fifth of the original one but the greatest violation of a linear constraint must be at
!! least MINCV = 0.2*DEL, in order to keep the interpolation points apart." WHY THE PREFERENCE?
!if (vnew >= 0.2_RP * vbig .or. (is_nan(vbig) .and. .not. is_nan(vnew))) then
!    bigcv = maximum(matprod(stmp, amat(:, trueloc(rstat == 1))) - rescon(trueloc(rstat == 1)))
!    cvtol = min(0.01_RP * sqrt(sstmp), TEN * norm(matprod(stmp, amat(:, iact(1:nact))), 'inf'))
!    ifeas = (bigcv <= cvtol)
!    !!ifeas = (bigcv < mincv)
!    !cvtol = ZERO
!    !!if (bigcv > ZERO .and. bigcv < 0.01_RP * sqrt(sstmp)) then
!    !if (bigcv > ZERO .and. bigcv <= 0.01_RP * sqrt(sstmp)) then
!    !    cvtol = maxval([ZERO, abs(matprod(stmp, amat(:, iact(1:nact))))])
!    !end if
!    !if (bigcv <= TEN * cvtol .or. bigcv >= mincv) then
!    if (ifeas .or. bigcv >= mincv) then
!        !ifeas = (bigcv <= TEN * cvtol)
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

end subroutine geostep


end module geometry_mod
