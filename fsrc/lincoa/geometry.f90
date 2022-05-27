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
! Last Modified: Friday, May 27, 2022 PM11:43:16
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: geostep


contains


subroutine geostep(iact, idz, knew, kopt, nact, amat, bmat, delbar, qfac, rescon, xopt, xpt, zmat, ifeas, s)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds S such that intends to improve the geometry of the interpolation set
! when XPT(:, KNEW) is changed to XOPT + S, where XOPT = XPT(:, KOPT).
!  S is chosen to provide a relatively large value of the modulus of
!       LFUNC(XOPT+S), subject to ||S|| .LE. DELBAR. A projected S is
!       calculated too, within the trust region, that does not alter the
!       residuals of the active constraints. The projected step is preferred
!       if its value of |LFUNC(XOPT+S)| is at least one fifth of the
!       original one but the greatest violation of a linear constraint must
!       be at least MINCV = 0.2*DELBAR in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes MINCV.
!     IFEAS will be set to TRUE or FALSE if XOPT+S is feasible or infeasible.
!
!     AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the same as the terms with these names
!     in subroutine LINCOB. KNEW is the index of the interpolation point that is going to be moved.
!     DELBAR is the current restriction on the length of S, which is never greater than the
!     current trust region radius DELTA.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TEN, TENTH, EPS, DEBUGGING
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
real(RP), intent(out) :: s(:)  ! S(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: rstat(size(amat, 2))
real(RP) :: cstrv
real(RP) :: cvtol
real(RP) :: dderiv(size(xpt, 2))
real(RP) :: den(size(xpt, 2))
real(RP) :: denabs
real(RP) :: distsq(size(xpt, 2))
real(RP) :: gg
real(RP) :: ghg
real(RP) :: glag(size(xpt, 1))
real(RP) :: grads(size(xpt, 1))
real(RP) :: gxpt(size(xpt, 2))
real(RP) :: pgrads(size(xpt, 1))
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: stplen(size(xpt, 2))
real(RP) :: tol
real(RP) :: vlagabs(size(xpt, 2))

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

!====================!
! Calculation starts !
!====================!

! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC. Set GLAG to the gradient of LFUNC at the trust region centre.
pqlag = omega_col(idz, zmat, knew)
glag = bmat(:, knew) + hess_mul(xopt, xpt, pqlag)

! Maximize |LFUNC| within the trust region on the lines through XOPT and other interpolation points,
! without considering the linear constraints. In the following, VLAGABS(K) is set to the maximum of
! |PHI_K(t)| subject to the trust-region constraint with PHI_K(t) = LFUNC((1-t)*XOPT + t*XPT(:, K)).
dderiv = matprod(glag, xpt) - inprod(glag, xopt) ! The derivatives PHI_K'(0).
distsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
! For each K /= KNEW, |PHI_K(t)| is maximized by STPLEN(K), the maximum being VLAGABS(K). Note that
! PHI_K(t) is a quadratic function with PHI_K'(0) = DDERIV(K) and PHI_K(0) = 0 = PHI_K(1).
stplen = -delbar / sqrt(distsq)
vlagabs = abs(stplen * (ONE - stplen) * dderiv)
! The maximization of |PHI_K(t)| is as follows. Note that PHI_K(t) is a quadratic function with
! PHI_K'(0) = DDERIV(K), PHI_K(0) = 0, and PHI_K(1) = 1.
if (dderiv(knew) * (dderiv(knew) - ONE) < 0) then
    stplen(knew) = -stplen(knew)
end if
vlagabs(knew) = abs(stplen(knew) * dderiv(knew)) + stplen(knew)**2 * abs(dderiv(knew) - ONE)
! It does not make sense to consider "straight line through XOPT and XPT(:, KOPT)". Thus we set
! VLAGABS(KOPT) to -1 so that KOPT is skipped in the comparison below.
vlagabs(kopt) = -ONE

! Find K so that VLAGABS(K) is maximized. We define K in a way slightly different from Powell's
! code, which sets K to MAXLOC(VLAGABS) by comparing the entries of VLAGABS sequentially.
! 1. If VLAGABS contains only NaN, which can happen, Powell's code leaves K uninitialized.
! 2. If VLAGABS(KNEW) = MINVAL(VLAG) = VLAG(K) with K < KNEW, Powell's code does not set K = KNEW.
k = knew
if (any(vlagabs > vlagabs(knew))) then
    k = maxloc(vlagabs, mask=(.not. is_nan(vlagabs)), dim=1)
    !!MATLAB: [~, k] = max(vlagabs, [], 'omitnan');
end if
s = stplen(k) * (xpt(:, k) - xopt)
den = calden(kopt, bmat, s, xpt, zmat, idz)
denabs = abs(den(knew))

! Replace S by a steepest ascent step from XOPT if the latter provides a larger value of DENABS.
gg = sum(glag**2)
if (gg > 0) then
    gxpt = matprod(glag, xpt)
    ghg = inprod(gxpt, pqlag * gxpt)  ! GHG = INPROD(G, HESS_MUL(G, XPT, PQLAG))
    if (ghg < 0) then
        grads = -(delbar / sqrt(gg)) * glag
    else  ! This ELSE covers the unlikely yet possible case where GHG is zero or even NaN.
        grads = (delbar / sqrt(gg)) * glag
    end if
    den = calden(kopt, bmat, grads, xpt, zmat, idz)
    if (abs(den(knew)) > denabs) then
        denabs = abs(den(knew))
        s = grads
    end if
end if

! Set IFEAS for the calculated S. RSTAT identifies the constraints that need evaluation. RSTAT(J)
! is -1, 0, or 1 respectively means constraint J is irrelevant, active, or inactive and relevant.
rstat = 1_IK
rstat(trueloc(abs(rescon) >= delbar)) = -1_IK
rstat(iact(1:nact)) = 0_IK
cstrv = maximum([ZERO, matprod(s, amat(:, trueloc(rstat >= 0))) - rescon(trueloc(rstat >= 0))])
ifeas = (cstrv <= 0)

! Return with the current S if NACT == 0 or NACT == N.
if (nact <= 0 .or. nact >= n) then
    return
end if

! When 0 < NACT < N, define PGRADS by maximizing |LFUNC| within the trust region from XOPT along the
! projection of GLAG to the column space of QFAC(:, NACT+1:N), which is the orthogonal complement of
! the space spanned by the active gradients. In precise arithmetic, moving along PGRADS does not
! change the values of the active constraints. This projected gradient step is preferred and will
! override S if it renders a denominator not too small and leads to good feasibility. This turns out
! critical for the performance of LINCOA. We first replace GLAG with the projected gradient, so that
! the calculation of PGRADS follows the same scheme as PGRAD.
glag = matprod(qfac(:, nact + 1:n), matprod(glag, qfac(:, nact + 1:n)))
!!MATLAB: glag = qfac(:, nact+1:n) * (glag' * qfac(:, nact+1:n))';
gg = sum(glag**2)
if (gg > 0) then
    gxpt = matprod(glag, xpt)
    ghg = inprod(gxpt, pqlag * gxpt)  ! GHG = INPROD(G, HESS_MUL(G, XPT, PQLAG))
    if (ghg < 0) then
        pgrads = -(delbar / sqrt(gg)) * glag
    else  ! This ELSE covers the unlikely yet possible case where GHG is zero or even NaN.
        pgrads = (delbar / sqrt(gg)) * glag
    end if
    den = calden(kopt, bmat, pgrads, xpt, zmat, idz)
    cstrv = maximum([ZERO, matprod(pgrads, amat(:, trueloc(rstat == 1))) - rescon(trueloc(rstat == 1))])
    !The purpose of CVTOL below is to provide a check on feasibility that includes a tolerance for
    !contributions from computer rounding errors. Note that CVTOL equals 0 in precise arithmetic.
    cvtol = min(0.01_RP * sqrt(sum(pgrads**2)), TEN * norm(matprod(pgrads, amat(:, iact(1:nact))), 'inf'))
    if (abs(den(knew)) > TENTH * denabs .and. cstrv <= cvtol) then
        s = pgrads
        ifeas = .true.  ! IFEAS = (CSTRV <= CVTOL)
    end if
end if

!====================!
!  Calculation ends  !
!====================!

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!N.B.: LINCOB tests whether xdiff > 0.1*RHO. Check whether it should be turned off.
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

end subroutine geostep


end module geometry_mod
