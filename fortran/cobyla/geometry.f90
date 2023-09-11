module geometry_cobyla_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Monday, August 07, 2023 AM03:54:43
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: assess_geo, setdrop_geo, setdrop_tr, geostep


contains


function assess_geo(delta, factor_alpha, factor_beta, sim, simi) result(adequate_geo)
!--------------------------------------------------------------------------------------------------!
! This function checks if an interpolation set has acceptable geometry as (14) of the COBYLA paper.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, ONE, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : isinv

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_beta
real(RP), intent(in) :: sim(:, :)   ! SIM(N, N+1)
real(RP), intent(in) :: simi(:, :)  ! SIMI(N, N)

! Outputs
logical :: adequate_geo

! Local variables
character(len=*), parameter :: srname = 'ASSESS_GEO'
integer(IK) :: n
real(RP) :: veta(size(sim, 1))
real(RP) :: vsig(size(sim, 1))
real(RP), parameter :: itol = TENTH

! Sizes
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(factor_alpha > 0 .and. factor_alpha < 1, '0 < FACTOR_ALPHA < 1', srname)
    call assert(factor_beta > 1, 'FACTOR_BETA > 1', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the values of sigma and eta.
! VETA(J) (1 <= J <= N) is the distance between vertices J and 0 (the best vertex) of the simplex.
! VSIG(J) is the distance from vertex J to its opposite face of the simplex. Thus VSIG <= VETA.
! N.B.: What about the distance from vertex N+1 to the its opposite face? Consider the simplex
! {V_{N+1}, V_{N+1} + L*e_1, ..., V_{N+1} + L*e_N}, where V_{N+1} is vertex N+1, namely the current
! "best" point, [e_1, ..., e_n] is an orthogonal matrix, and L is a constant in the order of DELTA.
! This simplex is optimal in the sense that the interpolation system has the minimal condition
! number, i.e., one. For this simplex, the distance from V_{N+1} to its opposite face is L/SQRT{N}.
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
adequate_geo = all(vsig >= factor_alpha * delta) .and. all(veta <= factor_beta * delta)

!====================!
!  Calculation ends  !
!====================!
end function assess_geo


function setdrop_tr(ximproved, d, delta, rho, sim, simi) result(jdrop)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds (the index) of a current interpolation point to be replaced with the
! trust-region trial point. See (19)--(22) of the COBYLA paper.
! N.B.:
! 1. If XIMPROVED == TRUE, then JDROP > 0 so that D is included into XPT. Otherwise, it is a bug.
! 2. COBYLA never sets JDROP = N+1.
! TODO: Check whether it improves the performance if JDROP = N+1 is allowed when XIMPROVED is TRUE.
! Note that UPDATEXFC should be revised accordingly.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, ONE, TENTH, DEBUGGING
use, non_intrinsic :: linalg_mod, only : matprod, isinv
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
logical, intent(in) :: ximproved
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: delta
real(RP), intent(in) :: rho
real(RP), intent(in) :: sim(:, :)   ! SIM(N, N+1)
real(RP), intent(in) :: simi(:, :)  ! SIMI(N, N)

! Outputs
integer(IK) :: jdrop

! Local variables
character(len=*), parameter :: srname = 'SETDROP_TR'
integer(IK) :: n
real(RP) :: distsq(size(sim, 2))
real(RP) :: weight(size(sim, 2))
real(RP) :: score(size(sim, 2))
real(RP) :: simid(size(simi, 1))
!real(RP) :: sigbar(size(sim, 1))
!real(RP) :: veta(size(sim, 1))
!real(RP) :: vsig(size(sim, 1))
real(RP), parameter :: itol = TENTH

! Sizes
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(delta >= rho .and. rho > 0, 'DELTA >= RHO > 0', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

!--------------------------------------------------------------------------------------------------!
! The following code is Powell's scheme for defining JDROP.
!--------------------------------------------------------------------------------------------------!
!! JDROP = 0 by default. It cannot be removed, as JDROP may not be set below in some cases (e.g.,
!! when XIMPROVED == FALSE, MAXVAL(ABS(SIMID)) <= 1, and MAXVAL(VETA) <= EDGMAX).
!jdrop = 0
!
!! SIMID(J) is the value of the J-th Lagrange function at D. It is the counterpart of VLAG in UOBYQA
!! and DEN in NEWUOA/BOBYQA/LINCOA, but it excludes the value of the (N+1)-th Lagrange function.
!simid = matprod(simi, d)
!if (any(abs(simid) > 1) .or. (ximproved .and. any(.not. is_nan(simid)))) then
!    jdrop = int(maxloc(abs(simid), mask=(.not. is_nan(simid)), dim=1), kind(jdrop))
!    !!MATLAB: [~, jdrop] = max(simid, [], 'omitnan');
!end if
!
!! VETA(J) is the distance from the J-th vertex of the simplex to the best vertex, taking the trial
!! point SIM(:, N+1) + D into account.
!if (ximproved) then
!    veta = sqrt(sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1))
!    !!MATLAB: veta = sqrt(sum((sim(:, 1:n) - d).^2));  % d should be a column! Implicit expansion
!else
!    veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
!end if
!
!! VSIG(J) (J=1, .., N) is the Euclidean distance from vertex J to the opposite face of the simplex.
!vsig = ONE / sqrt(sum(simi**2, dim=2))
!sigbar = abs(simid) * vsig
!
!! The following JDROP will overwrite the previous one if its premise holds.
!mask = (veta > factor_delta * delta .and. (sigbar >= factor_alpha * delta .or. sigbar >= vsig))
!if (any(mask)) then
!    jdrop = int(maxloc(veta, mask=mask, dim=1), kind(jdrop))
!    !!MATLAB: etamax = max(veta(mask)); jdrop = find(mask & ~(veta < etamax), 1, 'first');
!end if
!
!! Powell's code does not include the following instructions. With Powell's code, if SIMID consists
!! of only NaN, then JDROP can be 0 even when XIMPROVED == TRUE (i.e., D reduces the merit function).
!! With the following code, JDROP cannot be 0 when XIMPROVED == TRUE, unless VETA is all NaN, which
!! should not happen if X0 does not contain NaN, the trust-region/geometry steps never contain NaN,
!! and we exit once encountering an iterate containing Inf (due to overflow).
!if (ximproved .and. jdrop <= 0) then  ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
!    jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
!    !!MATLAB: [~, jdrop] = max(veta, [], 'omitnan');
!end if
!--------------------------------------------------------------------------------------------------!
! Powell's scheme ends here.
!--------------------------------------------------------------------------------------------------!


! The following definition of JDROP is inspired by SETDROP_TR in UOBYQA/NEWUOA/BOBYQA/LINCOA.
! It is simpler and works better than Powell's scheme. Note that we allow JDROP to be N+1 if
! IMPROVEX is TRUE, whereas Powell's code does not.
! See also (4.1) of Scheinberg-Toint-2010: Self-Correcting Geometry in Model-Based Algorithms for
! Derivative-Free Unconstrained Optimization, which refers to the strategy here as the "combined
! distance/poisedness criteria".

! DISTQ(J) is the square of the distance from the J-th vertex of the simplex to the "best" point so
! far, taking the trial point SIM(:, N+1) + D into account.
if (ximproved) then
    distsq(1:n) = sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1)
    !!MATLAB: distsq = sum((sim(:, 1:n) - d).^2);  % d should be a column! Implicit expansion
    distsq(n + 1) = sum(d**2)
else
    distsq(1:n) = sum(sim(:, 1:n)**2, dim=1)
    distsq(n + 1) = ZERO
end if

weight = max(ONE, distsq / max(rho, TENTH * delta)**2)  ! Similar to Powell's NEWUOA code
! Other possible definitions of WEIGHT.
! !weight = distsq  ! Similar to Powell's LINCOA code, but WRONG. See comments in LINCOA/geometry.f90.
! !weight = max(ONE, 25.0_RP * distsq / delta**2)  ! Similar to Powell's BOBYQA code, works well
! !weight = max(ONE, TEN * distsq / delta**2)
! !weight = max(ONE, 1.0E2_RP * distsq / delta**2)
! !weight = max(ONE, distsq / rho**2)  ! Similar to Powell's UOBYQA

! If 1 <= J <= N, SIMID(J) is the value of the J-th Lagrange function at D; the value of the
! (N+1)-th Lagrange function is 1 - SUM(SIMID). [SIMID, 1 - SUM(SIMID)] is the counterpart of
! VLAG in UOBYQA and DEN in NEWUOA/BOBYQA/LINCOA.
simid = matprod(simi, d)
score = weight * abs([simid, ONE - sum(simid)])

! If XIMPROVED = FALSE (D does not render a better X), set SCORE(N+1) = -1 to avoid JDROP = N+1.
if (.not. ximproved) then
    score(n + 1) = -ONE
end if

! The following IF works a bit better than `IF (ANY(SCORE > 1) .OR. ANY(SCORE > 0) .AND. XIMPROVED)`
! from Powell's UOBYQA and NEWUOA code.
if (any(score > 0)) then  ! Powell's BOBYQA and LINCOA code
    jdrop = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), kind(jdrop))
    !!MATLAB: [~, jdrop] = max(score, [], 'omitnan');
elseif (ximproved) then
    jdrop = int(maxloc(distsq, dim=1), kind(jdrop))
else
    jdrop = 0  ! We arrive here when XIMPROVED = FALSE and no entry of SCORE is positive.
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(jdrop >= 0 .and. jdrop <= n + 1, '0 <= JDROP <= N+1', srname)
    call assert(jdrop <= n .or. ximproved, 'JDROP >= 1 unless IMPROVEX = TRUE', srname)
    call assert(jdrop >= 1 .or. .not. ximproved, 'JDROP >= 1 unless IMPROVEX = FALSE', srname)
    ! JDROP >= 1 when XIMPROVED = TRUE unless NaN occurs in DISTSQ, which should not happen if the
    ! starting point does not contain NaN and the trust-region/geometry steps never contain NaN.
end if

end function setdrop_tr


function setdrop_geo(delta, factor_alpha, factor_beta, sim, simi) result(jdrop)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds (the index) of a current interpolation point to be replaced with a
! geometry-improving point. See (15)--(16) of the COBYLA paper.
! N.B.: COBYLA never sets JDROP = N + 1.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, ONE, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : isinv

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_beta
real(RP), intent(in) :: sim(:, :)   ! SIM(N, N+1)
real(RP), intent(in) :: simi(:, :)  ! SIMI(N, N)

! Outputs
integer(IK) :: jdrop

! Local variables
character(len=*), parameter :: srname = 'SETDROP_GEO'
integer(IK) :: n
real(RP) :: veta(size(sim, 1))
real(RP) :: vsig(size(sim, 1))
real(RP), parameter :: itol = TENTH

! Sizes
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
    call assert(factor_alpha > 0 .and. factor_alpha < 1, '0 < FACTOR_ALPHA < 1', srname)
    call assert(factor_beta > 1, 'FACTOR_BETA > 1', srname)
    call assert(.not. assess_geo(delta, factor_alpha, factor_beta, sim, simi), &
        & 'The geometry is not adequate when '//srname//' is called', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the values of sigma and eta.
! VSIG(J) (J=1, .., N) is the Euclidean distance from vertex J to the opposite face of the simplex.
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))

! Decide which vertex to drop from the simplex. It will be replaced with a new point to improve the
! acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
if (any(veta > factor_beta * delta)) then
    jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
    !!MATLAB: [~, jdrop] = max(veta, [], 'omitnan');
elseif (any(vsig < factor_alpha * delta)) then
    jdrop = int(minloc(vsig, mask=(.not. is_nan(vsig)), dim=1), kind(jdrop))
    !!MATLAB: [~, jdrop] = min(vsig, [], 'omitnan');
else
    ! We arrive here if VSIG and VETA are all NaN, which can happen due to NaN in SIM and SIMI,
    ! which should not happen unless there is a bug.
    jdrop = 0
end if

! Zaikun 20230202: What if we consider VETA and VSIG together? The following attempts do not work well.
! !jdrop = maxloc(sum(sim(:, 1:n)**2, dim=1) * sum(simi**2, dim=2)) ! Condition number
! !jdrop = maxloc(sum(sim(:, 1:n)**2, dim=1)**2 * sum(simi**2, dim=2)) ! Condition number times distance

!====================!
!  Calculation ends  !
!====================!

!Postconditions
if (DEBUGGING) then
    call assert(jdrop >= 1 .and. jdrop <= n, '1 <= JDROP <= N', srname)
end if
end function setdrop_geo


function geostep(jdrop, amat, bvec, conmat, cpen, cval, delta, fval, factor_gamma, simi) result(d)
!--------------------------------------------------------------------------------------------------!
! This function calculates a geometry step so that the geometry of the interpolation set is improved
! when SIM(:, JDRO_GEO) is replaced with SIM(:, N+1) + D. See (15)--(17) of the COBYLA paper.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, ONE, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm

implicit none

! Inputs
integer(IK), intent(in) :: jdrop
real(RP), intent(in) :: amat(:, :)
real(RP), intent(in) :: bvec(:)
real(RP), intent(in) :: conmat(:, :)    ! CONMAT(M, N+1)
real(RP), intent(in) :: cpen
real(RP), intent(in) :: cval(:)     ! CVAL(N+1)
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_gamma
real(RP), intent(in) :: fval(:)     ! FVAL(N+1)
real(RP), intent(in) :: simi(:, :)  ! SIMI(N, N)

! Outputs
real(RP) :: d(size(simi, 1))  ! D(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: m
integer(IK) :: m_lcon
integer(IK) :: n
real(RP) :: A(size(simi, 1), size(conmat, 1))
real(RP) :: cvnd
real(RP) :: cvpd
real(RP) :: g(size(simi, 1))
real(RP) :: vsigj

! Sizes
m_lcon = int(size(bvec), kind(m_lcon))
m = int(size(conmat, 1), kind(m))
n = int(size(simi, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(cpen > 0, 'CPEN > 0', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN/+Inf', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) == [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == NPT and CVAL does not contain negative NaN/+Inf', srname)
    call assert(jdrop >= 1 .and. jdrop <= n, '1 <= JDROP <= N', srname)
    call assert(factor_gamma > 0 .and. factor_gamma < 1, '0 < FACTOR_GAMMA < 1', srname)
end if

!====================!
! Calculation starts !
!====================!

! SIMI(JDROP, :) is a vector perpendicular to the face of the simplex to the opposite of vertex
! JDROP. Thus VSIGJ * SIMI(JDROP, :) is the unit vector in this direction.
vsigj = ONE / sqrt(sum(simi(jdrop, :)**2))

! Set D to the vector in the above-mentioned direction and with length FACTOR_GAMMA * DELTA. Since
! FACTOR_ALPHA < FACTOR_GAMMA < FACTOR_BETA, D improves the geometry of the simplex as per (14) of
! the COBYLA paper. This also explains why this subroutine does not replace DELTA with
! DELBAR = MAX(MIN(TENTH * SQRT(MAXVAL(DISTSQ)), HALF * DELTA), RHO) as in NEWUOA.
d = factor_gamma * delta * (vsigj * simi(jdrop, :))

! The code below chooses the direction of D according to an approximation of the merit function.
! See (17) of the COBYLA paper and  line 225 of Powell's cobylb.f.

! Calculate the coefficients of the linear approximations to the objective and constraint functions.
! N.B.: CONMAT and SIMI have been updated after the last trust-region step, but G and A have not.
! So we cannot pass G and A from outside.
g = matprod(fval(1:n) - fval(n + 1), simi)
A(:, 1:m_lcon) = amat
A(:, m_lcon + 1:m) = transpose(matprod(conmat(m_lcon + 1:m, 1:n) - spread(conmat(m_lcon + 1:m, n + 1), dim=2, ncopies=n), simi))
!!MATLAB: A(:, m_lcon+1:m) = simi'*(conmat(m_lcon+1:m, 1:n) - conmat(m_lcon+1:m, n+1))' % Implicit expansion for subtraction
! CVPD and CVND are the predicted constraint violation of D and -D by the linear models.
cvpd = maxval([ZERO, conmat(:, n + 1) + matprod(d, A)])
cvnd = maxval([ZERO, conmat(:, n + 1) - matprod(d, A)])
! Take -D if the linear models predict that its merit function value is lower.
if (-inprod(d, g) + cpen * cvnd < inprod(d, g) + cpen * cvpd) then
    d = -d
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    ! In theory, ||S|| == FACTOR_GAMMA*DELTA, which may be false due to rounding, but not too far.
    ! It is crucial to ensure that the geometry step is nonzero, which holds in theory.
    call assert(norm(d) > 0.9_RP * factor_gamma * delta .and. norm(d) <= 1.1_RP * factor_gamma * delta, &
        & '||D|| == FACTOR_GAMMA*DELTA', srname)
end if
end function geostep


end module geometry_cobyla_mod
