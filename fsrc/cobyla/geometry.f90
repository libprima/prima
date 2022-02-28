module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Monday, February 28, 2022 PM01:47:24
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: goodgeo, setdrop_geo, setdrop_tr, geostep


contains


function goodgeo(delta, factor_alpha, factor_beta, sim, simi) result(good_geo)
!--------------------------------------------------------------------------------------------------!
! This function checks whether an interpolation set has good geometry as (14) of the COBYLA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ONE, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : isinv

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_beta
real(RP), intent(in) :: sim(:, :)
real(RP), intent(in) :: simi(:, :)

! Outputs
logical :: good_geo

! Local variables
character(len=*), parameter :: srname = 'GOODGEO'
integer(IK) :: n
real(RP) :: pareta
real(RP) :: parsig
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
parsig = factor_alpha * delta
pareta = factor_beta * delta
! VETA(J) (1 <= J <= N) is the distance between vertices J and 0 (the best vertex) of the simplex.
! VSIG(J) is the distance from vertex J to the opposite face of the simplex. Thus VSIG <= VETA.
! But what about vertex N+1?
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
good_geo = all(vsig >= parsig) .and. all(veta <= pareta)

!====================!
!  Calculation ends  !
!====================!
end function goodgeo


function setdrop_tr(tr_success, d, delta, factor_alpha, factor_delta, sim, simi) result(jdrop)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds (the index) of a current interpolation point to be replaced by the
! trust-region trial point. See (19)--(22) of the COBYLA paper.
! N.B.:
! 1. If TR_SUCCESS == TRUE, then JDROP > 0 so that D is included into XPT. Otherwise, it is a bug.
! 2. COBYLA never sets JDROP = N + 1.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ONE, TENTH, DEBUGGING
use, non_intrinsic :: linalg_mod, only : matprod, isinv
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
logical, intent(in) :: tr_success
real(RP), intent(in) :: d(:)
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_delta
real(RP), intent(in) :: sim(:, :)
real(RP), intent(in) :: simi(:, :)

! Outputs
integer(IK) :: jdrop

! Local variables
character(len=*), parameter :: srname = 'SETDROP_TR'
integer(IK) :: n
real(RP) :: edgmax
real(RP) :: parsig
real(RP) :: sigbar(size(sim, 1))
real(RP) :: simid(size(sim, 1))
real(RP) :: veta(size(sim, 1))
real(RP) :: vsig(size(sim, 1))
real(RP), parameter :: itol = TENTH

! Sizes
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

! JDROP = 0 by default. It cannot be removed, as JDROP may not be set below in some cases (e.g.,
! when TR_SUCCESS == FALSE, MAXVAL(ABS(SIMID)) <= 1, and MAXVAL(VETA) <= EDGMAX).
jdrop = 0_IK

simid = matprod(simi, d)
if (any(abs(simid) > ONE) .or. (tr_success .and. any(.not. is_nan(simid)))) then
    jdrop = int(maxloc(abs(simid), mask=(.not. is_nan(simid)), dim=1), kind(jdrop))
end if

if (tr_success) then
    veta = sqrt(sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1))
else
    veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
end if
edgmax = factor_delta * delta
parsig = factor_alpha * delta
vsig = ONE / sqrt(sum(simi**2, dim=2))
sigbar = abs(simid) * vsig
! The following JDROP will overwrite the previous one if its premise holds.
if (any(veta > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig))) then
    jdrop = int(maxloc(veta, mask=(veta > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig)), &
        & dim=1), kind(jdrop))
end if

! Powell's code does not include the following instructions. With Powell's code, if SIMID consists
! of only NaN, then JDROP can be 0 even when TR_SUCCESS == TRUE (i.e., D reduces the merit function).
! With the following code, JDROP cannot be 0 when TR_SUCCESS == TRUE, unless VETA is all NaN, which
! should not happen if X0 does not contain NaN, the trust-region/geometry steps never contain NaN,
! and we exit once encountering an iterate containing Inf (due to overflow).
if (tr_success .and. jdrop <= 0) then  ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
    jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(jdrop >= 0 .and. jdrop <= n, '0 <= JDROP <= N', srname)
    call assert(jdrop >= 1 .or. .not. tr_success, 'JDROP >= 1 unless the trust-region step failed', srname)
end if

end function setdrop_tr


function setdrop_geo(delta, factor_alpha, factor_beta, sim, simi) result(jdrop)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds (the index) of a current interpolation point to be replaced by
! a geometry-improving point. See (15)--(16) of the COBYLA paper.
! N.B.: COBYLA never sets JDROP = N + 1.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ONE, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : isinv

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_beta
real(RP), intent(in) :: sim(:, :)
real(RP), intent(in) :: simi(:, :)

! Outputs
integer(IK) :: jdrop

! Local variables
character(len=*), parameter :: srname = 'SETDROP_GEO'
integer(IK) :: n
real(RP) :: pareta
real(RP) :: parsig
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
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the values of sigma and eta.
parsig = factor_alpha * delta
pareta = factor_beta * delta
! VSIG(J) (J=1, .., N) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))

! Decide which vertex to drop from the simplex. It will be replaced by a new point to improve
! acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
if (any(veta > pareta)) then
    jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
elseif (any(vsig < parsig)) then
    jdrop = int(minloc(vsig, mask=(.not. is_nan(vsig)), dim=1), kind(jdrop))
else
    ! We arrive here if VSIG and VETA are all NaN, which can happen due to NaN in SIM and SIMI,
    ! which should not happen unless there is a bug.
    jdrop = 0_IK
end if

!====================!
!  Calculation ends  !
!====================!

!Postconditions
if (DEBUGGING) then
    call assert(jdrop >= 1 .and. jdrop <= n, '1 <= JDROP <= N', srname)
end if
end function setdrop_geo


function geostep(jdrop, cpen, conmat, cval, delta, fval, factor_gamma, simi) result(d)
!--------------------------------------------------------------------------------------------------!
! This function calculates a geometry step so that the geometry of the interpolation set is improved
! when SIM(:, JDRO_GEO) is replaced by SIM(:, N+1) + D.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, ONE, TWO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf, is_neginf
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm

implicit none

! Inputs
integer(IK), intent(in) :: jdrop
real(RP), intent(in) :: conmat(:, :)
real(RP), intent(in) :: cpen
real(RP), intent(in) :: cval(:)
real(RP), intent(in) :: delta
real(RP), intent(in) :: factor_gamma
real(RP), intent(in) :: fval(:)
real(RP), intent(in) :: simi(:, :)

! Outputs
real(RP) :: d(size(simi, 1))

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: m
integer(IK) :: n
real(RP) :: A(size(simi, 1), size(conmat, 1) + 1)
real(RP) :: cvmaxn
real(RP) :: cvmaxp
real(RP) :: vsigj

! Sizes
m = int(size(conmat, 1), kind(m))
n = int(size(simi, 1), kind(m))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(cpen >= 0, 'CPEN >= 0', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN/+Inf', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) == [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_neginf(conmat)), 'CONMAT does not contain NaN/-Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == NPT and FVAL is not NaN/+Inf', srname)
    call assert(jdrop >= 1 .and. jdrop <= n, '1 <= JDROP <= N', srname)
end if

!====================!
! Calculation starts !
!====================!

! SIMI(JDROP, :) is a vector perpendicular to the face of the simplex to the opposite of vertex
! JDROP. Thus VSIGJ * SIMI(JDROP, :) is the unit vector in this direction.
vsigj = ONE / sqrt(sum(simi(jdrop, :)**2))

! Set D to the vector in the above-mentioned direction and with length FACTOR_GAMMA * DELTA. Since
! FACTOR_ALPHA < FACTOR_GAMMA < FACTOR_BETA, D improves the geometry of the simplex as per GOODGEO.
d = factor_gamma * delta * (vsigj * simi(jdrop, :))

! Calculate the coefficients of the linear approximations to the objective and constraint functions,
! placing minus the objective function gradient after the constraint gradients in the array A.
! When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result as
! Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
! finite-precision arithmetic.
A(:, 1:m) = transpose(matprod(conmat(:, 1:n) - spread(conmat(:, n + 1), dim=2, ncopies=n), simi))
A(:, m + 1) = matprod(fval(n + 1) - fval(1:n), simi)
cvmaxp = maxval([ZERO, -matprod(d, A(:, 1:m)) - conmat(:, n + 1)])
cvmaxn = maxval([ZERO, matprod(d, A(:, 1:m)) - conmat(:, n + 1)])
if (TWO * inprod(d, A(:, m + 1)) < cpen * (cvmaxp - cvmaxn)) then
    d = -d
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(norm(d) <= TWO * factor_gamma * delta, '|D| <= 2*FACTOR_GAMMA*DELTA', srname)
end if
end function geostep


end module geometry_mod
