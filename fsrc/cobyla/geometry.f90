module geometry_mod

implicit none
private
public :: goodgeo, geostep, setdrop_geo, setdrop_tr


contains

function goodgeo(factor_alpha, factor_beta, rho, sim, simi) result(good_geo)

use, non_intrinsic :: consts_mod, only : IK, RP, ONE, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop, verisize

implicit none

! Inputs
real(RP), intent(in) :: sim(:, :)
real(RP), intent(in) :: simi(:, :)
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_beta
real(RP), intent(in) :: rho

! Output
logical :: good_geo

! Local variables
integer(IK) :: n
real(RP) :: parsig
real(RP) :: pareta
real(RP) :: vsig(size(sim, 1))
real(RP) :: veta(size(sim, 1))
character(len=*), parameter :: srname = 'GOODGEO'

! Get and verify the sizes
n = size(sim, 1)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(SIM, 1) < 1')
    end if
    call verisize(sim, n, n + 1)
    call verisize(simi, n, n)
end if

! Calculate the values of sigma and eta.
parsig = factor_alpha * rho
pareta = factor_beta * rho
! VETA(J) (1 <= J <= N) is the distance between vertices J and 0 (the best vertex) of the simplex.
! VSIG(J) is the distance from vertex J to the opposite face of the simplex. Thus VSIG <= VETA.
! But what about vertex N+1?
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
good_geo = all(vsig >= parsig) .and. all(veta <= pareta)

end function goodgeo


function setdrop_tr(actrem, d, factor_alpha, factor_delta, rho, sim, simi) result(jdrop)
! This subroutine finds (the index) of a current interpolation point to be replaced by the
! trust-region trial point. See (19)--(21) of the COBYLA paper.
! N.B.:
! 1. If ACTREM > 0, then JDROP > 0 so that D is included into XPT. Otherwise, it is a bug.
! 2. COBYLA never sets JDROP = N + 1.

use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, ONE, DEBUGGING
use, non_intrinsic :: linalg_mod, only : matprod, inprod
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: debug_mod, only : errstop, verisize

implicit none

! Inputs
real(RP), intent(in) :: actrem
real(RP), intent(in) :: d(:)
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_delta
real(RP), intent(in) :: rho
real(RP), intent(in) :: sim(:, :)
real(RP), intent(in) :: simi(:, :)

! Output
integer(IK) :: jdrop

! Local variables
integer(IK) :: n
real(RP) :: veta(size(sim, 1))
real(RP) :: edgmax
real(RP) :: parsig
real(RP) :: sigbar(size(sim, 1))
real(RP) :: simid(size(sim, 1))
real(RP) :: vsig(size(sim, 1))
character(len=*), parameter :: srname = 'SETDROP_TR'

! Get and verify the sizes
n = size(sim, 1)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(SIM, 1) < 1')
    end if
    call verisize(d, n)
    call verisize(sim, n, n + 1)
    call verisize(simi, n, n)
end if

! JDROP = 0 by default. It cannot be removed, as JDROP may not be set below in some cases (e.g.,
! when ACTREM <= 0, MAXVAL(ABS(SIMID)) <= 1, and MAXVAL(VETA) <= EDGMAX).
jdrop = 0_IK

simid = matprod(simi, d)
if (any(abs(simid) > ONE) .or. (actrem > ZERO .and. any(.not. is_nan(simid)))) then
    jdrop = int(maxloc(abs(simid), mask=(.not. is_nan(simid)), dim=1), kind(jdrop))
end if

if (actrem > ZERO) then
    veta = sqrt(sum((sim(:, 1:n) - spread(d, dim=2, ncopies=n))**2, dim=1))
else
    veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
end if
edgmax = factor_delta * rho
parsig = factor_alpha * rho
vsig = ONE / sqrt(sum(simi**2, dim=2))
sigbar = abs(simid) * vsig
! The following JDROP will overwrite the previous one if its premise holds.
if (any(veta > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig))) then
    jdrop = int(maxloc(veta, mask=(veta > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig)), &
        & dim=1), kind(jdrop))
end if

! Powell's code does not include the following instructions. With Powell's code, if SIMID consists
! of only NaN, then JDROP can be 0 even when ACTREM > 0 (i.e., D reduces the merit function).
! With the following code, JDROP cannot be 0 when ACTREM > 0, unless VETA is all NaN, which should
! not happen if X0 does not contain NaN, the trust-region/geometry steps never contain NaN, and we
! exit once encountering an iterate containing Inf (due to overflow).
if (actrem > ZERO .and. jdrop <= 0) then  ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
    jdrop = int(maxloc(veta, mask=(.not. is_nan(veta)), dim=1), kind(jdrop))
end if

if (DEBUGGING) then
    if (actrem > ZERO .and. jdrop <= 0) then ! Write JDROP <= 0 instead of JDROP == 0 for robustness.
        ! This can happen only if NaN occurs in VETA, which should not happen if the starting point
        ! does not contain NaN and the trust-region/geometry steps never contain NaN.
        call errstop(srname, 'ACTREM > 0 but JDROP < 1')
    end if
end if

end function setdrop_tr


function setdrop_geo(factor_alpha, factor_beta, rho, sim, simi) result(jdrop)
! N.B.: COBYLA never sets JDROP = N + 1.

use, non_intrinsic :: consts_mod, only : IK, RP, ONE, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: debug_mod, only : errstop, verisize

implicit none

! Inputs
real(RP), intent(in) :: sim(:, :)
real(RP), intent(in) :: simi(:, :)
real(RP), intent(in) :: factor_alpha
real(RP), intent(in) :: factor_beta
real(RP), intent(in) :: rho

! Output
integer(IK) :: jdrop

! Local variables
integer(IK) :: n
real(RP) :: parsig
real(RP) :: pareta
real(RP) :: vsig(size(sim, 1))
real(RP) :: veta(size(sim, 1))
character(len=*), parameter :: srname = 'SETDROP_GEO'

! Get and verify the sizes.
n = size(sim, 1)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(SIM, 1) < 1')
    end if
    call verisize(sim, n, n + 1)
    call verisize(simi, n, n)
end if

! Calculate the values of sigma and eta.
parsig = factor_alpha * rho
pareta = factor_beta * rho
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
    ! We arrive here if VSIG and VETA are all NaN, which can happen due to NaN in SIM and SIMI.
    jdrop = 0_IK
end if

end function setdrop_geo


function geostep(jdrop, cpen, conmat, cval, fval, factor_gamma, rho, simi) result(d)

use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, ONE, TWO, DEBUGGING
use, non_intrinsic :: linalg_mod, only : matprod, inprod
use, non_intrinsic :: debug_mod, only : errstop, verisize

implicit none

! Inputs
integer(IK), intent(in) :: jdrop
real(RP), intent(in) :: simi(:, :)
real(RP), intent(in) :: factor_gamma
real(RP), intent(in) :: cpen
real(RP), intent(in) :: conmat(:, :)
real(RP), intent(in) :: cval(:)
real(RP), intent(in) :: fval(:)
real(RP), intent(in) :: rho

! Output
real(RP) :: d(size(simi, 1))

! Local variables
integer(IK) :: m
integer(IK) :: n
real(RP) :: cvmaxp
real(RP) :: cvmaxm
real(RP) :: vsig(size(simi, 1))
real(RP) :: A(size(simi, 1), size(conmat, 1) + 1)
character(len=*), parameter :: srname = 'GEOSTEP'

! Get and verify the sizes
m = size(conmat, 1)
n = size(simi, 1)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(SIMI) is invalid')
    end if
    call verisize(fval, n + 1)
    call verisize(conmat, m, n + 1)
    call verisize(cval, n + 1)
    call verisize(simi, n, n)
    if (jdrop < 1_IK .or. jdrop > n) then
        call errstop(srname, 'JDROP < 1 or JDROP > N')
    end if
end if

! VSIG(J) (J=1, .., N) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
vsig = ONE / sqrt(sum(simi**2, dim=2))
d = factor_gamma * rho * vsig(jdrop) * simi(jdrop, :)
! Calculate the coefficients of the linear approximations to the objective and constraint functions,
! placing minus the objective function gradient after the constraint gradients in the array A.
! When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result as
! Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
! finite-precision arithmetic.
A(:, 1:m) = transpose(matprod(conmat(:, 1:n) - spread(conmat(:, n + 1), dim=2, ncopies=n), simi))
A(:, m + 1) = matprod(fval(n + 1) - fval(1:n), simi)
cvmaxp = maxval([ZERO, -matprod(d, A(:, 1:m)) - conmat(:, n + 1)])
cvmaxm = maxval([ZERO, matprod(d, A(:, 1:m)) - conmat(:, n + 1)])
if (TWO * inprod(d, A(:, m + 1)) < cpen * (cvmaxp - cvmaxm)) then
    d = -d
end if

end function geostep


end module geometry_mod
