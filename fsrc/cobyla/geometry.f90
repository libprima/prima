module geometry_mod

contains

function goodgeo(factor_alpha, factor_beta, rho, sim, simi) result(good_geo)

use consts_mod, only : IK, RP, ONE, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize

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
character(len=SRNLEN) :: srname = 'GOODGEO'

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
! VSIG(J) (J=1, .., N) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
vsig = ONE / sqrt(sum(simi**2, dim=2))
veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
good_geo = all(vsig >= parsig) .and. all(veta <= pareta)

end function goodgeo


function setdrop_tr(actrem, d, factor_alpha, factor_delta, rho, sim, simi) result(jdrop)
! This subroutine finds (the index) of a current interpolation point to be replaced by the
! trust-region trial point. See (19)--(21) of the COBYLA paper.

use consts_mod, only : IK, RP, ZERO, ONE, DEBUGGING, SRNLEN
use lina_mod, only : matprod, inprod
use infnan_mod, only : is_nan
use debug_mod, only : errstop, verisize

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
real(RP) :: distx(size(sim, 1))
real(RP) :: edgmax
real(RP) :: parsig
real(RP) :: ratio
real(RP) :: sigbar(size(sim, 1))
real(RP) :: simid(size(sim, 1))
real(RP) :: vsig(size(sim, 1))
character(len=SRNLEN) :: srname = 'SETDROP_TR'

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

! JDROP = 0 by default. It cannot be removed, as JDROP is not set in all cases below.
jdrop = 0_IK

if (actrem <= ZERO) then
    ratio = ONE
else
    ratio = ZERO
end if
simid = matprod(simi, d)

if (maxval(abs(simid)) > ratio) then
    jdrop = int(maxloc(abs(simid), dim=1), kind(jdrop))
end if

if (actrem > ZERO) then
    distx = sqrt(sum((spread(d, dim=2, ncopies=n) - sim(:, 1:n))**2, dim=1))
else
    distx = sqrt(sum(sim(:, 1:n)**2, dim=1))
end if

edgmax = factor_delta * rho
parsig = factor_alpha * rho
vsig = ONE / sqrt(sum(simi**2, dim=2))
sigbar = abs(simid) * vsig

! The following JDROP will overwrite the previous one if its premise holds.
if (any(distx > edgmax .and. (sigbar >= parsig .or. sigbar >= vsig))) then
    jdrop = int(maxloc(distx, mask=(sigbar >= parsig .or. sigbar >= vsig), dim=1), kind(jdrop))
end if

! Powell's code does not include the following instructions. With Powell's code (i.e., the code
! above), THEORETICALLY, JDROP is positive if ACTREM > 0 (i.e., D reduces the merit function).
! However, JDROP may turn out 0 due to NaN in SIMID, SIGBAR, or DISTX. This may depend on the
! compiler and language. Here, we set explicitly JDROP = 0 in case NaN occurs in these arrays, which
! can happen in ill-conditioned problems, although rarely. Consequently, COBYLA will either take a
! geometry step or reduce RHO. (If SIMID contains NaN then so does SIGBAR. So we check only SIGBAR.)
if (is_nan(sum(abs(sigbar))) .or. is_nan(sum(distx))) then
    jdrop = 0_IK
end if

end function setdrop_tr


function setdrop_geo(factor_alpha, factor_beta, rho, sim, simi) result(jdrop)

use consts_mod, only : IK, RP, ONE, DEBUGGING, SRNLEN
use infnan_mod, only : is_nan
use debug_mod, only : errstop, verisize

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
character(len=SRNLEN) :: srname = 'SETDROP_GEO'

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
if (maxval(veta) > pareta) then
    jdrop = int(maxloc(veta, dim=1), kind(jdrop))
else
    jdrop = int(minloc(vsig, dim=1), kind(jdrop))
end if

! Set JDROPT = 0 in case of NaN in VSIG or VETA, which can happen due to NaN in SIM or SIMI.
! Powell's code does not include these instructions.
if (is_nan(sum(vsig)) .or. is_nan(sum(veta))) then
    jdrop = 0_IK
end if

end function setdrop_geo


function geostep(jdrop, cpen, datmat, factor_gamma, rho, simi) result(d)

use consts_mod, only : IK, RP, ZERO, ONE, TWO, DEBUGGING, SRNLEN
use lina_mod, only : matprod, inprod
use debug_mod, only : errstop, verisize

implicit none

! Inputs
integer(IK), intent(in) :: jdrop
real(RP), intent(in) :: simi(:, :)
real(RP), intent(in) :: factor_gamma
real(RP), intent(in) :: cpen
real(RP), intent(in) :: datmat(:, :)
real(RP), intent(in) :: rho

! Output
real(RP) :: d(size(simi, 1))

! Local variables
integer(IK) :: m
integer(IK) :: n
real(RP) :: cvmaxp
real(RP) :: cvmaxm
real(RP) :: vsig(size(simi, 1))
real(RP) :: A(size(simi, 1), size(datmat, 1) - 1)
character(len=SRNLEN) :: srname = 'GEOSTEP'

! Get and verify the sizes
m = size(datmat, 1) - 2
n = size(datmat, 2) - 1
if (DEBUGGING) then
    if (m < 0 .or. n < 1) then
        call errstop(srname, 'SIZE(DATMAT) is invalid')
    end if
    call verisize(simi, n, n)
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
A = transpose(matprod(datmat(1:m + 1, 1:n) - spread(datmat(1:m + 1, n + 1), dim=2, ncopies=n), simi))
A(:, m + 1) = -A(:, m + 1)
cvmaxp = maxval([ZERO, -matprod(d, A(:, 1:m)) - datmat(1:m, n + 1)])
cvmaxm = maxval([ZERO, matprod(d, A(:, 1:m)) - datmat(1:m, n + 1)])
if (TWO * inprod(d, A(:, m + 1)) < cpen * (cvmaxp - cvmaxm)) then
    d = -d
end if

end function geostep


end module geometry_mod
