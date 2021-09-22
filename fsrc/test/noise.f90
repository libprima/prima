module noise_mod
! This module provides functions that add noise to X or function values.
implicit none
private
public :: noisy, noisyx, noisyfun


contains


function noisy(x) result(noisy_x)

use, non_intrinsic :: consts_mod, only : RP, ONE, TWO
use, non_intrinsic :: param_mod, only : NOISE
use, non_intrinsic :: rand_mod, only : rand
implicit none

real(RP), intent(in) :: x
real(RP) :: noisy_x

noisy_x = x + NOISE * max(ONE, abs(x)) * (TWO * rand() - ONE)

end function noisy


function noisyx(x) result(noisy_x)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
use, non_intrinsic :: linalg_mod, only : norm
use, non_intrinsic :: param_mod, only : X_NOISE
use, non_intrinsic :: rand_mod, only : randn
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: noisy_x(size(x))

real(RP) :: r(size(x))

r = randn(int(size(r), IK))

noisy_x = x + X_NOISE * max(ONE, norm(x)) * (r / norm(r))

end function noisyx


function noisyfun(x, f) result(noify_f)
use, intrinsic :: ieee_arithmetic, only : IEEE_VALUE, IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF
use, non_intrinsic :: consts_mod, only : RP, ONE, TWO
use, non_intrinsic :: param_mod, only : F_NOISE
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(in) :: f
real(RP) :: noify_f

integer, allocatable :: seedsav(:)
integer :: seed
real(RP) :: r

call getseed(seedsav)
seed = ceiling(1.0E6_RP * (cos(1.0E6_RP * sin(1.0E6_RP * (abs(f) + ONE) * cos(1.0E6_RP * sum(abs(x)))))))
call setseed(max(abs(seed), 1))
call setseed(seedsav)
deallocate (seedsav)

r = TWO * rand() - ONE

noify_f = f + F_NOISE * max(ONE, abs(f)) * r

if (r > 0.75_RP) then
    noify_f = IEEE_VALUE(ONE, IEEE_POSITIVE_INF)
!elseif (r > 0.5_RP) then
    !noify_f = IEEE_VALUE(ONE, IEEE_QUIET_NAN)
elseif (r < -0.9995_RP) then
    noify_f = IEEE_VALUE(ONE, IEEE_NEGATIVE_INF)
end if
end function noisyfun


end module noise_mod
