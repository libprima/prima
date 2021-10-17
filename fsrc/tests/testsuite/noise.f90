module noise_mod
!--------------------------------------------------------------------------------------------------!
! This module provides functions that add noise to X or function values.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Sunday, October 17, 2021 PM07:50:57
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: noisy, noisyx, noisyfun

integer, parameter :: NOISE_TYPE_LEN = 64


contains


function noisy(x, noise_level) result(noisy_x)
!--------------------------------------------------------------------------------------------------!
! Return a noisy X, NOISE_LEVEL being the relative noise level.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, TWO
use, non_intrinsic :: param_mod, only : NOISE_DFT
use, non_intrinsic :: rand_mod, only : rand
implicit none

! Inputs
real(RP), intent(in) :: x
real(RP), optional, intent(in) :: noise_level

! Outputs
real(RP) :: noisy_x

! Local variables
real(RP) :: noise_level_loc

if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = NOISE_DFT
end if

noisy_x = x + noise_level_loc * max(abs(x), ONE) * (TWO * rand() - ONE)

end function noisy


function noisyx(x, noise_level) result(noisy_x)
!--------------------------------------------------------------------------------------------------!
! Return a noisy X, NOISE_LEVEL being the relative noise level.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
use, non_intrinsic :: linalg_mod, only : norm
use, non_intrinsic :: param_mod, only : XNOISE_DFT
use, non_intrinsic :: rand_mod, only : randn
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), optional, intent(in) :: noise_level

! Outputs
real(RP) :: noisy_x(size(x))

! Local variables
real(RP) :: r(size(x))
real(RP) :: noise_level_loc

if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = XNOISE_DFT
end if

r = randn(int(size(r), IK))

noisy_x = x + noise_level_loc * max(norm(x), ONE) * (r / norm(r))

end function noisyx


function noisyfun(x, f, noise_level, noise_type) result(noify_f)
!--------------------------------------------------------------------------------------------------!
! Return a noisy F, NOISE_LEVEL being the relative noise level, and NOISE_TYPE being 'uniform'
! or 'gaussian'.
! We implement the function in a way such that
! 1. the noise is uniquely determined by (X, F); multiple runs with the same (X, F) return the
! same result;
! 2. the noises for different (X, F) are (nearly) independent from each other.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, TWO, TEN, TENTH, EPS
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: param_mod, only : FNOISE_DFT
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand, randn
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: f
real(RP), optional, intent(in) :: noise_level
character(len=*), optional, intent(in) :: noise_type

! Outputs
real(RP) :: noify_f

! Local variables
character(len=NOISE_TYPE_LEN) :: noise_type_loc
integer, allocatable :: seedsav(:)
integer :: seed
real(RP) :: noise_level_loc
real(RP) :: r
real(RP) :: seedx
real(RP) :: seedf

! Read the optional inputs.
if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = FNOISE_DFT
end if
if (present(noise_type)) then
    noise_type_loc = noise_type
else
    noise_type_loc = 'gaussian'
end if

! Quick returns
! Without doing this, some compilers (e.g., Absoft 21.0) may have difficulties when evaluating COS.
if (.not. is_finite(f)) then
    noify_f = f
    return
end if

! Define the seed.
! Many compilers have difficulties in handling COS of huge variables. They may return invalid
! results (NaN or numbers with absolute values much larger than 1, observed on Absoft 21.0 and NAG
! Fortran 7.0) or even encounter segmentation faults (Absoft 21.0).
seedx = sum(cos(x * TEN**(-int(log10(abs(x) + EPS)))))
seedf = cos(f * TEN**(-int(log10(abs(f) + EPS))))
seed = ceiling(TENTH * real(huge(0), RP) * seedx * seedf)
seed = max(abs(seed), 1)

! Generate a random number R.
call getseed(seedsav)  ! Back up the current random seed in SEEDSAV.
call setseed(seed)  ! Set the random seed by SETSEED(SEED).
if (lower(trimstr(noise_type_loc)) == 'uniform') then
    r = TWO * rand() - ONE
else
    r = randn()
end if
call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
deallocate (seedsav)

! Define NOISY_F.
noify_f = f + noise_level_loc * max(abs(f), ONE) * r

! Inject faulty values into F.
r = rand()  ! Generate a random value that is independent of the NOISY_F above.
if (r > 0.8_RP) then
    noify_f = huge(0.0_RP)
!elseif (r > 0.7_RP) then
    !noify_f = IEEE_VALUE(0.0_RP, IEEE_QUIET_NAN)
elseif (r < 0.1_RP) then
    noify_f = -10.0_RP * max(abs(f), ONE)
end if
end function noisyfun


end module noise_mod
