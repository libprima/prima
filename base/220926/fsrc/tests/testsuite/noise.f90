module noise_mod
!--------------------------------------------------------------------------------------------------!
! This module provides functions that add noise to X or function values.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Thursday, January 06, 2022 PM12:27:02
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: pintrf_mod, only : OBJ, OBJCON
implicit none

private
public :: noisy
public :: noisy_calfun, noisy_calcfc
public :: orig_calfun, orig_calcfc

interface noisy
    module procedure noisy0, noisy1
end interface noisy

interface noisyfun
    module procedure noisyfun0, noisyfun1
end interface noisyfun

procedure(OBJ), pointer :: orig_calfun
procedure(OBJCON), pointer :: orig_calcfc
integer, parameter :: NOISE_TYPE_LEN = 64


contains


function noisy0(x, noise_level, noise_type, seed) result(noisy_x)
!--------------------------------------------------------------------------------------------------!
! Return a noisy X according to NOISE_LEVEL, NOISE_TYPE, and SEED.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, TWO
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: param_mod, only : NOISE_LEVEL_DFT, NOISE_TYPE_DFT
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand, randn
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Inputs
real(RP), intent(in) :: x
real(RP), intent(in), optional :: noise_level
character(len=*), intent(in), optional :: noise_type
integer, intent(in), optional :: seed

! Outputs
real(RP) :: noisy_x

! Local variables
character(len=NOISE_TYPE_LEN) :: noise_type_loc
integer, allocatable :: seedsav(:)
real(RP) :: noise_level_loc
real(RP) :: r

! Quick return
if (.not. is_finite(x)) then
    noisy_x = x
    return
end if

! Read the optional inputs.
if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = NOISE_LEVEL_DFT
end if
if (present(noise_type)) then
    noise_type_loc = noise_type
else
    noise_type_loc = NOISE_TYPE_DFT
end if

! Back up the random seed and then set the seed to the user-provided SEED if it is present.
if (present(seed)) then
    call getseed(seedsav)  ! Back up the current random seed in SEEDSAV.
    call setseed(seed)  ! Set the random seed by SETSEED(SEED).
end if

select case (lower(trimstr(noise_type_loc)))
case ('uniform')
    r = TWO * rand() - ONE
case ('gaussian')
    r = randn()
case default
    r = randn()
end select

noisy_x = x + noise_level_loc * max(abs(x), ONE) * r

! Recover the random seed if necessary.
if (present(seed)) then
    call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
!deallocate (seedsav)  ! SEEDSAV is deallocated automatically at return.
end if
end function noisy0


function noisy1(x, noise_level, noise_type, seed) result(noisy_x)
!--------------------------------------------------------------------------------------------------!
! Return a noisy X according to NOISE_LEVEL, NOISE_TYPE, and SEED.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
use, non_intrinsic :: linalg_mod, only : norm
use, non_intrinsic :: param_mod, only : NOISE_LEVEL_DFT, NOISE_TYPE_DFT
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand, randn
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in), optional :: noise_level
character(len=*), intent(in), optional :: noise_type
integer, intent(in), optional :: seed

! Outputs
real(RP) :: noisy_x(size(x))

! Local variables
character(len=NOISE_TYPE_LEN) :: noise_type_loc
integer(IK) :: i
integer, allocatable :: seedsav(:)
real(RP) :: noise_level_loc
real(RP) :: r(size(x))

! Read the optional inputs.
if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = NOISE_LEVEL_DFT
end if
if (present(noise_type)) then
    noise_type_loc = noise_type
else
    noise_type_loc = NOISE_TYPE_DFT
end if

! Back up the random seed and then set the seed to the user-provided SEED if it is present.
if (present(seed)) then
    call getseed(seedsav)  ! Back up the current random seed in SEEDSAV.
    call setseed(seed)  ! Set the random seed by SETSEED(SEED).
end if

select case (lower(trimstr(noise_type_loc)))
case ('sphere')
    r = randn(int(size(x), IK))
    r = r / norm(r)
    noisy_x = x + noise_level_loc * max(norm(x), ONE) * r
case ('ball')
    r = randn(int(size(x), IK))
    r = r / norm(r)
    noisy_x = x + noise_level_loc * max(norm(x), ONE) * rand() * r
case default
    noisy_x = [(noisy0(x(i), noise_level_loc, noise_type), i=1, int(size(x), IK))]
end select

! Recover the random seed if necessary.
if (present(seed)) then
    call setseed(seedsav)
!deallocate (seedsav)  ! SEEDSAV is deallocated automatically at return.
end if
end function noisy1


function noisyfun0(x, f, noise_level, noise_type) result(noisy_f)
!--------------------------------------------------------------------------------------------------!
! Return a noisy F according to NOISE_LEVEL, NOISE_TYPE, and SEED.
! We implement the function in a way such that
! 1. the noise is fully determined by (X, F); multiple runs with the same (X, F) return the same;
! 2. the noises for different (X, F) are (nearly) independent from each other.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, TEN, EPS
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: param_mod, only : NOISE_LEVEL_DFT, NOISE_TYPE_DFT
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: f
real(RP), intent(in), optional :: noise_level
character(len=*), intent(in), optional :: noise_type

! Outputs
real(RP) :: noisy_f

! Local variables
character(len=NOISE_TYPE_LEN) :: noise_type_loc
integer, allocatable :: seedsav(:)
integer :: seed
real(RP) :: noise_level_loc
real(RP) :: r
real(RP) :: seedx
real(RP) :: seedf

! Quick return
! Without doing this, some compilers (e.g., Absoft 21.0) may have difficulties when evaluating COS.
if (.not. is_finite(f)) then
    noisy_f = f
    return
end if

! Read the optional inputs.
if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = NOISE_LEVEL_DFT
end if
if (present(noise_type)) then
    noise_type_loc = noise_type
else
    noise_type_loc = NOISE_TYPE_DFT
end if

! Define the seed by scrambling X and F.
! Many compilers have difficulties in handling COS of huge variables. They may return invalid
! results (NaN or numbers with absolute values much larger than 1, observed on Absoft 21.0 and NAG
! Fortran 7.0) or even encounter segmentation faults (Absoft 21.0).
seedx = sum(cos(x * TEN**(-int(log10(abs(x) + EPS))))) / real(size(x), RP)
seedf = cos(f * TEN**(-int(log10(abs(f) + EPS))))
seed = ceiling(real(10**range(0), RP) * seedx * seedf)
seed = max(abs(seed), 1)

! Define NOISY_F.
noisy_f = noisy(f, noise_level_loc, noise_type_loc, seed)

! Inject faulty values into F.
call getseed(seedsav)  ! Back up the current random seed in SEEDSAV.
call setseed(seed)  ! Set the random seed by SETSEED(SEED).
r = rand()  ! Generate a random value that is independent of the NOISY_F above.
if (r > 0.8_RP) then
    noisy_f = rand() * huge(0.0_RP)  ! "Almost +Inf" values.
!elseif (r > 0.7_RP) then
    !noisy_f = IEEE_VALUE(0.0_RP, IEEE_QUIET_NAN)  ! NaN
elseif (r < 0.1_RP) then
    noisy_f = -10.0_RP * max(abs(f), ONE)  ! Discontinuously decreased values
end if
call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
!deallocate (seedsav)  ! SEEDSAV is deallocated automatically at return.
end function noisyfun0


function noisyfun1(x, f, noise_level, noise_type) result(noisy_f)
!--------------------------------------------------------------------------------------------------!
! Return a noisy F according to NOISE_LEVEL, NOISE_TYPE, and SEED.
! We implement the function in a way such that
! 1. the noise is fully determined by (X, F); multiple runs with the same (X, F) return the same;
! 2. the noises for different (X, F) are (nearly) independent from each other.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: param_mod, only : NOISE_LEVEL_DFT, NOISE_TYPE_DFT
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: f(:)
real(RP), intent(in), optional :: noise_level
character(len=*), intent(in), optional :: noise_type

! Outputs
real(RP) :: noisy_f(size(f))

! Local variables
character(len=NOISE_TYPE_LEN) :: noise_type_loc
integer(IK) :: i
real(RP) :: noise_level_loc

! Read the optional inputs.
if (present(noise_level)) then
    noise_level_loc = noise_level
else
    noise_level_loc = NOISE_LEVEL_DFT
end if
if (present(noise_type)) then
    noise_type_loc = noise_type
else
    noise_type_loc = NOISE_TYPE_DFT
end if

noisy_f = [(noisyfun0(x, f(i), noise_level_loc, noise_type_loc), i=1, int(size(f), IK))]
end function noisyfun1


subroutine noisy_calfun(x, f)
!--------------------------------------------------------------------------------------------------!
! Noisy version of ORIG_CALFUN.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: param_mod, only : NOISE_LEVEL_DFT, NOISE_TYPE_DFT
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
call orig_calfun(x, f)
f = noisyfun(x, f, noise_level=NOISE_LEVEL_DFT, noise_type=NOISE_TYPE_DFT)
end subroutine noisy_calfun


subroutine noisy_calcfc(x, f, constr)
!--------------------------------------------------------------------------------------------------!
! Noisy version of ORIG_CALCFC.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: param_mod, only : NOISE_LEVEL_DFT, NOISE_TYPE_DFT
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
call orig_calcfc(x, f, constr)
f = noisyfun(x, f, noise_level=NOISE_LEVEL_DFT, noise_type=NOISE_TYPE_DFT)
! N.B.: the constraints are CONSTR >= 0. In other words, we are "maximizing" CONSTR. Thus we impose
! noise to -CONSTR and then take the negative. Otherwise, NOISYFUN tends to introduce +Inf into
! CONSTR, which is not intended.
constr = -noisyfun(x, -constr, noise_level=NOISE_LEVEL_DFT, noise_type=NOISE_TYPE_DFT)
end subroutine noisy_calcfc


end module noise_mod
