module rand_mod
!--------------------------------------------------------------------------------------------------!
! This module provides procedures for setting the random seed and generating random numbers/arrays.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Mon 18 Aug 2025 10:16:31 AM CST
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getseed, setseed, rand, randn

interface setseed
! SETSEED sets the random seed by calling RANDOM_SEED.
! N.B.: with SEED being an integer, SETSEED(SEED) differs from SETSEED([SEED]).
    module procedure setseed0, setseed1
end interface setseed

interface rand
! RAND generates a random number/vector/matrix whose entries are iid sampled from U([0,1)).
! N.B.: Here, RAND(N) is an N-dimensional vector; in MATLAB, RAND(N) is an NxN matrix. DIFFERENT!
    module procedure rand0, rand1, rand2
end interface rand

interface randn
! RANDN generates a random number/vector/matrix whose entries are iid sampled from N(0,1). They are
! generated from U([0,1)) by the Box-Muller transform.
! N.B.: Here, RANDN(N) is an N-dimensional vector; in MATLAB, RANDN(N) is an NxN matrix. DIFFERENT!
    module procedure randn0, randn1, randn2
end interface randn


contains


function scramble_seed(n, seed) result(scrambled_seed)
!--------------------------------------------------------------------------------------------------!
! This function scrambles the input SEED using a linear congruential generator (LCG) algorithm.
! The output is a vector of length N, where each entry is a scrambled version of the input seed.
! The LCG parameters used are based on the MINSTD algorithm.
! -------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : INT64
implicit none

! Inputs
integer, intent(in) :: n
integer, intent(in) :: seed

! Outputs
integer :: scrambled_seed(n)

! LCG parameters (MINSTD)
integer(INT64), parameter :: modulus = 2147483647_INT64  ! 2^31 - 1
integer(INT64), parameter :: multiplier = 48271_INT64

! Local variables
integer :: i
integer :: j
integer(INT64) :: iseed

! Map seed to [1, modulus-1] in case it is larger (unlikely).
iseed = abs(modulo(int(seed, INT64), modulus))
if (iseed == 0) iseed = 1  ! Ensure non-zero seed

! Apply LCG to generate n scrambled seeds.
do i = 1, n
    do j = 1, 16  ! Loop to ensure good scrambling; without this, flang-new 20 produces poor results
        ! Compute the next ISEED. Since MULTIPLIER is a prime number, and the input ISEED < MODULUS,
        ! the output ISEED is non-zero.
        iseed = modulo(multiplier * iseed, modulus)  ! 64-bit safe multiplication, avoiding overflow
    end do
    ! Convert to the default integer and store.
    scrambled_seed(i) = int(iseed)
end do
end function scramble_seed


subroutine getseed(seed)
!--------------------------------------------------------------------------------------------------!
! This subroutine gets the current random seed, saving it in SEED.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: infos_mod, only : MEMORY_ALLOCATION_FAILS
implicit none

! Output
integer, allocatable, intent(out) :: seed(:)

! Local variable
character(len=*), parameter :: srname = 'GETSEED'
integer :: alloc_status
integer :: n  ! Should be a default INTEGER according to F2018.

call random_seed(size=n)

! The following line is unnecessary since F2003 as X is INTENT(OUT):
! !if (allocated(seed)) deallocate (seed)

! 1. The following allocation is NOT removable even in F2003.
! 2. Why not using SAFEALLOC? Because the kind of SEED is the default integer, while SAFEALLOC is
!    only implemented for INTEGER(IK), which may differ from the default integer.
allocate (seed(1:n), stat=alloc_status)
if (.not. (alloc_status == 0 .and. allocated(seed))) then
    call errstop(srname, 'Memory allocation fails', MEMORY_ALLOCATION_FAILS)
end if

call random_seed(get=seed)
end subroutine getseed


subroutine setseed0(seed)
!--------------------------------------------------------------------------------------------------!
! This procedure uses SEED to initialize the random seed. See SEED_TO_PUT for details.
! N.B.: We use exclusively the DEFAULT INTEGER in this procedure.
!--------------------------------------------------------------------------------------------------!
implicit none

integer, intent(in) :: seed
integer :: n  ! Should be a default INTEGER according to F2018.

call random_seed(size=n)
call random_seed(put=scramble_seed(n, seed))
end subroutine setseed0

subroutine setseed1(seed)
!--------------------------------------------------------------------------------------------------!
! This procedure sets the random seed by calling RANDOM_SEED.
! Let N be the size of the random seed on this platform. If SEED is present and SIZE(SEED) < N,
! then the random seed is to [SEED, 1, ..., 1]; if SEED is present and SIZE(SEED) >= N, then
! the random seed is set to SEED(1:N); if SEED is absent, then the random seed is set by invoking
! RANDOM_SEED with no argument, which sets the random seed with system-dependent random data.
! N.B.: Do NOT change the entries of SEED in any way, because this procedure may be used to restore
! random seeds that are saved before.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: infos_mod, only : MEMORY_ALLOCATION_FAILS
implicit none

integer, intent(in), optional :: seed(:)

character(len=*), parameter :: srname = 'SETSEED1'
integer :: alloc_status
integer :: n  ! Should be a default INTEGER according to F2018.
integer, allocatable :: seed_to_put(:)

if (present(seed)) then
    call random_seed(size=n)

    if (allocated(seed_to_put)) deallocate (seed_to_put)
    allocate (seed_to_put(1:n), stat=alloc_status)
    if (.not. (alloc_status == 0 .and. allocated(seed_to_put))) then
        call errstop(srname, 'Memory allocation fails', MEMORY_ALLOCATION_FAILS)
    end if

    seed_to_put = 1
    seed_to_put(1:min(size(seed), n)) = seed(1:min(size(seed), n))
    call random_seed(put=seed_to_put)
    deallocate (seed_to_put)
else
    call random_seed()
end if
end subroutine setseed1


function rand0() result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to a random number sampled from U([0, 1)).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO, ONE
implicit none

real(RP) :: x

call random_number(harvest=x)
! Zaikun 20250817: The following line takes the fractional part of K*X, with K being a reasonably
! large integer. If X follows U([0, 1)), then the fractional part of K*X follows U([0, 1)) as well.
! We do this because the X generated above seems not "random enough" with some compilers, e.g.,
! flang-new 20.1.8 (SCRAMBLE_SEED helps). We hope this transformation can improve the randomness.
x = x * real(10**min(range(0), range(x)), RP)
x = max(ZERO, min(ONE, x - real(floor(x), RP)))  ! MAX/MIN are for safety. Not needed in theory.
end function rand0

function rand1(n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an N-dimensional random vector with entries iid sampled from U([0, 1)).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE
implicit none

integer(IK), intent(in) :: n
real(RP) :: x(n)

call random_number(harvest=x)
! Zaikun 20250817: See the comment in RAND0 for the following two lines.
x = x * real(10**min(range(0), range(x)), RP)
x = max(ZERO, min(ONE, x - real(floor(x), RP)))  ! MAX/MIN are for safety. Not needed in theory.
end function

function rand2(m, n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an MxN random matrix with entries iid sampled from U([0, 1)).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE
implicit none

integer(IK), intent(in) :: m, n
real(RP) :: x(m, n)

call random_number(harvest=x)
! Zaikun 20250817: See the comment in RAND0 for the following two lines.
x = x * real(10**min(range(0), range(x)), RP)
x = max(ZERO, min(ONE, x - real(floor(x), RP)))  ! MAX/MIN are for safety. Not needed in theory.
end function


function randn0() result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to a random number sampled from N(0,1).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, TWO, PI
implicit none

real(RP) :: x

real(RP) :: u
real(RP) :: v

u = rand()
v = rand()
x = sqrt(-TWO * log(ONE - u)) * cos(TWO * PI * v)
end function randn0

function randn1(n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an N-dimensional random vector with entries iid sampled from N(0,1).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
implicit none

integer(IK), intent(in) :: n
real(RP) :: x(n)

real(RP) :: u(n)
real(RP) :: v(n)

u = rand(n)
v = rand(n)
x = sqrt(-TWO * log(ONE - u)) * cos(TWO * PI * v)
end function randn1

function randn2(m, n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an MxN random matrix with entries iid sampled from N(0,1).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
implicit none

integer(IK), intent(in) :: m, n
real(RP) :: x(m, n)

real(RP) :: u(m, n)
real(RP) :: v(m, n)

u = rand(m, n)
v = rand(m, n)
x = sqrt(-TWO * log(ONE - u)) * cos(TWO * PI * v)
end function randn2


end module rand_mod
