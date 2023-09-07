module rand_mod
!--------------------------------------------------------------------------------------------------!
! This module provides procedures for setting the random seed and generating random numbers/arrays.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Thursday, September 07, 2023 PM06:44:08
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
! N.B.: We use exclusively the DEFAULT INTEGER and the DOUBLE-PRECISION REAL in this procedure.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DP
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: infos_mod, only : MEMORY_ALLOCATION_FAILS
implicit none

integer, intent(in) :: seed

character(len=*), parameter :: srname = 'SETSEED0'
integer :: alloc_status
integer :: p
integer :: i
integer :: n  ! Should be a default INTEGER according to F2018.
integer, allocatable :: seed_to_put(:)
real(DP), allocatable :: cos_seed(:)

call random_seed(size=n)

if (allocated(seed_to_put)) deallocate (seed_to_put)
allocate (seed_to_put(1:n), stat=alloc_status)
if (.not. (alloc_status == 0 .and. allocated(seed_to_put))) then
    call errstop(srname, 'Memory allocation fails', MEMORY_ALLOCATION_FAILS)
end if

if (allocated(cos_seed)) deallocate (cos_seed)
allocate (cos_seed(1:n), stat=alloc_status)
if (.not. (alloc_status == 0 .and. allocated(cos_seed))) then
    call errstop(srname, 'Memory allocation fails', MEMORY_ALLOCATION_FAILS)
end if

! Some compilers cannot guarantee ABS(COS) <= 1 when the variable is huge. This may cause overflow.
! Note that 1.0_DP cannot be written as ONE, because KIND(ONE) = RP, which may not be DP.
cos_seed = min(max(cos(real([(i, i=seed - (n - 1), seed)], DP)), -1.0_DP), 1.0_DP)
seed_to_put = ceiling(0.9_DP * real(huge(0), DP) * cos_seed)
deallocate (cos_seed)
! P takes a `+1` at the end, so that it is guarantee to be positive.
p = int(real(huge(0), DP) / 1.0E2_DP) + 1
seed_to_put = modulo(seed_to_put, p) + 1

call random_seed(put=seed_to_put)
deallocate (seed_to_put)

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
use, non_intrinsic :: consts_mod, only : RP
implicit none

real(RP) :: x

call random_number(harvest=x)
end function rand0

function rand1(n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an N-dimensional random vector with entries iid sampled from U([0, 1)).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none

integer(IK), intent(in) :: n
real(RP) :: x(max(n, 0_IK))

call random_number(harvest=x)
end function

function rand2(m, n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an MxN random matrix with entries iid sampled from U([0, 1)).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none

integer(IK), intent(in) :: m, n
real(RP) :: x(max(m, 0_IK), max(n, 0_IK))

call random_number(harvest=x)
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

call random_number(harvest=u)
call random_number(harvest=v)
x = sqrt(-TWO * log(ONE - u)) * cos(TWO * PI * v)
end function randn0

function randn1(n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an N-dimensional random vector with entries iid sampled from N(0,1).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
implicit none

integer(IK), intent(in) :: n
real(RP) :: x(max(n, 0_IK))

real(RP) :: u(size(x))
real(RP) :: v(size(x))

call random_number(harvest=u)
call random_number(harvest=v)
x = sqrt(-TWO * log(ONE - u)) * cos(TWO * PI * v)
end function randn1

function randn2(m, n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to an MxN random matrix with entries iid sampled from N(0,1).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
implicit none

integer(IK), intent(in) :: m, n
real(RP) :: x(max(m, 0_IK), max(n, 0_IK))

real(RP) :: u(size(x, 1), size(x, 2))
real(RP) :: v(size(x, 1), size(x, 2))

call random_number(harvest=u)
call random_number(harvest=v)
x = sqrt(-TWO * log(ONE - u)) * cos(TWO * PI * v)
end function randn2


end module rand_mod
