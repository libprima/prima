! MEMORY_MOD is a module providing subroutines concerning memory management.
!
! Coded by Zaikun ZHANG in July 2020.
!
! Last Modified: Sunday, May 23, 2021 AM10:59:02


#include "ppf.h"

module memory_mod

implicit none
private
public :: safealloc, cstyle_sizeof

interface safealloc
    module procedure alloc_rvector, alloc_rmatrix
end interface safealloc

interface cstyle_sizeof
    module procedure size_of_sp, size_of_dp
#if __QP_AVAILABLE__ == 1
    module procedure size_of_qp
#endif
end interface cstyle_sizeof


contains


subroutine alloc_rvector(x, n)
! ALLOC_RVECTOR allocates the space for an allocatable REAL(RP)
! vector X, whose size is N after allocation.
use consts_mod, only : RP, IK, SRNLEN
use debug_mod, only : errstop
implicit none

! Input
integer(IK), intent(in) :: n

! Output
real(RP), allocatable, intent(out) :: x(:)

! Local variable
integer :: alloc_status
character(len=SRNLEN), parameter :: srname = 'ALLOC_RVECTOR'

! According to the Fortran 2003 standard, when a procedure is invoked,
! any allocated ALLOCATABLE object that is an actual argument associated
! with an INTENT(OUT) ALLOCATABLE dummy argument is deallocated. So it is
! unnecessary to write the following line in F2003 since X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
if (n >= 1) then
    x(1) = 0.0_RP
end if

end subroutine alloc_rvector


subroutine alloc_rmatrix(x, m, n)
! ALLOC_RMATRIX allocates the space for a REAL(RP) matrix X, whose
! size is (M, N) after allocation.
use consts_mod, only : RP, IK, SRNLEN
use debug_mod, only : errstop
implicit none

! Input
integer(IK), intent(in) :: m, n

! Output
real(RP), allocatable, intent(out) :: x(:, :)

! Local variable
integer :: alloc_status
character(len=SRNLEN), parameter :: srname = 'ALLOC_RMATRIX'

! Unnecessary to write the following line in F2003 since X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(m, n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
if (m * n >= 1) then
    x(1, 1) = 0.0_RP
end if

end subroutine alloc_rmatrix


pure function size_of_sp(x) result(y)
use consts_mod, only : SP, IK
implicit none
real(SP), intent(in) :: x
integer(IK) :: y
#if __USE_STORAGE_SIZE__ == 1
y = int(storage_size(x) / 8, kind(y))
#else
y = int(kind(x), kind(y)) ! Avoid complaint
y = int(4, kind(y))
#endif
end function size_of_sp


pure function size_of_dp(x) result(y)
use consts_mod, only : DP, IK
implicit none
real(DP), intent(in) :: x
integer(IK) :: y
#if __USE_STORAGE_SIZE__ == 1
y = int(storage_size(x) / 8, kind(y))
#else
y = int(kind(x), kind(y)) ! Avoid complaint
y = int(8, kind(y))
#endif
end function size_of_dp


#if __QP_AVAILABLE__ == 1

pure function size_of_qp(x) result(y)
use consts_mod, only : QP, IK
implicit none
real(QP), intent(in) :: x
integer(IK) :: y
#if __USE_STORAGE_SIZE__ == 1
y = int(storage_size(x) / 8, kind(y))
#else
y = int(kind(x), kind(y)) ! Avoid complaint
y = int(8, kind(y))
#endif
end function size_of_qp

#endif


end module memory_mod
