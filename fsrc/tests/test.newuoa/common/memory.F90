! MEMORY_MOD is a module providing subroutines concerning memory management.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Thursday, September 23, 2021 AM12:24:00


#include "ppf.h"

module memory_mod

implicit none
private
public :: safealloc, cstyle_sizeof

interface safealloc
    module procedure alloc_rvector, alloc_rmatrix, alloc_ivector, alloc_imatrix
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
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Input
integer(IK), intent(in) :: n

! Output
real(RP), allocatable, intent(out) :: x(:)

! Local variable
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RVECTOR'

! According to the Fortran 2003 standard, when a procedure is invoked,
! any allocated ALLOCATABLE object that is an actual argument associated
! with an INTENT(OUT) ALLOCATABLE dummy argument is deallocated. So it is
! unnecessary to write the following line since F2003 as X is INTENT(OUT):
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
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Input
integer(IK), intent(in) :: m, n

! Output
real(RP), allocatable, intent(out) :: x(:, :)

! Local variable
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RMATRIX'

! Unnecessary to write the following line since F2003 as X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(m, n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
if (m >= 1 .and. n >= 1) then
    x(1, 1) = 0.0_RP
end if

end subroutine alloc_rmatrix


subroutine alloc_ivector(x, n)
! ALLOC_IVECTOR allocates the space for an allocatable INTEGER(IK)
! vector X, whose size is N after allocation.
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Input
integer(IK), intent(in) :: n

! Output
integer(IK), allocatable, intent(out) :: x(:)

! Local variable
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_IVECTOR'

! According to the Fortran 2003 standard, when a procedure is invoked,
! any allocated ALLOCATABLE object that is an actual argument associated
! with an INTENT(OUT) ALLOCATABLE dummy argument is deallocated. So it is
! unnecessary to write the following line since F2003 as X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
if (n >= 1) then
    x(1) = 0_IK
end if

end subroutine alloc_ivector


subroutine alloc_imatrix(x, m, n)
! ALLOC_IMATRIX allocates the space for a INTEGER(IK) matrix X, whose
! size is (M, N) after allocation.
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Input
integer(IK), intent(in) :: m, n

! Output
integer(IK), allocatable, intent(out) :: x(:, :)

! Local variable
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_IMATRIX'

! Unnecessary to write the following line since F2003 as X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(m, n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
if (m >= 1 .and. n >= 1) then
    x(1, 1) = 0_IK
end if

end subroutine alloc_imatrix


pure function size_of_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP, IK
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
use, non_intrinsic :: consts_mod, only : DP, IK
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
use, non_intrinsic :: consts_mod, only : QP, IK
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
