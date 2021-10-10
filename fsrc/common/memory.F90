#include "ppf.h"

module memory_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning memory management.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Sunday, October 10, 2021 AM10:19:58
!--------------------------------------------------------------------------------------------------!

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
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(RP) vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, HUGENUM
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
real(RP), allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RVECTOR'

! According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
! object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
! deallocated. So it is unnecessary to write the following line since F2003 as X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
x = HUGENUM

end subroutine alloc_rvector


subroutine alloc_rmatrix(x, m, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(RP) matrix X, whose size is (M, N) after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, HUGENUM
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Inputs
integer(IK), intent(in) :: m, n

! Outputs
real(RP), allocatable, intent(out) :: x(:, :)

! Local variables
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
x = HUGENUM

end subroutine alloc_rmatrix


subroutine alloc_ivector(x, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable INTEGER(IK) vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
integer(IK), allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_IVECTOR'

! Unnecessary to write the following line since F2003 as X is INTENT(OUT):
!!if (allocated(x)) deallocate (x)

! Allocate memory for X
allocate (x(n), stat=alloc_status)
if (alloc_status /= 0) then
    call errstop(srname, 'Memory allocation fails.')
end if

! Use X; otherwise, compilers may complain.
x = huge(0_IK)

end subroutine alloc_ivector


subroutine alloc_imatrix(x, m, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for a INTEGER(IK) matrix X, whose size is (M, N) after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
implicit none

! Inputs
integer(IK), intent(in) :: m, n

! Outputs
integer(IK), allocatable, intent(out) :: x(:, :)

! Local variables
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
x = huge(0_IK)

end subroutine alloc_imatrix


pure function size_of_sp(x) result(y)
!--------------------------------------------------------------------------------------------------!
! Return the storage size of X in Bytes, X being a REAL(SP) scalar.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : SP, IK
implicit none
! Inputs
real(SP), intent(in) :: x
! Outputs
integer(IK) :: y

#if __USE_STORAGE_SIZE__ == 1
! We prefer STORAGE_SIZE to C_SIZEOF, because the former is intrinsic while the later requires the
! intrinsic module ISO_C_BINDING.
y = int(storage_size(x) / 8, kind(y))  ! Y = INT(C_SIZEOF(X), KIND(Y))
#else
y = int(kind(x), kind(y)) ! Avoid complaint
y = int(4, kind(y))  ! This is not portable
#endif
end function size_of_sp


pure function size_of_dp(x) result(y)
!--------------------------------------------------------------------------------------------------!
! Return the storage size of X in Bytes, X being a REAL(DP) scalar.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DP, IK
implicit none
! Inputs
real(DP), intent(in) :: x
! Outputs
integer(IK) :: y

#if __USE_STORAGE_SIZE__ == 1
y = int(storage_size(x) / 8, kind(y))
#else
y = int(kind(x), kind(y)) ! Avoid complaint
y = int(8, kind(y))  ! This is not portable
#endif
end function size_of_dp


#if __QP_AVAILABLE__ == 1

pure function size_of_qp(x) result(y)
!--------------------------------------------------------------------------------------------------!
! Return the storage size of X in Bytes, X being a REAL(QP) scalar.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : QP, IK
implicit none
! Inputs
real(QP), intent(in) :: x
! Outputs
integer(IK) :: y

#if __USE_STORAGE_SIZE__ == 1
y = int(storage_size(x) / 8, kind(y))
#else
y = int(kind(x), kind(y)) ! Avoid complaint
y = int(16, kind(y))  ! This is not portable
#endif
end function size_of_qp

#endif


end module memory_mod
