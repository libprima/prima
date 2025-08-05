#include "ppf.h"

module memory_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning memory management.
!
! In particular, the intrinsic ALLOCATE is wrapped into the procedure SAFEALLOC, which may be a
! controversial practice. We choose to do this because it has helped us a couple of times to locate
! bugs or problems in our code or even in compilers (e.g., Absoft). See the below for discussions:
! https://fortran-lang.discourse.group/t/best-practice-of-allocating-memory-in-fortran
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Sunday, November 13, 2022 PM02:05:50
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cstyle_sizeof
public :: safealloc

interface cstyle_sizeof
    module procedure size_of_sp, size_of_dp
#if __QP_AVAILABLE__ == 1
    module procedure size_of_qp
#endif
end interface cstyle_sizeof

interface safealloc
    module procedure alloc_lvector
    module procedure alloc_ivector, alloc_imatrix
    module procedure alloc_rvector_sp, alloc_rmatrix_sp
    module procedure alloc_rvector_dp, alloc_rmatrix_dp
#if __QP_AVAILABLE__ == 1
    module procedure alloc_rvector_qp, alloc_rmatrix_qp
#endif
end interface safealloc


contains


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


subroutine alloc_rvector_sp(x, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(SP) vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : SP, IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
real(SP), allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RVECTOR_SP'

! Preconditions (checked even not debugging)
call validate(n >= 0, 'N >= 0', srname)

! According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
! object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
! deallocated. So the following line is unnecessary since F2003 as X is INTENT(OUT):
! !if (allocated(x)) deallocate (x)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:n), stat=alloc_status)
x = -huge(x)  ! Costly if X is of a large size.
! N.B.: Do not write ALLOCATE (X(1:N), STAT=ALLOC_STATUS, SOURCE=-HUGE(X)), because
! 1. It is invalid to put X in the SOURCE specifier when it is being allocated;
! 2. Absoft does not support the SOURCE keyword as of 2022.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x) == n, 'SIZE(X) == N', srname)
end subroutine alloc_rvector_sp


subroutine alloc_rmatrix_sp(x, m, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(SP) matrix X, whose size is (M, N) after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : SP, IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: m, n

! Outputs
real(SP), allocatable, intent(out) :: x(:, :)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RMATRIX_SP'

! Preconditions (checked even not debugging)
call validate(m >= 0 .and. n >= 0, 'M >= 0, N >= 0', srname)

!if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:m, 1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x, 1) == m .and. size(x, 2) == n, 'SIZE(X) == [M, N]', srname)
end subroutine alloc_rmatrix_sp


subroutine alloc_rvector_dp(x, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(DP) vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DP, IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
real(DP), allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RVECTOR_DP'

! Preconditions (checked even not debugging)
call validate(n >= 0, 'N >= 0', srname)

! According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
! object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
! deallocated. So the following line is unnecessary since F2003 as X is INTENT(OUT):
! !if (allocated(x)) deallocate (x)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x) == n, 'SIZE(X) == N', srname)
end subroutine alloc_rvector_dp


subroutine alloc_rmatrix_dp(x, m, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(DP) matrix X, whose size is (M, N) after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DP, IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: m, n

! Outputs
real(DP), allocatable, intent(out) :: x(:, :)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RMATRIX_DP'

! Preconditions (checked even not debugging)
call validate(m >= 0 .and. n >= 0, 'M >= 0, N >= 0', srname)

!if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:m, 1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x, 1) == m .and. size(x, 2) == n, 'SIZE(X) == [M, N]', srname)
end subroutine alloc_rmatrix_dp


#if __QP_AVAILABLE__ == 1

subroutine alloc_rvector_qp(x, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(QP) vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : QP, IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
real(QP), allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RVECTOR_QP'

! Preconditions (checked even not debugging)
call validate(n >= 0, 'N >= 0', srname)

! According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
! object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
! deallocated. So the following line is unnecessary since F2003 as X is INTENT(OUT):
! !if (allocated(x)) deallocate (x)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x) == n, 'SIZE(X) == N', srname)
end subroutine alloc_rvector_qp


subroutine alloc_rmatrix_qp(x, m, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable REAL(QP) matrix X, whose size is (M, N) after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : QP, IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: m, n

! Outputs
real(QP), allocatable, intent(out) :: x(:, :)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_RMATRIX_QP'

! Preconditions (checked even not debugging)
call validate(m >= 0 .and. n >= 0, 'M >= 0, N >= 0', srname)

! !if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:m, 1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x, 1) == m .and. size(x, 2) == n, 'SIZE(X) == [M, N]', srname)
end subroutine alloc_rmatrix_qp

#endif


subroutine alloc_lvector(x, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable LOGICAL vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
logical, allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_LVECTOR'

! Preconditions (checked even not debugging)
call validate(n >= 0, 'N >= 0', srname)

! !if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
! Allocate memory for X. Initialize X to a compiler-independent value.
allocate (x(1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = .false.  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x) == n, 'SIZE(X) == N', srname)
end subroutine alloc_lvector


subroutine alloc_ivector(x, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for an allocatable INTEGER(IK) vector X, whose size is N after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
integer(IK), allocatable, intent(out) :: x(:)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_IVECTOR'

! Preconditions (checked even not debugging)
call validate(n >= 0, 'N >= 0', srname)

! !if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x) == n, 'SIZE(X) == N', srname)
end subroutine alloc_ivector


subroutine alloc_imatrix(x, m, n)
!--------------------------------------------------------------------------------------------------!
! Allocate space for a INTEGER(IK) matrix X, whose size is (M, N) after allocation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : validate
implicit none

! Inputs
integer(IK), intent(in) :: m, n

! Outputs
integer(IK), allocatable, intent(out) :: x(:, :)

! Local variables
integer :: alloc_status
character(len=*), parameter :: srname = 'ALLOC_IMATRIX'

! Preconditions (checked even not debugging)
call validate(m >= 0 .and. n >= 0, 'M >= 0, N >= 0', srname)

! !if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
! Allocate memory for X. Initialize X to a compiler-independent strange value.
allocate (x(1:m, 1:n), stat=alloc_status)  ! Absoft does not support the SOURCE keyword as of 2022.
x = -huge(x)  ! Costly if X is of a large size.

! Postconditions (checked even not debugging)
call validate(alloc_status == 0, 'Memory allocation succeeds (ALLOC_STATUS == 0)', srname)
call validate(allocated(x), 'X is allocated', srname)
call validate(size(x, 1) == m .and. size(x, 2) == n, 'SIZE(X) == [M, N]', srname)
end subroutine alloc_imatrix


end module memory_mod
