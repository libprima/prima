!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of memory.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 15-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! MEMORY_MOD is a module providing subroutines concerning memory management.


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

! Intermediate variable
      integer :: alloc_status
      character(len = SRNLEN), parameter :: srname = 'ALLOC_RVECTOR'

! Deallocate X if X is allocated. This is necessary.
      if (allocated(x)) then
          deallocate (x)
      end if

! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          call errstop(srname, 'Memory allocation fails.')
      end if

      x = 0.0_RP ! Set = 0; otherieise, compilers may complain.

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

! Intermediate variable
      integer :: alloc_status
      character(len = SRNLEN), parameter :: srname = 'ALLOC_RMATRIX'

! Deallocate X if X is allocated. This is necessary.
      if (allocated(x)) then
          deallocate (x)
      end if

! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          call errstop(srname, 'Memory allocation fails.')
      end if

      x = 0.0_RP ! Set X = 0; otherieise, compilers may complain.

      end subroutine alloc_rmatrix


      pure function size_of_sp(x) result(y)
      use consts_mod, only : SP, IK
      implicit none
      real(SP), intent(in) :: x
      integer(IK) :: y
#if __USE_STORAGE_SIZE__ == 1
          y = int(storage_size(x)/8, kind(y))
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
          y = int(storage_size(x)/8, kind(y))
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
          y = int(storage_size(x)/8, kind(y))
#else
          y = int(kind(x), kind(y)) ! Avoid complaint
          y = int(8, kind(y))
#endif
      end function size_of_qp

#endif


      end module memory_mod
