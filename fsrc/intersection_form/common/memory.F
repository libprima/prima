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
! on 09-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! MEMORY_MOD is a module providing subroutines concerning memory management.


      module memory_mod

      implicit none
      private
      public :: safealloc

      interface safealloc
          module procedure alloc_rvector, alloc_rmatrix
      end interface safealloc


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


      end module memory_mod