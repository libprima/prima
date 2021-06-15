!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of newuoa_exmp.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 15-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! This is an example to illustrate the usage of NEWUOA.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code.

!!!!!!!!!!!!!!!!!! THE MODULE THAT IMPLEMENTS CALFUN !!!!!!!!!!!!!!!!!!!
      module calfun_mod

      implicit none
      private
      public :: calfun

      contains

      subroutine calfun(x, f)
! The Chebyquad test problem (Fletcher, 1965)
      implicit none

      real(kind(0.0D0)), intent(in) :: x(:)
      real(kind(0.0D0)), intent(out) :: f
      real(kind(0.0D0)) :: y(size(x) + 1, size(x) + 1), tmp
      integer :: i, n

      n = size(x)

      y(1:n, 1) = 1.0D0
      y(1:n, 2) = 2.0D0 * x - 1.0D0
      do i = 2, n
          y(1:n, i + 1) = 2.0D0 * y(1:n, 2) * y(1:n, i) - y(1:n, i - 1)
      end do

      f = 0.0D0
      do i = 1, n + 1
          tmp = sum(y(1:n, i)) / real(n, kind(0.0D0))
          if (mod(i, 2) /= 0) then
              tmp = tmp + 1.0D0 / real(i * i - 2 * i, kind(0.0D0))
          end if
          f = f + tmp * tmp
      end do
      end subroutine calfun

      end module calfun_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!! THE MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!
      program newuoa_exmp

! The following line makes NEWUOA available.
!----------------------------------------------------------------------!
      use newuoa_mod, only : newuoa
!----------------------------------------------------------------------!

! The following line specifies which module provides CALFUN.
! If CALFUN is provided by an external subroutine instead of a module,
! then remove this line and uncomment the "external calfun" line below.
!----------------------------------------------------------------------!
      use calfun_mod, only : calfun
!----------------------------------------------------------------------!

      implicit none

      integer :: i, n, alloc_stat
      real(kind(0.0D0)) :: rhobeg, f
      real(kind(0.0D0)), allocatable :: x(:)

! If CALFUN is provided as an external subroutine, then remove the line
! of "use calfun_mod, only : calfun", and uncomment the following line:
!----------------------------------------------------------------------!
!external calfun
!----------------------------------------------------------------------!

      do n = 2, 10, 2
! Sets up the initial X for the Chebyquad problem.
          if (allocated(x)) deallocate (x)
          allocate (x(n), stat=alloc_stat)
          if (alloc_stat /= 0) print *, 'Memory allocation failed.'
          do i = 1, n
              x(i) = real(i, kind(0.0D0)) / real(n + 1, kind(0.0D0))
          end do

          rhobeg = 0.2D0 * x(1)

          print '(/1A, I2)', 'Result with N = ', n

! The following line illustrates how to call NEWUOA.
!------------------------------------------------------------------!
          call newuoa(calfun, x, f, rhobeg=rhobeg, iprint=2)
!------------------------------------------------------------------!
! In additon to the required arguments CALFUN, X, and F, the above
! illustration specifies also RHOBEG and IPRINT, which are optional.
! All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in NEWUOA. You can also ignore all
! the optional arguments and invoke NEWUOA by the following line.
!------------------------------------------------------------------!
! call newuoa(calfun, x, f)
!------------------------------------------------------------------!
      end do

      end program newuoa_exmp