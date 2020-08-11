!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of example.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 11-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine calfun (x, f)
!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
!     with NPT = 2N+1.
!
      implicit none
      real(kind(0.0D0)), intent(in) :: x(:)
      real(kind(0.0D0)), intent(out) :: f
      real(kind(0.0D0)) :: y(size(x) + 1, size(x) + 1), sum

      integer :: i, j, n, iw

      n = size(x)

      do j = 1, n
          y(1, j) = 1.0D0
          y(2, j) = 2.0D0*x(j) - 1.0D0
      end do
      do i = 2, n
          do j = 1, n
              y(i+1,j) = 2.0D0*y(2, j)*y(i, j) - y(i - 1, j)
          end do
      end do
      f=0.0D0
      iw=1
      do i = 1, n + 1
          sum = 0.0D0
          do j=1, n
              sum=sum+y(i,j)
          end do
          sum=sum/real(n, kind(0.0D0))
          if (iw .gt. 0) sum = sum + 1.0d0/real(i*i - 2*i, kind(0.0D0))
          iw = -iw
          f = f+sum*sum
      end do
      end subroutine calfun


      program example

      use newuoa_mod, only : newuoa
      implicit none

      integer :: i, iprint, n, alloc_stat
      real(kind(0.0D0)) :: rhobeg, f
      real(kind(0.0D0)), allocatable :: x(:)

!interface
!    subroutine calfun(x, f)
!    real(kind(0.0D0)), intent(in) :: x(:)
!    real(kind(0.0D0)), intent(out) :: f
!    end subroutine calfun
!end interface

      external calfun

      iprint = 2

      do n = 2, 8, 2
          if (allocated(x)) deallocate(x)
          allocate(x(n), stat = alloc_stat)
          if (alloc_stat /= 0) print *, 'Memory allocation failed.'
          do i = 1, n
              x(i) = real(i, kind(0.0D0))/real(n + 1, kind(0.0D0))
          end do

          rhobeg = 0.2D0*x(1)

          print '(//1X, 1A, I2)', 'Results with N = ', n

          call newuoa(calfun, x, f, rhobeg = rhobeg, iprint = iprint)
      end do

      end program example
