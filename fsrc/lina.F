! LINA is a module providing some basic linear algebra procedures.
#include "ppflags.h"

      module lina

      use consts, only : rp, one, two, zero
      implicit none

#if __USE_INTRINSIC_ALGEBRA__ == 0 
      interface matmul
          module procedure matmul12, matmul21, matmul22
      end interface matmul

      interface dot_product
          module procedure inprod
      end interface dot_product
#endif

#if __IMPROVE_POWELL_CODE__ == 0 || __FORTRAN_STANDARD__ < 2008
      interface hypot
          module procedure hypotenuse
      end interface hypot
#endif


      contains

#if __USE_INTRINSIC_ALGEBRA__ == 0 
      function matmul12(x, y) result(z)
      use warnerror, only : errmssg 
      implicit none
      real(kind = rp), dimension(:), intent(in) :: x
      real(kind = rp), dimension(:, :), intent(in) :: y
      real(kind = rp), dimension(size(y, 2)) :: z
      integer :: j
      if (size(x) /= size(y, 1)) then
          call errmssg('MATMUL12', 'SIZE(X) /= SIZE(Y, 1)')
          stop
      end if
      z = zero
      do j = 1, size(y, 2)
          z(j) = dot_product(x, y(:, j))
          !do i = 1, size(x)
          !    z(j) = z(j) + x(i)*y(i, j)
          !end do
      end do
      end function matmul12

      function matmul21(x, y) result(z)
      use warnerror, only : errmssg
      implicit none
      real(kind = rp), dimension(:, :), intent(in) :: x
      real(kind = rp), dimension(:), intent(in) :: y
      real(kind = rp), dimension(size(x, 1)) :: z
      integer :: i
      if (size(x, 2) /= size(y)) then
          call errmssg('MATMUL21', 'SIZE(X, 2) /= SIZE(Y)')
          stop 
      end if
      z = zero
      do i = 1, size(x, 2)
          z = z + x(:, i)*y(i)
      end do
      end function matmul21

      function matmul22(x, y) result(z)
      use warnerror, only : errmssg 
      implicit none
      real(kind = rp), dimension(:, :), intent(in) :: x, y
      real(kind = rp), dimension(size(x, 1), size(y, 2)) :: z
      integer :: i, j
      if (size(x, 2) /= size(y, 1)) then
          call errmssg('MATMUL22', 'SIZE(X, 2) /= SIZE(Y, 1)')
          stop 
      end if
      z = zero
      do j = 1, size(y, 2)
          do i = 1, size(x, 2)
              z(:, j) = z(:, j) + x(:, i)*y(i, j)
          end do
      end do
      end function matmul22

      function inprod(x, y) result(z)
      use warnerror, only : errmssg 
      implicit none
      real(kind = rp), dimension(:), intent(in) :: x, y
      real(kind = rp) :: z
      !z = sum(x*y)
      integer :: i
      if (size(x) /= size(y)) then
          call errmssg('INPROD', 'SIZE(X) /= SIZE(Y)')
          stop 
      end if
      z = zero
      do i = 1, size(x)
          z = z + x(i)*y(i)
      end do
      end function inprod
#endif

      function hypotenuse(x, y) result(r)
      ! HYPOTENUSE calculate sqrt(x^2+y^2) if __IMPROVE_POWELL_CODE__ 
      ! is 0 or the Fortran standard is lower than 2008.
#if __IMPROVE_POWELL_CODE__ != 0
      use consts, only : rp, one, two, zero
#endif
      implicit none
      real(rp), intent(in) :: x, y
      real(rp) :: r  
#if __IMPROVE_POWELL_CODE__ == 0
      r = sqrt(x**2 + y**2)  ! May encounter over/under flow
#else
      real(rp) :: ax, ay
      ax = abs(x)
      ay = abs(y)
      if (ax > ay) then  ! AX > AY ensures AX > 0.
          r = ax*sqrt(one + (ay/ax)**2)
      else if (ax < ay) then  ! AX < AY ensures AY > 0.
          r = ay*sqrt(one + (ax/ay)**2)
      else if (ax > zero) then
          r = ax*sqrt(two)
      else 
          r = zero
      end if
#endif
      end function hypotenuse

      function outprod(x, y) result(z)
      implicit none
      real(kind = rp), dimension(:), intent(in) :: x, y
      real(kind = rp), dimension(size(x), size(y)) :: z
      integer :: i
      do i = 1, size(y)
          z(:, i) = x*y(i)
      end do
      end function outprod


      end module lina
