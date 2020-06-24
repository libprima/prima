! LINA is a module providing some basic linear algebra procedures.

! Do we use the same linear algebra procedures as Powell did?
#ifdef __USE_POWELL_ALGEBRA__
#undef __USE_POWELL_ALGEBRA__
#endif
#define __USE_POWELL_ALGEBRA__ 1 


      module lina

      use consts, only : rp, one, two, zero
      implicit none

#if __USE_POWELL_ALGEBRA__ == 1 
      interface matmul
          module procedure matmul12, matmul21, matmul22
      end interface matmul

      interface dot_product
          module procedure inprod
      end interface dot_product
#endif


      contains

#if __USE_POWELL_ALGEBRA__ == 1 
      function matmul12(x, y)
      implicit none
      real(kind = rp), dimension(:), intent(in) :: x
      real(kind = rp), dimension(:, :), intent(in) :: y
      real(kind = rp), dimension(size(y, 2)) :: matmul12
      integer :: i, j
      if (size(x) /= size(y, 1)) then
          stop 'Error in MATMUL: RANK(X) = 1 and SIZE(X) /= SIZE(Y, 1).'
      end if
      matmul12 = zero
      do j = 1, size(y, 2)
          do i = 1, size(x)
              matmul12(j) = matmul12(j) + x(i)*y(i, j)
          end do
      end do
      end function matmul12

      function matmul21(x, y)
      implicit none
      real(kind = rp), dimension(:, :), intent(in) :: x
      real(kind = rp), dimension(:), intent(in) :: y
      real(kind = rp), dimension(size(x, 1)) :: matmul21 
      integer :: i
      if (size(x, 2) /= size(y)) then
          stop 'Error in MATMUL: RANK(Y) = 1 and SIZE(X, 2) /= SIZE(Y).'
      end if
      matmul21 = zero
      do i = 1, size(x, 2)
          matmul21 = matmul21 + x(:, i)*y(i) 
      end do
      end function matmul21

      function matmul22(x, y)
      implicit none
      real(kind = rp), dimension(:, :), intent(in) :: x, y
      real(kind = rp), dimension(size(x, 1), size(y, 2)) :: matmul22 
      integer :: i, j
      if (size(x, 2) /= size(y, 1)) then
          stop 'Error in MATMUL: SIZE(X, 2) /= SIZE(Y, 1).'
      end if
      matmul22 = zero
      do j = 1, size(y, 2)
          do i = 1, size(x, 2)
              matmul22(:, j) = matmul22(:, j) + x(:, i)*y(i, j) 
          end do
      end do
      end function matmul22

      function inprod(x, y)
      implicit none
      real(kind = rp), dimension(:), intent(in) :: x, y
      real(kind = rp) :: inprod 
      !inprod = sum(x*y)
      integer :: i
      if (size(x) /= size(y)) then
          stop 'Error in INPROD: SIZE(X) /= SIZE(Y).'
      end if
      inprod = zero
      do i = 1, size(x)
          inprod = inprod + x(i)*y(i)
      end do
      end function inprod
#endif

!      function outprod(x, y)
!      implicit none
!      real(kind = rp), dimension(:), intent(in) :: x, y
!      real(kind = rp), dimension(size(x), size(y)) :: outprod 
!      integer :: i
!      do i = 1, size(y)
!          outprod(:, i) = x*y(i)
!      end do
!      end function outprod

      subroutine givens(x, y, c, s, r)
      ! GIVENS calculates the parameters for the Givens rotation that
      ! [ c s] [x y]' = [r, 0]'
      ! [-s c] 
      ! Explicitly,
      ! c = x/sqrt(x^2+y^2)
      ! s = y/sqrt(x^2+y^2)
      ! r = sqrt(x^2+y^2)
      implicit none
      real(kind = rp), intent(in) :: x, y
      real(kind = rp), intent(out) :: c, s, r 
#if __USE_POWELL_ALGEBRA__ == 1
      ! Powell uses the following implementation.
      r = sqrt(x**2 + y**2)
      if (r > 0) then
          c = x/r
          s = y/r
      else
          c = sqrt(two)/two
          s = sqrt(two)/two
      end if
#else
      ! If calculated directly, r may suffer from over/underflow due to
      ! the squres. That's why we use the folloing scheme.  
      real(kind = rp) :: t, u
      if (abs(x) > abs(y)) then
          t = y/x
          u = sign(sqrt(one + t*t), x)
          ! SIGN is the sign copying funtion:
          ! SIGN(a, b) = (sign of b)*|a|, 
          ! where 'sign of b' is 1 if b > 0 and -1 if b < 0. 
          ! The behaviour of SIGN(., 0) depends on the compiler. 
          ! NOTE: MATLAB has also a SIGN function, but the definition is 
          ! SIGN(a) = 1 if a > 0, 0 if a = 0, and -1 if a < 0.
          c = one/u  
          s = c*t  
          r = x*u  
      else if (abs(x) < abs(y)) then
          t = x/y
          u = sign(sqrt(one + t*t), y)
          s = one/u
          c = s*t
          r = y*u
      else if (abs(x) > zero) then
          c = sign(sqrt(two)/two, x)
          s = sign(sqrt(two)/two, y)
          r = abs(x)*sqrt(two)
      else 
      ! This case may be combined with the last case. Here we implement
      ! it separately to avoid using SIGN(., 0), whose behaviour 
      ! depends on the compiler.
          c = sqrt(two)/two
          s = sqrt(two)/two
          r = zero
      end if
#endif
      end subroutine 

      end module lina
