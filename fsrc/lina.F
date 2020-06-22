! LINA is a module providing some basic linear algebra procedures.

! Do we use the intrinsic linear algebra functions?
#ifdef __USE_INTRINSIC_ALGEBRA__
#undef __USE_INTRINSIC_ALGEBRA__
#endif
#define __USE_INTRINSIC_ALGEBRA__ 0 


      module lina

      use consts, only : rp, zero
      implicit none

#if __USE_INTRINSIC_ALGEBRA__ == 0 
      interface matmul
          module procedure matmul12, matmul21, matmul22
      end interface matmul

      interface dot_product
          module procedure inprod
      end interface dot_product
#endif

      interface outprod
          module procedure outprod
      end interface outprod

      contains

#if __USE_INTRINSIC_ALGEBRA__ == 0 
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

      function outprod(x, y)
      implicit none
      real(kind = rp), dimension(:), intent(in) :: x, y
      real(kind = rp), dimension(size(x), size(y)) :: outprod 
      integer :: i
      do i = 1, size(y)
          outprod(:, i) = x*y(i)
      end do
      end function outprod

      end module lina
