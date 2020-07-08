! LINA is a module providing some basic linear algebra procedures.


#include "ppflags.h"

      module lina

      implicit none
      private
      public :: r1update, r2update, outprod

      interface r1update
          module procedure r1_with_scalar, r1_without_scalar
      end interface r1update

      interface r2update
          module procedure r2_with_scalar, r2_without_scalar
      end interface r2update

#if __USE_INTRINSIC_ALGEBRA__ == 0 
      public :: matmul, dot_product

      interface matmul
          module procedure matmul12, matmul21, matmul22
      end interface matmul

      interface dot_product
          module procedure inprod
      end interface dot_product
#endif


#if __IMPROVE_POWELL_CODE__ == 0 || __FORTRAN_STANDARD__ < 2008
      public :: hypot
      interface hypot
          module procedure hypotenuse
      end interface hypot
#endif


      contains

      subroutine r1_with_scalar(A, alpha, x, y)
      ! R1_WITH_SCALAR sets 
      ! A = A + ALPHA*( X*Y^T ), 
      ! where A is an MxN matrix, ALPHA is a real scalar, X is an
      ! M-dimenional vector, and Y is an N-dimenional vector.
      ! Note that we should calculate the increment as ALPHA*( X*Y^T )
      ! rather than ALPHA*X*Y^T. Because, when X = Y, ALPHA*( X*Y^T ) 
      ! is a symmetric matrix while ALPHA*X*Y^T is asymmetric due to 
      ! the non-associativity of floating-point multiplication.
      use consts, only : RP, DEBUG_MODE
      use warnerror, only : errmssg
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer :: j 

      if (DEBUG_MODE) then
          if (size(A, 1) /= size(y) .or. size(A, 2) /= size(y)) then
              call errmssg('R1UPDATE', 'SIZE(A) is invalid')
              stop
          end if
      end if

      do j = 1, size(y)
          A(:, j) = A(:, j) + alpha*( x*y(j) )
          ! Here we should calculate the increment as ALPHA*( X*Y(J) )
          ! rather than ALPHA*X*Y(J). See the comments above.
      end do
      end subroutine r1_with_scalar


      subroutine r1_without_scalar(A, x, y)
      ! R1_WITHOUT_SCALAR sets
      ! A = A + X*Y^T, 
      ! where A is an MxN matrix, X is an M-dimenional vector, and Y is
      ! an N-dimenional vector.
      use consts, only : RP, DEBUG_MODE
      use warnerror, only : errmssg
      implicit none
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer :: j 

      if (DEBUG_MODE) then
          if (size(A, 1) /= size(y) .or. size(A, 2) /= size(y)) then
              call errmssg('R1UPDATE', 'SIZE(A) is invalid')
              stop
          end if
      end if

      do j = 1, size(y)
          A(:, j) = A(:, j) + x*y(j)
      end do
      end subroutine r1_without_scalar

      subroutine r2_with_scalar(A, alpha, x, y, beta, u, v)
      ! R2_WITH_SCALAR sets 
      ! A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ), 
      ! where A is an MxN matrix, ALPHA and BETA are real scalars, X and 
      ! U are M-dimenional vectors, Y and V are N-dimenional vectors.
      ! Note:
      ! 1. We should calculate the increment as 
      ! ALPHA*( X*Y^T ) + BETA*( U*V^T )
      ! rather than 
      ! ALPHA*X*Y^T + BETA*U*V^T
      ! to ensure the symmetry when X = Y and U = V.
      ! 2. We should set A to 
      ! A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) )
      ! rather than 
      ! A + ALPHA*( X*Y^T ) + BETA*( U*V^T )
      ! to ensure the semmetry when ALPHA = BETA, X = V, and Y = U.
      use consts, only : RP, DEBUG_MODE
      use warnerror, only : errmssg
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: beta 
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(in) :: u(:)  ! U(SIZE(X))
      real(RP), intent(in) :: v(:)  ! V(SIZE(Y))
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer :: j 

      if (DEBUG_MODE) then
          if (size(u) /= size(x)) then
              call errmssg('R2_WITH_SCALAR', 'SIZE(U) /= SIZE(X)')
              stop
          end if
          if (size(v) /= size(y)) then
              call errmssg('R2_WITH_SCALAR', 'SIZE(V) /= SIZE(Y)')
              stop
          end if
          if (size(A, 1) /= size(y) .or. size(A, 2) /= size(y)) then
              call errmssg('R2_WITH_SCALAR', 'SIZE(A) is invalid')
              stop
          end if
      end if

      do j = 1, size(y)
          A(:, j) = A(:, j) + ( alpha*( x*y(j) ) + beta*(u*v(j)) )
      end do
      end subroutine r2_with_scalar

      subroutine r2_without_scalar(A, x, y, u, v)
      ! R2_WITHOUT_SCALAR sets 
      ! A = A + ( X*Y^T + U*V^T ), 
      ! where A is an MxN matrix, X and U are M-dimenional vectors, 
      ! Y and V are N-dimenional vectors.
      ! Note that we should set A to 
      ! A + ( X*Y^T +  U*V^T )
      ! rather than 
      ! A + X*Y^T +  U*V^T 
      ! to ensure the semmetry when X = V and Y = U.
      use consts, only : RP, DEBUG_MODE
      use warnerror, only : errmssg
      implicit none
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(in) :: u(:)  ! U(SIZE(X))
      real(RP), intent(in) :: v(:)  ! V(SIZE(Y))
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer :: j 

      if (DEBUG_MODE) then
          if (size(u) /= size(x)) then
              call errmssg('R2_WITH_SCALAR', 'SIZE(U) /= SIZE(X)')
              stop
          end if
          if (size(v) /= size(y)) then
              call errmssg('R2_WITH_SCALAR', 'SIZE(V) /= SIZE(Y)')
              stop
          end if
          if (size(A, 1) /= size(y) .or. size(A, 2) /= size(y)) then
              call errmssg('R2_WITH_SCALAR', 'SIZE(A) is invalid')
              stop
          end if
      end if

      do j = 1, size(y)
          A(:, j) = A(:, j) + ( x*y(j) + u*v(j) )
      end do
      end subroutine r2_without_scalar


#if __USE_INTRINSIC_ALGEBRA__ == 0 

      function matmul12(x, y) result(z)
      use consts, only : RP, DEBUG_MODE
      use warnerror, only : errmssg 
      implicit none
      real(RP), dimension(:), intent(in) :: x
      real(RP), dimension(:, :), intent(in) :: y
      real(RP), dimension(size(y, 2)) :: z
      integer :: j
      
      if (DEBUG_MODE) then
          if (size(x) /= size(y, 1)) then
              call errmssg('MATMUL12', 'SIZE(X) /= SIZE(Y, 1)')
              stop
          end if
      end if

      do j = 1, size(y, 2)
          z(j) = dot_product(x, y(:, j))
      end do
      end function matmul12

      function matmul21(x, y) result(z)
      use consts, only : RP, ZERO, DEBUG_MODE
      use warnerror, only : errmssg
      implicit none
      real(RP), dimension(:, :), intent(in) :: x
      real(RP), dimension(:), intent(in) :: y
      real(RP), dimension(size(x, 1)) :: z
      integer :: j 

      if (DEBUG_MODE) then
          if (size(x, 2) /= size(y)) then
              call errmssg('MATMUL21', 'SIZE(X, 2) /= SIZE(Y)')
              stop 
          end if
      end if

      z = ZERO
      do j = 1, size(x, 2)
          z = z + x(:, j)*y(j)
      end do
      end function matmul21

      function matmul22(x, y) result(z)
      use consts, only : RP, ZERO, DEBUG_MODE
      use warnerror, only : errmssg 
      implicit none
      real(RP), dimension(:, :), intent(in) :: x, y
      real(RP), dimension(size(x, 1), size(y, 2)) :: z
      integer :: i, j

      if (DEBUG_MODE) then
          if (size(x, 2) /= size(y, 1)) then
              call errmssg('MATMUL22', 'SIZE(X, 2) /= SIZE(Y, 1)')
              stop 
          end if
      end if

      z = ZERO
      do j = 1, size(y, 2)
          do i = 1, size(x, 2)
              z(:, j) = z(:, j) + x(:, i)*y(i, j)
          end do
      end do
      end function matmul22

      function inprod(x, y) result(z)
      use consts, only : RP, ZERO, DEBUG_MODE
      use warnerror, only : errmssg 
      implicit none
      real(RP), dimension(:), intent(in) :: x, y
      real(RP) :: z
      integer :: i

      if (DEBUG_MODE) then
          if (size(x) /= size(y)) then
              call errmssg('INPROD', 'SIZE(X) /= SIZE(Y)')
              stop 
          end if
      end if

      !z = sum(x*y)
      z = ZERO
      do i = 1, size(x)
          z = z + x(i)*y(i)
      end do
      end function inprod

#endif

      function hypotenuse(x, y) result(r)
      ! HYPOTENUSE calculate sqrt(x^2+y^2) if __IMPROVE_POWELL_CODE__ 
      ! is 0 or the Fortran standard is lower than 2008.

#if __IMPROVE_POWELL_CODE__ == 0
      use consts, only : RP
#else
      use consts, only : RP, ONE, TWO, ZERO
#endif

      implicit none
      real(RP), intent(in) :: x, y
      real(RP) :: r  

#if __IMPROVE_POWELL_CODE__ == 0
      r = sqrt(x**2 + y**2)  ! May encounter over/under flow
#else
      real(RP) :: ax, ay
      ax = abs(x)
      ay = abs(y)
      if (ax > ay) then  ! AX > AY ensures AX > 0.
          r = ax*sqrt(ONE + (ay/ax)**2)
      else if (ax < ay) then  ! AX < AY ensures AY > 0.
          r = ay*sqrt(ONE + (ax/ay)**2)
      else if (ax > ZERO) then
          r = ax*sqrt(TWO)
      else 
          r = ZERO
      end if
#endif
      end function hypotenuse

      function outprod(x, y) result(z)
      use consts, only : RP
      implicit none
      real(RP), dimension(:), intent(in) :: x, y
      real(RP), dimension(size(x), size(y)) :: z
      integer :: i
      do i = 1, size(y)
          z(:, i) = x*y(i)
      end do
      end function outprod


      end module lina
