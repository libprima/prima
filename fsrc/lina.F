! LINA is a module providing some basic linear algebra procedures.
! To improve the performance of these procedures, especially matmul, one
! can customize their implementations, according to the resources 
! (hardware, e.g., cache, and libraries, e.g., BLAS) available and the
! sizes of the matrices/vectors. 


#include "ppf.h"

      module lina_mod

      implicit none
      private
      public :: r1update, r2update
      public :: symmetrize


      interface r1update
          module procedure r1_sym, r1
      end interface r1update

      interface r2update
          module procedure r2_sym, r2
      end interface r2update

! When interfaced with MATLAB, the intrinsic matmul and dot_product seem
! not as efficient as the implementations below (mostly by loops). This
! may well depend on the machine (e.g., cache size), compiler, compiling
! options, and MATLAB version.
#if __USE_INTRINSIC_ALGEBRA__ == 0 
      public :: matmul, dot_product

      interface matmul
          module procedure matmul12, matmul21, matmul22
      end interface matmul

      interface dot_product
          module procedure inprod
      end interface dot_product
#endif


! Check whether HYPOT is available as an intrinsic.
! F2008 provides the intrinsic HYPOT; MATLAB, Python, Julia, R, etc.,
! contain built-in functions corresponding to HYPOT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef __HYPOT_AVAILABLE__
#undef __HYPOT_AVAILABLE__
#endif
#define __HYPOT_AVAILABLE__ 0 
! HYPOT is included in F2008.
#if __FORTRAN_STANDARD__ >= 2008 
#undef __HYPOT_AVAILABLE__
#define __HYPOT_AVAILABLE__ 1
#endif
! HYPOT is included in gfortran since 4.4.7 
#ifdef __GFORTRAN__
#if __GNUC__ >= 5
#undef __HYPOT_AVAILABLE__
#define __HYPOT_AVAILABLE__ 1
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if __IMPROVE_POWELL_CODE__ == 0 || __HYPOT_AVAILABLE__ == 0
      public :: hypot
      interface hypot
          module procedure hypotenuse
      end interface hypot
#endif


      contains

      subroutine r1_sym(A, alpha, x)
      ! R1_SYM sets 
      ! A = A + ALPHA*( X*X^T ), 
      ! where A is an NxN matrix, ALPHA is a scalar, and X is an
      ! N-dimenional vector.
      use consts_mod, only : RP, IK, ONE, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))
      integer(IK) :: n, j 
      character(len = SRNLEN), parameter :: srname ='R1_SYM'
      ! Be careful with initialization! 
      ! In Fortran >=90, the initialization in the declaration implies
      ! the "save" attribute. If the variable is not a parameter, it may
      ! casue unwanted behavior.

      if (DEBUG_MODE) then
          if (size(A, 1) /= size(x) .or. size(A, 2) /= size(x)) then
              call errstop(srname, 'SIZE(A) is invalid')
          end if
      end if

      n = size(x)

#if __IMPROVE_POWELL_CODE__ == 0
      ! Only update the LOWER TRIANGULAR part of A.
      ! Both of the following cases are invoked in NEWUOA.
      if (abs(alpha - ONE) > ZERO) then
          do j = 1, n
              A(j : n, j) = A(j : n, j) + alpha*x(j : n)*x(j)
          end do
      else
          do j = 1, n
              A(j : n, j) = A(j : n, j) + x(j : n)*x(j)
          end do
      end if
      ! Set A(UPPER_TRIANGULAR) = A(LOWER_TRIANGULAR).
      call symmetrize(A) 
#else
      ! Here we calculate the increment as ALPHA*( X*X^T ), which 
      ! is guaranteed symmetric even in floating-point arithmetic.
      ! Both of the following cases are invoked in NEWUOA.
      if (abs(alpha - ONE) > ZERO) then
          do j = 1, n
              A(:, j) = A(:, j) + alpha*( x*x(j) )
          end do
      else
          do j = 1, n
              A(:, j) = A(:, j) + x*x(j)
          end do
      end if
#endif          
      if (DEBUG_MODE) then
          call verisym(A, ZERO)
      end if
      end subroutine r1_sym

      subroutine r1(A, alpha, x, y)
      ! R1 sets 
      ! A = A + ALPHA*( X*Y^T ), 
      ! where A is an MxN matrix, ALPHA is a real scalar, X is an
      ! M-dimenional vector, and Y is an N-dimenional vector.
      use consts_mod, only : RP, IK, ONE, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer(IK) :: j 
      character(len = SRNLEN), parameter :: srname = 'R1'
      ! Be careful with initialization! 
      ! In Fortran >=90, the initialization in the declaration implies
      ! the "save" attribute. If the variable is not a parameter, it may
      ! casue unwanted behavior.

      if (DEBUG_MODE) then
          if (size(A, 1) /= size(x) .or. size(A, 2) /= size(y)) then
              call errstop(srname, 'SIZE(A) is invalid')
          end if
      end if

#if __IMPROVE_POWELL_CODE__ == 0
      ! Both of the following cases are invoked in NEWUOA.
      if (abs(alpha - ONE) > ZERO) then
          do j = 1, size(y)
              A(:, j) = A(:, j) + alpha*x*y(j)
          end do
      else
          do j = 1, size(y)
              A(:, j) = A(:, j) + x*y(j)
          end do
      end if
#else
      ! Here we calculate the increment as ALPHA*( X*Y^T ) rather than 
      ! ALPHA*X*Y^T to ensure that, when X = Y, the increment is 
      ! symmetric even in floating-point arithmetic. 
      ! Both of the following cases are invoked in NEWUOA.
      if (abs(alpha - ONE) > ZERO) then
          do j = 1, size(y)
              A(:, j) = A(:, j) + alpha*( x*y(j) )
          end do
      else
          do j = 1, size(y)
              A(:, j) = A(:, j) + x*y(j)
          end do
      end if
#endif
      end subroutine r1


      subroutine r2_sym(A, alpha, x, y)
      ! R2_SYM sets 
      ! A = A + ALPHA*( X*Y^T + Y*X^T ), 
      ! where A is an NxN matrix, X and Y are N-dimenional vectors, and
      ! alpha is a sclalar.
      use consts_mod, only : RP, IK, ONE, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))
      integer(IK) :: n, j 
      character(len = SRNLEN), parameter :: srname = 'R2_SYM'

      if (DEBUG_MODE) then
          if (size(x) /= size(y)) then
              call errstop(srname, 'SIZE(X) /= SIZE(Y)')
          end if
          if (size(A, 1) /= size(x) .or. size(A, 2) /= size(x)) then
              call errstop(srname, 'SIZE(A) is invalid')
          end if
      end if

      n = size(x)

#if __IMPROVE_POWELL_CODE__ == 0
      ! Only update the LOWER TRIANGULAR part of A.
      if (abs(alpha - ONE) > ZERO) then  ! NOT invoked in NEWUOA.
          do j = 1, n
              A(j : n, j) = A(j : n, j) + alpha*x(j : n)*y(j) +         &
     &         alpha*y(j : n)*x(j)
          end do
      else  ! This case is invoked in NEWUOA.
          do j = 1, n
              A(j : n, j) = A(j : n, j) + x(j : n)*y(j) + y(j : n)*x(j)
          end do
      end if
      ! Set A(UPPER_TRIANGULAR) = A(LOWER_TRIANGULAR).
      call symmetrize(A)
#else
      ! Here we calculate the increment as ALPHA*( X*Y^T + Y*X^T), which
      ! is guaranteed symmetric even in floating-point arithmetic. 
      if (abs(alpha - ONE) > ZERO) then  ! NOT invoked in NEWUOA.
          do j = 1, n
              A(:, j) = A(:, j) + alpha*( x*y(j) + y*x(j) )
          end do
      else  ! This case is invoked in NEWUOA.
          do j = 1, n
              A(:, j) = A(:, j) + ( x*y(j) + y*x(j) )
          end do
      end if
#endif
      if (DEBUG_MODE) then
          call verisym(A, ZERO)
      end if
      end subroutine r2_sym

      subroutine r2(A, alpha, x, y, beta, u, v)
      ! R2 sets 
      ! A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ), 
      ! where A is an MxN matrix, ALPHA and BETA are real scalars, X and 
      ! U are M-dimenional vectors, Y and V are N-dimenional vectors.
      use consts_mod, only : RP, IK, ONE, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: beta 
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(in) :: u(:)  ! U(SIZE(X))
      real(RP), intent(in) :: v(:)  ! V(SIZE(Y))
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer(IK) :: j 
      character(len = SRNLEN), parameter :: srname = 'R2'

      if (DEBUG_MODE) then
          if (size(u) /= size(x)) then
              call errstop(srname, 'SIZE(U) /= SIZE(X)')
          end if
          if (size(v) /= size(y)) then
              call errstop(srname, 'SIZE(V) /= SIZE(Y)')
          end if
          if (size(A, 1) /= size(x) .or. size(A, 2) /= size(y)) then
              call errstop(srname, 'SIZE(A) is invalid')
          end if
      end if

#if __IMPROVE_POWELL_CODE__ == 0 
      if (abs(alpha - ONE) > ZERO .or. abs(beta - ONE) > ZERO) then
          ! This case is NOT invoked in NEWUOA.
          do j = 1, size(y)
              A(:, j) = A(:, j) + alpha*x*y(j) + beta*u*v(j)
          end do
      else  ! This case is invoked in NEWUOA.
          do j = 1, size(y)
              A(:, j) = A(:, j) + x*y(j) + u*v(j)
          end do
      end if
#else
      if (abs(alpha - beta) > ZERO) then  ! NOT invoked in NEWUOA.
          do j = 1, size(y)
              A(:, j) = A(:, j) + (alpha*(x*y(j)) + beta*(u*v(j)))
          end do
      else if (abs(alpha - ONE) > ZERO) then  ! NOT invoked in NEWUOA.
      ! Here update A to A + ALPHA*( X*Y^T + U*V^T ), a semmetric update
      ! when X = V and Y = U, even in floating-point arithmetic.
          do j = 1, size(y)
              A(:, j) = A(:, j) + alpha*( x*y(j)  + u*v(j) )
          end do
      else  ! This case is invoked in NEWUOA.
      ! Here update A to A + ( X*Y^T + U*V^T ), a semmetric update 
      ! when X = V and Y = U, even in floating-point arithmetic.
          do j = 1, size(y)
              A(:, j) = A(:, j) + ( x*y(j) + u*v(j) )
          end do
      end if
#endif
      end subroutine r2


#if __USE_INTRINSIC_ALGEBRA__ == 0 

      function matmul12(x, y) result(z)
      use consts_mod, only : RP, IK, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), dimension(:), intent(in) :: x
      real(RP), dimension(:, :), intent(in) :: y
      real(RP), dimension(size(y, 2)) :: z
      integer(IK) :: j
      character(len = SRNLEN), parameter :: srname = 'MATMUL12' 
      
      if (DEBUG_MODE) then
          if (size(x) /= size(y, 1)) then
              call errstop(srname, 'SIZE(X) /= SIZE(Y, 1)')
          end if
      end if

      do j = 1, size(y, 2)
          ! When interfaced with MATLAB, the following seems more
          ! efficient than a loop, which is strange because inprod
          ! itself is implemented by a loop. This may well depend on the
          ! machine (e.g., cache size), compiler, compiling options,
          ! and MATLAB version.
          z(j) = inprod(x, y(:, j))
          ! Here we can also write dot_prod(x, y(:, j)), it will still
          ! call inprod according to the generic interface in lina_mod.
      end do
      end function matmul12

      function matmul21(x, y) result(z)
      use consts_mod, only : RP, IK, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), dimension(:, :), intent(in) :: x
      real(RP), dimension(:), intent(in) :: y
      real(RP), dimension(size(x, 1)) :: z
      integer(IK) :: j 
      character(len = SRNLEN), parameter :: srname = 'MATMUL21' 

      if (DEBUG_MODE) then
          if (size(x, 2) /= size(y)) then
              call errstop(srname, 'SIZE(X, 2) /= SIZE(Y)')
          end if
      end if

      z = ZERO
      do j = 1, size(x, 2)
          z = z + x(:, j)*y(j)
      end do
      end function matmul21

      function matmul22(x, y) result(z)
      use consts_mod, only : RP, IK, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), dimension(:, :), intent(in) :: x, y
      real(RP), dimension(size(x, 1), size(y, 2)) :: z
      integer(IK) :: i, j
      character(len = SRNLEN), parameter :: srname = 'MATMUL22'

      if (DEBUG_MODE) then
          if (size(x, 2) /= size(y, 1)) then
              call errstop(srname, 'SIZE(X, 2) /= SIZE(Y, 1)')
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
      use consts_mod, only : RP, IK, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop 
      implicit none
      real(RP), dimension(:), intent(in) :: x, y
      real(RP) :: z
      integer(IK) :: i
      character(len = SRNLEN), parameter :: srname = 'INPROD'

      if (DEBUG_MODE) then
          if (size(x) /= size(y)) then
              call errstop(srname, 'SIZE(X) /= SIZE(Y)')
          end if
      end if

      !z = sum(x*y)  
      ! Using sum seems not as efficient as a loop when interfaced with 
      ! MATLAB, but this may well depend on the machine (e.g., cache 
      ! size), compiler, compiling options, and MATLAB version.  
      z = ZERO
      do i = 1, size(x)
          z = z + x(i)*y(i)
      end do
      end function inprod

#endif

      pure elemental function hypotenuse(x, y) result(r)
      ! HYPOTENUSE calculate sqrt(x^2+y^2) if __IMPROVE_POWELL_CODE__ 
      ! is 0 or if the compiler does not provide HYPOT as an intrinsic.
      ! F2008 includes the intrinsic HYPOT; MATLAB, Python, Julia, R,
      ! etc., contain built-in functions corresponding to HYPOT.

#if __IMPROVE_POWELL_CODE__ == 0
      use consts_mod, only : RP
#else
      use consts_mod, only : RP, ONE, TWO, ZERO
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

      !subroutine grot(Z, i, j, k)
      !end subroutine grot
      
      subroutine symmetrize(A)
#if __IMPROVE_POWELL_CODE__ == 0
      use consts_mod, only : RP, IK, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), intent(inout) :: A(:, :)
      integer(IK) :: j
#else
      use consts_mod, only : RP, HALF, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      implicit none
      real(RP), intent(inout) :: A(:, :)
#endif

      character(len = SRNLEN), parameter :: srname = 'SYMMETRIZE'

      if (DEBUG_MODE) then
          if (size(A, 1) /= size(A, 2)) then
              call errstop(srname, 'A is not square')
          end if
      end if

#if __IMPROVE_POWELL_CODE__ == 0
      ! A is symmetrized by setting 
      ! A(UPPER_TRIANGULAR) = A(LOWER_TRIANGULAR).
      do j = 1, size(A, 1) 
          A(1 : j - 1, j) = A(j, 1 : j - 1)  
      end do
#else
      A = A + transpose(A)
      A = A*HALF
#endif
      if (DEBUG_MODE) then
          call verisym(A, ZERO)
      end if
      end subroutine symmetrize

      subroutine verisym(A, tol)
      ! VERISYM verifies whether a matrix A is symmetric up to TOL.
      use consts_mod, only : RP, ONE, ZERO, SRNLEN
      use warnerror_mod, only : errstop
      implicit none

      real(RP), intent(in) :: A(:, :)
      real(RP), intent(in) :: tol

      character(len = SRNLEN), parameter :: srname = 'VERISYM'

      if (size(A, 1) /= size(A, 2)) then
          call errstop(srname, 'A is not square')
      end if
      if (tol > ZERO) then
          if (maxval(abs((A - transpose(A)))) >                         &
     &     tol*max(maxval(abs(A)), ONE)) then
               call errstop(srname, 'A is not symmetric up to TOL')
          end if
      else
          if (maxval(abs((A - transpose(A)))) > ZERO) then
               call errstop(srname, 'A is not symmetric up to TOL')
          end if
      end if
      end subroutine verisym

      end module lina_mod
