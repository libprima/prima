! LINA is a module providing some basic linear algebra procedures.
! To improve the performance of these procedures, especially matmul, one
! can customize their implementations, according to the resources 
! (hardware, e.g., cache, and libraries, e.g., BLAS) available and the
! sizes of the matrices/vectors. 


#include "ppf.h"

      module lina_mod

      implicit none
      private
      public :: r1update, r2update, grota
      public :: symmetrize
      public :: verisize
      public :: xpy_dot_z, xdy_plus_a, Ax_plus_y, xA_plus_y


      interface r1update
          module procedure r1_sym, r1
      end interface r1update

      interface r2update
          module procedure r2_sym, r2
      end interface r2update

      interface verisize
          module procedure verisize_real_1, verisize_real_2
          module procedure verisize_int_1, verisize_int_2
      end interface verisize

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


      contains

      subroutine r1_sym(A, alpha, x)
      ! R1_SYM sets 
      ! A = A + ALPHA*( X*X^T ), 
      ! where A is an NxN matrix, ALPHA is a scalar, and X is an
      ! N-dimenional vector.
      use consts_mod, only : RP, IK, ONE, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))
      integer(IK) :: n, j 

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname ='R1_SYM'
      ! Be careful with initialization! 
      ! In Fortran >=90, the initialization in the declaration implies
      ! the "save" attribute. If the variable is not a parameter, it may
      ! casue unwanted behavior.
      if (size(A, 1) /= size(x) .or. size(A, 2) /= size(x)) then
          call errstop(srname, 'SIZE(A) is invalid')
      end if
#endif

      n = int(size(x), kind(n))

#if __USE_POWELL_ALGEBRA__ == 1
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
      ! Set A(UPPER_TRI) = A(LOWER_TRI).
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

#if __DEBUG_MODE__ == 1
      call verisym(A, ZERO)
#endif
      end subroutine r1_sym

      subroutine r1(A, alpha, x, y)
      ! R1 sets 
      ! A = A + ALPHA*( X*Y^T ), 
      ! where A is an MxN matrix, ALPHA is a real scalar, X is an
      ! M-dimenional vector, and Y is an N-dimenional vector.
      use consts_mod, only : RP, IK, ONE, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer(IK) :: n, j 

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'R1'
      ! Be careful with initialization! 
      ! In Fortran >=90, the initialization in the declaration implies
      ! the "save" attribute. If the variable is not a parameter, it may
      ! casue unwanted behavior.
      if (size(A, 1) /= size(x) .or. size(A, 2) /= size(y)) then
          call errstop(srname, 'SIZE(A) is invalid')
      end if
#endif

      n = int(size(y), kind(n))

#if __USE_POWELL_ALGEBRA__ == 1
      ! Both of the following cases are invoked in NEWUOA.
      if (abs(alpha - ONE) > ZERO) then
          do j = 1, n
              A(:, j) = A(:, j) + alpha*x*y(j)
          end do
      else
          do j = 1, n
              A(:, j) = A(:, j) + x*y(j)
          end do
      end if
#else
      ! Here we calculate the increment as ALPHA*( X*Y^T ) rather than 
      ! ALPHA*X*Y^T to ensure that, when X = Y, the increment is 
      ! symmetric even in floating-point arithmetic. 
      ! Both of the following cases are invoked in NEWUOA.
      if (abs(alpha - ONE) > ZERO) then
          do j = 1, n
              A(:, j) = A(:, j) + alpha*( x*y(j) )
          end do
      else
          do j = 1, n
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
      use consts_mod, only : RP, IK, ONE, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))
      integer(IK) :: n, j 

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'R2_SYM'
      if (size(x) /= size(y)) then
          call errstop(srname, 'SIZE(X) /= SIZE(Y)')
      end if
      if (size(A, 1) /= size(x) .or. size(A, 2) /= size(x)) then
          call errstop(srname, 'SIZE(A) is invalid')
      end if
#endif

      n = int(size(x), kind(n))

#if __USE_POWELL_ALGEBRA__ == 1
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
      ! Set A(UPPER_TRI) = A(LOWER_TRI).
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

#if __DEBUG_MODE__ == 1
      call verisym(A, ZERO)
#endif
      end subroutine r2_sym

      subroutine r2(A, alpha, x, y, beta, u, v)
      ! R2 sets 
      ! A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ), 
      ! where A is an MxN matrix, ALPHA and BETA are real scalars, X and 
      ! U are M-dimenional vectors, Y and V are N-dimenional vectors.
      use consts_mod, only : RP, IK, ONE, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: alpha
      real(RP), intent(in) :: beta 
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(in) :: u(:)  ! U(SIZE(X))
      real(RP), intent(in) :: v(:)  ! V(SIZE(Y))
      real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))
      integer(IK) :: n, j 

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'R2'
      if (size(u) /= size(x)) then
          call errstop(srname, 'SIZE(U) /= SIZE(X)')
      end if
      if (size(v) /= size(y)) then
          call errstop(srname, 'SIZE(V) /= SIZE(Y)')
      end if
      if (size(A, 1) /= size(x) .or. size(A, 2) /= size(y)) then
          call errstop(srname, 'SIZE(A) is invalid')
      end if
#endif

      n = int(size(y), kind(n))

#if __USE_POWELL_ALGEBRA__ == 1 
      if (abs(alpha - ONE) > ZERO .or. abs(beta - ONE) > ZERO) then
          ! This case is NOT invoked in NEWUOA.
          do j = 1, n
              A(:, j) = A(:, j) + alpha*x*y(j) + beta*u*v(j)
          end do
      else  ! This case is invoked in NEWUOA.
          do j = 1, n
              A(:, j) = A(:, j) + x*y(j) + u*v(j)
          end do
      end if
#else
      if (abs(alpha - beta) > ZERO) then  ! NOT invoked in NEWUOA.
          do j = 1, n
              A(:, j) = A(:, j) + (alpha*(x*y(j)) + beta*(u*v(j)))
          end do
      else if (abs(alpha - ONE) > ZERO) then  ! NOT invoked in NEWUOA.
      ! Here update A to A + ALPHA*( X*Y^T + U*V^T ), a semmetric update
      ! when X = V and Y = U, even in floating-point arithmetic.
          do j = 1, n 
              A(:, j) = A(:, j) + alpha*( x*y(j)  + u*v(j) )
          end do
      else  ! This case is invoked in NEWUOA.
      ! Here update A to A + ( X*Y^T + U*V^T ), a semmetric update 
      ! when X = V and Y = U, even in floating-point arithmetic.
          do j = 1, n
              A(:, j) = A(:, j) + ( x*y(j) + u*v(j) )
          end do
      end if
#endif
      end subroutine r2


#if __USE_INTRINSIC_ALGEBRA__ == 0 

      function matmul12(x, y) result(z)
      use consts_mod, only : RP, IK
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:, :)
      real(RP) :: z(size(y, 2))
      integer(IK) :: j
      
#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'MATMUL12' 
      if (size(x) /= size(y, 1)) then
          call errstop(srname, 'SIZE(X) /= SIZE(Y, 1)')
      end if
#endif

      do j = 1, int(size(y, 2), kind(j))
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
      use consts_mod, only : RP, IK, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: x(:, :)
      real(RP), intent(in) :: y(:)
      real(RP) :: z(size(x, 1))
      integer(IK) :: j 

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'MATMUL21' 
      if (size(x, 2) /= size(y)) then
          call errstop(srname, 'SIZE(X, 2) /= SIZE(Y)')
      end if
#endif

      z = ZERO
      do j = 1, int(size(x, 2), kind(j))
          z = z + x(:, j)*y(j)
      end do
      end function matmul21

      function matmul22(x, y) result(z)
      use consts_mod, only : RP, IK, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif
      implicit none
      real(RP), intent(in) :: x(:, :)
      real(RP), intent(in) :: y(:, :)
      real(RP) :: z(size(x, 1), size(y, 2))
      integer(IK) :: i, j

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'MATMUL22'
      if (size(x, 2) /= size(y, 1)) then
          call errstop(srname, 'SIZE(X, 2) /= SIZE(Y, 1)')
      end if
#endif

      z = ZERO
      do j = 1, int(size(y, 2), kind(j))
          do i = 1, int(size(x, 2), kind(i))
              z(:, j) = z(:, j) + x(:, i)*y(i, j)
          end do
      end do
      end function matmul22

      function inprod(x, y) result(z)
      use consts_mod, only : RP, IK, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop 
#endif
      implicit none
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP) :: z
      integer(IK) :: i

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'INPROD'
      if (size(x) /= size(y)) then
          call errstop(srname, 'SIZE(X) /= SIZE(Y)')
      end if
#endif

      !z = sum(x*y)  
      ! Using sum seems not as efficient as a loop when interfaced with 
      ! MATLAB, but this may well depend on the machine (e.g., cache 
      ! size), compiler, compiling options, and MATLAB version.  
      z = ZERO
      do i = 1, int(size(x), kind(i))
          z = z + x(i)*y(i)
      end do
      end function inprod

#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check whether HYPOT is available as an intrinsic.
! F2008 provides the intrinsic HYPOT; MATLAB, Python, Julia, R, etc.,
! contain built-in functions corresponding to HYPOT.
#ifdef __HYPOT_AVAILABLE__
#undef __HYPOT_AVAILABLE__
#endif
#define __HYPOT_AVAILABLE__ 0 
! HYPOT is included in F2008.
#if __FORTRAN_STANDARD__ >= 2008 
#undef __HYPOT_AVAILABLE__
#define __HYPOT_AVAILABLE__ 1
#endif
!! HYPOT is included in gfortran since 4.4.7 
!#ifdef __GFORTRAN__
!#if __GNUC__ >= 5
!#undef __HYPOT_AVAILABLE__
!#define __HYPOT_AVAILABLE__ 1
!#endif
!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine grota(A, i, j, k)
      ! GROTA sets A = A*G, 
      ! where G is the Givens rotation that 
      !    [A(K, 1), ..., A(K, I), ..., A(K, J), ..., A(K, N)] * G
      !  = [A(K, 1), ...,    R   , ...,    0   , ..., A(K, N)]
      ! with R = SQRT(A(K, I)^2 + A(K, J)^2).
      ! Indeed, the [(I, I), (I, J); (J, I), (J, J)] block of G is 
      ! [C, -S; S, C], where C = A(K, I)/R, S = A(K, J)/R, and the other
      ! entries are 1 (diagonal) or 0. Therefore, A*G is identical to A
      ! except for columns I and J:
      ! new colum I =  C*A(:, I) + S*A(:, J)
      ! new colum J = -S*A(:, I) + C*A(:, J)

      use consts_mod, only : IK, RP, ZERO
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN
      use warnerror_mod, only : errstop
#endif      
      implicit none
      integer(IK), intent(in) :: i
      integer(IK), intent(in) :: j
      integer(IK), intent(in) :: k
      real(RP), intent(inout) :: A(:, :)

      real(RP) :: c, s, r, Atemp(size(A, 1))

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'GROTA'
      if (i == j .or. min(i, j) < 1 .or. max(i, j) > size(A, 2)) then
          call errstop(srname, 'I or J is invalid')
      end if
      if (k < 1 .or. k > size(A, 1)) then
          call errstop(srname, 'K is invalid')
      end if
#endif

      if (abs(A(k, j)) >  ZERO) then
          c = A(k, i)
          s = A(k, j)
#if __USE_POWELL_ALGEBRA__ == 1 || __HYPOT_AVAILABLE__ == 0 
          r = sqrt(c**2 + s**2)
#else
          r = hypot(c, s) 
#endif
          c = c/r
          s = s/r
          Atemp = A(:, i)
          A(:, i) =  c*Atemp + s*A(:, j)
          A(:, j) = -s*Atemp + c*A(:, j) 
          A(k, j) = ZERO
#if __USE_POWELL_ALGEBRA__ == 0
          A(k, i) = r  
      else if (A(k, i) < ZERO) then
          ! In this case, C = -1, S = 0. Thus it is reasonbale to set
          ! the following. Doing this will ensure the continuity of this
          ! subroutine as an operator.
              A(:, i) = -A(:, i)  
              A(:, j) = -A(:, j)  
#endif
      end if
      end subroutine grota
      
      subroutine symmetrize(A)
#if __USE_POWELL_ALGEBRA__ == 1
      use consts_mod, only : RP, IK
#else
      use consts_mod, only : RP, HALF
#endif

#if __DEBUG_MODE__ == 1
      use consts_mod, only : ZERO, SRNLEN
      use warnerror_mod, only : errstop
#endif

      implicit none
      real(RP), intent(inout) :: A(:, :)

#if __USE_POWELL_ALGEBRA__ == 1
      integer(IK) :: j
#endif

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'SYMMETRIZE'
      if (size(A, 1) /= size(A, 2)) then
          call errstop(srname, 'A is not square')
      end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
      ! A is symmetrized by setting A(UPPER_TRI) = A(LOWER_TRI).
      do j = 1, int(size(A, 1), kind(j)) 
          A(1 : j - 1, j) = A(j, 1 : j - 1)
      end do
#else
      A = A + transpose(A)
      A = A*HALF
#endif

#if __DEBUG_MODE__ == 1
      call verisym(A, ZERO)
#endif
      end subroutine symmetrize

#if __DEBUG_MODE__ == 1
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
#endif

      subroutine verisize_real_1(x, n)
      use consts_mod, only : RP, IK, SRNLEN
      use warnerror_mod, only : errstop
      real(RP), intent(in) :: x(:)
      integer(IK), intent(in) :: n

      character(len = SRNLEN), parameter :: srname = 'VERISIZE_REAL_1'

      if (size(x) /= n) then
          call errstop(srname, 'SIZE(X) /= N')
      end if
      end subroutine verisize_real_1

      subroutine verisize_real_2(x, m, n)
      use consts_mod, only : RP, IK, SRNLEN
      use warnerror_mod, only : errstop
      real(RP), intent(in) :: x(:, :)
      integer(IK), intent(in) :: m
      integer(IK), intent(in) :: n

      character(len = SRNLEN), parameter :: srname = 'VERISIZE_REAL_2'

      if (size(x, 1) /= m) then
          call errstop(srname, 'SIZE(X, 1) /= M')
      end if
      if (size(x, 2) /= n) then
          call errstop(srname, 'SIZE(X, 2) /= N')
      end if
      end subroutine verisize_real_2

      subroutine verisize_int_1(x, n)
      use consts_mod, only : IK, SRNLEN
      use warnerror_mod, only : errstop
      integer(IK), intent(in) :: x(:)
      integer(IK), intent(in) :: n

      character(len = SRNLEN), parameter :: srname = 'VERISIZE_INT_1'

      if (size(x) /= n) then
          call errstop(srname, 'SIZE(X) /= N')
      end if
      end subroutine verisize_int_1

      subroutine verisize_int_2(x, m, n)
      use consts_mod, only : IK, SRNLEN
      use warnerror_mod, only : errstop
      integer(IK), intent(in) :: x(:, :)
      integer(IK), intent(in) :: m
      integer(IK), intent(in) :: n

      character(len = SRNLEN), parameter :: srname = 'VERISIZE_INT_2'

      if (size(x, 1) /= m) then
          call errstop(srname, 'SIZE(X, 1) /= M')
      end if
      if (size(x, 2) /= n) then
          call errstop(srname, 'SIZE(X, 2) /= N')
      end if
      end subroutine verisize_int_2


      function xpy_dot_z (x, y, z) result(t)

      use consts_mod, only : RP
#if __USE_POWELL_ALGEBRA__ == 1
      use consts_mod, only : IK, ZERO
#endif
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN 
      use warnerror_mod, only : errstop
#endif
      implicit none

      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(in) :: z(:) 
      real(RP) :: t
#if __USE_POWELL_ALGEBRA__ == 1
      integer(IK) :: i
#endif

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'XPY_DOT_Z'
      if (size(x) /= size(z) .or. size(y) /= size(z)) then
          call errstop(srname, 'SIZE(X) or SIZE(Y) /= SIZE(Z)')
      end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
      t = ZERO
      do i = 1, size(z)
          t = t + x(i)*z(i) + y(i)*z(i)
      end do
#else
      t = dot_product(x + y, z)
#endif
      end function xpy_dot_z


      function xdy_plus_a(x, y, a) result(t)

      use consts_mod, only : RP
#if __USE_POWELL_ALGEBRA__ == 1
      use consts_mod, only : IK
#endif
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN 
      use warnerror_mod, only : errstop
#endif
      implicit none

      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP), intent(in) :: a 
      real(RP) :: t
#if __USE_POWELL_ALGEBRA__ == 1
      integer(IK) :: i
#endif

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'XDY_PLUS_A'
      if (size(x) /= size(y)) then
          call errstop(srname, 'SIZE(X) /= SIZE(Y)')
      end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
      t = a
      do i = 1, size(x)
          t = t + x(i)*y(i)
      end do
#else
      t = dot_product(x, y) + a
#endif
      end function xdy_plus_a


      function Ax_plus_y(A, x, y) result(z)

      use consts_mod, only : RP
#if __USE_POWELL_ALGEBRA__ == 1
      use consts_mod, only : IK
#endif
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN 
      use warnerror_mod, only : errstop
#endif
      implicit none

      real(RP), intent(in) :: A(:, :)
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP) :: z(size(y)) 
#if __USE_POWELL_ALGEBRA__ == 1
      integer(IK) :: j 
#endif

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'AX_PLUS_Y'
      if (size(x) /= size(A, 2) .or. size(y) /= size(A, 1)) then
          call errstop(srname, 'SIZE(A) /= (SIZE(Y), SIZE(X))')
      end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
      z = y
      do j = 1, size(A, 2)
          z = z + A(:, j)*x(j)
      end do
#else
      z = matmul(A, x) + y
#endif
      end function Ax_plus_y


      function xA_plus_y(A, x, y) result(z)

      use consts_mod, only : RP
#if __USE_POWELL_ALGEBRA__ == 1
      use consts_mod, only : IK
#endif
#if __DEBUG_MODE__ == 1
      use consts_mod, only : SRNLEN 
      use warnerror_mod, only : errstop
#endif
      implicit none

      real(RP), intent(in) :: A(:, :)
      real(RP), intent(in) :: x(:)
      real(RP), intent(in) :: y(:)
      real(RP) :: z(size(y)) 
#if __USE_POWELL_ALGEBRA__ == 1
      integer(IK) :: i 
#endif

#if __DEBUG_MODE__ == 1
      character(len = SRNLEN), parameter :: srname = 'XA_PLUS_Y'
      if (size(x) /= size(A, 1) .or. size(y) /= size(A, 2)) then
          call errstop(srname, 'SIZE(A) /= (SIZE(X), SIZE(Y))')
      end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
      z = y
      do i = 1, size(A, 2)
          z = z + x(i)*A(i, :)
      end do
#else
      z = matmul(x, A) + y
#endif
      end function xA_plus_y

      end module lina_mod
