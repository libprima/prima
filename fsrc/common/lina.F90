! LINA is a module providing some basic linear algebra procedures. To improve the performance of
! these procedures, especially matprod, one can customize their implementations according to the
! resources (hardware, e.g., cache, and libraries, e.g., BLAS) available and the sizes of the
! matrices/vectors.

! N.B.
! When implementing the code by MATLAB, Python, ..., note the following.
! 1. We should follow the implementation with __USE_POWELL_ALGEBRA__ == 0, which uses matrix/vector
! operations instead of loops.
! 2. We should not implement the code in subroutines/functions but write it inline, because the
! code is short using matrix/vector operations, and because the overhead of subroutine/function
! calling can be high in these languages. Here we implement them as subroutines/functions in order
! to align with Powell's original code, which cannot be translated directly to matrix/vector
! operations that produce the same results in floating-point arithmetic.

! Coded by Zaikun ZHANG in July 2020.
!
! Last Modified: Friday, July 30, 2021 AM01:33:19


#include "ppf.h"

module lina_mod

implicit none
private
! Mathematically, inprod = dot_product, matprod = matmul
public :: inprod, matprod, outprod
public :: r1update, r2update, symmetrize
public :: xpy_dot_z, xdy_plus_a, Ax_plus_y, xA_plus_y
public :: eye
public :: planerot
public :: calquad
public :: hessmul
public :: isminor

interface matprod
! N.B.:
! 1. When __USE_INTRINSIC_ALGEBRA__ = 0, matprod22(x, y) may differ from matmul(x, y) due to
! finite-precision arithmetic. This means that the implementation of matmul is not a naive triple
! loop. For the moment (2021-07-04), the difference has not been observed on matprod12 and
! matprod21. This depends on the platform.
! 2. When interfaced with MATLAB, the intrinsic matmul and dot_product seem not as efficient as the
! implementations below (mostly by loops). This may depend on the machine (e.g., cache size),
! compiler, compiling options, and MATLAB version.
    module procedure matprod12, matprod21, matprod22
end interface matprod

interface r1update
    module procedure r1_sym, r1
end interface r1update

interface r2update
    module procedure r2_sym, r2
end interface r2update


contains


subroutine r1_sym(A, alpha, x)
! R1_SYM sets
! A = A + ALPHA*( X*X^T ),
! where A is an NxN matrix, ALPHA is a scalar, and X is an N-dimenional vector.
use consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : ZERO, SRNLEN
use debug_mod, only : errstop, verisym
#endif

implicit none
real(RP), intent(in) :: alpha
real(RP), intent(in) :: x(:)
real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: n, j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'R1_SYM'
! Be careful with initialization!
! In Fortran >=90, the initialization in the declaration implies the "save" attribute.
! If the variable is not a parameter, it may casue unwanted behavior.
if (size(A, 1) /= size(x) .or. size(A, 2) /= size(x)) then
    call errstop(srname, 'SIZE(A) is invalid')
end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
n = int(size(x), kind(n))
! Only update the LOWER TRIANGULAR part of A. Both of the following cases are invoked in NEWUOA.
do j = 1, n
    A(j:n, j) = A(j:n, j) + alpha * x(j:n) * x(j)
end do
call symmetrize(A)  ! Set A(UPPER_TRI) by COPYING A(LOWER_TRI).
#else
! For some reason, A + alpha*outprod(x,x), A + (outprod(alpha*x, x) + outprod(x, alpha*x))/2,
! A + symmetrize(x, alpha*x), or A = A + sign(alpha) * outprod(sqrt(|alpha|) * x, sqrt(|alpha|) * x)
! does not work as well as the following lines for NEWUOA, where SYMMETRIZE should set A(UPPER_TRI)
! by COPYING A(LOWER_TRI) rather than set A = (A'+A)/2. This is essentially the same as Powell's
! code, although it computes A(UPPER_TRI) unnecessarily. When X is rather small or large, calculating
! OUTPROD(X, X) can be a bad idea, even though it guarantees symmetry in finite-precision arithmetic.
A = A + outprod(alpha * x, x)
call symmetrize(A)
#endif

#if __DEBUGGING__ == 1
call verisym(A, ZERO)
#endif
end subroutine r1_sym


subroutine r1(A, alpha, x, y)
! R1 sets
! A = A + ALPHA*( X*Y^T ),
! where A is an MxN matrix, ALPHA is a real scalar, X is an M-dimenional vector, and Y is an
! N-dimenional vector.
use consts_mod, only : RP

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: alpha
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'R1'
! Be careful with initialization!
! In Fortran >=90, the initialization in the declaration implies the "save" attribute. If the
! variable is not a parameter, it may casue unwanted behavior.
if (size(A, 1) /= size(x) .or. size(A, 2) /= size(y)) then
    call errstop(srname, 'SIZE(A) is invalid')
end if
#endif

A = A + outprod(alpha * x, y)
!A = A + alpha * outprod(x, y)
end subroutine r1


subroutine r2_sym(A, alpha, x, y)
! R2_SYM sets
! A = A + ALPHA*( X*Y^T + Y*X^T ),
! where A is an NxN matrix, X and Y are N-dimenional vectors, and alpha is a scalar.
use consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : ZERO, SRNLEN
use debug_mod, only : errstop, verisym
#endif

implicit none
real(RP), intent(in) :: alpha
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: n, j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'R2_SYM'
if (size(x) /= size(y)) then
    call errstop(srname, 'SIZE(X) /= SIZE(Y)')
end if
if (size(A, 1) /= size(x) .or. size(A, 2) /= size(x)) then
    call errstop(srname, 'SIZE(A) is invalid')
end if
#endif


#if __USE_POWELL_ALGEBRA__ == 1
n = int(size(x), kind(n))
do j = 1, n
    A(j:n, j) = A(j:n, j) + alpha * x(j:n) * y(j) + alpha * y(j:n) * x(j)
end do
! Set A(UPPER_TRI) by copying A(LOWER_TRI).
call symmetrize(A)
#else
! For some reason, A = A + ALPHA * (OUTPROD(X, Y) + OUTPROD(Y, X)) does not work as well as the
! following lines for NEWUOA, where SYMMETRIZE should set A(UPPER_TRI) by copying A(LOWER_TRI),
! although ALPHA*( X*Y^T + Y*X^T) is guaranteed symmetric even in floating-point arithmetic. These
! lines are essentially the same as Powell's code, although it calculates A(UPPER_TRI) unnecessarily.
A = A + outprod(alpha * x, y) + outprod(alpha * y, x)
call symmetrize(A)
#endif

#if __DEBUGGING__ == 1
call verisym(A, ZERO)
#endif
end subroutine r2_sym


subroutine r2(A, alpha, x, y, beta, u, v)
! R2 sets
! A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ),
! where A is an MxN matrix, ALPHA and BETA are real scalars, X and U are M-dimenional vectors,
! Y and V are N-dimenional vectors.
use consts_mod, only : RP

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: alpha
real(RP), intent(in) :: beta
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), intent(in) :: u(:)          ! U(SIZE(X))
real(RP), intent(in) :: v(:)          ! V(SIZE(Y))
real(RP), intent(inout) :: A(:, :)    ! A(SIZE(X), SIZE(Y))


#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'R2'
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


A = A + outprod(alpha * x, y) + outprod(beta * u, v)
!A = A + (alpha * outprod(x, y) + beta * outprod(u, v))
end subroutine r2


function matprod12(x, y) result(z)

#if __USE_INTRINSIC_ALGEBRA__ == 1
use consts_mod, only : RP
#else
use consts_mod, only : RP, IK
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:, :)
real(RP) :: z(size(y, 2))

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'MATPROD12'
if (size(x) /= size(y, 1)) then
    call errstop(srname, 'SIZE(X) /= SIZE(Y, 1)')
end if
#endif

#if __USE_INTRINSIC_ALGEBRA__ == 1
z = matmul(x, y)
#else
do j = 1, int(size(y, 2), kind(j))
    ! When interfaced with MATLAB, the following seems more efficient than a loop, which is strange
    ! since inprod itself is implemented by a loop. This may depend on the machine (e.g., cache
    ! size), compiler, compiling options, and MATLAB version.
    z(j) = inprod(x, y(:, j))
end do
#endif
end function matprod12


function matprod21(x, y) result(z)

#if __USE_INTRINSIC_ALGEBRA__ == 1
use consts_mod, only : RP
#else
use consts_mod, only : RP, IK, ZERO
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:, :)
real(RP), intent(in) :: y(:)
real(RP) :: z(size(x, 1))

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'MATPROD21'
if (size(x, 2) /= size(y)) then
    call errstop(srname, 'SIZE(X, 2) /= SIZE(Y)')
end if
#endif

#if __USE_INTRINSIC_ALGEBRA__ == 1
z = matmul(x, y)
#else
z = ZERO
do j = 1, int(size(x, 2), kind(j))
    z = z + x(:, j) * y(j)
end do
#endif
end function matprod21


function matprod22(x, y) result(z)
! N.B.: When __USE_INTRINSIC_ALGEBRA__ = 0, matprod22(x, y) may differ from matmul(x, y) due to
! finite-precision arithmetic. This means that the implementation of matmul is not a naive triple
! loop. Of course, the implementation depends on the platform.

#if __USE_INTRINSIC_ALGEBRA__ == 1
use consts_mod, only : RP
#else
use consts_mod, only : RP, IK, ZERO
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:, :)
real(RP), intent(in) :: y(:, :)
real(RP) :: z(size(x, 1), size(y, 2))

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: i, j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'MATPROD22'
if (size(x, 2) /= size(y, 1)) then
    call errstop(srname, 'SIZE(X, 2) /= SIZE(Y, 1)')
end if
#endif

#if __USE_INTRINSIC_ALGEBRA__ == 1
z = matmul(x, y)
#else
z = ZERO
do j = 1, int(size(y, 2), kind(j))
    do i = 1, int(size(x, 2), kind(i))
        z(:, j) = z(:, j) + x(:, i) * y(i, j)
    end do
end do
#endif
end function matprod22


function inprod(x, y) result(z)

#if __USE_INTRINSIC_ALGEBRA__ == 1
use consts_mod, only : RP
#else
use consts_mod, only : RP, IK, ZERO
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP) :: z

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: i
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'INPROD'
if (size(x) /= size(y)) then
    call errstop(srname, 'SIZE(X) /= SIZE(Y)')
end if
#endif

#if __USE_INTRINSIC_ALGEBRA__ == 1
z = dot_product(x, y)
#else
!z = sum(x*y)
! Using sum seems not as efficient as a loop when interfaced with MATLAB, but this may depend
! on the machine (e.g., cache size), compiler, compiling options, and MATLAB version.
z = ZERO
do i = 1, int(size(x), kind(i))
    z = z + x(i) * y(i)
end do
#endif
end function inprod

function outprod(x, y) result(z)
! OUTPROD calculates the outer product of X and Y, i.e., Z = X*Y^T, regarding both X and Y as columns
use consts_mod, only : RP, IK
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), dimension(size(x), size(y)) :: z

integer(IK) :: i
do i = 1, int(size(y), kind(i))
    z(:, i) = x * y(i)
end do
end function outprod

function eye(m, n) result(x)
! EYE is a function similar to the MATLAB function with the same name.
use consts_mod, only : RP, IK, ZERO, ONE
use memory_mod, only : safealloc
implicit none
integer(IK), intent(in) :: m
integer(IK), intent(in), optional :: n
real(RP), allocatable :: x(:, :)
integer(IK) :: i

if (present(n)) then
    call safealloc(x, max(m, 0_IK), max(n, 0_IK))
else
    call safealloc(x, max(m, 0_IK), max(m, 0_IK))
end if
if (size(x, 1) * size(x, 2) > 0) then
    x = ZERO
    do i = 1, int(min(size(x, 1), size(x, 2)), kind(i))
        x(i, i) = ONE
    end do
end if
end function eye


function planerot(x) result(G)
! This function returns a 2-by-2 orthogonal matrix G so that y = G*x has y(2) = 0 (same as in MATLAB).
use consts_mod, only : RP, ZERO

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:)
real(RP) :: G(2, 2)

real(RP) :: c, s, r

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'PLANEROT'
if (size(x) /= 2) then
    call errstop(srname, 'SIZE(X) /= 2')
end if
#endif

c = x(1)
s = x(2)
if (abs(s) > ZERO) then
    r = sqrt(c**2 + s**2)
    c = c / r
    s = s / r
    G = reshape([c, s, -s, c], [2, 2])
#if __USE_POWELL_ALGEBRA__ == 0
! In this case, it is reasonable to set G = -eye(2,2) so that G is continuous with respect to X.
elseif (c < ZERO) then
    G = -eye(2, 2)
#endif
else
    G = eye(2, 2)
end if
end function planerot

!!subroutine grota(A, i, j, k)
!!! GROTA sets A = A*G,
!!! where G is the Givens rotation that
!!!    [A(K, 1), ..., A(K, I), ..., A(K, J), ..., A(K, N)] * G
!!!  = [A(K, 1), ..., R   , ..., 0   , ..., A(K, N)]
!!! with R = SQRT(A(K, I)^2 + A(K, J)^2).
!!! Indeed, the [(I, I), (I, J); (J, I), (J, J)] block of G is [C, -S; S, C], where C = A(K, I)/R,
!!! S = A(K, J)/R, and the other entries are 1 (diagonal) or 0. Therefore, A*G is identical to A
!!! except for columns I and J:
!!! new column I =  C*A(:, I) + S*A(:, J)
!!! new column J = -S*A(:, I) + C*A(:, J)
!
!!use consts_mod, only : IK, RP, ZERO
!
!!#if __DEBUGGING__ == 1
!!use consts_mod, only : SRNLEN
!!use debug_mod, only : errstop
!!#endif
!
!!implicit none
!!integer(IK), intent(in) :: i
!!integer(IK), intent(in) :: j
!!integer(IK), intent(in) :: k
!!real(RP), intent(inout) :: A(:, :)
!
!!#if __DEBUGGING__ == 1
!!character(len=SRNLEN), parameter :: srname = 'GROTA'
!!if (i == j .or. min(i, j) < 1 .or. max(i, j) > size(A, 2)) then
!!    call errstop(srname, 'I or J is invalid')
!!end if
!!if (k < 1 .or. k > size(A, 1)) then
!!    call errstop(srname, 'K is invalid')
!!end if
!!#endif
!
!!if (abs(A(k, j)) > ZERO) then
!!    A(:, [i, j]) = matprod(A(:, [i, j]), planerot(A(k, [i, j])))
!!#if __USE_POWELL_ALGEBRA__ == 0
!!    A(k, i) = sqrt(A(k, i)**2 + A(k, j)**2)  ! A(K, I) = R in precise arithmetic
!!else if (A(k, i) < ZERO) then
!!    ! In this case, C = -1, S = 0. Thus it is reasonable to set the following. Doing this will
!!    ! ensure the continuity of this subroutine as an operator.
!!    A(:, i) = -A(:, i)
!!    A(:, j) = -A(:, j)
!!#endif
!!end if
!!end subroutine grota


subroutine symmetrize(A)
! SYMMETRIZE(A) symmetrizes A.
! Here, we assume that A is a matrix that is symmetric in precise arithmetic, and its asymmetry
! comes only from rounding errors.

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : RP, IK
#else
use consts_mod, only : RP, HALF
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : ZERO, SRNLEN
use debug_mod, only : errstop, verisym
#endif

implicit none
real(RP), intent(inout) :: A(:, :)

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'SYMMETRIZE'
if (size(A, 1) /= size(A, 2)) then
    call errstop(srname, 'A is not square')
end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
! A is symmetrized by setting A(UPPER_TRI) = A(LOWER_TRI).
do j = 1, int(size(A, 1), kind(j))
    A(1:j - 1, j) = A(j, 1:j - 1)
end do
#else
A = A + transpose(A)
A = A * HALF
#endif

#if __DEBUGGING__ == 1
call verisym(A, ZERO)
#endif
end subroutine symmetrize


function xpy_dot_z(x, y, z) result(t)
! t = (x + y)'*z

use consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK, ZERO
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), intent(in) :: z(:)
real(RP) :: t

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: i
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'XPY_DOT_Z'
if (size(x) /= size(z) .or. size(y) /= size(z)) then
    call errstop(srname, 'SIZE(X) or SIZE(Y) /= SIZE(Z)')
end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
t = ZERO
do i = 1, int(size(z), kind(i))
    t = t + x(i) * z(i) + y(i) * z(i)
end do
#else
! Very strangely, none of the following works as well as the above loop for NEWUOA. They are all
! equivalent mathematically but different numerically.
!t = inprod(x, z) + inprod(y, z)
!t = inprod(x + y, z)
t = sum(x * z + y * z)
#endif
end function xpy_dot_z


function xdy_plus_a(x, y, a) result(t)
! t = x'*y + a

use consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), intent(in) :: a
real(RP) :: t

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: i
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'XDY_PLUS_A'
if (size(x) /= size(y)) then
    call errstop(srname, 'SIZE(X) /= SIZE(Y)')
end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
t = a
do i = 1, int(size(x), kind(i))
    t = t + x(i) * y(i)
end do
#else
t = inprod(x, y) + a
#endif
end function xdy_plus_a


function Ax_plus_y(A, x, y) result(z)
! z = A*x + y (imagine x, y, and z as columns)

use consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif
implicit none

real(RP), intent(in) :: A(:, :)
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP) :: z(size(y))

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'AX_PLUS_Y'
if (size(x) /= size(A, 2) .or. size(y) /= size(A, 1)) then
    call errstop(srname, 'SIZE(A) /= (SIZE(Y), SIZE(X))')
end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
z = y
do j = 1, int(size(A, 2), kind(j))
    z = z + A(:, j) * x(j)
end do
#else
! The following line works a bit better than Powell's loop for NEWUOA.
z = matprod(A, x) + y
#endif
end function Ax_plus_y


function xA_plus_y(x, A, y) result(z)
! z = x*A + y (imagine x, y, z as rows).
! Note that Fortran does not distinguish rows and columns.
! If implemented in MATLAB, pay attention to the sizes of the vectors/matrices; if x and y are
! columns, then z = (x'*A)' + y (do not calculate A'*x + y, which involves transposing A).

use consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : SRNLEN
use debug_mod, only : errstop
#endif
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(in) :: A(:, :)
real(RP), intent(in) :: y(:)
real(RP) :: z(size(y))

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: i
#endif

#if __DEBUGGING__ == 1
character(len=SRNLEN), parameter :: srname = 'XA_PLUS_Y'
if (size(x) /= size(A, 1) .or. size(y) /= size(A, 2)) then
    call errstop(srname, 'SIZE(A) /= (SIZE(X), SIZE(Y))')
end if
#endif

#if __USE_POWELL_ALGEBRA__ == 1
z = y
do i = 1, int(size(A, 1), kind(i))
    z = z + x(i) * A(i, :)
end do
#else
! The following does not work as well as Powell's loop for NEWUOA.
z = matprod(x, A) + y
! Note that Fortran does not distinguish rows and columns.
! If implemented in MATLAB, pay attention to the sizes of the vectors/matrices; if x and y are
! columns, then z = (x'*A)' + y (do not calculate A'*x+y, which involves transposing A).
#endif
end function xA_plus_y


function calquad(d, gq, hq, pq, x, xpt) result(qdiff)
! CALQUAD calculates
! QDIFF = Q(X + D) - Q(X)
! with Q being the quadratic function defined via (GQ, HQ, PQ) by
! Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .

use consts_mod, only : RP, HALF

#if __USE_POWELL_ALGEBRA__ == 1
use consts_mod, only : IK, ZERO
#endif

#if __DEBUGGING__ == 1
use consts_mod, only : IK, SRNLEN
use debug_mod, only : errstop, verisize
#endif

implicit none

! Inputs
real(RP), intent(in) :: d(:)      ! D(N)
real(RP), intent(in) :: gq(:)     ! GQ(N)
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in) :: x(:)      ! X(N)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Output
real(RP) :: qdiff

! Local variable
real(RP) :: s(size(x))

#if __USE_POWELL_ALGEBRA__ == 1
real(RP) :: w(size(pq)), t
integer(IK) :: i, ih, j
#endif

#if __DEBUGGING__ == 1
integer(IK) :: n, npt
character(len=SRNLEN), parameter :: srname = 'CALQUAD'
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
if (n == 0 .or. npt < n + 2) then
    call errstop(srname, 'SIZE(XPT) is invalid')
end if
call verisize(d, n)
call verisize(x, n)
call verisize(gq, n)
call verisize(hq, n, n)
call verisize(pq, npt)
#endif

#if __USE_POWELL_ALGEBRA__ == 1
s = x + d
! First order term and explicit second order term
qdiff = ZERO
ih = 0
do j = 1, int(size(d), kind(j))
    qdiff = qdiff + d(j) * gq(j)
    do i = 1, j
        ih = int(ih + 1, kind(ih))
        t = d(i) * s(j) + d(j) * x(i)
        if (i == j) then
            t = HALF * t
        end if
        qdiff = qdiff + t * hq(i, j)
    end do
end do
!-----------------------------------------------------------!
! Powell's original code calculates W externally and pass it
! as an input to CALQUAD. Here W is calculated internally.
w = matprod(d, xpt) !---------------------------------------!
w = w * (HALF * w + matprod(x, xpt)) !----------------------!
!-----------------------------------------------------------!
! Implicit second order term
do i = 1, int(size(pq), kind(i))
    qdiff = qdiff + pq(i) * w(i)
end do
#else
! The order of calculation seems quite important. The following order seems to work well.
! 1st order term
qdiff = inprod(d, gq)
s = HALF * d + x  ! Different from the above version.
! implicit 2nd-order term
qdiff = qdiff + sum(pq * (matprod(s, xpt) * matprod(d, xpt)))
! explicit 2nd-order term
qdiff = qdiff + inprod(s, matprod(hq, d))
! The following implementations do not work as well as the above one.
!qdiff = qdiff + inprod(d, matprod(hq, s))
!qdiff = qdiff + sum(hq * outprod(s, d))
!qdiff = qdiff + HALF*(inprod(d, matprod(hq, s)) + inprod(s, matprod(hq, d)))
#endif
end function calquad


function hessmul(hq, pq, xpt, y) result(hy)
! This subroutine calculates HESSIAN*Y, where HESSIAN consists of an explicit part HQ and an
! implicit part PQ in Powell's way: HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).

use consts_mod, only : RP

#if __DEBUGGING__ == 1
use consts_mod, only : IK, SRNLEN
use debug_mod, only : errstop, verisize
#endif

implicit none

! Inputs
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in) :: y(:)      ! X(N)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Output
real(RP) :: hy(size(hq, 1))

#if __DEBUGGING__ == 1
integer(IK) :: n, npt
character(len=SRNLEN), parameter :: srname = 'HESSMUL'
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
if (n == 0 .or. npt < n + 2) then
    call errstop(srname, 'SIZE(XPT) is invalid')
end if
call verisize(y, n)
call verisize(hq, n, n)
call verisize(pq, npt)
#endif

hy = matprod(hq, y) + matprod(xpt, pq * matprod(y, xpt))
end function hessmul


function isminor(x, ref) result(is_minor)
! This subroutine tests whether X is minor compared to REF. It is used by Powell in, e.g., COBYLA.
use consts_mod, only : RP, TENTH, TWO
implicit none

real(RP), intent(in) :: x
real(RP), intent(in) :: ref
logical :: is_minor
real(RP) :: refa, refb

refa = abs(ref) + TENTH * abs(x)
refb = abs(ref) + TWO * TENTH * abs(x)
is_minor = abs(ref) >= refa .or. refa >= refb
end function


end module lina_mod
