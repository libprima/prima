#include "ppf.h"

module linalg_mod
!--------------------------------------------------------------------------------------------------
! This module provides some basic linear algebra procedures. To improve the performance of
! these procedures, especially matprod, one can customize their implementations according to the
! resources (hardware, e.g., cache, and libraries, e.g., BLAS) available and the sizes of the
! matrices/vectors.
!
! N.B.
! When implementing the code by MATLAB, Python, ..., note the following.
! 1. We should follow the implementation with __USE_POWELL_ALGEBRA__ == 0, which uses matrix/vector
! operations instead of loops.
! 2. We should not implement the code in subroutines/functions but write it inline, because the
! code is short using matrix/vector operations, and because the overhead of subroutine/function
! calling can be high in these languages. Here we implement them as subroutines/functions in order
! to align with Powell's original code, which cannot be translated directly to matrix/vector
! operations that produce the same results in floating-point arithmetic.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Tuesday, November 02, 2021 PM11:02:16
!--------------------------------------------------------------------------------------------------

implicit none
private
! Mathematically, inprod = dot_product, matprod = matmul
public :: inprod, matprod, outprod
public :: r1update, r2update, symmetrize
public :: xpy_dot_z, xdy_plus_a, Ax_plus_y, xA_plus_y
public :: eye
public :: inv
public :: planerot
public :: calquad, errquad, hessmul
public :: errh
public :: isminor
public :: issymmetric
public :: norm
public :: sort

interface matprod
! N.B.:
! 1. When __USE_INTRINSIC_ALGEBRA__ = 0, matprod22(x, y) may differ from matmul(x, y) due to
! finite-precision arithmetic. This means that the implementation of matmul is not a naive triple
! loop. The difference has been observed on matprod22 and matprod12. The second case occurred on
! Oct. 11, 2021 in the trust-region subproblem solver of COBYLA, and it took enormous time to find
! out that Powell's code and the modernized code behaved differently due to matmul and matprod12
! when calculating RESMAX (in Powell's code) and CSTRV (in the modernized code) when stage 2 starts.
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

interface eye
    module procedure eye1, eye2
end interface eye

interface isminor
    module procedure isminor0, isminor1
end interface isminor

interface sort
    module procedure sort_i1, sort_i2
end interface sort


contains


subroutine r1_sym(A, alpha, x)
! R1_SYM sets
! A = A + ALPHA*( X*X^T ),
! where A is an NxN matrix, ALPHA is a scalar, and X is an N-dimenional vector.
use, non_intrinsic :: consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop, assert
#endif

implicit none
real(RP), intent(in) :: alpha
real(RP), intent(in) :: x(:)
real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(X))

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: n, j
#endif

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'R1_SYM'
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
call assert(issymmetric(A), 'A is symmetric', srname)
#endif
end subroutine r1_sym


subroutine r1(A, alpha, x, y)
! R1 sets
! A = A + ALPHA*( X*Y^T ),
! where A is an MxN matrix, ALPHA is a real scalar, X is an M-dimenional vector, and Y is an
! N-dimenional vector.
use, non_intrinsic :: consts_mod, only : RP

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: alpha
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP), intent(inout) :: A(:, :)  ! A(SIZE(X), SIZE(Y))

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'R1'
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
use, non_intrinsic :: consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop, assert
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
character(len=*), parameter :: srname = 'R2_SYM'
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
call assert(issymmetric(A), 'A is symmetric', srname)
#endif
end subroutine r2_sym


subroutine r2(A, alpha, x, y, beta, u, v)
! R2 sets
! A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ),
! where A is an MxN matrix, ALPHA and BETA are real scalars, X and U are M-dimenional vectors,
! Y and V are N-dimenional vectors.
use, non_intrinsic :: consts_mod, only : RP

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
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
character(len=*), parameter :: srname = 'R2'
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
use, non_intrinsic :: consts_mod, only : RP
#else
use, non_intrinsic :: consts_mod, only : RP, IK
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:, :)
real(RP) :: z(size(y, 2))

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'MATPROD12'
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
use, non_intrinsic :: consts_mod, only : RP
#else
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:, :)
real(RP), intent(in) :: y(:)
real(RP) :: z(size(x, 1))

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'MATPROD21'
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
use, non_intrinsic :: consts_mod, only : RP
#else
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:, :)
real(RP), intent(in) :: y(:, :)
real(RP) :: z(size(x, 1), size(y, 2))

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: i, j
#endif

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'MATPROD22'
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
use, non_intrinsic :: consts_mod, only : RP
#else
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
#endif

implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
real(RP) :: z

#if __USE_INTRINSIC_ALGEBRA__ == 0
integer(IK) :: i
#endif

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'INPROD'
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
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
! Input
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)
! Output
real(RP) :: z(size(x), size(y))

integer(IK) :: i
do i = 1, int(size(y), kind(i))
    z(:, i) = x * y(i)
end do
end function outprod


pure function eye1(n) result(x)
! EYE1 is the univariate case of EYE, a function similar to the MATLAB function with the same name.
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE
implicit none
! Input
integer(IK), intent(in) :: n
! Output
real(RP) :: x(max(n, 0_IK), max(n, 0_IK))
! Local variable
integer(IK) :: i
if (size(x, 1) * size(x, 2) > 0) then
    x = ZERO
    do i = 1, int(min(size(x, 1), size(x, 2)), kind(i))
        x(i, i) = ONE
    end do
end if
end function eye1


pure function eye2(m, n) result(x)
! EYE2 is the bivariate case of EYE, a function similar to the MATLAB function with the same name.
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE
implicit none
! Inputs
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n
! Output
real(RP) :: x(max(m, 0_IK), max(n, 0_IK))
! Local variable
integer(IK) :: i
if (size(x, 1) * size(x, 2) > 0) then
    x = ZERO
    do i = 1, int(min(size(x, 1), size(x, 2)), kind(i))
        x(i, i) = ONE
    end do
end if
end function eye2


function inv(A, matrix_type) result(B)
! This function calculates the inverse of a given matrix of special type. The matrix is ASSUMED
! TO BE SMALL AND INVERTIBLE. The function is implemented NAIVELY. It is NOT coded for general
! purposes but only for the usage in this project.
! Only the following matrix types (MATRIX_TYPE) are supported, which is sufficient for this project:
! ltri/LTRI: lower triangular
! utri/UTRI: upper triangular
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, DEBUGGING
use, non_intrinsic :: debug_mod, only : verisize, errstop
implicit none

! Input
real(RP), intent(in) :: A(:, :)
character(len=*), intent(in) :: matrix_type

! Output
real(RP) :: B(size(A, 1), size(A, 1))

! Local variables
integer(IK) :: i
integer(IK) :: n
character(len=*), parameter :: srname = 'INV'

n = int(size(A, 1), kind(n))
if (DEBUGGING) then
    call verisize(A, n, n)
    if (.not. (matrix_type == 'ltri' .or. matrix_type == 'LTRI' .or. matrix_type == 'utri' .or. matrix_type == 'UTRI')) then
        call errstop(srname, 'Unknown matrix type: '//matrix_type)
    end if
end if

if (matrix_type == 'ltri' .or. matrix_type == 'LTRI') then
    B = ZERO
    do i = 1, n
        B(i, i) = ONE / A(i, i)
        B(i, 1:i - 1) = -matprod(A(i, 1:i - 1) / A(i, i), B(1:i - 1, 1:i - 1))
    end do
else if (matrix_type == 'utri' .or. matrix_type == 'UTRI') then
    B = ZERO
    do i = 1, n
        B(i, i) = ONE / A(i, i)
        B(1:i - 1, i) = -matprod(B(1:i - 1, 1:i - 1), A(1:i - 1, i) / A(i, i))
    end do
end if

end function inv


function planerot(x) result(G)
! As in MATLAB, PLANEROT(X) returns a 2x2 Givens matrix G for X in R^2 so that Y = G*X has Y(2) = 0.
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : verisize
#endif

implicit none
real(RP), intent(in) :: x(:)
real(RP) :: G(2, 2)

real(RP) :: c, s, r, scaling

#if __DEBUGGING__ == 1
call verisize(x, 2_IK)
#endif

if (abs(x(2)) > ZERO) then
    r = sqrt(sum(x**2))
    if (r > ZERO .and. r <= HUGENUM) then
        c = x(1) / r
        s = x(2) / r
    else
        scaling = maxval(abs(x))
        r = sqrt(sum((x / scaling)**2))
        c = (x(1) / scaling) / r
        s = (x(2) / scaling) / r
    end if
    G = reshape([c, -s, s, c], [2, 2])
elseif (x(1) < ZERO) then
    ! Setting G = -EYE(2,2) in this case ensures the continuity of G with respect to X except at 0.
    G = -eye(2_IK)
else
    G = eye(2_IK)
end if
end function planerot


! A stabilized Givens rotation that avoids over/underflow and keeps continuity (see wikipedia)
!function [c, s, r] = givens_rotation(a, b)
!    if b == 0;
!        c = sign(a);  % Continuity
!        if (c == 0);
!            c = 1.0; % Unlike other languages, MatLab's sign function returns 0 on input 0.
!        end;
!        s = 0;
!        r = abs(a);
!    elseif a == 0;
!        c = 0;
!        s = sign(b);
!        r = abs(b);
!    elseif abs(a) > abs(b);
!        t = b / a;
!        u = sign(a) * sqrt(1 + t * t);
!        c = 1 / u;
!        s = c * t;
!        r = a * u;
!    else
!        t = a / b;
!        u = sign(b) * sqrt(1 + t * t);
!        s = 1 / u;
!        c = s * t;
!        r = b * u;
!    end
!end


subroutine symmetrize(A)
! SYMMETRIZE(A) symmetrizes A.
! Here, we assume that A is a matrix that is symmetric in precise arithmetic, and its asymmetry
! comes only from rounding errors.

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : RP, IK
#else
use, non_intrinsic :: consts_mod, only : RP, HALF
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop, assert
#endif

implicit none
real(RP), intent(inout) :: A(:, :)

#if __USE_POWELL_ALGEBRA__ == 1
integer(IK) :: j
#endif

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'SYMMETRIZE'
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
call assert(issymmetric(A), 'A is symmetrized', srname)
#endif
end subroutine symmetrize


function xpy_dot_z(x, y, z) result(t)
! t = (x + y)'*z

use, non_intrinsic :: consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK, ZERO
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
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
character(len=*), parameter :: srname = 'XPY_DOT_Z'
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

use, non_intrinsic :: consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
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
character(len=*), parameter :: srname = 'XDY_PLUS_A'
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

use, non_intrinsic :: consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
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
character(len=*), parameter :: srname = 'AX_PLUS_Y'
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

use, non_intrinsic :: consts_mod, only : RP

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
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
character(len=*), parameter :: srname = 'XA_PLUS_Y'
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


function calquad(d, gq, hq, pq, x, xpt) result(qred)
! This function evaluates QRED = Q(X) - Q(X + D) with Q being the quadratic function defined via
! (GQ, HQ, PQ) by
! Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .

use, non_intrinsic :: consts_mod, only : RP, HALF

#if __USE_POWELL_ALGEBRA__ == 1
use, non_intrinsic :: consts_mod, only : IK, ZERO
#endif

#if __DEBUGGING__ == 1
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop, verisize
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
real(RP) :: qred

! Local variable
real(RP) :: s(size(x))

#if __USE_POWELL_ALGEBRA__ == 1
real(RP) :: w(size(pq)), t
integer(IK) :: i, ih, j
#endif

#if __DEBUGGING__ == 1
integer(IK) :: n, npt
character(len=*), parameter :: srname = 'CALQUAD'
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
if (n < 1 .or. npt < n + 2) then
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
qred = ZERO
ih = 0_IK
do j = 1, int(size(d), kind(j))
    qred = qred - d(j) * gq(j)
    do i = 1, j
        ih = int(ih + 1, kind(ih))
        t = d(i) * s(j) + d(j) * x(i)
        if (i == j) then
            t = HALF * t
        end if
        qred = qred - t * hq(i, j)
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
    qred = qred - pq(i) * w(i)
end do
#else
! The order of calculation seems quite important. The following order seems to work well.
! 1st order term
qred = -inprod(d, gq)
s = HALF * d + x  ! Different from the above version.
! implicit 2nd-order term
qred = qred - sum(pq * (matprod(s, xpt) * matprod(d, xpt)))
! explicit 2nd-order term
qred = qred - inprod(s, matprod(hq, d))
! The following implementations do not work as well as the above one.
!qred = qred - inprod(d, matprod(hq, s))
!qred = qred - sum(hq * outprod(s, d))
!qred = qred - HALF*(inprod(d, matprod(hq, s)) + inprod(s, matprod(hq, d)))
#endif
end function calquad


function errquad(gq, hq, pq, xpt, fval) result(err)
!--------------------------------------------------------------------------------------------------!
! This function calculates the maximal relative error of Q in interpolating FVAL on XPT.
! Here, Q is the quadratic function defined via (GQ, HQ, PQ) by
! Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
implicit none

! Inputs
real(RP), intent(in) :: gq(:)     ! GQ(N)
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)
real(RP), intent(in) :: fval(:)   ! FVAL(NPT)

! Outputs
real(RP) :: err

! Local variables
character(len=*), parameter :: srname = 'ERRQUAD'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: fmq(size(xpt, 2))
real(RP) :: qval(size(xpt, 2))
real(RP) :: zeros(size(xpt, 1))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt, 'SIZE(FVAL) = NPT', srname)
    call assert(.not. any(is_nan(fval) .or. is_posinf(fval)), 'FVAL is not NaN or +Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

zeros = ZERO
qval = [(-calquad(xpt(:, k), gq, hq, pq, zeros, xpt), k=1, npt)]
if (.not. all(is_finite(qval))) then
    err = HUGENUM
else
    fmq = fval - qval
    err = (maxval(fmq) - minval(fmq)) !/ max(ONE, maxval(abs(fval)))
end if

!====================!
!  Calculation ends  !
!====================!
end function errquad


function hessmul(hq, pq, xpt, y) result(hy)
! This function calculates HESSIAN*Y, where HESSIAN consists of an explicit part HQ and an
! implicit part PQ in Powell's way: HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).

use, non_intrinsic :: consts_mod, only : RP

#if __DEBUGGING__ == 1
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop, verisize
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
character(len=*), parameter :: srname = 'HESSMUL'
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
if (n < 1 .or. npt < n + 2) then
    call errstop(srname, 'SIZE(XPT) is invalid')
end if
call verisize(y, n)
call verisize(hq, n, n)
call verisize(pq, npt)
#endif

hy = matprod(hq, y) + matprod(xpt, pq * matprod(y, xpt))
end function hessmul


function errh(idz, bmat, zmat, xpt) result(err)
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: bmat(:, :)
real(RP), intent(in) :: zmat(:, :)
real(RP), intent(in) :: xpt(:, :)

! Outputs
real(RP) :: err

! Local variables
character(len=*), parameter :: srname = 'ERRH'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: A(size(xpt, 2), size(xpt, 2))
real(RP) :: e(3, 3)
real(RP) :: maxabs
real(RP) :: Omega(size(xpt, 2), size(xpt, 2))
real(RP) :: U(size(xpt, 2), size(xpt, 2))
real(RP) :: V(size(xpt, 1), size(xpt, 2))
real(RP) :: r(size(xpt, 2))
real(RP) :: s(size(xpt, 1))
real(RP) :: t(size(xpt, 2))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

! N.B.: The (NPT+1)th column (and row) of H is not saved. It is [r; t(1), s]  below.
! W = [A, ONES(NPT, 1), XPT^T; ONES(1, NPT), ZERO, ZEROS(1, N); XPT, ZEROS(N, 1), ZEROS(N, N)]
! H = [Omega, r, BMAT(:, 1:NPT)^T; r^T, t(1), s^T, BMAT(:, 1:NPT), s, BMAT(:, NPT+1:NPT+N)]

A = HALF * matprod(transpose(xpt), xpt)**2
Omega = -matprod(zmat(:, 1:idz - 1), transpose(zmat(:, 1:idz - 1))) + &
    & matprod(zmat(:, idz:npt - n - 1), transpose(zmat(:, idz:npt - n - 1)))
maxabs = maxval([ONE, maxval(abs(A)), maxval(abs(Omega)), maxval(abs(bmat))])
U = eye(npt) - matprod(A, Omega) - matprod(transpose(xpt), bmat(:, 1:npt))
V = -matprod(bmat(:, 1:npt), A) - matprod(bmat(:, npt + 1:npt + n), xpt)
r = sum(U, dim=1) / real(npt, RP)
s = sum(V, dim=2) / real(npt, RP)
t = -matprod(A, r) - matprod(s, xpt)
e(1, 1) = maxval(maxval(U, dim=1) - minval(U, dim=1))
e(1, 2) = maxval(t) - minval(t)
e(1, 3) = maxval(maxval(V, dim=2) - minval(V, dim=2))
e(2, 1) = maxval(abs(sum(Omega, dim=1)))
e(2, 2) = abs(sum(r) - ONE)
e(2, 3) = maxval(abs(sum(bmat(:, 1:npt), dim=2)))
e(3, 1) = maxval(abs(matprod(xpt, Omega)))
e(3, 2) = maxval(abs(matprod(xpt, r)))
e(3, 3) = maxval(abs(matprod(xpt, transpose(bmat(:, 1:npt))) - eye(n)))
err = maxval(e) / maxabs
end function


pure function isminor0(x, ref) result(is_minor)
! This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
use, non_intrinsic :: consts_mod, only : RP, TENTH, TWO
implicit none

real(RP), intent(in) :: x
real(RP), intent(in) :: ref
logical :: is_minor
real(RP) :: refa, refb

refa = abs(ref) + TENTH * abs(x)
refb = abs(ref) + TWO * TENTH * abs(x)
is_minor = abs(ref) >= refa .or. refa >= refb
end function isminor0


function isminor1(x, ref) result(is_minor)
! This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
use, non_intrinsic :: consts_mod, only : IK, RP
#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : errstop
#endif
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(in) :: ref(:)
logical :: is_minor(size(x))
integer(IK) :: i

#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'ISMINOR1'
if (size(x) /= size(ref)) then
    call errstop(srname, 'SIZE(X) /= SIZE(REF)')
end if
#endif

is_minor = [(isminor0(x(i), ref(i)), i=1, int(size(x), IK))]
end function isminor1

pure function issymmetric(A, tol) result(is_symmetric)
! This function tests whether A is symmetric up to TOL.
use, non_intrinsic :: consts_mod, only : RP, ONE, ZERO
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_symmetric

is_symmetric = .true.
if (size(A, 1) /= size(A, 2)) then
    is_symmetric = .false.
elseif (.not. all(is_nan(A) .eqv. is_nan(transpose(A)))) then
    is_symmetric = .false.
elseif (.not. present(tol) .and. any(abs(A - transpose(A)) > ZERO)) then
    is_symmetric = .false.
elseif (present(tol)) then
    ! Do not merge the next line with the last, as Fortran may not evaluate the logical expression
    ! in the short-circuit way.
    if (any(abs(A - transpose(A)) > abs(tol) * max(maxval(abs(A)), ONE))) then
        is_symmetric = .false.
    end if
end if

end function issymmetric


pure function norm(x, p) result(y)
! This function calculates the P-norm of a vector X.
use, non_intrinsic :: consts_mod, only : RP, ONE, ZERO
use, non_intrinsic :: infnan_mod, only : is_finite, is_posinf
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in), optional :: p

! Outputs
real(RP) :: y

! Local variables
real(RP) :: scaling

if (size(x) == 0) then
    y = ZERO
elseif (.not. all(is_finite(x))) then
    ! If X contains NaN, then Y is NaN. Otherwise, Y is Inf when X contains +/-Inf.
    y = sum(abs(x))
elseif (.not. any(abs(x) > ZERO)) then
    ! The following is incorrect without checking the last case, as X may be all NaN.
    y = ZERO
else
    if (.not. present(p)) then
        y = sqrt(sum(x**2))
    else if (is_posinf(p)) then
        y = maxval(abs(x))
    else
        scaling = max(tiny(1.0_RP), sqrt(maxval(abs(x)) * minval(abs(x), mask=(abs(x) > ZERO))))
        y = scaling * sum(abs(x / scaling)**p)**(ONE / p)
    end if
end if

end function norm

function sort_i1(x, direction) result(y)
! This function sorts X according to DIRECTION, which should be 'ascend' (default) or 'descend'.
use, non_intrinsic :: consts_mod, only : IK
#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : assert
#endif
implicit none

! Inputs
integer(IK), intent(in) :: x(:)
character(len=*), optional, intent(in) :: direction

! Outputs
integer(IK) :: y(size(x))

! Local variables
#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'SORT_I1'
#endif
integer(IK) :: i
integer(IK) :: n
integer(IK) :: newn
logical :: ascending

!====================!
! Calculation starts !
!====================!

ascending = .true.
if (present(direction)) then
    if (direction == 'descend' .or. direction == 'DESCEND') then
        ascending = .false.
    end if
end if

y = x
n = int(size(y), kind(n))
do while (n > 1)  ! Bubble sort.
    newn = 0
    do i = 2, n
        if ((y(i - 1) > y(i) .and. ascending) .or. (y(i - 1) < y(i) .and. .not. ascending)) then
            y([i - 1_IK, i]) = y([i, i - 1_IK])
            newn = i
        end if
    end do
    n = newn
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
#if __DEBUGGING__ == 1
if (ascending) then
    call assert(all(y(1:n - 1) <= y(2:n)), 'Y is ascending', srname)
else
    call assert(all(y(1:n - 1) >= y(2:n)), 'Y is descending', srname)
end if
#endif
end function sort_i1

function sort_i2(x, dim, direction) result(y)
! This function sorts a matrix X according to DIM (1 or 2) and DIRECTION ('ascend' or 'descend').
#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : assert
#endif
use, non_intrinsic :: consts_mod, only : IK
implicit none

! Inputs
integer(IK), intent(in) :: x(:, :)
integer, optional, intent(in) :: dim
character(len=*), optional, intent(in) :: direction

! Outputs
integer(IK) :: y(size(x, 1), size(x, 2))

! Local variables
#if __DEBUGGING__ == 1
character(len=*), parameter :: srname = 'SORT_I2'
integer(IK) :: n
#endif
character(len=100) :: direction_loc
integer :: dim_loc
integer(IK) :: i

!====================!
! Calculation starts !
!====================!

dim_loc = 1
if (present(dim)) then
    dim_loc = dim
end if

direction_loc = 'ascend'
if (present(direction)) then
    direction_loc = direction
end if

y = x
if (dim_loc == 1) then
    do i = 1, int(size(x, 2), IK)
        y(:, i) = sort_i1(y(:, i), trim(direction_loc))
    end do
else
    do i = 1, int(size(x, 1), IK)
        y(i, :) = sort_i1(y(i, :), trim(direction_loc))
    end do
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
#if __DEBUGGING__ == 1
if (dim_loc == 1) then
    n = int(size(y, 1), kind(n))
    if (trim(direction_loc) == 'ascend' .or. trim(direction_loc) == 'ASCEND') then
        call assert(all(y(1:n - 1, :) <= y(2:n, :)), 'Y is ascending along dimension 1', srname)
    else
        call assert(all(y(1:n - 1, :) >= y(2:n, :)), 'Y is descending along dimension 1', srname)
    end if
else
    n = int(size(y, 2), kind(n))
    if (trim(direction_loc) == 'ascend' .or. trim(direction_loc) == 'ASCEND') then
        call assert(all(y(:, 1:n - 1) <= y(:, 2:n)), 'Y is ascending along dimension 2', srname)
    else
        call assert(all(y(:, 1:n - 1) >= y(:, 2:n)), 'Y is descending along dimension 2', srname)
    end if
end if
#endif
end function sort_i2


end module linalg_mod
