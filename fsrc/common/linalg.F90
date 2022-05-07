#include "ppf.h"

module linalg_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some basic linear algebra procedures. To improve the performance of
! these procedures, especially MATPROD, one can customize their implementations according to the
! resources (hardware, e.g., cache, and libraries, e.g., BLAS) available and the sizes of the
! matrices/vectors.
!
! N.B.: When implementing the code by MATLAB, Python, ..., note the following.
! 1. Most of the procedures here are intrinsic to the languages or available in standard libraries.
! If available, they should not be implemented from scratch like we do here.
! 2. For the procedures that are not available, it may be better to code them inline instead of as
! external functions, because the code is usually short using matrix/vector operations, and because
! the overhead of function calling can be high in these languages.
! 3. In Fortran, we implement the procedures as subroutines/functions here for several reasons.
! 3.1.) Most of the procedures are not available intrinsically in Fortran.
! 3.2.) We want to start with an implementation that is verifiably faithful to Powell's original code.
! To achieve such faithfulness, it is not always possible to use the matrix/vector operations that
! are intrinsically available in Fortran, the most noticeable examples being DOT_PRODUCT (INPROD) and
! MATMUL (MATPROD). Powell implemented all matrix/vector operations by loops, which may not be the
! case for intrinsic procedures such as MATMUL and DOT_PRODUCT. Different implementations lead to
! slightly different results due to rounding, and hence the verification of faithfulness will fail.
! 3.3.) As of 20220507, with some compilers, the performance of Fortran's intrinsic matrix/vector
! procedures may not be as performant as naive loops, let alone highly optimized libraries like
! BLAS. Concentrating all the linear algebra procedures at one place like here, we can optimize them
! in a relatively easy way when necessary.
!
! TODO: To avoid stack overflows, functions that return a potentially large array (e.g., MATPROD)
! should declare the array as ALLOCATABLE rather than automatic.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Saturday, May 07, 2022 PM06:24:17
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: inprod, matprod, outprod  ! Mathematically, INPROD = DOT_PRODUCT, MATPROD = MATMUL
public :: r1update, r2update, symmetrize
public :: Ax_plus_y
public :: eye
public :: diag
public :: hypotenuse, planerot
public :: project
public :: inv, isinv
public :: solve
public :: lsqr
public :: isminor
public :: issymmetric, isorth, istriu
public :: sort
public :: trueloc, falseloc
public :: minimum, maximum
public :: norm
public :: linspace
public :: hessenberg
public :: int

interface matprod
! N.B.:
! 1. When __USE_INTRINSIC_ALGEBRA__ = 0, matprod22(x, y) may differ from matmul(x, y) in
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

interface project
    module procedure project1, project2
end interface

interface lsqr
    module procedure lsqr_Rdiag, lsqr_Rfull
end interface

interface isminor
    module procedure isminor0, isminor1
end interface isminor

interface sort
    module procedure sort_i1, sort_i2
end interface sort

interface minimum
    module procedure minimum1, minimum2
end interface minimum

interface maximum
    module procedure maximum1, maximum2
end interface maximum

interface norm
    module procedure p_norm, named_norm
end interface norm

interface linspace
    module procedure linspace_r, linspace_i
end interface linspace

interface hessenberg
    module procedure hessenberg_hhd_trid, hessenberg_full
end interface hessenberg

interface int
    module procedure logical_to_int
end interface int


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
! If the variable is not a parameter, it may cause unwanted behavior.
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
call symmetrize(A)  ! Copy A(LOWER_TRI) to A(UPPER_TRI).
#else
! For some reason, A + alpha*outprod(x,x), A + (outprod(alpha*x, x) + outprod(x, alpha*x))/2,
! A + symmetrize(x, alpha*x), or A = A + sign(alpha) * outprod(sqrt(|alpha|) * x, sqrt(|alpha|) * x)
! does not work as well as the following lines for NEWUOA, where SYMMETRIZE should copy A(LOWER_TRI)
! to A(UPPER_TRI) rather than set A = (A'+A)/2. This is essentially the same as Powell's code,
! although it computes A(UPPER_TRI) unnecessarily. When X is rather small or large, calculating
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
! variable is not a parameter, it may cause unwanted behavior.
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
! Copy A(LOWER_TRI) to A(UPPER_TRI).
call symmetrize(A)
#else
! For some reason, A = A + ALPHA * (OUTPROD(X, Y) + OUTPROD(Y, X)) does not work as well as the
! following lines for NEWUOA, where SYMMETRIZE should copy A(LOWER_TRI) to A(UPPER_TRI),
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
!--------------------------------------------------------------------------------------------------!
! OUTPROD calculates the outer product of X and Y, i.e., Z = X*Y^T, regarding X and Y as columns.
!--------------------------------------------------------------------------------------------------!
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
!--------------------------------------------------------------------------------------------------!
! EYE1 is the univariate case of EYE, a function similar to the MATLAB function with the same name.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE
implicit none
! Inputs
integer(IK), intent(in) :: n
! Outputs
real(RP) :: x(max(n, 0_IK), max(n, 0_IK))
! Local variables
integer(IK) :: i
if (size(x, 1) * size(x, 2) > 0) then
    x = ZERO
    do i = 1, int(min(size(x, 1), size(x, 2)), kind(i))
        x(i, i) = ONE
    end do
end if
end function eye1


pure function eye2(m, n) result(x)
!--------------------------------------------------------------------------------------------------!
! EYE2 is the bivariate case of EYE, a function similar to the MATLAB function with the same name.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE
implicit none
! Inputs
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n
! Outputs
real(RP) :: x(max(m, 0_IK), max(n, 0_IK))
! Local variables
integer(IK) :: i
if (size(x, 1) * size(x, 2) > 0) then
    x = ZERO
    do i = 1, int(min(size(x, 1), size(x, 2)), kind(i))
        x(i, i) = ONE
    end do
end if
end function eye2


function solve(A, b) result(x)
!--------------------------------------------------------------------------------------------------!
! This function solves the linear system A*X = B. We assume that A is a square matrix that is small
! and invertible, and B is a vector of length SIZE(A, 1). The implementation is NAIVE.
! TODO: Better to implement it into several subfunctions: triu, tril, and general square.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in) :: b(:)

! Outputs
real(RP) :: x(size(A, 2))

! Local variables
character(len=*), parameter :: srname = 'SOLVE'
integer(IK) :: P(size(A, 1))
integer(IK) :: i
integer(IK) :: n
real(RP) :: Q(size(A, 1), size(A, 1))
real(RP) :: R(size(A, 1), size(A, 2))
real(RP) :: tol

! Sizes
n = int(size(A, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(A, 1) == size(A, 2), 'A is square', srname)
    call assert(size(b) == size(A, 1), 'SIZE(B) == SIZE(A, 1)', srname)
end if

!====================!
! Calculation starts !
!====================!

if (n <= 0) then  ! Of course, N < 0 should never happen.
    return
end if

if (istril(A)) then
    x(1) = b(1) / A(1, 1)  ! Ensure N >= 1!
    ! Indeed, the last line should be merged to the following loop, but some compilers (particularly
    ! Flang 7.0.1) are buggy concerning the array sections here when I == 1.
    ! See https://github.com/flang-compiler/flang/issues/1238
    do i = 2, n
        x(i) = (b(i) - inprod(A(i, 1:i - 1), x(1:i - 1))) / A(i, i)
    end do
elseif (istriu(A)) then  ! This case is invoked in LINCOA.
    x(n) = b(n) / A(n, n)  ! Ensure N >= 1!
    ! Indeed, the last line should be merged to the following loop, but some compilers (particularly
    ! Flang 7.0.1) are buggy concerning the array sections here when I == N.
    ! See https://github.com/flang-compiler/flang/issues/1238
    do i = n - 1_IK, 1, -1
        x(i) = (b(i) - inprod(A(i, i + 1:n), x(i + 1:n))) / A(i, i)
    end do
else
    ! This is NOT the best algorithm for linear systems, but since the QR subroutine is available ...
    call qr(A, Q, R, P)
    x = matprod(b, Q)
    x(n) = x(n) / R(n, n)  ! Ensure N >= 1!
    do i = n - 1_IK, 1, -1
        x(i) = (x(i) - inprod(R(i, i + 1:n), x(i + 1:n))) / R(i, i)
    end do
    x(P) = x
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(x) == size(A, 2), 'SIZE(X) == SIZE(A, 2)', srname)
    tol = max(1.0E-8_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n + 1_IK, RP)))
    call assert(norm(matprod(A, x) - b) <= tol * maxval([ONE, norm(b), norm(x)]), 'A*X == B', srname)
end if

end function solve


function inv(A) result(B)
!--------------------------------------------------------------------------------------------------!
! This function calculates the inverse of a matrix A, which is ASSUMED TO BE SMALL AND INVERTIBLE.
! The function is implemented NAIVELY. It is NOT coded for general purposes but only for the usage
! in this project. Indeed, only the lower triangular case is used.
! TODO: extend this function to calculate the pseudo inverse of any matrix of full rank. Better to
! implement it into several subfunctions: triu with M >= N, tril with M <= N; general with M >= N,
! general with M <= N, etc.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)

! Outputs
real(RP) :: B(size(A, 1), size(A, 1))

! Local variables
character(len=*), parameter :: srname = 'INV'
integer(IK) :: P(size(A, 1))
integer(IK) :: InvP(size(A, 1))
integer(IK) :: i
integer(IK) :: n
real(RP) :: Q(size(A, 1), size(A, 1))
real(RP) :: R(size(A, 1), size(A, 1))
real(RP) :: tol

! Sizes
n = int(size(A, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(A, 1) == size(A, 2), 'A is square', srname)
end if

!====================!
! Calculation starts !
!====================!

if (n <= 0) then  ! Of course, N < 0 should never happen.
    return
end if

if (istril(A)) then
    ! This case is invoked in COBYLA.
    R = transpose(A)  ! Take transpose to work on columns.
    B = ZERO
    do i = 1, n
        B(i, i) = ONE / R(i, i)
        B(1:i - 1, i) = -matprod(B(1:i - 1, 1:i - 1), R(1:i - 1, i) / R(i, i))
    end do
    B = transpose(B)
elseif (istriu(A)) then
    B = ZERO
    do i = 1, n
        B(i, i) = ONE / A(i, i)
        B(1:i - 1, i) = -matprod(B(1:i - 1, 1:i - 1), A(1:i - 1, i) / A(i, i))
    end do
else
    ! This is NOT the best algorithm for the inverse, but since the QR subroutine is available ...
    call qr(A, Q, R, P)
    R = transpose(R)  ! Take transpose to work on columns.
    do i = n, 1, -1
        B(:, i) = (Q(:, i) - matprod(B(:, i + 1:n), R(i + 1:n, i))) / R(i, i)
    end do
    InvP(P) = linspace(1_IK, n, n)  ! The inverse permutation
    B = transpose(B(:, InvP))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(B, 1) == n .and. size(B, 2) == n, 'SIZE(B) == [N, N]', srname)
    call assert(istril(B) .or. .not. istril(A), 'If A is lower triangular, then so is B', srname)
    call assert(istriu(B) .or. .not. istriu(A), 'If A is upper triangular, then so is B', srname)
    tol = max(1.0E-8_RP, min(1.0E-1_RP, 1.0E10_RP * EPS * real(n + 1_IK, RP)))
    call assert(isinv(A, B, tol), 'B = A^{-1}', srname)
end if

end function inv


function isinv(A, B, tol) result(is_inv)
use, non_intrinsic :: consts_mod, only : RP, IK, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in) :: B(:, :)
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_inv

! Local variables
character(len=*), parameter :: srname = 'ISINV'
real(RP) :: tol_loc
integer(IK) :: n

! Sizes
n = int(size(A, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(A, 1) == size(A, 2), 'A is suqare', srname)
    call assert(size(B, 1) == size(B, 2), 'B is suqare', srname)
    call assert(size(A, 1) == size(B, 1), 'SIZE(A) == SIZE(B)', srname)
end if

if (present(tol)) then
    tol_loc = tol
else
    tol_loc = min(1.0E-3_RP, 1.0E2_RP * EPS * real(max(size(A, 1), size(A, 2)), RP))
end if
tol_loc = maxval([tol_loc, tol_loc * maxval(abs(A)), tol_loc * maxval(abs(B))])

is_inv = all(abs(matprod(A, B) - eye(n)) <= tol_loc) .or. all(abs(matprod(B, A) - eye(n)) <= tol_loc)
end function isinv


subroutine qr(A, Q, R, P)
!--------------------------------------------------------------------------------------------------!
! This subroutine calculates the QR factorization of A, possibly with column pivoting, so that
! A = Q*R (if no pivoting) or A(:, P) = Q*R (if pivoting), where the columns of Q are orthonormal,
! and R is upper triangular.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)

! Outputs
real(RP), intent(out), optional :: Q(:, :)
real(RP), intent(out), optional :: R(:, :)
integer(IK), intent(out), optional :: P(:)

! Local variables
character(len=*), parameter :: srname = 'QR'
logical :: pivote
integer(IK) :: i
integer(IK) :: j
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
real(RP) :: G(2, 2)
real(RP) :: Q_loc(size(A, 1), size(A, 1))
real(RP) :: T(size(A, 2), size(A, 1))
real(RP) :: tol

if (.not. (present(Q) .or. present(R) .or. present(R))) then
    return
end if

! Sizes
m = int(size(A, 1), kind(m))
n = int(size(A, 2), kind(n))

! Preconditions
if (DEBUGGING) then
    if (present(Q)) then
        call assert(size(Q, 1) == m .and. (size(Q, 2) == m .or. size(Q, 2) == min(m, n)), &
            & 'SIZE(Q) == [M, N] .or. SIZE(Q) == [M, MIN(M, N)]', srname)
    end if
    if (present(R)) then
        call assert((size(R, 1) == m .or. size(R, 1) == min(m, n)) .and. size(R, 2) == n, &
            & 'SIZE(R) == [M, N] .or. SIZE(R) == [MIN(M, N), N]', srname)
    end if
    if (present(Q) .and. present(R)) then
        call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    end if
    if (present(P)) then
        call assert(size(P) == n, 'SIZE(P) == N', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

pivote = (present(P))
Q_loc = eye(m)
T = transpose(A)  ! T is the transpose of R. We consider T in order to work on columns.
if (pivote) then
    P = linspace(1_IK, n, n)
end if

do j = 1, n
    if (pivote) then
        k = int(maxloc(sum(T(j:n, j:m)**2, dim=2), dim=1), kind(k))
        if (k > 1 .and. k <= n - j + 1) then
            k = k + j - 1_IK
            P([j, k]) = P([k, j])
            T([j, k], :) = T([k, j], :)
        end if
    end if
    do i = m, j + 1_IK, -1_IK
        G = transpose(planerot(T(j, [j, i])))
        T(j, [j, i]) = [hypotenuse(T(j, j), T(j, i)), ZERO]
        !T(j, [j, i]) = [sqrt(T(j, j)**2 + T(j, i)**2), ZERO]
        T(j + 1:n, [j, i]) = matprod(T(j + 1:n, [j, i]), G)
        Q_loc(:, [j, i]) = matprod(Q_loc(:, [j, i]), G)
    end do
end do

if (present(Q)) then
    Q = Q_loc(:, 1:size(Q, 2))
end if
if (present(R)) then
    R = transpose(T(:, 1:size(R, 1)))
end if

! Postconditions
if (DEBUGGING) then
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E4_RP * EPS * real(max(m, n) + 1_IK, RP)))
    call assert(isorth(Q_loc, tol), 'The columns of Q are orthonormal', srname)
    call assert(istril(T, tol), 'R is upper triangular', srname)
    if (pivote) then
        call assert(all(abs(matprod(Q_loc, transpose(T)) - A(:, P)) <= &
                        max(tol, tol * maxval(abs(A)))), 'A(:, P) == Q*R', srname)
        do j = 1, min(m, n) - 1_IK
            ! The following test cannot be passed on ill-conditioned problems.
            !call assert(abs(T(j, j)) + max(tol, tol * abs(T(j, j))) >= &
            !    & abs(T(j + 1, j + 1)), '|R(J, J)| >= |R(J + 1, J + 1)|', srname)
            call assert(all(T(j, j)**2 + max(tol, tol * T(j, j)**2) >= &
                & sum(T(j + 1:n, j:min(m, n))**2, dim=2)), &
                & 'R(J, J)^2 >= SUM(R(J : MIN(M, N), J + 1 : N).^2', srname)
        end do
    else
        call assert(all(abs(matprod(Q_loc, transpose(T)) - A) <= max(tol, tol * maxval(abs(A)))), &
            & 'A == Q*R', srname)
    end if
end if
end subroutine qr


function lsqr_Rdiag(A, b, Q, Rdiag) result(x)
!--------------------------------------------------------------------------------------------------!
! This function solves the linear least squares problem min ||A*x - b||_2 by the QR factorization.
! This function is used in COBYLA, where,
! 1. Q is supplied externally (called Z);
! 2. RDIAG (the diagonal of R) is supplied externally (called ZDOTA);
! 3. A HAS FULL COLUMN RANK;
! 4. It seems that b (CGRAD and DNEW) is in the column space of A (not sure yet).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)  ! A(M, N)
real(RP), intent(in) :: b(:)  ! B(M)
real(RP), intent(in), optional :: Q(:, :)  ! Q(M, :), SIZE(Q, 2) = M or MIN(M, N)
real(RP), intent(in), optional :: Rdiag(:) ! Rdiag(MIN(M, N))

! Outputs
real(RP) :: x(size(A, 2))

! Local variables
character(len=*), parameter :: srname = 'LSQR_RDIAG'
logical :: pivote
integer(IK) :: i
integer(IK) :: j
integer(IK) :: m
integer(IK) :: n
integer(IK) :: P(size(A, 2))
integer(IK) :: rank
real(RP) :: Q_loc(size(A, 1), min(size(A, 1), size(A, 2)))
real(RP) :: Rdiag_loc(min(size(A, 1), size(A, 2)))
real(RP) :: tol
real(RP) :: y(size(b))
real(RP) :: yq
real(RP) :: yqa

! Sizes
m = int(size(A, 1), kind(m))
n = int(size(A, 2), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(b) == m, 'SIZE(B) == M', srname)
    if (present(Q)) then
        call assert(size(Q, 1) == m .and. (size(Q, 2) == m .or. size(Q, 2) == min(m, n)), &
            & 'SIZE(Q) == [M, N] .or. SIZE(Q) == [M, MIN(M, N)]', srname)
        tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E6_RP * EPS * real(max(m, n) + 1_IK, RP)))
        call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    end if
    if (present(Rdiag)) then
        call assert(size(Rdiag) == min(m, n), 'SIZE(R) == MIN(M, N)', srname)
        call assert(present(Q), 'Rdiag is present only if Q is present', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

if (n <= 0) then  ! Of course, N < 0 should never happen.
    return
end if

if (present(Q)) then
    Q_loc = Q(:, 1:size(Q_loc, 2))
    if (present(Rdiag)) then
        Rdiag_loc = Rdiag
    else
        Rdiag_loc = [(inprod(Q_loc(:, i), A(:, i)), i=1, min(m, n))]
        !!MATLAB: Rdiag_loc = sum(Q_loc(:, 1:min(m,n)) .* A(:, 1:min(m,n)), 1);  % Row vector
    end if
    rank = min(m, n)
    pivote = .false.
end if

!if (.not. present(Q) .or. any(abs(Rdiag_loc) <= ZERO)) then  ! This is more reasonable
if (.not. present(Q)) then
    call qr(A, Q=Q_loc, P=P)
    Rdiag_loc = [(inprod(Q_loc(:, i), A(:, P(i))), i=1, min(m, n))]
    !!MATLAB: Rdiag_loc = sum(Q_loc(:, 1:min(m,n)) .* A(:, P(1:min(m,n))), 1);  % Row vector
    rank = maxval([0_IK, trueloc(abs(Rdiag_loc) > 0)])
    pivote = .true.
end if

x = ZERO
y = b  ! Local copy of B; B is INTENT(IN) and should not be modified.

do i = rank, 1, -1
    if (pivote) then
        j = P(i)
    else
        j = i
    end if
    ! The following IF comes from Powell. It forces X(J) = 0 if deviations from this value can be
    ! attributed to computer rounding errors. This is a favorable choice in the context of COBYLA.
    yq = inprod(y, Q_loc(:, i))
    yqa = inprod(abs(y), abs(Q_loc(:, i)))
    if (isminor(yq, yqa)) then
        x(j) = ZERO
    else
        x(j) = yq / Rdiag_loc(i)
        y = y - x(j) * A(:, j)
    end if
end do

!====================!
!  Calculation ends  !
!====================!

!! Postconditions
!if (DEBUGGING) then
!    ! The following test cannot be passed.
!    !call assert(norm(matprod(b - matprod(A, x), A)) <= max(tol, tol * norm(matprod(b, A))), &
!    !    & 'A*X is the projection of B to the column space of A', srname)
!end if
end function lsqr_Rdiag

function lsqr_Rfull(b, Q, R) result(x)
!--------------------------------------------------------------------------------------------------!
! This function solves the linear least squares problem min ||A*x - b||_2 by the QR factorization.
! This function is used in LINCOA, where,
! 1. The economy-size QR factorization is supplied externally (Q is called QFAC and R is called RFAC);
! 3. R is non-singular.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: b(:)  ! B(M)
real(RP), intent(in) :: Q(:, :)  ! Q(M, N)
real(RP), intent(in) :: R(:, :)  ! R(N, N)

! Outputs
real(RP) :: x(size(R, 2))

! Local variables
character(len=*), parameter :: srname = 'LSQR_RFULL'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol

! Sizes
m = int(size(Q, 1), kind(m))
n = int(size(R, 2), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= n .and. n >= 0, 'M >= N >= 0', srname)
    call assert(size(b) == m, 'SIZE(B) == M', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) == n, 'SIZE(Q) == [M, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E6_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(size(R, 1) == n .and. size(R, 2) == n, 'SIZE(R) == [N, N]', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
end if

!====================!
! Calculation starts !
!====================!

if (n <= 0) then  ! Of course, N < 0 should never happen.
    return
end if

x = matprod(b, Q)
do i = n, 1, -1
    do j = i + 1_IK, n
        x(i) = x(i) - R(i, j) * x(j)
    end do
    x(i) = x(i) / R(i, i)
end do
!--------------------------------------------------------------------------------------------------!
! The following is equivalent to the above, yet the above version works slightly better in LINCOA.
!!do i = n, 1_IK, -1_IK
!!    x(i) = (inprod(Q(:, i), b) - inprod(R(i, i + 1:n), x(i + 1:n))) / R(i, i)
!!end do
!--------------------------------------------------------------------------------------------------!

!====================!
!  Calculation ends  !
!====================!

end function lsqr_Rfull


function diag(A, k) result(D)
!--------------------------------------------------------------------------------------------------!
! This function takes the K-th diagonal of the matrix A, K = 0 (default) corresponding to the main
! diagonal, K > 0 above the main diagonal, and K < 0 below the main diagonal. When |K| exceeds the
! number of rows or columns in A, then the function returns an empty rank-1 array.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: memory_mod, only : safealloc
implicit none
! Inputs
real(RP), intent(in) :: A(:, :)
integer(IK), intent(in), optional :: k
! Outputs
real(RP), allocatable :: D(:)
! Local variables
integer(IK) :: dlen
integer(IK) :: i
integer(IK) :: k_loc

if (present(k)) then
    k_loc = k
else
    k_loc = 0_IK
end if

dlen = int(min(size(A, 1), size(A, 2)), IK)
dlen = max(0_IK, dlen - abs(k_loc))  ! We allow |K| to exceed the number of rows/columns in A.
call safealloc(D, dlen)
if (k_loc >= 0) then
    D = [(A(i, i + k_loc), i=1, dlen)]
else
    D = [(A(i - k_loc, i), i=1, dlen)]
end if
end function diag


function isbanded(A, lwidth, uwidth, tol) result(is_banded)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the matrix A banded within the bandwidth specified by LWIDTH and
! UWIDTH up to the tolerance TOL.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none
! Inputs
real(RP), intent(in) :: A(:, :)
integer(IK), intent(in) :: lwidth
integer(IK), intent(in) :: uwidth
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_banded

! Local variables
character(len=*), parameter :: srname = 'ISBANDED'
integer(IK) :: i
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol_loc

if (DEBUGGING) then
    call assert(lwidth >= 0 .and. uwidth >= 0, 'LWIDTH >= 0 .and. UWIDTH >= 0', srname)
    if (present(tol)) then
        call assert(tol >= 0, 'TOL >= 0', srname)
    end if
end if

tol_loc = ZERO
if (present(tol)) then
    tol_loc = max(tol, tol * maxval(abs(A)))
end if
if (is_nan(tol_loc)) then
    tol_loc = ZERO
end if

m = int(size(A, 1), kind(m))
n = int(size(A, 2), kind(n))

is_banded = .true.
do i = 1, n
    is_banded = (all(abs(A(i + lwidth + 1:m, i)) <= tol_loc) .and. all(abs(A(1:i - uwidth - 1, i)) <= tol_loc))
    if (.not. is_banded) then
        exit
    end if
end do
end function isbanded


function istril(A, tol) result(is_tril)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the matrix A is lower triangular up to the tolerance TOL.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
implicit none
! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_tril

! Local variables
integer(IK) :: width
real(RP) :: tol_loc

if (present(tol)) then
    tol_loc = tol
else
    tol_loc = ZERO
end if
width = int(max(0, size(A, 1) - 1), kind(width))
is_tril = isbanded(A, width, 0_IK, tol_loc)
end function istril


function istriu(A, tol) result(is_triu)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the matrix A is upper triangular up to the tolerance TOL.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
implicit none
! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_triu

! Local variables
integer(IK) :: width
real(RP) :: tol_loc

if (present(tol)) then
    tol_loc = tol
else
    tol_loc = ZERO
end if
width = int(max(0, size(A, 2) - 1), kind(width))
is_triu = isbanded(A, 0_IK, width, tol_loc)
end function istriu


function isorth(A, tol) result(is_orth)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the matrix A has orthonormal columns up to the tolerance TOL.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none
! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_orth

! Local variables
integer(IK) :: n

n = int(size(A, 2), kind(n))

if (n > size(A, 1)) then
    is_orth = .false.
else if (is_nan(sum(abs(A)))) then
    is_orth = .false.
else
    if (present(tol)) then
        is_orth = all(abs(matprod(transpose(A), A) - eye(n)) <= max(tol, tol * maxval(abs(A))))
    else
        is_orth = all(abs(matprod(transpose(A), A) - eye(n)) <= 0)
    end if
end if
end function isorth


function project1(x, v) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns the projection of X to SPAN(V).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_inf, is_finite, is_nan
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: v(:)

! Outputs
real(RP) :: y(size(x))

! Local variables
character(len=*), parameter :: srname = 'PROJECT1'
real(RP) :: u(size(v))
real(RP) :: scaling
real(RP) :: tol

! Preconditions
if (DEBUGGING) then
    call assert(size(x) == size(v), 'SIZE(X) == SIZE(V)', srname)
end if

!====================!
! Calculation starts !
!====================!

if (all(abs(x) <= ZERO) .or. all(abs(v) <= ZERO)) then
    y = ZERO
elseif (any(is_nan(x)) .or. any(is_nan(v))) then
    y = sum(x) + sum(v)  ! Set Y to NaN
elseif (any(is_inf(v))) then
    u = ZERO
    u(trueloc(is_inf(v))) = sign(ONE, v)
    !!MATLAB: u = 0; u(isinf(v)) = sign(v)
    u = u / norm(u)
    !y = inprod(x, u) * u
    scaling = maxval(abs(x))  ! The scaling seems to reduce the rounding error.
    y = scaling * inprod(x / scaling, u) * u
else
    u = v / norm(v)
    !y = inprod(x, u) * u
    scaling = maxval(abs(x))  ! The scaling seems to reduce the rounding error.
    y = scaling * inprod(x / scaling, u) * u
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    if (is_finite(norm(x)) .and. is_finite(norm(v))) then
        tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E6_RP * EPS))
        call assert(norm(y) <= (ONE + tol) * norm(x), 'NORM(Y) <= NORM(X)', srname)
        call assert(norm(x - y) <= (ONE + tol) * norm(x), 'NORM(X - Y) <= NORM(X)', srname)
        ! The following test may not be passed.
        call assert(abs(inprod(x - y, v)) <= max(tol, tol * max(norm(x - y) * norm(v), abs(inprod(x, v)))), &
           & 'X - Y is orthogonal to V', srname)
    end if
end if
end function project1


function project2(x, V) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns the projection of X to RANGE(V).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_inf, is_finite, is_nan
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: V(:, :)

! Outputs
real(RP) :: y(size(x))

! Local variables
character(len=*), parameter :: srname = 'PROJECT2'
real(RP) :: U(size(V, 1), min(size(V, 1), size(V, 2)))
real(RP) :: V_loc(size(V, 1), size(V, 2))
real(RP) :: tol

! Preconditions
if (DEBUGGING) then
    call assert(size(x) == size(V, 1), 'SIZE(X) == SIZE(V, 1)', srname)
end if

!====================!
! Calculation starts !
!====================!

if (size(V, 2) == 1) then
    y = project1(x, V(:, 1))
elseif (all(abs(x) <= ZERO) .or. all(abs(V) <= ZERO)) then
    y = ZERO
elseif (any(is_nan(x)) .or. any(is_nan(V))) then
    y = sum(x) + sum(V)  ! Set Y to NaN
elseif (any(is_inf(V))) then
    where (is_inf(V))
        V_loc = sign(ONE, V)
    elsewhere
        V_loc = ZERO
    end where
    !!MATLAB: V_loc = 0; V_loc(isinf(V)) = sign(V);
    call qr(V_loc, Q=U)
    y = matprod(U, matprod(x, U))
else
    call qr(V, Q=U)
    y = matprod(U, matprod(x, U))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    if (is_finite(norm(x)) .and. is_finite(sum(V**2))) then
        tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E6_RP * EPS))
        call assert(norm(y) <= (ONE + tol) * norm(x), 'NORM(Y) <= NORM(X)', srname)
        call assert(norm(x - y) <= (ONE + tol) * norm(x), 'NORM(X - Y) <= NORM(X)', srname)
        ! The following test may not be passed.
        call assert(norm(matprod(x - y, V)) <= max(tol, tol * max(norm(x - y) * sqrt(sum(V**2)), &
            & norm(matprod(x, V)))), 'X - Y is orthogonal to V', srname)
    end if
end if
end function project2


function hypotenuse(x1, x2) result(r)
! HYPOTENUSE(X1, X2) returns SQRT(X1^2 + X2^2), handling over/underflow.
use, non_intrinsic :: consts_mod, only : RP, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
implicit none

! Inputs
real(RP), intent(in) :: x1
real(RP), intent(in) :: x2

! Outputs
real(RP) :: r

! Local variables
character(len=*), parameter :: srname = 'HYPOTENUSE'
real(RP) :: y(2)

!====================!
! Calculation starts !
!====================!

if (.not. is_finite(x1)) then
    r = abs(x1)
elseif (.not. is_finite(x2)) then
    r = abs(x2)
else
    y = abs([x1, x2])
    y = [minval(y), maxval(y)]
    !if (y(1) > sqrt(REALMIN) .and. y(2) < sqrt(HUGENUM / 2.1_RP)) then
    !    r = sqrt(sum(y**2))
    !elseif (y(2) > 0) then
    !    r = y(2) * sqrt((y(1) / y(2))**2 + ONE)
    !else
    !    r = ZERO
    !end if
    ! Scaling seems to improve the precision in general.
    if (y(2) > 0) then
        r = y(2) * sqrt((y(1) / y(2))**2 + ONE)
    else
        r = ZERO
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    if (is_nan(x1) .or. is_nan(x2)) then
        call assert(is_nan(r), 'R is NaN if X1 or X2 is NaN', srname)
    else
        call assert(r >= abs(x1) .and. r >= abs(x2) .and. r <= abs(x1) + abs(x2), &
            & 'MAX{ABS(X1), ABS(X2)} <= R <= ABS(X1) + ABS(X2)', srname)
    end if
end if

end function hypotenuse


function planerot(x) result(G)
!--------------------------------------------------------------------------------------------------!
! As in MATLAB, PLANEROT(X) returns a 2x2 Givens matrix G for X in R^2 so that Y = G*X has Y(2) = 0.
! Roughly speaking, using a MATLAB-style formulation of matrices,
! G = [X(1)/R, X(2)/R; -X(2)/R, X(1)/R] with R = SQRT(X(1)^2+X(2)^2), and G*X = [R; 0].
! 0. We need to take care of the possibilities of R = 0, Inf, NaN, and over/underflow.
! 1. The G defined above is continuous with respect to X except at 0. Following this definition,
! G = [sign(X(1)), 0; 0, sign(X(1))] if X(2) = 0, G = [0, sign(X(2)); -sign(X(2)), 0] if X(2) = 0.
! Yet some implementations ignore the signs, leading to discontinuity and numerical instability.
! 2. Difference from MATLAB: if X contains NaN or consists of only Inf, MATLAB returns a NaN matrix,
! but we return an identity matrix or a matrix of +/-SQRT(2). We intend to keep G always orthogonal.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO, ONE, REALMIN, EPS, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_inf
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP) :: G(2, 2)

! Local variables
character(len=*), parameter :: srname = 'PLANEROT'
real(RP) :: c
real(RP) :: s
real(RP) :: r
real(RP) :: t
real(RP) :: u
real(RP) :: tol

! Preconditions
if (DEBUGGING) then
    call assert(size(x) == 2, 'SIZE(X) == 2', srname)
end if

!====================!
! Calculation starts !
!====================!

! Define C = X(1) / R and S = X(2) / R with R = HYPOT(X(1), X(2)). Handle Inf/NaN, over/underflow.
if (any(is_nan(x))) then
    ! In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
    c = ONE
    s = ZERO
elseif (all(is_inf(x))) then
    ! In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
    c = sign(1 / sqrt(2.0_RP), x(1))
    s = sign(1 / sqrt(2.0_RP), x(2))
elseif (abs(x(1)) <= 0 .and. abs(x(2)) <= 0) then  ! X(1) == 0 == X(2).
    c = ONE
    s = ZERO
elseif (abs(x(2)) <= EPS * abs(x(1))) then
    ! N.B.:
    ! 0. With <= instead of <, this case covers X(1) == 0 == X(2), which is treated above separately
    ! to avoid the confusing SIGN(., 0) (see 1).
    ! 1. SIGN(A, 0) = ABS(A) in Fortran but sign(0) = 0 in MATLAB, Python, Julia, and R!!!
    ! 2. Taking SIGN(X(1)) into account ensures the continuity of G with respect to X except at 0.
    c = sign(ONE, x(1))  !!MATLAB: c = sign(x(1))
    s = ZERO
elseif (abs(x(1)) <= EPS * abs(x(2))) then
    ! N.B.: SIGN(A, X) = ABS(A) * sign of X /= A * sign of X !!! Therefore, it is WRONG to define G
    ! as SIGN(RESHAPE([ZERO, -ONE, ONE, ZERO], [2, 2]), X(2)). This mistake was committed on
    ! 20211206 and took a whole day to debug! NEVER use SIGN on arrays unless you are really sure.
    c = ZERO
    s = sign(ONE, x(2))  !!MATLAB: s = sign(x(2))
else
    ! Here is the normal case. It implements the Givens rotation in a stable & continuous way as in:
    ! Bindel, D., Demmel, J., Kahan, W., and Marques, O. (2002). On computing Givens rotations
    ! reliably and efficiently. ACM Transactions on Mathematical Software (TOMS), 28(2), 206-238.
    ! N.B.: 1. Modern compilers compute SQRT(REALMIN) and SQRT(HUGENUM/2.1) at compilation time.
    ! 2. The direct calculation without involving T and U seems to work better; use it if possible.
    if (all(abs(x) > sqrt(REALMIN) .and. abs(x) < sqrt(HUGENUM / 2.1_RP))) then
        ! Do NOT use HYPOTENUSE here; the best implementation for one may be suboptimal for the other
        r = sqrt(sum(x**2))
        c = x(1) / r
        s = x(2) / r
    elseif (abs(x(1)) > abs(x(2))) then
        t = x(2) / x(1)
        u = sign(sqrt(ONE + t**2), x(1))  !!MATLAB: u = sign(x(1))*sqrt(ONE + t**2)
        c = ONE / u
        s = t / u
    else
        t = x(1) / x(2)
        u = sign(sqrt(ONE + t**2), x(2))  !!MATLAB: u = sign(x(2))*sqrt(ONE + t**2)
        c = t / u
        s = ONE / u
    end if
end if

G = reshape([c, -s, s, c], [2, 2])  !!MATLAB: G = [c, s; -s, c]

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(G, 1) == 2 .and. size(G, 2) == 2, 'SIZE(G) == [2, 2]', srname)
    call assert(all(is_finite(G)), 'G is finite', srname)
    call assert(abs(G(1, 1) - G(2, 2)) + abs(G(1, 2) + G(2, 1)) <= 0, &
        & 'G(1,1) == G(2,2), G(1,2) = -G(2,1)', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E6_RP * EPS))
    call assert(isorth(G, tol), 'G is orthonormal', srname)
    if (all(is_finite(x) .and. abs(x) < sqrt(HUGENUM / 2.1_RP))) then
        r = sqrt(sum(x**2))
        call assert(norm(matprod(G, x) - [r, ZERO]) <= max(tol, tol * r), 'G*x = [|x|, 0]', srname)
    end if
end if

end function planerot


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
! A is symmetrized by copying A(LOWER_TRI) to A(UPPER_TRI).
! N.B.: The following assumes that A(LOWER_TRI) has been properly defined.
do j = 1, int(size(A, 1), kind(j))
    A(1:j - 1, j) = A(j, 1:j - 1)
end do
#else
! N.B.: The following assumes that both A(UPPER_TRI) and A(LOWER_TRI) have been properly defined.
A = HALF * (A + transpose(A))
#endif

#if __DEBUGGING__ == 1
call assert(issymmetric(A), 'A is symmetrized', srname)
#endif
end subroutine symmetrize


function Ax_plus_y(A, x, y) result(z)
!--------------------------------------------------------------------------------------------------!
! z = A*x + y (imagine x, y, and z as columns)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)

! Outputs
real(RP) :: z(size(y))

! Local variables
character(len=*), parameter :: srname = 'AX_PLUS_Y'
integer(IK) :: j

! Preconditions
if (DEBUGGING) then
    call assert(size(x) == size(A, 2) .and. size(y) == size(A, 1), 'SIZE(A) == [SIZE(Y), SIZE(X)]', &
        & srname)
end if

!====================!
! Calculation starts !
!====================!

!--------------------------------------------------------------------------------------------------!
! In BIGLAG of NEWUOA, the following loop works numerically better than Z = MATPROD(A, X) + Y. Why?
!--------------------------------------------------------------------------------------------------!
z = y
do j = 1, int(size(A, 2), kind(j))
    z = z + A(:, j) * x(j)
end do

!====================!
!  Calculation ends  !
!====================!

end function Ax_plus_y


pure function isminor0(x, ref) result(is_minor)
! This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
! In precise arithmetic, ISMINOR(X, REF) is TRUE if and only if X == 0; in floating-point
! arithmetic, ISMINOR(X, REF) is true if X is zero or its nonzero value can be attributed to
! computer rounding errors according to REF.
! Larger SENSITIVITY means the function is more strict/precise.
use, non_intrinsic :: consts_mod, only : RP, TENTH, TWO
implicit none

real(RP), intent(in) :: x
real(RP), intent(in) :: ref
logical :: is_minor
real(RP), parameter :: sensitivity = TENTH
real(RP) :: refa
real(RP) :: refb

refa = abs(ref) + sensitivity * abs(x)
refb = abs(ref) + TWO * sensitivity * abs(x)
is_minor = (abs(ref) >= refa .or. refa >= refb)
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

function issymmetric(A, tol) result(is_symmetric)
!--------------------------------------------------------------------------------------------------!
! This function tests whether A is symmetric up to TOL.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, SYMTOL_DFT, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)
real(RP), intent(in), optional :: tol

! Outputs
logical :: is_symmetric

! Local variables
character(len=*), parameter :: srname = 'ISSYMMETRIC'
real(RP) :: tol_loc

! Preconditions
if (DEBUGGING) then
    if (present(tol)) then
        call assert(tol >= 0, 'TOL >= 0', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

tol_loc = SYMTOL_DFT
if (present(tol)) then
    tol_loc = tol
end if

! In Fortran, the following instructions cannot be written as the following Boolean expression:
!!IS_SYMMETRIC = (SIZE(A, 1)==SIZE(A, 2) .AND. &
!!    & ALL(IS_NAN(A) .EQV. IS_NAN(TRANSPOSE(A))) .AND. &
!!    & .NOT. ANY(ABS(A - TRANSPOSE(A)) > TOL_LOC * MAX(MAXVAL(ABS(A)), ONE)))
! This is because Fortran may not perform short-circuit evaluation of this expression. If A is not
! square, then IS_NAN(A) .EQV. IS_NAN(TRANSPOSE(A)) and A - TRANSPOSE(A) are invalid.
! In addition, since Inf - Inf is NaN, we cannot replace ANY(ABS(A - TRANSPOSE(A)) > TOL_LOC ...)
! with .NOT. ALL(ABS(A - TRANSPOSE(A)) <= TOL_LOC ...).
is_symmetric = .true.
if (size(A, 1) /= size(A, 2)) then
    is_symmetric = .false.
elseif (.not. all(is_nan(A) .eqv. is_nan(transpose(A)))) then
    is_symmetric = .false.
elseif (any(abs(A - transpose(A)) > tol_loc * max(maxval(abs(A)), ONE))) then
    is_symmetric = .false.
end if

!====================!
!  Calculation ends  !
!====================!

end function issymmetric


function p_norm(x, p) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates the P-norm of a vector X.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, TWO, ZERO
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_posinf
implicit none

! Inputs
real(RP), intent(in) :: x(:)
real(RP), intent(in), optional :: p

! Outputs
real(RP) :: y

! Local variables
character(len=*), parameter :: srname = 'P_NORM'
real(RP) :: scaling
real(RP) :: p_loc

if (present(p)) then
    p_loc = p
    call assert(p >= 0, 'P >= 0', srname)
else
    p_loc = TWO
end if

if (size(x) == 0) then
    y = ZERO
elseif (p_loc <= 0) then
    y = real(count(abs(x) > 0), kind(y))
elseif (.not. all(is_finite(x))) then
    ! If X contains NaN, then Y is NaN. Otherwise, Y is Inf when X contains +/-Inf.
    y = sum(abs(x))
elseif (.not. any(abs(x) > 0)) then
    ! The following is incorrect without checking the last case, as X may be all NaN.
    y = ZERO
else
    if (is_posinf(p_loc)) then
        ! If SIZE(X) = 0, then MAXVAL(ABS(X)) = -HUGE(X); since we have handled such a case in the
        ! above, it is OK to write Y = MAXVAL(ABS(X)) below, but we append 0 for robustness.
        y = maxval([abs(x), ZERO])
    elseif (.not. present(p) .or. abs(p_loc - TWO) <= 0) then
        !y = sqrt(sum(x**2))
        scaling = maxval(abs(x))  ! The scaling seems to reduce the rounding error.
        y = scaling * sqrt(sum((x / scaling)**2))
    else
        !scaling = max(REALMIN, sqrt(maxval(abs(x)) * minval(abs(x), mask=(abs(x) > ZERO))))
        scaling = maxval(abs(x))  ! The scaling seems to reduce the rounding error.
        y = scaling * sum(abs(x / scaling)**p_loc)**(ONE / p_loc)
    end if
end if
end function p_norm

function named_norm(x, nname) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates named norms of a vector X.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ZERO
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Inputs
real(RP), intent(in) :: x(:)
character(len=*), intent(in) :: nname

! Outputs
real(RP) :: y

! Local variables
character(len=*), parameter :: srname = 'NAMED_NORM'

if (size(x) == 0) then
    y = ZERO
elseif (.not. all(is_finite(x))) then
    ! If X contains NaN, then Y is NaN. Otherwise, Y is Inf when X contains +/-Inf.
    y = sum(abs(x))
elseif (.not. any(abs(x) > 0)) then
    ! The following is incorrect without checking the last case, as X may be all NaN.
    y = ZERO
else
    select case (lower(trimstr(nname)))
    case ('fro')
        y = p_norm(x)  ! 2-norm, which is the default case of P_NORM.
    case ('inf')
        ! If SIZE(X) = 0, then MAXVAL(ABS(X)) = -HUGE(X); since we have handled such a case in the
        ! above, it is OK to write Y = MAXVAL(ABS(X)) below, but we append 0 for robustness.
        y = maxval([abs(x), ZERO])
    case default
        y = -huge(x)  ! Without this, compiles complain that Y may be uninitialized.
        call errstop(srname, 'Unknown name of norm: '//trimstr(nname))
    end select
end if
end function named_norm


function sort_i1(x, direction) result(y)
! This function sorts X according to DIRECTION, which should be 'ascend' (default) or 'descend'.
use, non_intrinsic :: consts_mod, only : IK
#if __DEBUGGING__ == 1
use, non_intrinsic :: debug_mod, only : assert
#endif
implicit none

! Inputs
integer(IK), intent(in) :: x(:)
character(len=*), intent(in), optional :: direction

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
integer, intent(in), optional :: dim
character(len=*), intent(in), optional :: direction

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


pure elemental function logical_to_int(x) result(y)
!--------------------------------------------------------------------------------------------------!
! LOGICAL_TO_INT(.TRUE.) = 1, LOGICAL_TO_INT(.FALSE.) = 0
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
implicit none
logical, intent(in) :: x
integer(IK) :: y
y = merge(tsource=1_IK, fsource=0_IK, mask=x)
end function logical_to_int


function trueloc(x) result(loc)
!--------------------------------------------------------------------------------------------------!
! Similar to the `find` function in MATLAB, TRUELOC returns the indices where X is true.
! The motivation for this function is the fact that Fortran does not support logical indexing. See,
! for example, https://fortran-lang.discourse.group/t/indexing-arrays-by-an-array-of-logicals
! 1. MATLAB, Python, Julia, and R support logical indexing, so that the Fortran code Y(TRUELOC(X))
! can simply be translated to Y(X).
! 2. If the return of TRUELOC is NOT used for indexing, its analogs in other languages are:
! MATLAB -- find, Python -- numpy.argwhere, Julia -- findall, R -- which.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc
implicit none
logical, intent(in) :: x(:)
integer(IK), allocatable :: loc(:)  ! INTEGER(IK) :: LOC(COUNT(X)) does not work with Absoft 22.0
character(len=*), parameter :: srname = 'TRUELOC'
integer(IK) :: n
call safealloc(loc, int(count(x), IK))  ! Removable in F03.
n = int(size(x), IK)
loc = pack(linspace(1_IK, n, n), mask=x)
if (DEBUGGING) then
    call assert(all(loc >= 1 .and. loc <= n), '1 <= LOC <= N', srname)
    call assert(size(loc) == count(x), 'SIZE(LOC) == COUNT(X)', srname)
    call assert(all(x(loc)), 'X(LOC) is all TRUE', srname)
end if
end function trueloc


function falseloc(x) result(loc)
!--------------------------------------------------------------------------------------------------!
! FALSELOC = TRUELOC(.NOT. X)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc
implicit none
logical, intent(in) :: x(:)
integer(IK), allocatable :: loc(:)  ! INTEGER(IK) :: LOC(COUNT(.NOT.X)) does not work with Absoft 22.0
character(len=*), parameter :: srname = 'FALSELOC'
call safealloc(loc, int(count(.not. x), IK))  ! Removable in F03.
loc = trueloc(.not. x)
if (DEBUGGING) then
    call assert(all(loc >= 1 .and. loc <= size(x)), '1 <= LOC <= N', srname)
    call assert(size(loc) == size(x) - count(x), 'SIZE(LOC) == SIZE(X) - COUNT(X)', srname)
    call assert(all(.not. x(loc)), 'X(LOC) is all FALSE', srname)
end if
end function falseloc


function minimum1(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns NaN if X contains NaN; otherwise, it returns MINVAL(X). Vector version.
! F2018 does not specify MINVAL(X) when X contains NaN, which motivates this function. The behavior
! of this function is the same as the following functions in various languages:
! MATLAB: min(x, [], 'includenan')
! Python: numpy.min(x)
! Julia: minimum(x)
! R: min(x)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP) :: y

! Local variables
character(len=*), parameter :: srname = 'MINIMUM1'
real(RP) :: nan_test

!y = merge(tsource=sum(x), fsource=minval(x), mask=any(is_nan(x)))
nan_test = sum(abs(x))  ! 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
y = merge(tsource=nan_test, fsource=minval(x), mask=is_nan(nan_test))

if (DEBUGGING) then
    call assert(.not. any(x < y), 'No entry of X is smaller than Y', srname)
    call assert((.not. is_nan(y)) .or. any(is_nan(x)), 'Y is not NaN unless X contains NaN', srname)
    call assert(is_nan(y) .or. .not. any(is_nan(x)), 'Y is NaN if X contains NaN', srname)
end if
end function minimum1

function minimum2(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns NaN if X contains NaN; otherwise, it returns MINVAL(X). Matrix version.
! F2018 does not specify MINVAL(X) when X contains NaN, which motivates this function. The behavior
! of this function is the same as the following functions in various languages:
! MATLAB: min(x, [], 'all', 'includenan')
! Python: numpy.min(x)
! Julia: minimum(x)
! R: min(x)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: x(:, :)

! Outputs
real(RP) :: y

! Local variables
character(len=*), parameter :: srname = 'MINIMUM2'
real(RP) :: nan_test

!y = merge(tsource=sum(x), fsource=minval(x), mask=any(is_nan(x)))
nan_test = sum(abs(x))  ! 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
y = merge(tsource=nan_test, fsource=minval(x), mask=is_nan(nan_test))

if (DEBUGGING) then
    call assert(.not. any(x < y), 'No entry of X is smaller than Y', srname)
    call assert((.not. is_nan(y)) .or. any(is_nan(x)), 'Y is not NaN unless X contains NaN', srname)
    call assert(is_nan(y) .or. .not. any(is_nan(x)), 'Y is NaN if X contains NaN', srname)
end if
end function minimum2


function maximum1(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns NaN if X contains NaN; otherwise, it returns MAXVAL(X). Vector version.
! F2018 does not specify MAXVAL(X) when X contains NaN, which motivates this function. The behavior
! of this function is the same as the following functions in various languages:
! MATLAB: max(x, [], 'includenan')
! Python: numpy.max(x)
! Julia: maximum(x)
! R: max(x)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP) :: y

! Local variables
character(len=*), parameter :: srname = 'MAXIMUM1'
real(RP) :: nan_test

!y = merge(tsource=sum(x), fsource=maxval(x), mask=any(is_nan(x)))
nan_test = sum(abs(x))  ! 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
y = merge(tsource=nan_test, fsource=maxval(x), mask=is_nan(nan_test))

if (DEBUGGING) then
    call assert(.not. any(x > y), 'No entry of X is larger than Y', srname)
    call assert((.not. is_nan(y)) .or. any(is_nan(x)), 'Y is not NaN unless X contains NaN', srname)
    call assert(is_nan(y) .or. .not. any(is_nan(x)), 'Y is NaN if X contains NaN', srname)
end if
end function maximum1

function maximum2(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns NaN if X contains NaN; otherwise, it returns MAXVAL(X). Matrix version.
! F2018 does not specify MAXVAL(X) when X contains NaN, which motivates this function. The behavior
! of this function is the same as the following functions in various languages:
! MATLAB: max(x, [],  'all', 'includenan')
! Python: numpy.max(x)
! Julia: maximum(x)
! R: max(x)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
real(RP), intent(in) :: x(:, :)

! Outputs
real(RP) :: y

! Local variables
character(len=*), parameter :: srname = 'MAXIMUM2'
real(RP) :: nan_test

!y = merge(tsource=sum(x), fsource=maxval(x), mask=any(is_nan(x)))
nan_test = sum(abs(x))  ! 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
y = merge(tsource=nan_test, fsource=maxval(x), mask=is_nan(nan_test))

if (DEBUGGING) then
    call assert(.not. any(x > y), 'No entry of X is larger than Y', srname)
    call assert((.not. is_nan(y)) .or. any(is_nan(x)), 'Y is not NaN unless X contains NaN', srname)
    call assert(is_nan(y) .or. .not. any(is_nan(x)), 'Y is NaN if X contains NaN', srname)
end if
end function maximum2


function linspace_r(xstart, xstop, n) result(x)
!--------------------------------------------------------------------------------------------------!
! Similar to the function `linspace` in MATLAB and Python, this function generates N evenly spaced
! numbers, the space between the consecutive points being (XSTOP-XSTART)/(N-1).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
real(RP), intent(in) :: xstart
real(RP), intent(in) :: xstop
integer(IK), intent(in) :: n

! Outputs
real(RP), allocatable :: x(:)

! Local variables
character(len=*), parameter :: srname = 'LINSPACE_R'
integer(IK) :: i
integer(IK) :: nm
real(RP) :: xunit

!====================!
! Calculation starts !
!====================!

call safealloc(x, max(n, 0_IK))

if (n <= 0) then  ! Quick return when N <= 0.
    return
end if

nm = n - 1_IK

if (n == 1 .or. (xstart <= xstop .and. xstop <= xstart)) then
    x = xstop
elseif (abs(xstart) <= abs(xstop) .and. abs(xstop) <= abs(xstart)) then
    xunit = xstop / real(nm, RP)
    x = xunit * real([(i, i=-nm, nm, 2_IK)], RP)
    if (modulo(nm, 2_IK) == 0) then
        x(1_IK + nm / 2_IK) = ZERO
    end if
else
    xunit = (xstop - xstart) / real(nm, RP)
    x = xstart + xunit * real([(i, i=0, nm)], RP)
end if

if (n >= 1) then  ! Indeed, N < 1 cannot happen due to the quick return when N <= 0.
    x(1) = xstart
    x(n) = xstop
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(x) == max(n, 0_IK), 'SIZE(X) == MAX(N, 0)', srname)
end if
end function linspace_r

function linspace_i(xstart, xstop, n) result(x)
!--------------------------------------------------------------------------------------------------!
! This function returns INT(LINSPACE_R(REAL(XSTART, RP), REAL(XSTOP, RP), N), IK).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in) :: xstart
integer(IK), intent(in) :: xstop
integer(IK), intent(in) :: n

! Outputs
integer(IK), allocatable :: x(:)

! Local variables
character(len=*), parameter :: srname = 'LINSPACE_I'

!====================!
! Calculation starts !
!====================!

call safealloc(x, max(n, 0_IK))
x = int(linspace_r(real(xstart, RP), real(xstop, RP), n), IK)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(x) == max(n, 0_IK), 'SIZE(X) == MAX(N, 0)', srname)
end if
end function linspace_i


subroutine hessenberg_hhd_trid(A, tdiag, tsubdiag)
!--------------------------------------------------------------------------------------------------!
! This subroutine applies Householder transformations to obtain a tridiagonal matrix that is similar
! a SYMMETRIC matrix A. The tridiagonal matrix is the Hessenberg form of A; its diagonal will be
! stored in TDIAD, and the sub diagonal in TSUBDIAG. The Householder vectors will be stored in the
! lower triangular part of A.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! In-outputs
real(RP), intent(inout) :: A(:, :)

! Outputs
real(RP), intent(out) :: tdiag(:)
real(RP), intent(out) :: tsubdiag(:)

! Local variables
character(len=*), parameter :: srname = 'HESSENBERG_HHD_TRID'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: k
integer(IK) :: n
real(RP) :: Asubd
real(RP) :: colsq
real(RP) :: w(size(A, 1))
real(RP) :: wz
real(RP) :: z(size(A, 1))

! Sizes
n = int(size(A, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    ! Even though we only need the lower triangular part of A, we assume that something is wrong if
    ! this subroutine is invoked with a non-symmetric matrix A.
    call assert(issymmetric(A), 'A is symmetric', srname)
    call assert(size(tdiag) == n, 'SIZE(TDIAG) == N', srname)
    call assert(size(tsubdiag) == max(0_IK, n - 1_IK), 'SIZE(TDIAG) == MAX(0, N-1)', srname)
end if

!====================!
! Calculation starts !
!====================!

if (n <= 0) then  ! Quick return when N <= 0. Of course, N < 0 is impossible.
    return
end if

tdiag = diag(A)

do k = 1, n - 1
    colsq = sum(A(k + 2:n, k)**2)
    if (colsq <= 0) then
        tsubdiag(k) = A(k + 1, k)
        A(k + 1, k) = ZERO
        cycle
    end if

    Asubd = A(k + 1, k)
    tsubdiag(k) = sign(sqrt(colsq + Asubd**2), Asubd)
    A(k + 1, k) = -colsq / (Asubd + tsubdiag(k))
    w(k + 1:n) = sqrt(TWO / (colsq + A(k + 1, k)**2)) * A(k + 1:n, k)
    A(k + 1:n, k) = w(k + 1:n)
    z(k + 1:n) = tdiag(k + 1:n) * w(k + 1:n)

    do j = k + 1, n - 1
        z(j + 1:n) = z(j + 1:n) + A(j + 1:n, j) * w(j)
        do i = j + 1, n
            z(j) = z(j) + A(i, j) * w(i)
        end do
    end do

    wz = inprod(w(k + 1:n), z(k + 1:n))

    tdiag(k + 1:n) = tdiag(k + 1:n) + w(k + 1:n) * (wz * w(k + 1:n) - TWO * z(k + 1:n))
    do j = k + 1, n
        A(j + 1:n, j) = A(j + 1:n, j) - w(j + 1:n) * z(j) - w(j) * (z(j + 1:n) - wz * w(j + 1:n))
    end do
end do

!====================!
!  Calculation ends  !
!====================!

end subroutine hessenberg_hhd_trid

subroutine hessenberg_full(A, H, Q)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds a Hessenberg matrix H (all entries below the subdiagonal are 0) such that
! H = Q^T*A*Q, where Q is a orthogonal matrix that may also be returned. A will stay unchanged.
! As of 20220507, we ONLY handle the case with a symmetric A by invoking HESSENBERG_HHD_TRID.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)

! Outputs
real(RP), intent(out) :: H(:, :)
real(RP), intent(out), optional :: Q(:, :)

! Local variables
character(len=*), parameter :: srname = 'HESSENBERG_FULL'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: n
real(RP) :: tdiag(size(A, 1))
real(RP) :: tsubdiag(size(A, 1) - 1)

! Debugging variables
real(RP) :: tol

! Sizes
n = int(size(A, 1), kind(n))

!====================!
! Calculation starts !
!====================!

! Preconditions
if (DEBUGGING) then
    call assert(size(A, 1) == size(A, 2), 'A is square', srname)
    call assert(size(H, 1) == n .and. size(H, 2) == n, 'SIZE(H) == [N, N]', srname)
    if (present(Q)) then
        call assert(size(Q, 1) == n .and. size(Q, 2) == n, 'SIZE(Q) == [N, N]', srname)
    end if
end if

if (n <= 0) then  ! Quick return when N <= 0. Of course, N < 0 is impossible.
    return
end if

H = A
call hessenberg_hhd_trid(H, tdiag, tsubdiag)

if (present(Q)) then
    Q = eye(n, n)
    do j = n - 1, 1, -1
        do i = 1, n
            Q(j + 1:n, i) = Q(j + 1:n, i) - inprod(Q(j + 1:n, i), H(j + 1:n, j)) * H(j + 1:n, j)
        end do
    end do
end if

H = ZERO
do j = 1, n
    H(j, j) = tdiag(j)
    if (j > 1) then
        H(j - 1, j) = tsubdiag(j - 1)
    end if
    if (j < n) then
        H(j + 1, j) = tsubdiag(j)
    end if
end do

!====================!
!  Calculation ends  !
!====================!

if (DEBUGGING) then
    call assert(size(H, 1) == n .and. size(H, 2) == n, 'SIZE(H) == [N, N]', srname)
    call assert(isbanded(H, 1_IK, n - 1_IK), 'H is a Hessenberg matrix', srname)
    call assert(issymmetric(H) .or. .not. issymmetric(A), 'H is symmetric if so is A', srname)
    if (present(Q)) then
        call assert(size(Q, 1) == n .and. size(Q, 2) == n, 'SIZE(Q) == [N, N]', srname)
        tol = max(1.0E-8_RP, min(1.0E-1_RP, 1.0E10_RP * EPS * real(n, RP)))
        call assert(isorth(Q, tol), 'Q is orthogonal', srname)
        call assert(all(abs(matprod(Q, H) - matprod(A, Q)) <= tol * maxval(abs(A))), 'Q*H = Q*Q', srname)
    end if
end if

end subroutine hessenberg_full


end module linalg_mod
