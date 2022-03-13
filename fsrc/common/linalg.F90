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
! 2. Most of the subroutines/functions here are better coded inline, because the code is short
! using matrix/vector operations, and because the overhead of subroutine/function calling can be
! high in these languages. Here we implement them as subroutines/functions in order to align with
! Powell's original code, which cannot be translated directly to matrix/vector operations that
! produce the same results in floating-point arithmetic.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Sunday, March 13, 2022 PM11:50:49
!--------------------------------------------------------------------------------------------------

implicit none
private
! Mathematically, inprod = dot_product, matprod = matmul
public :: inprod, matprod, outprod
public :: r1update, r2update, symmetrize
public :: Ax_plus_y
public :: eye
public :: hypotenuse, planerot, lsqr
public :: project
public :: inv, isinv
public :: qradd, qrexc
public :: calquad, errquad, hess_mul
public :: omega_col, omega_mul, omega_inprod
public :: errh
public :: isminor
public :: issymmetric, isorth, istriu
public :: norm
public :: sort
public :: int
public :: trueloc, falseloc
public :: minimum, maximum

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

interface project
    module procedure project1, project2
end interface

interface qradd
    module procedure qradd_Rdiag, qradd_Rfull
end interface

interface qrexc
    module procedure qrexc_Rdiag, qrexc_Rfull
end interface

interface isminor
    module procedure isminor0, isminor1
end interface isminor

interface sort
    module procedure sort_i1, sort_i2
end interface sort

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
integer(IK) :: PI(size(A, 1))
integer(IK) :: i
integer(IK) :: n
real(RP) :: Q(size(A, 1), size(A, 1))
real(RP) :: R(size(A, 1), size(A, 1))
real(RP) :: tol

! Sizes
n = int(size(A, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(A, 1) == size(A, 2), 'A is squre', srname)
end if

!====================!
! Calculation starts !
!====================!

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
    ! This is NOT the best algorithm for INV, but since the QR subroutine is available ...
    call qr(A, Q, R, P)
    R = transpose(R)  ! Take transpose to work on columns.
    do i = n, 1, -1
        B(:, i) = (Q(:, i) - matprod(B(:, i + 1:n), R(i + 1:n, i))) / R(i, i)
    end do
    PI(P) = [(i, i=1, n)]  ! The inverse permutation
    B = transpose(B(:, PI))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(B, 1) == n .and. size(B, 2) == n, 'SIZE(B) == [N, N]', srname)
    call assert(istril(B) .or. .not. istril(A), 'If A is lower triangular, then so is B', srname)
    call assert(istriu(B) .or. .not. istriu(A), 'If A is upper triangular, then so is B', srname)
    tol = max(1.0E-8_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n + 1_IK, RP)))
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
    P = [(j, j=1, n)]
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


function lsqr(A, b, Q, Rdiag) result(x)
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
real(RP), intent(in) :: A(:, :)
real(RP), intent(in) :: b(:)
real(RP), intent(in), optional :: Q(:, :)
real(RP), intent(in), optional :: Rdiag(:)

! Outputs
real(RP) :: x(size(A, 2))

! Local variables
character(len=*), parameter :: srname = 'LSQR'
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

if (present(Q)) then
    Q_loc = Q(:, 1:size(Q_loc, 2))
    if (present(Rdiag)) then
        Rdiag_loc = Rdiag
    else
        Rdiag_loc = [(inprod(Q_loc(:, i), A(:, i)), i=1, min(m, n))]
    end if
    rank = min(m, n)
    pivote = .false.
end if

!if (.not. present(Q) .or. any(abs(Rdiag_loc) <= ZERO)) then  ! This is more reasonable
if (.not. present(Q)) then
    call qr(A, Q=Q_loc, P=P)
    Rdiag_loc = [(inprod(Q_loc(:, i), A(:, P(i))), i=1, min(m, n))]
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
end function lsqr


function diag(A) result(D)
!--------------------------------------------------------------------------------------------------!
! This function takes the main diagonal of the matrix A.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
! Inputs
real(RP), intent(in) :: A(:, :)
! Outputs
real(RP) :: D(min(size(A, 1), size(A, 2)))
! Local variables
integer(IK) :: i
D = [(A(i, i), i=1, int(size(D), IK))]
end function diag


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
integer(IK) :: i
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol_loc

if (present(tol)) then
    tol_loc = max(tol, tol * maxval(abs(A)))
else
    tol_loc = ZERO
end if
m = int(size(A, 1), kind(m))
n = int(size(A, 2), kind(n))
is_tril = .true.
do i = 1, min(m, n)
    if (any(abs(A(1:i - 1, i)) > tol_loc)) then
        is_tril = .false.
        exit
    end if
end do
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
integer(IK) :: i
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol_loc

if (present(tol)) then
    tol_loc = max(tol, tol * maxval(abs(A)))
else
    tol_loc = ZERO
end if
m = int(size(A, 1), kind(m))
n = int(size(A, 2), kind(n))
is_triu = .true.
do i = 1, min(m, n)
    if (any(abs(A(i + 1:m, i)) > tol_loc)) then
        is_triu = .false.
        exit
    end if
end do
end function istriu


function isorth(A, tol) result(is_orth)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the matrix A has orthonormal columns up to the tolerance TOL.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
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
    where (is_inf(v))
        u = sign(ONE, v)
    elsewhere
        u = ZERO
    end where
    u = u / norm(u)
!    y = inprod(x, u) * u
    scaling = maxval(abs(x))  ! The scaling seems to reduce the rounding error.
    y = scaling * inprod(x / scaling, u) * u
else
    u = v / norm(v)
!    y = inprod(x, u) * u
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
    where (.not. is_inf(V))
        V_loc = ZERO
    elsewhere
        V_loc = sign(ONE, V)
    end where
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
! However, we need to take care of the possibilities of R = 0, Inf, NaN, and over/underflow.
! The G defined in this way is continuous with respect to X except at 0. Following this definition,
! G = [sign(X(1)), 0; 0, sign(X(1))] if X(2) = 0, G = [0, sign(X(2)); -sign(X(2)), 0] if X(2) = 0.
! Yet some implementations ignore the signs, leading to discontinuity and numerical instability.
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
    ! In this case, MATLAB sets G to NaN(2, 2). We refrain from doing this to keep G orthogonal.
    c = ONE
    s = ZERO
elseif (all(is_inf(x))) then
    ! In this case, MATLAB sets G to NaN(2, 2). We refrain from doing this to keep G orthogonal.
    c = sign(1 / sqrt(2.0_RP), x(1))
    s = sign(1 / sqrt(2.0_RP), x(2))
elseif (abs(x(1)) <= 0 .and. abs(x(2)) <= 0) then  ! X(1) == 0 == X(2).
    c = ONE
    s = ZERO
elseif (abs(x(2)) <= EPS * abs(x(1))) then
    ! N.B.:
    ! 0. With <= instead of <, this case covers X(1)==0==X(2), which is treated above separately to
    ! avoid the confusing SIGN(.,0) (see 1).
    ! 1. SIGN(A, 0) = ABS(A) in Fortran but sign(0) = 0 in MATLAB, Python, Julia, and R!!!
    ! 2. Taking SIGN(X(1)) into account ensures the continuity of G with respect to X except at 0.
    c = sign(ONE, x(1))  ! MATLAB: c = sign(x(1))
    s = ZERO
elseif (abs(x(1)) <= EPS * abs(x(2))) then
    ! N.B.: SIGN(A, X) = ABS(A) * sign of X /= A * sign of X !!! Therefore, it is WRONG to define G
    ! as SIGN(RESHAPE([ZERO, -ONE, ONE, ZERO], [2, 2]), X(2)). ! This mistake was committed on on
    ! 20211206, taking a whole day to debug! NEVER use SIGN on arrays unless you are really sure.
    c = ZERO
    s = sign(ONE, x(2))  ! MATLAB: s = sign(x(2))
else
    ! The following is a stable and continuous implementation of the Givens rotation. It follows
    ! Bindel, D., Demmel, J., Kahan, W., & Marques, O. (2002). On computing Givens rotations reliably
    ! and efficiently. ACM Transactions on Mathematical Software (TOMS), 28(2), 206-238.
    ! 1. Modern compilers compute SQRT(REALMIN) and SQRT(HUGENUM/2.1) at compilation time.
    ! 2. The direct calculation without involving T and U seems to work better; use it if possible.
    if (minval(abs(x)) > sqrt(REALMIN) .and. maxval(abs(x)) < sqrt(HUGENUM / 2.1_RP)) then
        ! Do NOT use HYPOTENUSE here; the best implementation for one may not be the best for the other
        r = sqrt(sum(x**2))
        c = x(1) / r
        s = x(2) / r
    elseif (abs(x(1)) > abs(x(2))) then
        t = x(2) / x(1)
        u = sign(sqrt(ONE + t**2), x(1))  ! MATLAB: u = sign(x(1))*sqrt(ONE + t**2)
        c = ONE / u
        s = t / u
    else
        t = x(1) / x(2)
        u = sign(sqrt(ONE + t**2), x(2))  ! MATLAB: u = sign(x(2))*sqrt(ONE + t**2)
        c = t / u
        s = ONE / u
    end if
end if

G = reshape([c, -s, s, c], [2, 2])  ! MATLAB: G = [c, s; -s, c]

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
    if (maxval(abs(x)) < sqrt(HUGENUM / 2.1_RP)) then
        r = sqrt(sum(x**2))
        call assert(norm(matprod(G, x) - [r, ZERO]) <= max(tol, tol * r), 'G*x = [|x|, 0]', srname)
    end if
end if

end function planerot


subroutine qradd_Rdiag(c, Q, Rdiag, n)  ! Used in COBYLA
!--------------------------------------------------------------------------------------------------!
! This subroutine updates the QR factorization of an MxN matrix A of full column rank, attempting to
! add a new column C is to this matrix as the LAST column while maintaining the full-rankness.
! Case 1. If C is not in range(A) (theoretically, it implies N < M), then the new matrix is [A, C];
! Case 2. If C is in range(A), then the new matrix is [A(:, 1:N-1), C].
! N.B.:
! 0. Instead of R, this subroutine updates RDIAG, which is diag(R), with a size at most M and at
! least MIN(M, N+1). The number is MIN(M, N+1) rather than MIN(M, N) as N may be augmented by 1 in
! the subroutine.
! 1. With the two cases specified as above, this function does not need A as an input.
! 2. The subroutine changes only Q(:, NSAVE+1:M) (NSAVE is the original value of N)
! and R(:, N) (N takes the updated value).
! 3. Indeed, when C is in range(A), Powell wrote in comments that "set IOUT to the index of the
! constraint (here, column of A -- Zaikun) to be deleted, but branch if no suitable index can be
! found". The idea is to replace a column of A by C so that the new matrix still has full rank
! (such a column must exist unless C = 0). But his code essentially sets IOUT = N always. Maybe he
! found this worked well enough in practice. Meanwhile, Powell's code includes a snippet that can
! never be reached, which was probably intended to deal with the case with IOUT =/= N.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
implicit none

! Inputs
real(RP), intent(in) :: c(:)  ! C(M)

! In-outputs
integer(IK), intent(inout) :: n
real(RP), intent(inout) :: Q(:, :)  ! Q(M, M)
real(RP), intent(inout) :: Rdiag(:)  ! MIN(M, N+1) <= SIZE(Rdiag) <= M

! Local variables
character(len=*), parameter :: srname = 'QRADD_RDIAG'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: nsave
real(RP) :: cq(size(Q, 2))
real(RP) :: cqa(size(Q, 2))
real(RP) :: G(2, 2)
real(RP) :: Qsave(size(Q, 1), n)  ! Debugging only
real(RP) :: Rdsave(n)  ! Debugging only
real(RP) :: tol

! Sizes
m = int(size(Q, 2), kind(m))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 0 .and. n <= m, '0 <= N <= M', srname)  ! N = 0 is possible.
    call assert(size(c) == m, 'SIZE(C) == M', srname)
    call assert(size(Rdiag) >= min(m, n + 1_IK) .and. size(Rdiag) <= m, 'MIN(M, N+1) <= SIZE(Rdiag) <= M', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) == m, 'SIZE(Q) == [M, M]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E6_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(Q, tol), 'The columns of Q are orthonormal', srname)  !! Costly!
    Qsave = Q(:, 1:n)  ! For debugging only
    Rdsave = Rdiag(1:n)  ! For debugging only
end if

!====================!
! Calculation starts !
!====================!

nsave = n  ! Needed for debugging (only).

! As in Powell's COBYLA, CQ is set to 0 at the positions where CQ is negligible according to ISMINOR.
! It may not be the best choice when the subroutine is used elsewhere, e.g., LINCOA. Tests needed.
cq = matprod(c, Q)
cqa = matprod(abs(c), abs(Q))
where (isminor(cq, cqa))  ! Code in MATLAB: CQ(ISMINOR(CQ, CQA)) = ZERO
    cq = ZERO
end where

! Update Q so that the columns of Q(:, N+2:M) are orthogonal to C. This is done by applying a 2D
! Givens rotation to Q(:, [K, K+1]) from the right to zero C'*Q(:, K+1) out for K = N+1, ..., M-1.
! Nothing will be done if N >= M-1.
do k = m - 1_IK, n + 1_IK, -1
    if (abs(cq(k + 1)) > 0) then
        ! Powell wrote CQ(K+1) /= 0 instead of ABS(CQ(K+1)) > 0. The two differ if CQ(K+1) is NaN.
        ! If we apply the rotation below when CQ(K+1) = 0, then CQ(K) will get updated to |CQ(K)|.
        G = planerot(cq([k, k + 1_IK]))
        Q(:, [k, k + 1_IK]) = matprod(Q(:, [k, k + 1_IK]), transpose(G))
        cq(k) = hypotenuse(cq(k), cq(k + 1))
        !cq(k) = sqrt(cq(k)**2 + cq(k + 1)**2)
    end if
end do

! Augment N by 1 if C is not in range(A).
! The two IFs cannot be merged as Fortran may evaluate CQ(N+1) even if N>=M, leading to a SEGFAULT.
if (n < m) then
    if (abs(cq(n + 1)) > 0) then  ! C is not in range(A).
        ! Powell wrote CQ(N+1) /= 0 instead of ABS(CQ(N+1)) > 0. The two differ if CQ(N+1) is NaN.
        n = n + 1_IK
    end if
end if

! Update RDIAG so that RDIAG(N) = CQ(N) = INPROD(C, Q(:, N)). Note that N may have been augmented.
if (n >= 1 .and. n <= m) then  ! Indeed, N > M should not happen unless the input is wrong.
    Rdiag(n) = cq(n)  ! Indeed, RDIAG(N) = INPROD(C, Q(:, N))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(n >= nsave .and. n <= min(nsave + 1_IK, m), 'NSAV <= N <= MIN(NSAV + 1, M)', srname)
    call assert(size(Rdiag) >= n .and. size(Rdiag) <= m, 'N <= SIZE(Rdiag) <= M', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) == m, 'SIZE(Q) == [M, M]', srname)
    call assert(isorth(Q, tol), 'The columns of Q are orthonormal', srname)  !! Costly!
    call assert(all(abs(Q(:, 1:nsave) - Qsave(:, 1:nsave)) <= 0), 'Q(:, 1:NSAVE) is unchanged', srname)
    call assert(all(abs(Rdiag(1:n - 1) - Rdsave(1:n - 1)) <= 0), 'Rdiag(1:N-1) is unchanged', srname)
    if (n < m .and. is_finite(norm(c))) then
        call assert(norm(matprod(c, Q(:, n + 1:m))) <= max(tol, tol * norm(c)), 'C^T*Q(:, N+1:M) == 0', srname)
    end if
    if (n >= 1) then  ! N = 0 is possible.
        call assert(abs(inprod(c, Q(:, n)) - Rdiag(n)) <= max(tol, tol * inprod(abs(c), abs(Q(:, n)))) &
            & .or. .not. is_finite(Rdiag(n)), 'C^T*Q(:, N) == Rdiag(N)', srname)
    end if
end if
end subroutine qradd_Rdiag


subroutine qrexc_Rdiag(A, Q, Rdiag, i)  ! Used in COBYLA
!--------------------------------------------------------------------------------------------------!
! This subroutine updates the QR factorization for an MxN matrix A = Q*R so that the updated Q and
! R form a QR factorization of [A_1, ..., A_{I-1}, A_{I+1}, ..., A_N, A_I], which is the matrix
! obtained by rearranging columns [I, I+1, ..., N] of A to [I+1, ..., N, I]. Here, A is ASSUMED TO
! BE OF FULL COLUMN RANK, Q is a matrix whose columns are orthogonal, and R, which is not present,
! is an upper triangular matrix whose diagonal entries are nonzero. Q and R need not to be square.
! N.B.:
! 0. Instead of R, this subroutine updates RDIAG, which is diag(R), the size being N.
! 1. With L = SIZE(Q, 2) = SIZE(R, 1), we have M >= L >= N. Most often, L = M or N.
! 2. The subroutine changes only Q(:, I:N) and Rdiag(I:N).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)

! In-outputs
real(RP), intent(inout) :: Q(:, :)  ! Q(M, :), N <= SIZE(Q, 2) <= M
real(RP), intent(inout) :: Rdiag(:)  ! Rdiag(N)
integer(IK), intent(in) :: i

! Local variables
character(len=*), parameter :: srname = 'QREXC_RDIAG'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
real(RP) :: Anew(size(A, 1), size(A, 2))
real(RP) :: G(2, 2)
real(RP) :: Qsave(size(Q, 1), size(Q, 2))  ! Debugging only
real(RP) :: QtAnew(size(Q, 2), size(A, 2))  ! Debugging only
real(RP) :: Rdsave(i)  ! Debugging only
real(RP) :: tol

! Sizes
m = int(size(A, 1), kind(m))
n = int(size(A, 2), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. n <= m, '1 <= N <= M', srname)
    call assert(i >= 1 .and. i <= n, '1 <= i <= N', srname)
    call assert(size(Rdiag) == n, 'SIZE(Rdiag) == N', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) >= n .and. size(Q, 2) <= m, &
        & 'SIZE(Q, 1) == M, N <= SIZE(Q, 2) <= M', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(Q, tol), 'The columns of Q are orthonormal', srname)  !! Costly!
    Qsave = Q  ! For debugging only.
    Rdsave = Rdiag(1:i) ! For debugging only.
end if

!====================!
! Calculation starts !
!====================!

if (i <= 0 .or. i >= n) then
    ! Only I == N is really needed, as 1 <= I <= N unless the input is wrong.
    return
end if

! Let R be the upper triangular matrix in the QR factorization, namely R = Q^T*A.
! For each K, find the Givens rotation G with G*R([K, K+1], :) = [HYPT, 0], and update Q(:, [K,K+1])
! to Q(:, [K, K+1])*G^T. Then R = Q^T*A is an upper triangular matrix as long as A(:, [K, K+1]) is
! updated to A(:, [K+1, K]). Indeed, this new upper triangular matrix can be obtained by first
! updating R([K, K+1], :) to G*R([K, K+1], :) and then exchanging its columns K and K+1; at the same
! time, entries K and K+1 of R's diagonal RDIAG become [HYPT, -(RDIAG(K+1) / HYPT) * RDIAG(K)].
! After this is done for each K = 1, ..., N-1, we obtain the QR factorization of the matrix that
! rearranges columns [I, I+1, ..., N] of A as [I+1, ..., N, I].
! Powell's code, however, is slightly different: before everything, he first exchanged columns K and
! K+1 of Q (as well as rows K and K+1 of R). This makes sure that the entires of the update RDIAG
! are all positive if it is the case for the original RDIAG.
do k = i, n - 1_IK
    !hypt = hypotenuse(Rdiag(k + 1), inprod(Q(:, k), A(:, k + 1)))
    !hypt = sqrt(Rdiag(k + 1)**2 + inprod(Q(:, k), A(:, k + 1))**2)
    G = planerot([Rdiag(k + 1), inprod(Q(:, k), A(:, k + 1))])
    Q(:, [k, k + 1_IK]) = matprod(Q(:, [k + 1_IK, k]), transpose(G))

    ! Powell's code updates RDIAG in the following way.
    !----------------------------------------------------------------!
    !!Rdiag([k, k + 1_IK]) = [hypt, (Rdiag(k + 1) / hypt) * Rdiag(k)]!
    !----------------------------------------------------------------!
    ! Note that RDIAG(N) inherits all rounding in RDIAG(I:N-1) and Q(:, I:N-1) and hence contain
    ! significant errors. Thus we may modify Powell's code to set only RDIAG(K) = HYPT here and then
    ! calculate RDIAG(N) by an inner product after the loop. Nevertheless, we simply calculate RDIAG
    ! from scratch we do below.
end do

! Calculate RDIAG(I:N) from scratch.
Rdiag(i:n - 1) = [(inprod(Q(:, k), A(:, k + 1)), k=i, n - 1_IK)]
Rdiag(n) = inprod(Q(:, n), A(:, i))  ! Calculate RDIAG(N) from scratch. See the comments above.

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(Rdiag) == n, 'SIZE(Rdiag) == N', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) >= n .and. size(Q, 2) <= m, &
        & 'SIZE(Q, 1) == M, N <= SIZE(Q, 2) <= M', srname)
    call assert(isorth(Q, tol), 'The columns of Q are orthonormal', srname)  !! Costly!
    call assert(all(abs(Q(:, 1:i - 1) - Qsave(:, 1:i - 1)) <= 0) .and. &
        & all(abs(Q(:, n + 1:) - Qsave(:, n + 1:)) <= 0), 'Q is unchanged except Q(:, I:N)', srname)
    call assert(all(abs(Rdiag(1:i - 1) - Rdsave(1:i - 1)) <= 0), 'Rdiag(1:I-1) is unchanged', srname)
    Anew = reshape([A(:, 1:i - 1), A(:, i + 1:n), A(:, i)], shape(A))
    QtAnew = matprod(transpose(Q), Anew)
    call assert(istriu(QtAnew, tol), 'Q^T*Anew is upper triangular', srname)
    ! The following test may fail if RDIAG is not calculated from scratch.
    call assert(norm(diag(QtAnew) - Rdiag) <= max(tol, tol * norm([(inprod(abs(Q(:, k)), &
        & abs(Anew(:, k))), k=1, n)])), 'Rdiag == diag(Q^T*Anew)', srname)
end if
end subroutine qrexc_Rdiag


subroutine qradd_Rfull(c, Q, R, n)  ! Used in LINCOA
!--------------------------------------------------------------------------------------------------!
! This subroutine updates the QR factorization of an MxN matrix A = Q*R(:, 1:N) when a new column C
! is appended to this matrix A as the LAST column.
! N.B.:
! 0. Different from QRADD_RDIAG, it seems that QRADD_RFULL does not try to maintain that A is of
! full column rank after the update; QRADD_RFULL always append C to A, and always increase N by 1,
! but QRADD_RDIAG does so only if C is not in the column space of A.
! 1. At entry, Q is a MxM orthonormal matrix, and R is a MxL upper triangular matrix with N < L <= M.
! 2. The subroutine changes only Q(:, N+1:M) and R(:, N+1), where N takes the original value.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: c(:)

! In-outputs
integer(IK), intent(inout) :: n
real(RP), intent(inout) :: Q(:, :)  ! Q(M, M)
real(RP), intent(inout) :: R(:, :)  ! R(M, :), N+1 <= SIZE(R, 2) <= M

! Local variables
character(len=*), parameter :: srname = 'QRADD_RFULL'
integer(IK) :: k
integer(IK) :: m
real(RP) :: cq(size(Q, 2))
real(RP) :: G(2, 2)
real(RP) :: Qsave(size(Q, 1), n)
real(RP) :: Rsave(size(R, 1), n)
real(RP) :: tol

! Sizes
m = int(size(Q, 1), kind(m))

if (DEBUGGING) then
    call assert(n >= 0 .and. n <= m - 1, '0 <= N <= M - 1', srname)
    call assert(size(c) == m, 'SIZE(C) == M', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) == m, 'SIZE(Q) = [M, M]', srname)
    call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    call assert(size(R, 2) >= n + 1 .and. size(R, 2) <= m, 'N+1 <= SIZE(R, 2) <= M', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    Qsave = Q(:, 1:n)  ! For debugging only.
    Rsave = R(:, 1:n)  ! For debugging only.
end if

cq = matprod(c, Q)

! Update Q so that the columns of Q(:, N+2:M) are orthogonal to C. This is done by applying a 2D
! Givens rotation to Q(:, [K, K+1]) from the right to zero C'*Q(:, K+1) out for K = N+1, ..., M-1.
! Nothing will be done if N >= M-1.
do k = m - 1_IK, n + 1_IK, -1
    if (abs(cq(k + 1)) > 0) then  ! Powell: IF (ABS(CQ(K + 1)) > 1.0D-20 * ABS(CQ(K))) THEN
        G = planerot(cq([k, k + 1_IK]))
        Q(:, [k, k + 1_IK]) = matprod(Q(:, [k, k + 1_IK]), transpose(G))
        cq(k) = sqrt(cq(k)**2 + cq(k + 1)**2)
    end if
end do

R(1:n, n + 1) = matprod(c, Q(:, 1:n))

! Maintain the positiveness of the diagonal entries of R.
if (cq(n + 1) < 0) then
    Q(:, n + 1) = -Q(:, n + 1)
end if
R(n + 1, n + 1) = abs(cq(n + 1))

n = n + 1_IK

if (DEBUGGING) then
    call assert(n >= 1 .and. n <= m, '1 <= N <= M', srname)
    call assert(size(Q, 1) == m .and. size(Q, 2) == m, 'SIZE(Q) = [M, M]', srname)
    call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    call assert(size(R, 2) >= n .and. size(R, 2) <= m, 'N <= SIZE(R, 2) <= M', srname)
    call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    !call assert(.not. any(abs(Q(:, 1:n - 1) - Qsave(:, 1:n - 1)) > 0), 'Q(:, 1:N-1) is unchanged', srname)
    !call assert(.not. any(abs(R(:, 1:n - 1) - Rsave(:, 1:n - 1)) > 0), 'R(:, 1:N-1) is unchanged', srname)
    ! If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead.
    call assert(all(abs(Q(:, 1:n - 1) - Qsave(:, 1:n - 1)) <= 0), 'Q(:, 1:N-1) is unchanged', srname)
    call assert(all(abs(R(:, 1:n - 1) - Rsave(:, 1:n - 1)) <= 0), 'R(:, 1:N-1) is unchanged', srname)
end if
end subroutine qradd_Rfull


subroutine qrexc_Rfull(Q, R, i)  ! Used in LINCOA
!--------------------------------------------------------------------------------------------------!
! This subroutine updates the QR factorization for an MxN matrix A = Q*R so that the updated Q and
! R form a QR factorization of [A_1, ..., A_{I-1}, A_{I+1}, ..., A_N, A_I], which is the matrix
! obtained by rearranging columns [I, I+1, ..., N] of A to [I+1, ..., N, I]. At entry, A = Q*R
! is ASSUMED TO BE OF FULL COLUMN RANK, Q is a matrix whose columns are orthogonal, and R is an
! upper triangular matrix whose diagonal entries are all nonzero. Q and R need not to be square.
! N.B.:
! 1. With L = SIZE(Q, 2) = SIZE(R, 1), we have M >= L >= N. Most often, L = M or N.
! 2. The subroutine changes only Q(:, I:N) and R(:, I:N).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: i

! In-outputs
real(RP), intent(inout) :: Q(:, :)  ! Q(M, :), SIZE(Q, 2) <= M
real(RP), intent(inout) :: R(:, :)  ! R(:, N), SIZE(R, 1) >= N

! Local variables
character(len=*), parameter :: srname = 'QREXC_RFULL'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
real(RP) :: G(2, 2)
real(RP) :: hypt
real(RP) :: Qsave(size(Q, 1), size(Q, 2))  ! Debugging only
real(RP) :: Rsave(size(R, 1), i)  ! Debugging only
real(RP) :: tol

! Sizes
m = int(size(Q, 1), kind(m))
n = int(size(R, 2), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. n <= m, '1 <= N <= M', srname)
    call assert(i >= 1 .and. i <= n, '1 <= I <= N', srname)
    call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    call assert(size(Q, 2) >= n .and. size(Q, 2) <= m, 'N <= SIZE(Q, 2) <= M', srname)
    call assert(size(R, 1) >= n .and. size(R, 1) <= m, 'N <= SIZE(R, 1) <= M', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    Qsave = Q  ! For debugging only.
    Rsave = R(:, 1:i)  ! For debugging only.
end if

!====================!
! Calculation starts !
!====================!

if (i <= 0 .or. i >= n) then
    ! Only I == N is really needed, as 1 <= I <= N unless the input is wrong.
    return
end if

! For each K, find the Givens rotation G with G*R([K, K+1], K+1) = [HYPT, 0]. Then make two updates.
! First, update Q(:, [K, K+1]) to Q(:, [K, K+1])*G^T, and R([K, K+1], :) to G*R[K+1, K], :), which
! keeps Q*R unchanged and maintains orthogonality of Q's columns. Second, exchange columns K and K+1
! of R. Then R becomes upper triangular, and the new product Q*R exchanges columns K and K+1 of
! the original one. After this is done for each K = 1, ..., N-1, we obtain the QR factorization of
! the matrix that rearranges columns [I, I+1, ..., N] of A as [I+1, ..., N, I].
! Powell's code, however, is slightly different: before everything, he first exchanged columns K and
! K+1 of Q as well as rows K and K+1 of R. This makes sure that the diagonal entries of the updated
! R are all positive if it is the case for the original R.
do k = i, n - 1_IK
    G = planerot(R([k + 1_IK, k], k + 1))  ! G = [c, -s; s, c]
    hypt = sqrt(R(k, k + 1)**2 + R(k + 1, k + 1)**2)  ! HYPT must be calculated before R is updated
    !!hypt = G(1, 1) * R(k + 1, k + 1) + G(1, 2) * R(k, k + 1)  ! Does not perform well on 20220312
    !!hypt = hypotenuse(R(k + 1, k + 1), R(k, k + 1))  ! Does not perform well on 20220312

    ! Update Q(:, [K, K+1]).
    Q(:, [k, k + 1_IK]) = matprod(Q(:, [k + 1_IK, k]), transpose(G))

    ! Update R([K, K+1], :).
    R([k, k + 1_IK], k:n) = matprod(G, R([k + 1_IK, k], k:n))
    R(1:k + 1, [k, k + 1_IK]) = R(1:k + 1, [k + 1_IK, k])
    ! N.B.: The above two lines implement the following while noting that R is upper triangular.
    !!R([k, k + 1_IK], :) = matprod(G, R([k + 1_IK, k], :))  ! No need for R([K, K+1], 1:K-1) = 0
    !!R(:, [k, k + 1_IK]) = R(:, [k + 1_IK, k])  ! No need for R(K+2:, [K, K+1]) = 0

    ! Revise R([K, K+1], K). Changes nothing in theory but seems good for the practical performance.
    R([k, k + 1_IK], k) = [hypt, ZERO]

    !----------------------------------------------------------------------------------------------!
    ! The following code performs the update without exchanging columns K and K+1 of Q or rows K and
    ! K+1 of R beforehand. If the diagonal entries of the original R are positive, then all the
    ! updated ones become negative.
    !
    !G = planerot(R([k, k + 1_IK], k + 1))
    !hypt = sqrt(R(k, k + 1)**2 + R(k + 1, k + 1)**2)
    !!hypt = G(1, 1) * R(k, k + 1) + G(1, 2) * R(k+1, k + 1)  ! Does not perform well on 20220312
    !!hypt = hypotenuse(R(k, k + 1), R(k + 1, k + 1))  ! Does not perform well on 20220312
    !
    !Q(:, [k, k + 1_IK]) = matprod(Q(:, [k, k + 1_IK]), transpose(G))
    !
    !R([k, k + 1_IK], k:n) = matprod(G, R([k, k + 1_IK], k:n))
    !R(1:k + 1, [k, k + 1_IK]) = R(1:k + 1, [k + 1_IK, k])
    !R([k, k + 1_IK], k) = [hypt, ZERO]
    !----------------------------------------------------------------------------------------------!
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    call assert(size(Q, 2) >= n .and. size(Q, 2) <= m, 'N <= SIZE(Q, 2) <= M', srname)
    call assert(size(R, 1) >= n .and. size(R, 1) <= m, 'N <= SIZE(R, 1) <= M', srname)
    call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    !call assert(.not. any(abs(Q(:, 1:i - 1) - Qsave(:, 1:i - 1)) > 0) .and. &
    !    & .not. any(abs(Q(:, n + 1:) - Qsave(:, n + 1:)) > 0), 'Q is unchanged except Q(:, I:N)', srname)
    !call assert(.not. any(abs(R(:, 1:i - 1) - Rsave(:, 1:i - 1)) > 0), 'R(:, 1:I-1) is unchanged', srname)
    ! If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead.
    call assert(all(abs(Q(:, 1:i - 1) - Qsave(:, 1:i - 1)) <= 0) .and. &
        & all(abs(Q(:, n + 1:) - Qsave(:, n + 1:)) <= 0), 'Q is unchanged except Q(:, I:N)', srname)
    call assert(all(abs(R(:, 1:i - 1) - Rsave(:, 1:i - 1)) <= 0), 'R(:, 1:I-1) is unchanged', srname)
end if

end subroutine qrexc_Rfull


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


function calquad(d, gq, hq, pq, x, xpt) result(qred)
!--------------------------------------------------------------------------------------------------!
! This function evaluates QRED = Q(X) - Q(X + D) with Q being the quadratic function defined via
! (GQ, HQ, PQ) by
! Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
!--------------------------------------------------------------------------------------------------!

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
    call assert(size(gq) == n, 'SIZE(GQ) == N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(fval) == npt, 'SIZE(FVAL) == NPT', srname)
    call assert(.not. any(is_nan(fval) .or. is_posinf(fval)), 'FVAL is not NaN/+Inf', srname)
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
!
end function errquad


function hess_mul(hq, pq, xpt, x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates HESSIAN*X, where HESSIAN consists of an explicit part HQ and an
! implicit part PQ in Powell's way: HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
implicit none

! Inputs
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in) :: x(:)      ! X(N)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Outputs
real(RP) :: y(size(hq, 1))

! Local variables
character(len=*), parameter :: srname = 'HESSMUL'
integer(IK) :: j
integer(IK) :: n
integer(IK) :: npt

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(x) == n, 'SIZE(Y) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

!--------------------------------------------------------------------------------!
!----------! y = matprod(hq, x) + matprod(xpt, pq * matprod(x, xpt)) !-----------!
! The following loop works numerically better than the last line (but why?).
!--------------------------------------------------------------------------------!
y = matprod(xpt, pq * matprod(x, xpt))
do j = 1, n
    y = y + hq(:, j) * x(j)
end do

!====================!
!  Calculation ends  !
!====================!

end function hess_mul


function errh(idz, bmat, zmat, xpt) result(err)
!--------------------------------------------------------------------------------------------------!
! This function calculates the error in H as the inverse of W. See (3.12) of the NEWUOA paper.
! N.B.: The (NPT+1)th column (row) of H is not contained in [BMAT, ZMAT]. It is [r; t(1); s] below.
! In the complete form, the W and H in (3.12) of the NEWUOA paper are as follows.
! W = [A, ONES(NPT, 1), XPT^T; ONES(1, NPT), ZERO, ZEROS(1, N); XPT, ZEROS(N, 1), ZEROS(N, N)]
! H = [Omega, r, BMAT(:, 1:NPT)^T; r^T, t(1), s^T, BMAT(:, 1:NPT), s, BMAT(:, NPT+1:NPT+N)]
!--------------------------------------------------------------------------------------------------!
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
err = maxval(e) / (maxabs * real(n + npt, RP))

!====================!
!  Calculation ends  !
!====================!

end function errh


function omega_col(idz, zmat, k) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates Y = column K of OMEGA, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: k
real(RP), intent(in) :: zmat(:, :)

! Outputs
real(RP) :: y(size(zmat, 1))

! Local variables
character(len=*), parameter :: srname = 'OMEGA_COL'
real(RP) :: zk(size(zmat, 2))

! Preconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(k >= 1 .and. idz <= size(zmat, 1), '1 <= K <= SIZE(ZMAT, 1)', srname)
end if

!====================!
! Calculation starts !
!====================!

zk = zmat(k, :)
zk(1:idz - 1) = -zk(1:idz - 1)
y = matprod(zmat, zk)

!====================!
!  Calculation ends  !
!====================!

end function omega_col


function omega_mul(idz, zmat, x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates Y = OMEGA*X, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: zmat(:, :)
real(RP), intent(in) :: x(:)

! Outputs
real(RP) :: y(size(zmat, 1))

! Local variables
character(len=*), parameter :: srname = 'OMEGA_MUL'
real(RP) :: xz(size(zmat, 2))

! Preconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(x) == size(zmat, 1), 'SIZE(X) == SIZE(ZMAT, 1)', srname)
end if

!====================!
! Calculation starts !
!====================!

xz = matprod(x, zmat)
xz(1:idz - 1) = -xz(1:idz - 1)
y = matprod(zmat, xz)

!====================!
!  Calculation ends  !
!====================!

end function omega_mul


function omega_inprod(idz, zmat, x, y) result(p)
!--------------------------------------------------------------------------------------------------!
! This function calculates P = X^T*OMEGA*Y, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: zmat(:, :)
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)

! Outputs
real(RP) :: p

! Local variables
character(len=*), parameter :: srname = 'OMEGA_INPROD'
real(RP) :: xz(size(zmat, 2))
real(RP) :: yz(size(zmat, 2))

! Preconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(x) == size(zmat, 1), 'SIZE(X) == SIZE(ZMAT, 1)', srname)
    call assert(size(y) == size(zmat, 1), 'SIZE(Y) == SIZE(ZMAT, 1)', srname)
end if

!====================!
! Calculation starts !
!====================!

xz = matprod(x, zmat)
xz(1:idz - 1) = -xz(1:idz - 1)
yz = matprod(y, zmat)
p = inprod(xz, yz)

!====================!
!  Calculation ends  !
!====================!

end function omega_inprod


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


function norm(x, p) result(y)
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
character(len=*), parameter :: srname = 'NORM'
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
elseif (.not. any(abs(x) > ZERO)) then
    ! The following is incorrect without checking the last case, as X may be all NaN.
    y = ZERO
else
    if (is_posinf(p_loc)) then
        y = maxval(abs(x))
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
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: memory_mod, only : safealloc
implicit none
logical, intent(in) :: x(:)
integer(IK), allocatable :: loc(:)  ! INTEGER(IK) :: LOC(COUNT(X)) does not work with Absoft 22.0
integer(IK) :: i
call safealloc(loc, int(count(x), IK))  ! Removable in F03.
loc = pack([(i, i=1_IK, int(size(x), IK))], mask=x)
end function trueloc


function falseloc(x) result(loc)
!--------------------------------------------------------------------------------------------------!
! FALSELOC = TRUELOC(.NOT. X)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: memory_mod, only : safealloc
implicit none
logical, intent(in) :: x(:)
integer(IK), allocatable :: loc(:)  ! INTEGER(IK) :: LOC(COUNT(.NOT.X)) does not work with Absoft 22.0
call safealloc(loc, int(count(.not. x), IK))  ! Removable in F03.
loc = trueloc(.not. x)
end function falseloc


function minimum(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns NaN if X contains NaN; otherwise, it returns MINVAL(X).
! F2018 does not specify MINVAL(X) when X contain NaN, which motivates this function.
! Regarding NaN, the behavior of MINIMUM is the same as the following functions in various languages:
! MATLAB: min(x, [], 'includenan')
! Python: numpy.min(x)
! Julia: minimum(x)
! R: min(x)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none
real(RP), intent(in) :: x(:)
real(RP) :: y
y = merge(tsource=sum(x), fsource=minval(x), mask=any(is_nan(x)))
end function minimum


function maximum(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns NaN if X contains NaN; otherwise, it returns MAXVAL(X).
! F2018 does not specify MAXVAL(X) when X contain NaN, which motivates this function.
! Regarding NaN, the behavior of MAXIMUM is the same as the following functions in various languages:
! MATLAB: max(x, [], 'includenan')
! Python: numpy.max(x)
! Julia: maximum(x)
! R: max(x)
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none
real(RP), intent(in) :: x(:)
real(RP) :: y
y = merge(tsource=sum(x), fsource=maxval(x), mask=any(is_nan(x)))
end function maximum


end module linalg_mod
