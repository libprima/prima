module powalg_mod
!--------------------------------------------------------------------------------------------------
! This module provides some Powell-style linear algebra procedures.
!
! TODO: To avoid stack overflows, functions that return a potentially large array should declare
! the array as ALLOCATABLE rather than automatic.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Sunday, April 17, 2022 PM03:14:21
!--------------------------------------------------------------------------------------------------

implicit none
private
public :: qradd, qrexc
public :: calquad, errquad
public :: hess_mul
public :: omega_col, omega_mul, omega_inprod
public :: updateh, errh
public :: calvlag, calbeta

interface qradd
    module procedure qradd_Rdiag, qradd_Rfull
end interface

interface qrexc
    module procedure qrexc_Rdiag, qrexc_Rfull
end interface

interface calquad
    module procedure calquad_gq, calquad_gopt
end interface calquad


contains


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
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm, planerot, hypotenuse, isorth, isminor, trueloc
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
!------------------------------------------------------------!
real(RP) :: Qsave(size(Q, 1), n)  ! Debugging only
real(RP) :: Rdsave(n)  ! Debugging only
real(RP) :: tol  ! Debugging only
!------------------------------------------------------------!

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

! As in Powell's COBYLA, CQ is set to 0 at the positions with CQ being negligible as per ISMINOR.
! This may not be the best choice if the subroutine is used in other contexts, e.g., LINCOA.
cq = matprod(c, Q)
cqa = matprod(abs(c), abs(Q))
cq(trueloc(isminor(cq, cqa))) = ZERO  !!MATLAB: cq(isminor(cq, cqa)) = zero

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


subroutine qradd_Rfull(c, Q, R, n)  ! Used in LINCOA
!--------------------------------------------------------------------------------------------------!
! This subroutine updates the QR factorization of an MxN matrix A = Q*R(:, 1:N) when a new column C
! is appended to this matrix A as the LAST column.
! N.B.:
! 0. Different from QRADD_RDIAG, it seems that QRADD_RFULL does not try to maintain that A is of
! full column rank after the update; QRADD_RFULL always append C to A, and always increase N by 1,
! but QRADD_RDIAG does so only if C is not in the column space of A.
! 1. At entry, Q is a MxM orthonormal matrix, and R is a MxL upper triangular matrix with N < L <= M.
! 2. The subroutine changes only Q(:, N+1:M) and R(:, N+1) with N taking the original value.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, planerot, isorth, istriu
implicit none

! Inputs
real(RP), intent(in) :: c(:)  ! C(M)

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
!------------------------------------------------------------!
real(RP) :: Anew(size(Q, 1), n + 1)  ! Debugging only
real(RP) :: Qsave(size(Q, 1), n)  ! Debugging only
real(RP) :: Rsave(size(R, 1), n)  ! Debugging only
real(RP) :: tol  ! Debugging only
!------------------------------------------------------------!

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
    Anew = reshape([matprod(Q, R(:, 1:n)), c], shape(Anew))
    Qsave = Q(:, 1:n)  ! For debugging only.
    Rsave = R(:, 1:n)  ! For debugging only.
end if

cq = matprod(c, Q)

! Update Q so that the columns of Q(:, N+2:M) are orthogonal to C. This is done by applying a 2D
! Givens rotation to Q(:, [K, K+1]) from the right to zero C'*Q(:, K+1) out for K = N+1, ..., M-1.
! Nothing will be done if N >= M-1.
do k = m - 1_IK, n + 1_IK, -1
    if (abs(cq(k + 1)) > 0) then  ! Powell: IF (ABS(CQ(K + 1)) > 1.0D-20 * ABS(CQ(K))) THEN
        G = planerot(cq([k, k + 1_IK]))  ! G = [c, -s; s, c]. It improves the performance of LINCOA
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

    !!call assert(.not. any(abs(Q(:, 1:n - 1) - Qsave(:, 1:n - 1)) > 0), 'Q(:, 1:N-1) is unchanged', srname)
    !!call assert(.not. any(abs(R(:, 1:n - 1) - Rsave(:, 1:n - 1)) > 0), 'R(:, 1:N-1) is unchanged', srname)
    ! If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead of the last two.
    call assert(all(abs(Q(:, 1:n - 1) - Qsave(:, 1:n - 1)) <= 0), 'Q(:, 1:N-1) is unchanged', srname)
    call assert(all(abs(R(:, 1:n - 1) - Rsave(:, 1:n - 1)) <= 0), 'R(:, 1:N-1) is unchanged', srname)

    ! The following test may fail.
    call assert(all(abs(Anew - matprod(Q, R(:, 1:n))) <= max(tol, tol * maxval(abs(Anew)))), 'Anew = Q*R', srname)
end if
end subroutine qradd_Rfull


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
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm, planerot, isorth, istriu, diag
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)  ! A(M, N)

! In-outputs
real(RP), intent(inout) :: Q(:, :)  ! Q(M, :), N <= SIZE(Q, 2) <= M
real(RP), intent(inout) :: Rdiag(:)  ! Rdiag(N)
integer(IK), intent(in) :: i

! Local variables
character(len=*), parameter :: srname = 'QREXC_RDIAG'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
real(RP) :: G(2, 2)
!------------------------------------------------------------!
real(RP) :: Anew(size(A, 1), size(A, 2))  ! Debugging only
real(RP) :: Qsave(size(Q, 1), size(Q, 2))  ! Debugging only
real(RP) :: QtAnew(size(Q, 2), size(A, 2))  ! Debugging only
real(RP) :: Rdsave(i)  ! Debugging only
real(RP) :: tol  ! Debugging only
!------------------------------------------------------------!

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
!!MATLAB: Rdiag(i:n-1) = sum(Q(:, i:n-1) .* A(:, i+1:n), 1);  % Row vector
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

    Qsave(:, i:n) = Q(:, i:n)
    call assert(all(abs(Q - Qsave) <= 0), 'Q is unchanged except Q(:, I:N)', srname)
    call assert(all(abs(Rdiag(1:i - 1) - Rdsave(1:i - 1)) <= 0), 'Rdiag(1:I-1) is unchanged', srname)

    Anew = reshape([A(:, 1:i - 1), A(:, i + 1:n), A(:, i)], shape(Anew))
    QtAnew = matprod(transpose(Q), Anew)
    call assert(istriu(QtAnew, tol), 'Q^T*Anew is upper triangular', srname)
    ! The following test may fail if RDIAG is not calculated from scratch.
    call assert(norm(diag(QtAnew) - Rdiag) <= max(tol, tol * norm([(inprod(abs(Q(:, k)), &
        & abs(Anew(:, k))), k=1, n)])), 'Rdiag == diag(Q^T*Anew)', srname)
    !!MATLAB: norm(diag(QtAnew) - Rdiag) <= max(tol, tol * norm(sum(abs(Q(:, 1:n)) .* abs(Anew), 1)))
end if
end subroutine qrexc_Rdiag


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
use, non_intrinsic :: linalg_mod, only : matprod, planerot, isorth, istriu
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
!------------------------------------------------------------!
real(RP) :: Anew(size(Q, 1), size(R, 2))  ! Debugging only
real(RP) :: Qsave(size(Q, 1), size(Q, 2))  ! Debugging only
real(RP) :: Rsave(size(R, 1), i)  ! Debugging only
real(RP) :: tol  ! Debugging only
!------------------------------------------------------------!

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
    Anew = matprod(Q, R)
    Anew = reshape([Anew(:, 1:i - 1), Anew(:, i + 1:n), Anew(:, i)], shape(Anew))
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
! keeps Q*R unchanged and maintains the orthogonality of Q's columns. Second, exchange columns K and
! K+1 of R. Then R becomes upper triangular, and the new product Q*R exchanges columns K and K+1 of
! the original one. After this is done for each K = 1, ..., N-1, we obtain the QR factorization of
! the matrix that rearranges columns [I, I+1, ..., N] of A as [I+1, ..., N, I].
! Powell's code, however, is slightly different: before everything, he first exchanged columns K and
! K+1 of Q as well as rows K and K+1 of R. This makes sure that the diagonal entries of the updated
! R are all positive if it is the case for the original R.
do k = i, n - 1_IK
    G = planerot(R([k + 1_IK, k], k + 1))  ! G = [c, -s; s, c]. It improves the performance of LINCOA
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

    Qsave(:, i:n) = Q(:, i:n)
    !!call assert(.not. any(abs(Q - Qsave) > 0), 'Q is unchanged except Q(:, I:N)', srname)
    !!call assert(.not. any(abs(R(:, 1:i - 1) - Rsave(:, 1:i - 1)) > 0), 'R(:, 1:I-1) is unchanged', srname)
    ! If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead of the last two.
    call assert(all(abs(Q - Qsave) <= 0), 'Q is unchanged except Q(:, I:N)', srname)
    call assert(all(abs(R(:, 1:i - 1) - Rsave(:, 1:i - 1)) <= 0), 'R(:, 1:I-1) is unchanged', srname)

    ! The following test may fail.
    call assert(all(abs(Anew - matprod(Q, R)) <= max(tol, tol * maxval(abs(Anew)))), 'Anew = Q*R', srname)
end if

end subroutine qrexc_Rfull


function calquad_gq(d, gq, hq, pq, x, xpt) result(qred)
!--------------------------------------------------------------------------------------------------!
! This function evaluates QRED = Q(XBASE + X) - Q(XBSE + X + D) with Q being the quadratic function
! defined via (GQ, HQ, PQ) by
! Q(XBASE + S) = <S, GQ> + 0.5*<S, HESSIAN*S>,
! with HESSIAN consisting of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
! N.B.: GQ is the gradient of Q at XBASE; the value of XBASE is irrelevant in the calculation.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, issymmetric
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
character(len=*), parameter :: srname = 'CALQUAD_GQ'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: n
integer(IK) :: npt
real(RP) :: s(size(x))
real(RP) :: t
real(RP) :: w(size(pq))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
end if

!====================!
! Calculation starts !
!====================!

s = x + d

! First-order term and explicit second order term
qred = ZERO
do j = 1, n
    qred = qred - d(j) * gq(j)
    do i = 1, j
        t = d(i) * s(j) + d(j) * x(i)
        if (i == j) then
            t = HALF * t
        end if
        qred = qred - t * hq(i, j)
    end do
end do

! Implicit second-order term
w = matprod(d, xpt)
w = w * (HALF * w + matprod(x, xpt))
do i = 1, npt
    qred = qred - pq(i) * w(i)
end do

!--------------------------------------------------------------------------------------------------!
! The following is a loop-free implementation, which should be applied in MATLAB/Python/R/Julia.
!--------------------------------------------------------------------------------------------------!
!! The order of calculation seems quite important. The following order seems to work well.
!! First-order term
!qred = -inprod(d, gq)
!s = HALF * d + x  ! Different from the above version.
!! Implicit second-order term
!qred = qred - sum(pq * (matprod(s, xpt) * matprod(d, xpt)))
!! Explicit second-order term
!qred = qred - inprod(s, matprod(hq, d))
!! In Fortran, the following implementations do not work as well as the above line.
!!qred = qred - inprod(d, matprod(hq, s))
!!qred = qred - HALF*(inprod(d, matprod(hq, s)) + inprod(s, matprod(hq, d)))
!--------------------------------------------------------------------------------------------------!

!====================!
!  Calculation ends  !
!====================!

end function calquad_gq


function calquad_gopt(d, gopt, hq, pq, xpt) result(qred)
!--------------------------------------------------------------------------------------------------!
! This function evaluates QRED = Q(XOPT) - Q(XOPT+D) = -<D, GOPT> - 0.5*<D, HESSIAN*D>,
! with HESSIAN consisting of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, issymmetric
implicit none

! Inputs
real(RP), intent(in) :: d(:)      ! D(N)
real(RP), intent(in) :: gopt(:)   ! GOPT(N)
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Output
real(RP) :: qred

! Local variable
character(len=*), parameter :: srname = 'CALQUAD_GOPT'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: w(size(pq))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
end if

!====================!
! Calculation starts !
!====================!

!!--------------------------------------------------------------------------------------------------!
!! The following is Powell's scheme in LINCOA.
!! First-order term and explicit second-order term
!qred = ZERO
!do j = 1, n
!    qred = qred - d(j) * gopt(j)
!    do i = 1, j
!        t = d(i) * d(j)
!        if (i == j) then
!            t = HALF * t
!        end if
!        qred = qred - t * hq(i, j)
!    end do
!end do

!! Implicit second-order term
!w = matprod(d, xpt)
!do i = 1, npt
!    qred = qred - HALF * pq(i) * w(i) * w(i)  ! In BOBYQA, it is QRED - HALF * PQ(I) * W(I)**2.
!end do
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! The following is a loop-free implementation, which should be applied in MATLAB/Python/R/Julia.
!--------------------------------------------------------------------------------------------------!
w = matprod(d, xpt)
qred = -inprod(d, gopt + HALF * matprod(hq, d)) - HALF * inprod(w, pq * w)
!!MATLAB: qred = -d'*(gopt + 0.5*hq*d) - 0.5*w'*(pq*w);
!--------------------------------------------------------------------------------------------------!

!====================!
!  Calculation ends  !
!====================!

end function calquad_gopt


function errquad(gq, hq, pq, xpt, fval) result(err)
!--------------------------------------------------------------------------------------------------!
! This function calculates the maximal relative error of Q in interpolating FVAL on XPT.
! Here, Q is the quadratic function defined via (GQ, HQ, PQ) by
! Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
! with HESSIAN consisting of an explicit part HQ and an implicit part PQ in Powell's way:
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: linalg_mod, only : issymmetric
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
!!MATLAB: qval = cellfun(@(x) -calquad(x, gq, hq, pq, zeros, xpt), num2cell(xpt, 1));  % Row vector
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


function hess_mul(x, xpt, pq, hq) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates HESSIAN*X, with HESSIAN consisting of an explicit part HQ (0 if absent)
! and an implicit part PQ in Powell's way: HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, issymmetric
implicit none

! Inputs
real(RP), intent(in) :: x(:)      ! X(N)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in), optional :: hq(:, :)  ! HQ(N, N)

! Outputs
real(RP) :: y(size(x))

! Local variables
character(len=*), parameter :: srname = 'HESS_MUL'
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
    call assert(size(x) == n, 'SIZE(Y) == N', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    if (present(hq)) then
        call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

!--------------------------------------------------------------------------------!
!----------! y = matprod(hq, x) + matprod(xpt, pq * matprod(x, xpt)) !-----------!
! The following loop works numerically better than the last line (but why?).
!--------------------------------------------------------------------------------!
y = matprod(xpt, pq * matprod(x, xpt))
if (present(hq)) then
    do j = 1, n
        y = y + hq(:, j) * x(j)
    end do
end if

!====================!
!  Calculation ends  !
!====================!

end function hess_mul


function omega_col(idz, zmat, k) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates Y = column K of OMEGA. As Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
! OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod
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
! This function calculates Y = OMEGA*X. As Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
! OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod
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
! This function calculates P = X^T*OMEGA*Y. As Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
! OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, inprod
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


function errh(idz, bmat, zmat, xpt) result(err)
!--------------------------------------------------------------------------------------------------!
! This function calculates the error in H as the inverse of W. See (3.12) of the NEWUOA paper.
! N.B.: The (NPT+1)th column (row) of H is not contained in [BMAT, ZMAT]. It is [r; t(1); s] below.
! In the complete form, using MATLAB-style notation, the W and H in the NEWUOA paper are as follows.
! W = [A, ONES(NPT, 1), XPT^T; ONES(1, NPT), ZERO, ZEROS(1, N); XPT, ZEROS(N, 1), ZEROS(N, N)]
! H = [Omega, r, BMAT(:, 1:NPT)^T; r^T, t(1), s^T, BMAT(:, 1:NPT), s, BMAT(:, NPT+1:NPT+N)]
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, eye, issymmetric
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


subroutine updateh(knew, kopt, idz, d, xpt, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine updates arrays BMAT and ZMAT together with IDZ, in order to replace the
! interpolation point XPT(:, KNEW) by XNEW = XPT(:, KOPT) + D. See Section 4 of the NEWUOA paper.
! [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
! the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
! ZMAT holds a factorization of the leading NPT*NPT submatrix OMEGA of H, the factorization being
! OMEGA = ZMAT*Diag(S)*ZMAT^T with S(1:IDZ-1)= -1 and S(IDZ : NPT-N-1) = +1; BMAT holds the last N
! ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th column of
! H are not stored as they are unnecessary for the calculation.
!
! N.B.: The most natural (not necessarily the best numerically) version of UPDATEH should work based
! on [KNEW, XNEW - XPT(:,KNEW)] instead of [KNEW, KOPT, D], because the update is independent of KOPT
! in theory. UPDATEH needs KOPT only for calculating VLAG and BETA, which can also be found by
!!VLAG = CALVLAG(KNEW, BMAT, XNEW - XPT(:, KNEW), XPT, ZMAT, IDZ)
!!BETA = CALBETA(KNEW, BMAT, XNEW - XPT(:, KNEW), XPT, ZMAT, IDZ)
! Theoretically (but not numerically), they should return the same VLAG and BETA as the calls below.
! However, as observed on 20220412, such an implementation can lead to significant errors in H!
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: GROT(2, 2), V1(N), V2(N), VLAG(NPT+N), W(NPT+N)
! Size of local arrays: REAL(RP)*(4+4*N+2*NPT)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: info_mod, only : INFO_DFT, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : matprod, planerot, symmetrize, issymmetric, outprod!, r2update
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: d(:)  ! D(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
integer(IK), intent(inout) :: idz
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
integer(IK), intent(out), optional :: info

! Local variables
character(len=*), parameter :: srname = 'UPDATEH'
integer(IK) :: j
integer(IK) :: ja
integer(IK) :: jb
integer(IK) :: jl
integer(IK) :: n
integer(IK) :: npt
real(RP) :: alpha
real(RP) :: beta
real(RP) :: denom
real(RP) :: grot(2, 2)
real(RP) :: scala
real(RP) :: scalb
real(RP) :: sqrtdn
real(RP) :: tau
real(RP) :: temp
real(RP) :: tempa
real(RP) :: tempb
real(RP) :: v1(size(bmat, 1))
real(RP) :: v2(size(bmat, 1))
real(RP) :: vlag(size(bmat, 2))
real(RP) :: w(size(bmat, 2))
real(RP) :: ztest

! Debugging variables
real(RP) :: beta_test
real(RP) :: tol
real(RP), allocatable :: vlag_test(:)
!real(RP), allocatable :: xpt_test(:, :)

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)

    ! Theoretically, CALVLAG and CALBETA should be independent of the reference point XPT(:, KOPT).
    ! So we test the following. By the implementation of CALVLAG and CALBETA, we are indeed testing
    ! H*[w(X_KOPT) - w(X_KNEW)] = e_KOPT - e_KNEW. Thus H = W^{-1} is also tested to some extend.
    if (knew >= 1) then
        tol = 1.0E-2_RP  ! W and H are quite ill-conditioned, so we do not test a high precision.
        call safealloc(vlag_test, npt + n)
        vlag_test = calvlag(knew, bmat, d + (xpt(:, kopt) - xpt(:, knew)), xpt, zmat, idz)
        call wassert(maxval(abs(vlag_test - calvlag(kopt, bmat, d, xpt, zmat, idz))) <= &
            & tol * maxval([ONE, abs(vlag_test)]) .or. RP == kind(0.0), 'VLAG_TEST == VLAG', srname)
        deallocate (vlag_test)
        beta_test = calbeta(knew, bmat, d + (xpt(:, kopt) - xpt(:, knew)), xpt, zmat, idz)
        call wassert(abs(beta_test - calbeta(kopt, bmat, d, xpt, zmat, idz)) <= &
            & tol * max(ONE, abs(beta_test)) .or. RP == kind(0.0), 'BETA_TEST == BETA', srname)
    end if

    ! The following is too expensive to check.
    !call wassert(errh(idz, bmat, zmat, xpt) <= tol .or. RP == kind(0.0), &
    !    & 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
end if

!====================!
! Calculation starts !
!====================!

if (present(info)) then
    info = INFO_DFT
end if

! We must not do anything if KNEW is 0. This can only happen sometimes after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    return
end if

! Calculate VLAG and BETA according to D.
! VLAG contains the components of the vector H*w of the updating formula (4.11) in the NEWUOA paper,
! and BETA holds the value of the parameter that has this name.
! N.B.: Powell's original comments mention that VLAG is "the vector THETA*WCHECK + e_b of the
! updating formula (6.11)", which does not match the published version of the NEWUOA paper.
vlag = calvlag(kopt, bmat, d, xpt, zmat, idz)
beta = calbeta(kopt, bmat, d, xpt, zmat, idz)

! Apply Givens rotations to put zeros in the KNEW-th row of ZMAT and set JL. After this,
! ZMAT(KNEW, :) contains at most two nonzero entries ZMAT(KNEW, 1) and ZMAT(KNEW, JL), one
! corresponding to all the columns of ZMAT that has a coefficient -1 in the factorization of
! OMEGA (if any), and the other corresponding to all the columns with +1. In specific,
! 1. If IDZ = 1 (all coefficients are +1 for the columns of ZMAT in the factorization of OMEGA ) or
! NPT - N (all the coefficients are -1), then JL = 1, and ZMAT(KNEW, 1) is L2-norm of ZMAT(KNEW, :);
! 2. If 2 <= IDZ <= NPT - N -1, then JL = IDZ, and ZMAT(KNEW, 1) is L2-norm of ZMAT(KNEW, 1 : IDZ-1),
! while ZMAT(KNEW, JL) is L2 norm of ZMAT(KNEW, IDZ : NPT-N-1).
! See (4.15)--(4.17) of the NEWUOA paper and the elaboration around them.
ztest = 1.0E-20_RP * maxval(abs(zmat))  ! Taken from BOBYQA. It is implicitly zero in NEWUOA/LINCOA.
jl = 1_IK  ! In the loop below, if 2 <= J < IDZ, then JL = 1; if IDZ < J <= NPT-N-1, then JL = IDZ.
do j = 2_IK, npt - n - 1_IK
    if (j == idz) then
        jl = idz  ! Do nothing but changing JL from 1 to IDZ. It occurs at most once along the loop.
        cycle
    end if

    ! Powell's condition in NEWUOA/LINCOA for the IF ... THEN below: IF (ZMAT(KNEW, J) /= 0) THEN
    ! A possible alternative: IF (ABS(ZMAT(KNEW, J)) > 1.0E-20_RP * ABS(ZMAT(KNEW, JL))) THEN
    if (abs(zmat(knew, j)) > ztest) then
        ! Multiply a Givens rotation to ZMAT from the right so that ZMAT(KNEW, [JL,J]) becomes [*,0].
        grot = planerot(zmat(knew, [jl, j]))  !!MATLAB: grot = planerot(zmat(knew, [jl, j])')
        zmat(:, [jl, j]) = matprod(zmat(:, [jl, j]), transpose(grot))
    end if
    zmat(knew, j) = ZERO
end do

! Set W(1:NPT) to the first NPT components of the KNEW-th column of H.
! Here, Powell's code calculates this column by the ZMAT that has been updated by the rotation as
! above. Theoretically, W(1:NPT) can also be calculated before the rotation by calling OMEGA_COL,
! which is tempting to do because the rotation may introduce some rounding errors to ZMAT. However,
! according to a test on 20220411, this alternative does not improve the performance. We can also
! calculate VLAG and BETA using the updated ZMAT here, which does not lead to improvements either.
if (idz <= 1) then  ! IDZ <= 0 is impossible unless there is a bug.
    !tempa = zmat(knew,1)  ! Not needed anymore. Retained for the comments below.
    w(1:npt) = zmat(knew, 1) * zmat(:, 1)
else
    !tempa = -zmat(knew,1)  ! Not needed anymore. Retained for the comments below.
    w(1:npt) = -zmat(knew, 1) * zmat(:, 1)
end if
if (jl > 1) then  ! In this case, 1 < JL == IDZ < NPT - N.
    w(1:npt) = w(1:npt) + zmat(knew, jl) * zmat(:, jl)
end if

! Calculate the parameters of the updating formula (4.18)--(4.20) in the NEWUOA paper.
alpha = w(knew)
tau = vlag(knew)
denom = alpha * beta + tau**2
if (abs(denom) <= 0 .or. is_nan(denom)) then
    ! 1. Up to here, only ZMAT is rotated, which does not change H in precise arithmetic, so it is
    ! fine if we do not revert ZMAT to the un-updated version.
    ! 2. After UPDATEH returns, ideally, the algorithm should do something to rectify the damaging
    ! rounding. However, nothing is done in the current (20220412) version of NEWUOA/LINCOA, while
    ! Powell's version of LINCOA terminates immediately. Note that H is not updated at all here.
    if (present(info)) then
        info = DAMAGING_ROUNDING
    end if
    return
end if
sqrtdn = sqrt(abs(denom))
! After the following line, VLAG = H*w - e_t in the NEWUOA paper.
vlag(knew) = vlag(knew) - ONE

if (jl == 1) then
    ! Complete the updating of ZMAT when there is only 1 nonzero in ZMAT(KNEW, :) after the rotation.
    ! This is the normal case, as IDZ = 1 in precise arithmetic; it also covers the rare case that
    ! IDZ = NPT-N, meaning that OMEGA = -ZMAT*ZMAT^T. See (4.18) of the NEWUOA paper for details.
    ! Note that (4.18) updates Z_{NPT-N-1}, but the code here updates ZMAT(:, 1). Correspondingly,
    ! we implicitly update S_1 to SIGN(DENOM)*S_1 according to (4.18). If IDZ = NPT-N before the
    ! update, then IDZ is reduced by 1, and we need to switch ZMAT(:, 1) and ZMAT(:, IDZ) to maintain
    ! that S_J = -1 iff 1 <= J < IDZ, which is done after the END IF together with another case.

    !---------------------------------------------------------------------------------------------!
    ! Up to now, TEMPA = ZMAT(KNEW, 1) if IDZ = 1 and TEMPA = -ZMAT(KNEW, 1) if IDZ >= 2. However,
    ! according to (4.18) of the NEWUOA paper, TEMPB should always be ZMAT(KNEW, 1)/SQRTDN
    ! regardless of IDZ. Therefore, the following definition of TEMPB is inconsistent with (4.18).
    ! This is probably a BUG. See also Lemma 4 and (5.13) of Powell's paper "On updating the inverse
    ! of a KKT matrix". However, the inconsistency is hardly observable in practice, because JL = 1
    ! implies IDZ = 1 in precise arithmetic.
    !--------------------------------------------!
    !!tempb = tempa/sqrtdn
    !!tempa = tau/sqrtdn
    !--------------------------------------------!
    ! Here is the corrected version (only TEMPB is changed).
    tempa = tau / sqrtdn
    tempb = zmat(knew, 1) / sqrtdn
    !---------------------------------------------------------------------------------------------!

    ! The following line updates ZMAT(:, 1) according to (4.18) of the NEWUOA paper.
    zmat(:, 1) = tempa * zmat(:, 1) - tempb * vlag(1:npt)

    !---------------------------------------------------------------------------------------------!
    ! Zaikun 20220411: The update of IDZ is decoupled from the update of ZMAT, located after END IF.
    !---------------------------------------------------------------------------------------------!
    ! The following six lines from Powell's NEWUOA code are obviously problematic --- SQRTDN is
    ! always nonnegative. According to (4.18) of the NEWUOA paper, "SQRTDN < 0" and "SQRTDN >= 0"
    ! below should be both revised to "DENOM < 0". See also the corresponding part of the LINCOA
    ! code. Note that the NEWUOA paper uses SIGMA to denote DENOM. Check also Lemma 4 and (5.13) of
    ! Powell's paper "On updating the inverse of a KKT matrix". Note that the BOBYQA code does not
    ! have this part, as it does not have IDZ at all.
    !--------------------------------------------!
    !!if (idz == 1 .and. sqrtdn < 0) then
    !!    idz = 2
    !!end if
    !!if (idz >= 2 .and. sqrtdn >= 0) then
    !!    reduce_idz = .true.
    !!end if
    !--------------------------------------------!
    !! This is the corrected version, copied from LINCOA.
    !if (denom < 0) then
    !    if (idz == 1) then
    !        idz = 2_IK
    !    else
    !        reduce_idz = .true.
    !    end if
    !end if
    !---------------------------------------------------------------------------------------------!
else
    ! Complete the updating of ZMAT in the alternative case: ZMAT(KNEW, :) has 2 nonzeros. See (4.19)
    ! and (4.20) of the NEWUOA paper.
    ! First, set JA and JB so that ZMAT(: [JA, JB]) corresponds to [Z_1, Z_2] in (4.19) when BETA>=0,
    ! and corresponds to [Z2, Z1] in (4.20) when BETA<0. In this way, the update of ZMAT(:, [JA, JB])
    ! follows the same scheme regardless of BETA. Indeed, since S_1 = 1 and S_2 = -1 in (4.19)-(4.20)
    ! as elaborated above the equations, ZMAT(:, [1, JL]) always correspond to [Z_2, Z_1].
    if (beta >= 0) then  ! ZMAT(:, [JA, JB]) corresponds to [Z_1, Z_2] in (4.19)
        ja = jl
        jb = 1_IK
    else  ! ZMAT(:, [JA, JB]) corresponds to [Z_2, Z_1] in (4.20)
        ja = 1_IK
        jb = jl
    end if
    ! Now update ZMAT(:, [ja, jb]) according to (4.19)--(4.20) of the NEWUOA paper.
    temp = zmat(knew, jb) / denom
    tempa = temp * beta
    tempb = temp * tau
    temp = zmat(knew, ja)
    scala = ONE / sqrt(abs(beta) * temp**2 + tau**2)  ! 1/SQRT(ZETA) in (4.19)-(4.20) of NEWUOA paper
    scalb = scala * sqrtdn
    zmat(:, ja) = scala * (tau * zmat(:, ja) - temp * vlag(1:npt))
    zmat(:, jb) = scalb * (zmat(:, jb) - tempa * w(1:npt) - tempb * vlag(1:npt))

    !----------------------------------------------------------------------------------------------!
    ! Zaikun 20220411: The update of IDZ is decoupled from the update of ZMAT, located after END IF.
    !----------------------------------------------------------------------------------------------!
    !! If and only if DENOM < 0, IDZ will be revised according to the sign of BETA.
    !! See (4.19)--(4.20) of the NEWUOA paper.
    !if (denom < 0) then
    !    if (beta < 0) then
    !        idz = int(idz + 1, kind(idz))
    !    else
    !        reduce_idz = .true.
    !    end if
    !end if
    !----------------------------------------------------------------------------------------------!
end if

!--------------------------------------------------------------------------------------------------!
! Zaikun 20220411: The update of IDZ is decoupled from the update of ZMAT, located right below.
!--------------------------------------------------------------------------------------------------!
!! IDZ is reduced in the following case. Then exchange ZMAT(:, 1) and ZMAT(:, IDZ).
!if (reduce_idz) then
!    idz = int(idz - 1, kind(idz))
!    if (idz > 1) then
!        zmat(:, [1_IK, idz]) = zmat(:, [idz, 1_IK])
!    end if
!end if
!--------------------------------------------------------------------------------------------------!

! According to (4.18) and (4.19)--(4.20) of the NEWUOA paper, the coefficients {S_J} need update iff
! DENOM < 0, in which case one of the S_J will flip the sign when multiplied by SIGN(DENOM), leading
! to an increase of IDZ (if S_J flipped from 1 to -1) or a decrease (if S_J flipped from -1 to 1).
if (denom < 0) then
    if (idz == 1 .or. (idz < npt - n .and. beta < 0)) then  ! (4.18), (4.20) of the NEWUOA paper
        idz = idz + 1_IK
    elseif (idz == npt - n .or. (idz > 1 .and. beta >= 0)) then  ! (4.18), (4.19) of the NEWUOA paper
        idz = idz - 1_IK
        ! Exchange ZMAT(:, 1) and ZMAT(:. IDZ) if IDZ > 1. Why? No matter whether the update is
        ! given by (4.18) (IDZ = NPT-N) or (4.19) (1 < IDZ < NPT-N and BETA >= 0), we have S_1 = +1
        ! and S_{IDZ} = -1 at this moment (unless IDZ = 1). Thus we need to exchange ZMAT(:, 1) with
        ! ZMAT(:, IDZ) and implicitly S_1 with S_IDZ to maintain that S_J = -1 iff 1 <= J < IDZ.
        ! Note that, in the case of (4.18), ZMAT(:, 1) (and implicitly S_1) rather than
        ! ZMAT(:, NPT-N-1) was updated by the code above.
        if (idz > 1) then
            zmat(:, [1_IK, idz]) = zmat(:, [idz, 1_IK])
        end if
    end if
end if

! Finally, update the matrix BMAT. It implements the last N rows of (4.11) in the NEWUOA paper.
w(npt + 1:npt + n) = bmat(:, knew)
v1 = (alpha * vlag(npt + 1:npt + n) - tau * w(npt + 1:npt + n)) / denom
v2 = (-beta * w(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) / denom
bmat = bmat + outprod(v1, vlag) + outprod(v2, w) !!call r2update(bmat, ONE, v1, vlag, ONE, v2, w)
! Numerically, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be symmetric.
call symmetrize(bmat(:, npt + 1:npt + n))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    ! The following is too expensive to check.
    !call safealloc(xpt_test, n, npt)
    !xpt_test = xpt
    !xpt_test(:, knew) = xpt(:, kopt) + d
    !call wassert(errh(idz, bmat, zmat, xpt_test) <= tol .or. RP == kind(0.0), &
    !    & 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
    !deallocate (xpt_test)
end if
end subroutine updateh


!--------------------------------------------------------------------------------------------------!
! CALVLAG and CALBETA are subroutine that calculate VLAG and BETA for a given step D. Both VLAG and
! BETA are critical for the updating procedure of H, which is detailed formula (4.11) of the NEWUOA
! paper. See (4.12) for the definition of BETA, and VLAG is indeed H*w without the (NPT+1)the entry;
! (4.25)--(4.26) formulate the actual calculating scheme of VLAG and BETA.
! In languages like MATLAB/Python/Julia/R, CALVLAG and CALBETA should be implemented into one
! single function, as they share most of the calculation. We separate them in Fortran (at the
! expense of repeating some calculation) because Fortran functions can only have one output.
!
! Explanation on the matrix H in (3.12) and w(X) in (6.3) of the NEWUOA paper and WCHECK in the code:
! 0. As defined in (6.3) of the paper w(X)(K) = 0.5*[(X-X_0)^T*(X_K-X_0)]^2 for K = 1, ..., NPT,
! w(X)(NPT+1) = 1, and w(X)(NPT+2:NPT+N+1) = X - X_0. As in (3.12) of the paper, w(X_K) is the K-th
! column in the coefficient matrix W of the KKT system for the interpolation problem
! Minimize ||nabla^2 Q||_F s.t. Q(X_K) = Y_K, K = 1, ..., NPT.
! This is why w(X) is ubiquitous in the code.
! 1. By (3.9) of the paper, the solution to the above interpolation problem is Q(X) = Y^T*H*w(X),
! with H = W^{-1}, and Y being the vector [Y_1; ...; Y_NPT; 0; ...; 0] with N trailing zeros. In
! particular, the K-th Lagrange function of this interpolation problem is e_K^T*H*w(X), namely the
! K-th entry of the vector H*w(X). This is why H*w(X) appears as the vector VLAG in the code.
! As a consequence, SUM(VLAG(1:NPT)) = 1 in theory.
! 2. As above, H can provide us interpolants and the Lagrange functions. Thus the code maintains H.
! 3. Since the matrix H is W^{-1} as defined in (3.12) of the paper, we have H*w(X_K) = e_K for
! any K in {1, ..., NPT}.
! 4. When the interpolation set is updated by replacing X_K with X, W is correspondingly updated by
! changing the K-th column from w(X_t) to w(X). This is why the update of H = W^{-1} must involve
! H*[w(X) - w(X_t)] = H*w(X) - e_t.
! 5. As explained above, the vector H*w(X) is essential to the algorithm. The quantity w(X)^T*H*w(X)
! is also needed in the update of H (particularly by BETA). However, they can be tricky to calculate,
! because much cancellation can happen when X_0 is far away from the interpolation set, as explained
! in (7.8)--(7.10) of the paper and the discussions around. To overcome the difficulty, note that
! H*w(X) = H*[w(X) - w(X_K)] + H*w(X_K) = H*[w(X) - w(X_K)] + e_K, and
! w(x)^T*H*w(X) = [w(X) - w(X_K)]^T*H*[w(X) - w(X_K)] + 2*w(X)(K) - w(X_K)(K).
! The dependence of w(X)-w(X_K) on X_0 is weaker, which reduces (but does not resolve) the difficulty.
! In theory, these formulas are invariant with respect to K. In the code, this means that
! CALVLAG(K, BMAT, X - XPT(:, K), ZMAT, IDZ) and CALBETA(K, BMAT, X - XPT(:, K), ZMAT, IDZ)
! are invariant with respect to K in {1, ..., NPT}. Powell's code uses always K = KOPT.
! 6. Since the (NPT+1)-th entry of w(X) - w(X_K) is 0, the above formulas do not require the
! (NPT+1)-th column of H, which is not stored in the code.
! 7. WCHECK contains the first NPT entries of w-v for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, with w = w(X) and v = w(X_KOPT); it is also hat{w} in (6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
! 8. Assume that the |D| ~ DELTA, |XPT| ~ |XOPT|, and DELTA < |XOPT|. Then WCHECK is of the order
! DELTA*|XOPT|^3, which is can be huge at the beginning of the algorithm and quickly become tiny.
!--------------------------------------------------------------------------------------------------!

function calvlag(kopt, bmat, d, xpt, zmat, idz) result(vlag)
!--------------------------------------------------------------------------------------------------!
! This function calculates VLAG = H*w for a given step D. See (4.25) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: VLAG(NPT+N), WCHECK(NPT), XOPT(N)
! Size of local arrays: REAL(RP)*(2*NPT+2*N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, issymmetric

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
integer(IK), intent(in), optional :: idz  ! Absent in BOBYQA, being equivalent to IDZ = 1

! Outputs
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))    ! VLAG(NPT + N)

! Local variables
character(len=*), parameter :: srname = 'CALVLAG'
integer(IK) :: idz_loc
integer(IK) :: n
integer(IK) :: npt
real(RP) :: tol  ! For debugging only
real(RP) :: wcheck(size(zmat, 1))
real(RP) :: xopt(size(xpt, 1))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Read IDZ, which is not present in BOBYQA, being equivalent to IDZ = 1.
idz_loc = 1_IK
if (present(idz)) then
    idz_loc = idz
end if

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz_loc >= 1 .and. idz_loc <= size(zmat, 2) + 1, '1 <= ID <= SIZE(ZMAT, 2) + 1', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)

    !--------------------------------------------------------------------------------------!
    ! Disable the test for the moment, as it cannot pass in BOBYQA.!!!!!!!!!!!!!!!!!!!!!!!!!
    call wassert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    !--------------------------------------------------------------------------------------!

    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

! Set WCHECK to the first NPT entries of (w-v) for w and v in (4.10) and (4.24) of the NEWUOA paper.
wcheck = matprod(d, xpt)
wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

! The following two lines set VLAG to H*(w-v).
vlag(1:npt) = omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt))
vlag(npt + 1:npt + n) = matprod(bmat, [wcheck, d])
! The following line is equivalent to the above one, but handles WCHECK and D separately.
!!vlag(npt + 1:npt + n) = matprod(bmat(:, 1:npt), wcheck) + matprod(bmat(:, npt + 1:npt + n), d)

! The following line sets VLAG(KOPT) to the correct value.
vlag(kopt) = vlag(kopt) + ONE

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(vlag) == npt + n, 'SIZE(VLAG) == NPT + N', srname)
    tol = max(1.0E-8_RP, min(1.0E-1_RP, 1.0E10_RP * EPS * real(npt + n, RP)))
    call wassert(abs(sum(vlag(1:npt)) - ONE) / real(npt, RP) <= tol .or. RP == kind(0.0), &
        & 'SUM(VLAG(1:NPT)) == 1', srname)
end if

end function calvlag


function calbeta(kopt, bmat, d, xpt, zmat, idz) result(beta)
!--------------------------------------------------------------------------------------------------!
! This function calculates BETA for a given step D. See (4.12) and (4.26) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: BW(N), BD(N), WCHECK(NPT), XOPT(N)
! Size of local arrays: REAL(RP)*(3*NPT+4*N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod, issymmetric

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
integer(IK), intent(in), optional :: idz  ! Absent in BOBYQA, being equivalent to IDZ = 1

! Outputs
real(RP) :: beta

! Local variables
character(len=*), parameter :: srname = 'CALBETA'
integer(IK) :: idz_loc
integer(IK) :: n
integer(IK) :: npt
real(RP) :: dsq
real(RP) :: dvlag
real(RP) :: dxopt
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))
real(RP) :: wcheck(size(zmat, 1))
real(RP) :: wmv(size(xpt, 1) + size(xpt, 2))
real(RP) :: wvlag
real(RP) :: xopt(size(xpt, 1))
real(RP) :: xoptsq

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Read IDZ, which is absent from BOBYQA, being equivalent to IDZ = 1.
idz_loc = 1_IK
if (present(idz)) then
    idz_loc = idz
end if

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz_loc >= 1 .and. idz_loc <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)

    !--------------------------------------------------------------------------------------!
    ! Disable the test for the moment, as it cannot pass in BOBYQA.!!!!!!!!!!!!!!!!!!!!!!!!!
    call wassert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    !--------------------------------------------------------------------------------------!

    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

! Set WCHECK to the first NPT entries of (w-v) for w and v in (4.10) and (4.24) of the NEWUOA paper.
wcheck = matprod(d, xpt)
wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

! WMV is the vector (w-v) for w and v in (4.10) and (4.24) of the NEWUOA paper.
wmv = [wcheck, d]
! The following two lines set VLAG to H*(w-v).
vlag(1:npt) = omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt))
vlag(npt + 1:npt + n) = matprod(bmat, wmv)
! The following line is equivalent to the above one, but handles WCHECK and D separately.
!!vlag(npt + 1:npt + n) = matprod(bmat(:, 1:npt), wcheck) + matprod(bmat(:, npt + 1:npt + n), d)

! BETA = HALF*|XOPT + D|^4 - (W-V)'*H*(W-V) - [XOPT'*(X+XOPT)]^2 + HALF*|XOPT|^4. See equations
! (4.12) and (4.26) of the NEWUOA paper.
dxopt = inprod(d, xopt)
dsq = inprod(d, d)
xoptsq = inprod(xopt, xopt)
dvlag = inprod(d, vlag(npt + 1:npt + n))
wvlag = inprod(wcheck, vlag(1:npt))
beta = dxopt**2 + dsq * (xoptsq + dxopt + dxopt + HALF * dsq) - dvlag - wvlag
!---------------------------------------------------------------------------------------------------!
! The last line is equivalent to either of the following lines, but performs better numerically.
!!beta = dxopt**2 + dsq * (xoptsq + dxopt + dxopt + HALF * dsq) - inprod(vlag, wmv) ! Not good
!!beta = dxopt**2 + dsq * (xoptsq + dxopt + dxopt + HALF * dsq) - wvlag - dvlag  ! Bad
!---------------------------------------------------------------------------------------------------!

! N.B.:
! 1. Mathematically, the following two quantities are equal:
! DXOPT**2 + DSQ * (XOPTSQ + DXOPT + DXOPT + HALF * DSQ) ,
! HALF * (INPROD(X, X)**2 + INPROD(XOPT, XOPT)**2) - INPROD(X, XOPT)**2 with X = XOPT + D.
! However, the first (by Powell) is a better numerical scheme. According to the first formulation,
! this quantity is in the oder of |D|^2*|XOPT|^2 if |XOPT| >> |D|, which is normally the case.
! However, each term in the second formulation has an order of |XOPT|^4. Thus much cancellation will
! occur in the second formulation. In addition, the first formulation contracts the rounding error
! in (XOPTSQ + DXOPT + DXOPT + HALF * DSQ) by a factor of |D|^2, which is typically small.
! 2. We can evaluate INPROD(VLAG, WMV) as INPROD(VLAG(1:NPT), WCHECK) + INPROD(VLAG(NPT+1:NPT+N),D)
! if it is desirable to handle WCHECK and D separately due to their significantly different magnitudes.

! The following line sets VLAG(KOPT) to the correct value if we intend to output VLAG.
!!vlag(kopt) = vlag(kopt) + ONE

!====================!
!  Calculation ends  !
!====================!

end function calbeta


end module powalg_mod
