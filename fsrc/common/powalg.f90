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
! Last Modified: Friday, April 08, 2022 AM11:07:41
!--------------------------------------------------------------------------------------------------

implicit none
private
public :: calquad, errquad
public :: hess_mul
public :: errh
public :: omega_col, omega_mul, omega_inprod
public :: qradd, qrexc

interface calquad
    module procedure calquad_gq, calquad_gopt
end interface calquad

interface qradd
    module procedure qradd_Rdiag, qradd_Rfull
end interface

interface qrexc
    module procedure qrexc_Rdiag, qrexc_Rfull
end interface


contains


function calquad_gq(d, gq, hq, pq, x, xpt) result(qred)
!--------------------------------------------------------------------------------------------------!
! This function evaluates QRED = Q(XBASE + X) - Q(XBSE + X + D) with Q being the quadratic function
! defined via (GQ, HQ, PQ) by
! Q(XBASE + S) = <S, GQ> + 0.5*<S, HESSIAN*S>,
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
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
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
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
! N.B.: INPROD(W, PQ*W) is mathematically equal to INPROD(PQ, W**2). However, the latter may be a
! bad formulation for computation, because the entires of W are in the tiny order of O(DELTA^2).
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
! where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
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


function hess_mul(hq, pq, xpt, x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates HESSIAN*X, where HESSIAN consists of an explicit part HQ and an
! implicit part PQ in Powell's way: HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, issymmetric
implicit none

! Inputs
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)     ! PQ(NPT)
real(RP), intent(in) :: x(:)      ! X(N)
real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Outputs
real(RP) :: y(size(hq, 1))

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


function omega_col(idz, zmat, k) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates Y = column K of OMEGA, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
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
! This function calculates Y = OMEGA*X, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
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
! This function calculates P = X^T*OMEGA*Y, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
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
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm, planerot, hypotenuse, isorth, isminor
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


end module powalg_mod
