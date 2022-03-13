module getact_mod
!--------------------------------------------------------------------------------------------------!
! This module provides the GETACT subroutine of LINCOA, which picks the current active set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Sunday, March 13, 2022 PM07:21:55
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, dd, dw, vlam)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, validate
use, non_intrinsic :: linalg_mod, only : inprod, eye, istriu

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: snorm

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)

! Outputs
real(RP), intent(out) :: dd
real(RP), intent(out) :: dw(:)  ! DW(N)  ; better name?
real(RP), intent(out) :: vlam(:)  ! VLAM(N)

! Local variables
character(len=*), parameter :: srname = 'GETACT'
real(RP) :: w(size(g))
real(RP) :: ctol, ddsav, dnorm, summ, tdel, temp, test, tinynum,   &
&        violmx, vmult
integer(IK) :: i, ic, iflag, j, k, l

integer(IK) :: m
integer(IK) :: n

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(g), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)

    !----------------------------------------------------------------------------------------------!
    !tol == ???
    !call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    !call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    !----------------------------------------------------------------------------------------------!

    call assert(size(dw) == n, 'SIZE(DW) == N', srname)
    call assert(size(vlam) == n, 'SIZE(VLAM) == N', srname)
end if


!     N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms
!       with these names in SUBROUTINE LINCOB. The current values must be
!       set on entry. NACT, IACT, QFAC and RFAC are kept up to date when
!       GETACT changes the current active set.
!     SNORM, RESNEW, RESACT, G and DW are the same as the terms with these
!       names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are
!       also kept up to date.
!     VLAM and W are used for working space, the vector VLAM being reserved
!       for the Lagrange multipliers of the calculation. Their lengths must
!       be at least N.
!     The main purpose of GETACT is to pick the current active set. It is
!       defined by the property that the projection of -G into the space
!       orthogonal to the active constraint normals is as large as possible,
!       subject to this projected steepest descent direction moving no closer
!       to the boundary of every constraint whose current residual is at most
!       0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are
!       all appropriate to this choice of active set.
!     Occasionally this projected direction is ZERO, and then the final value
!       of W(1) is set to ZERO. Otherwise, the direction itself is returned
!       in DW, and W(1) is set to the square of the length of the direction.
!
!     Set some constants and a temporary VLAM.
!
tinynum = real(tiny(0.0), RP)
tdel = 0.2_RP * snorm
ddsav = TWO * inprod(g, g)
vlam = ZERO
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
if (nact == 0) then
    qfac = eye(n)
    goto 100
end if
!
!     Remove any constraints from the initial active set whose residuals
!       exceed TDEL.
!
iflag = 1
ic = nact
40 if (resact(ic) > tdel) goto 800
50 ic = ic - 1
if (ic > 0) goto 40
!
!     Remove any constraints from the initial active set whose Lagrange
!       multipliers are nonnegative, and set the surviving multipliers.
!
iflag = 2
60 if (nact == 0) goto 100
ic = nact
70 temp = ZERO
do i = 1, n
    temp = temp + qfac(i, ic) * g(i)
end do
if (ic < nact) then
    do j = ic + 1, nact
        temp = temp - rfac(ic, j) * vlam(j)
    end do
end if
if (temp >= ZERO) goto 800
vlam(ic) = temp / rfac(ic, ic)
ic = ic - 1_IK
if (ic > 0) goto 70
!
!     Set the new search direction D. Terminate if the 2-norm of D is ZERO
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
100 if (nact == n) goto 290
!!if (nact < 0) return !??? See next line.
do j = nact + 1, n  ! Here we have to ensure NACT >= 0; is it guaranteed in theory?
    w(j) = ZERO
    do i = 1, n
        w(j) = w(j) + qfac(i, j) * g(i)
    end do
end do
dd = ZERO
do i = 1, n
    dw(i) = ZERO
    do j = nact + 1, n  ! Here we have to ensure NACT >= 0; is it guaranteed in theory?
        dw(i) = dw(i) - w(j) * qfac(i, j)
    end do
    dd = dd + dw(i)**2
end do
if (dd >= ddsav) goto 290
if (dd == ZERO) goto 300
ddsav = dd
dnorm = sqrt(dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
l = 0
if (m > 0) then
    test = dnorm / snorm
    violmx = ZERO
    do j = 1, m
        if (resnew(j) > ZERO .and. resnew(j) <= tdel) then
            summ = ZERO
            do i = 1, n
                summ = summ + amat(i, j) * dw(i)
            end do
            if (summ > test * resnew(j)) then
                if (summ > violmx) then
                    l = j
                    violmx = summ
                end if
            end if
        end if
    end do
    ctol = ZERO
    temp = 0.01_RP * dnorm
    if (violmx > ZERO .and. violmx < temp) then
        if (nact > 0) then
            do k = 1, nact
                j = iact(k)
                summ = ZERO
                do i = 1, n
                    summ = summ + dw(i) * amat(i, j)
                end do
                ctol = max(ctol, abs(summ))
            end do
        end if
    end if
end if
w(1) = ONE
if (l == 0) goto 300
if (violmx <= TEN * ctol) goto 300
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
!--------------------------------------------------------------------------------------------------!
! Zaikun 20220305: If NACT >= N, then NACTP >= N+1, RFAC(IDIAT+J) and QFAC(I, NACTP) will be
! invalid. Is it guaranteed that NACT < N in theory? Probably yes because of line number 100, where
! NACT == N leads to return.
if (nact >= n) goto 300
!--------------------------------------------------------------------------------------------------!
call validate(nact < n, 'NACT < N', srname)
call qradd(amat(:, l), qfac, rfac, nact)

iact(nact) = l
resact(nact) = resnew(l)
vlam(nact) = ZERO
resnew(l) = ZERO
!
!     Set the components of the vector VMU in W.
!
220 continue

w(nact) = ONE / rfac(nact, nact)**2
if (nact > 1) then
    do i = nact - 1, 1, -1
        summ = ZERO
        do j = i + 1, nact
            summ = summ - rfac(i, j) * w(j)
        end do
        w(i) = summ / rfac(i, i)
    end do
end if
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
vmult = violmx
ic = 0
j = 1
250 continue
if (j < nact) then
    if (vlam(j) >= vmult * w(j)) then
        ic = j
        vmult = vlam(j) / w(j)
    end if
    j = j + 1
    goto 250
end if
! Very strangely, if we (mistakenly) change 'J = N, 1, -1' to 'J = N-1, 1, -1' in
! "Apply Givens rotations to the last (N-NACT) columns of QFAC", then the following lines
! encounter a SEGFAULT when this subroutine is called with NACT = 0 and we arrive here with NACT = IC = 0.
do j = 1, nact
    vlam(j) = vlam(j) - vmult * w(j)
end do
if (ic > 0) vlam(ic) = ZERO
violmx = max(violmx - vmult, ZERO)
if (ic == 0) violmx = ZERO
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
iflag = 3
!--------------------------------------------------------------------------------------------------!
! Zaikun 2021 July, 20220305:
! If NACT <= 0, then IC <= 0, and hence memory errors will occur when accessing VLAM(IC),
! IACT(IC), RESNEW(IACT(IC)). Is NACT >= 1 ensured theoretically? What about NACT <= N?
if (nact <= 0) goto 300  ! What about DD and W(1)???
!--------------------------------------------------------------------------------------------------!
ic = nact
270 continue
if (vlam(ic) < ZERO) goto 280
resnew(iact(ic)) = max(resact(ic), tinynum)
goto 800
280 ic = ic - 1
if (ic > 0) goto 270
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
if (violmx > ZERO) goto 220
if (nact < n) goto 100
290 dd = ZERO
!300 w(1) = dd
300 continue


!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    !----------------------------------------------------------------------------------------------!
    !tol == ???
    !call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    !call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    !----------------------------------------------------------------------------------------------!

    call assert(size(dw) == n, 'SIZE(DW) == N', srname)
    call assert(size(vlam) == n, 'SIZE(VLAM) == N', srname)
end if
return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
800 continue


call validate(ic <= nact, 'IC <= NACT', srname)  ! When IC > NACT, the following lines are invalid.
resnew(iact(ic)) = max(resact(ic), tinynum)

!call qrexc(qfac, rfac(:, 1:nact), ic)
! It suffices to pass only the first NACT columns of QFAC and the first NACT rows of RFAC as follows.
call qrexc(qfac(:, 1:nact), rfac(1:nact, 1:nact), ic)

iact(ic:nact) = [iact(ic + 1:nact), iact(ic)]
resact(ic:nact) = [resact(ic + 1:nact), resact(ic)]
vlam(ic:nact) = [vlam(ic + 1:nact), vlam(ic)]
nact = nact - 1
if (iflag == 1) then
    goto 50
elseif (iflag == 2) then
    goto 60
elseif (iflag == 3) then
    goto 280
end if
end subroutine getact


subroutine qradd(c, Q, R, n)
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : istriu, inprod, matprod, planerot
implicit none

integer(IK), intent(inout) :: n
real(RP), intent(in) :: c(:)
real(RP), intent(inout) :: Q(:, :)
real(RP), intent(inout) :: R(:, :)


character(len=*), parameter :: srname = 'QRADD'
integer(IK) :: k, m
real(RP) :: G(2, 2)
real(RP) :: cq(size(Q, 2))

m = int(size(Q, 1), kind(m))

if (DEBUGGING) then
    call assert(n >= 0 .and. n <= m - 1, '0 <= N <= M - 1', srname)
    call assert(size(c) == m, 'SIZE(C) == M', srname)
    call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    call assert(size(Q, 2) >= n .and. size(Q, 2) <= m, 'N <= SIZE(Q, 2) <= M', srname)
    call assert(size(R, 1) >= n .and. size(R, 1) <= m, 'N <= SIZE(R, 1) <= M', srname)
    !tol == ???
    !call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    !Qsave = Q  ! For debugging only.
    !Rsave = R(:, 1:i - 1)  ! For debugging only.
end if

cq = matprod(c, Q)

! Update Q so that the columns of Q(:, N+2:M) are orthogonal to C. This is done by applying a 2D
! Givens rotation to Q(:, [K, K+1]) from the right to zero C'*Q(:, K+1) out for K = N+1, ..., M-1.
! Nothing will be done if N >= M-1.
do k = m - 1_IK, n + 1_IK, -1
    !if (abs(cq(k + 1)) > 1.0D-20 * abs(cq(k))) then  ! Powell's version
    if (abs(cq(k + 1)) > 0) then
        G = planerot(cq([k, k + 1_IK]))
        Q(:, [k, k + 1_IK]) = matprod(Q(:, [k, k + 1_IK]), transpose(G))
        cq(k) = sqrt(cq(k)**2 + cq(k + 1)**2)
    end if
end do

R(1:n, n + 1) = matprod(c, Q(:, 1:n))

if (cq(n + 1) < 0) then
    Q(:, n + 1) = -Q(:, n + 1)
end if
R(n + 1, n + 1) = abs(cq(n + 1))

n = n + 1_IK

end subroutine qradd


subroutine qrexc(Q, R, i)
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
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, planerot, istriu

! Inputs
integer(IK), intent(in) :: i

! In-outputs
real(RP), intent(inout) :: Q(:, :)  ! Q(M, :), SIZE(Q, 2) <= M
real(RP), intent(inout) :: R(:, :)  ! R(:, N), SIZE(R, 1) >= N

! Local variables
character(len=*), parameter :: srname = 'QREXC'
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
real(RP) :: G(2, 2)
real(RP) :: hypt
real(RP) :: Qsave(size(Q, 1), size(Q, 2))  ! Debugging only
real(RP) :: Rsave(size(R, 1), max(i - 1_IK, 0_IK))  ! I >= 1 if the input is correct; debugging only

! Sizes
m = size(Q, 1)
n = size(R, 2)

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. n <= m, '1 <= N <= M', srname)
    call assert(i >= 1 .and. i <= n, '1 <= I <= N', srname)
    call assert(size(Q, 2) == size(R, 1), 'SIZE(Q, 2) == SIZE(R, 1)', srname)
    call assert(size(Q, 2) >= n .and. size(Q, 2) <= m, 'N <= SIZE(Q, 2) <= M', srname)
    call assert(size(R, 1) >= n .and. size(R, 1) <= m, 'N <= SIZE(R, 1) <= M', srname)
    !tol == ???
    !call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    Qsave = Q  ! For debugging only.
    Rsave = R(:, 1:i - 1)  ! For debugging only.
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
    !tol == ???
    !call assert(isorth(Q, tol), 'The columns of Q are orthogonal', srname)
    call assert(istriu(R), 'R is upper triangular', srname)
    call assert(.not. any(abs(Q(:, 1:i - 1) - Qsave(:, 1:i - 1)) > 0) .and. &
        & .not. any(abs(Q(:, n + 1:) - Qsave(:, n + 1:)) > 0), 'Q is unchanged except Q(:, I:N)', srname)
    call assert(.not. any(abs(R(:, 1:i - 1) - Rsave) > 0), 'R(:, 1:N-1) is unchanged', srname)
    ! If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead.
    !call assert(all(abs(Q(:, 1:i - 1) - Qsave(:, 1:i - 1)) <= 0) .and. &
    !    & all(abs(Q(:, n + 1:) - Qsave(:, n + 1:)) <= 0), 'Q is unchanged except Q(:, I:N)', srname)
    !call assert(all(abs(R(:, 1:i - 1) - Rsave) <= 0), 'R(:, 1:N-1) is unchanged', srname)
end if

end subroutine qrexc

end module getact_mod
