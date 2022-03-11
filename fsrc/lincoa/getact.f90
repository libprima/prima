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
! Last Modified: Saturday, March 12, 2022 AM12:54:20
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, dd, dw, vlam)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
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
real(RP) :: cosv, ctol, cval, ddsav, dnorm, rdiag,   &
&        sinv, sprod, summ, sval, tdel, temp, test, tinynum,   &
&        violmx, vmult
integer(IK) :: i, ic, idiag, iflag, j, jc, jcp, jdiag,  &
&           k, l, nactp, jw

integer(IK) :: m
integer(IK) :: n

real(RP) :: rfacv(size(g) * (size(g) + 1) / 2)


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
nactp = nact + 1
rdiag = ZERO

do j = n, 1, -1
!do j = n - 1, 1, -1  ! Why does this lead to SEGFAULT in DEGENLPA when calling profile('lincoa')?
    sprod = ZERO
    do i = 1, n
        sprod = sprod + qfac(i, j) * amat(i, l)
    end do
    if (j <= nact) then
        rfac(j, nact + 1) = sprod
    else
        !if (abs(rdiag) <= 1.0D-20 * abs(sprod)) then
        if (j == n .or. abs(rdiag) <= 1.0D-20 * abs(sprod)) then
            ! Zaikun 20220304: what if j = n ? How to ensure that we will not go to the else? Is it
            ! because RDIAG = 0 when J = n? What if SPROD = NaN?
            ! Why do we observe out of boundary error in stest_i8_r4_d1_tst?
            ! Why does stest encounter NaN in X when calling evaluate?
            rdiag = sprod
        else
            temp = sqrt(sprod * sprod + rdiag * rdiag)
            cosv = sprod / temp
            sinv = rdiag / temp
            rdiag = temp
            do i = 1, n
                !if (j >= n) cycle
                temp = cosv * qfac(i, j) + sinv * qfac(i, j + 1)
                qfac(i, j + 1) = -sinv * qfac(i, j) + cosv * qfac(i, j + 1)
                qfac(i, j) = temp
            end do
        end if
    end if
end do

if (rdiag < ZERO) then
    do i = 1, n
        qfac(i, nactp) = -qfac(i, nactp)
    end do
end if
rfac(nact + 1, nact + 1) = abs(rdiag)

nact = nactp
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

!-------------------!
rfacv = m2v(rfac)
!-------------------!

resnew(iact(ic)) = max(resact(ic), tinynum)
jc = ic
810 if (jc < nact) then
    jcp = jc + 1
    idiag = jc * jcp / 2
    jw = idiag + jcp
    temp = sqrt(rfacv(jw - 1)**2 + rfacv(jw)**2)
    cval = rfacv(jw) / temp
    sval = rfacv(jw - 1) / temp
    rfacv(jw - 1) = sval * rfacv(idiag)
    rfacv(jw) = cval * rfacv(idiag)
    rfacv(idiag) = temp
    if (jcp < nact) then
        do j = jcp + 1, nact
            temp = sval * rfacv(jw + jc) + cval * rfacv(jw + jcp)
            rfacv(jw + jcp) = cval * rfacv(jw + jc) - sval * rfacv(jw + jcp)
            rfacv(jw + jc) = temp
            jw = jw + j
        end do
    end if
    jdiag = idiag - jc
    do i = 1, n
        if (i < jc) then
            temp = rfacv(idiag + i)
            rfacv(idiag + i) = rfacv(jdiag + i)
            rfacv(jdiag + i) = temp
        end if
        temp = sval * qfac(i, jc) + cval * qfac(i, jcp)
        qfac(i, jcp) = cval * qfac(i, jc) - sval * qfac(i, jcp)
        qfac(i, jc) = temp
    end do
    iact(jc) = iact(jcp)
    resact(jc) = resact(jcp)
    vlam(jc) = vlam(jcp)
    jc = jcp
    goto 810
end if

!-------------------!
rfac = v2m(rfacv)
!-------------------!

!--------------------------------------------------------------------------------------------------!
! Zaikun 20220305: without the next line, SEGFAULT may occur below line number 100. Better choice? GOTO 100?
if (nact <= 0) goto 300  ! What about DD and W(1)???
!--------------------------------------------------------------------------------------------------!
nact = nact - 1
if (iflag == 1) then
    goto 50
elseif (iflag == 2) then
    goto 60
elseif (iflag == 3) then
    goto 280
end if
end subroutine getact


function v2m(rfacv) result(rfacm)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
use, non_intrinsic :: debug_mod, only : assert
implicit none
real(RP), intent(in) :: rfacv(:)
real(RP) :: rfacm((floor(sqrt(real(8 * size(rfacv) + 1))) - 1) / 2, (floor(sqrt(real(8 * size(rfacv) + 1))) - 1) / 2)
integer(IK) :: n, i, j, ir

n = int(size(rfacm, 1), kind(n))
call assert(size(rfacv) == n * (n + 1) / 2, 'SIZE(RFACV) = N*(N+1)/2', 'v2m')

rfacm = ZERO
ir = 0_IK
do j = 1_IK, n
    do i = 1_IK, j
        ir = ir + 1_IK
        rfacm(i, j) = rfacv(ir)
    end do
end do

end function v2m

function m2v(rfacm) result(rfacv)
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
real(RP), intent(in) :: rfacm(:, :)
real(RP) :: rfacv((size(rfacm, 1) * (size(rfacm, 1) + 1)) / 2)
integer(IK) :: n, i, j, ir

n = int(size(rfacm, 1), kind(n))
ir = 0_IK
do j = 1_IK, n
    do i = 1_IK, j
        ir = ir + 1_IK
        rfacv(ir) = rfacm(i, j)
    end do
end do

end function m2v

end module getact_mod
