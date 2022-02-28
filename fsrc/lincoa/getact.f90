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
! Last Modified: Monday, February 28, 2022 PM04:14:54
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, dd, dw, vlam)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : inprod

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
real(RP), intent(inout) :: rfac(:)  ! RFAC(N, N)

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
integer(IK) :: i, ic, idiag, iflag, j, jc, jcp, jdiag, jw,   &
&           k, l, nactp

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
    !call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    !call assert(istriu(rfac, tol), 'RFAC is upper triangular', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(size(rfac) == n * (n + 1_IK) / 2_IK, 'RFAC is n-by-n and upper triangular', srname)
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
tdel = 0.2D0 * snorm
ddsav = inprod(g, g) + inprod(g, g)
vlam = ZERO
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
if (nact == 0) then
    qfac = ZERO
    do i = 1, n
        qfac(i, i) = ONE
    end do
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
idiag = (ic * ic + ic) / 2
if (ic < nact) then
    jw = idiag + ic
    do j = ic + 1, nact
        temp = temp - rfac(jw) * vlam(j)
        jw = jw + j
    end do
end if
if (temp >= ZERO) goto 800
vlam(ic) = temp / rfac(idiag)
ic = ic - 1_IK
if (ic > 0) goto 70
!
!     Set the new search direction D. Terminate if the 2-norm of D is ZERO
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
100 if (nact == n) goto 290
do j = nact + 1, n
    w(j) = ZERO
    do i = 1, n
        w(j) = w(j) + qfac(i, j) * g(i)
    end do
end do
dd = ZERO
do i = 1, n
    dw(i) = ZERO
    do j = nact + 1, n
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
    temp = 0.01D0 * dnorm
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
if (violmx <= 10.0D0 * ctol) goto 300
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ONEs required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
nactp = nact + 1
idiag = (nactp * nactp - nactp) / 2
rdiag = ZERO
do j = n, 1, -1
    sprod = ZERO
    do i = 1, n
        sprod = sprod + qfac(i, j) * amat(i, l)
    end do
    if (j <= nact) then
        rfac(idiag + j) = sprod
    else
        if (abs(rdiag) <= 1.0D-20 * abs(sprod)) then
            rdiag = sprod
        else
            temp = sqrt(sprod * sprod + rdiag * rdiag)
            cosv = sprod / temp
            sinv = rdiag / temp
            rdiag = temp
            do i = 1, n
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
rfac(idiag + nactp) = abs(rdiag)
nact = nactp
iact(nact) = l
resact(nact) = resnew(l)
vlam(nact) = ZERO
resnew(l) = ZERO
!
!     Set the components of the vector VMU in W.
!
220 w(nact) = ONE / rfac((nact * nact + nact) / 2)**2
if (nact > 1) then
    do i = nact - 1, 1, -1
        idiag = (i * i + i) / 2
        jw = idiag + i
        summ = ZERO
        do j = i + 1, nact
            summ = summ - rfac(jw) * w(j)
            jw = jw + j
        end do
        w(i) = summ / rfac(idiag)
    end do
end if
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
vmult = violmx
ic = 0
j = 1
250 if (j < nact) then
    if (vlam(j) >= vmult * w(j)) then
        ic = j
        vmult = vlam(j) / w(j)
    end if
    j = j + 1
    goto 250
end if
do j = 1, nact
    vlam(j) = vlam(j) - vmult * w(j)
end do
if (ic > 0) vlam(ic) = ZERO
violmx = max(violmx - vmult, ZERO)
if (ic == 0) violmx = ZERO
!
!     Reduce the active set if necessary, so that all compONEnts of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
iflag = 3
ic = nact
!!!! If NACT=0, then IC = 0, and hence IACT(IC) is undefined, which leads to memory error when
!RESNEW(IACT(IC)) is accessed.
270 if (vlam(ic) < ZERO) goto 280
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
300 return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by ONE.
!
800 resnew(iact(ic)) = max(resact(ic), tinynum)
jc = ic
810 if (jc < nact) then
    jcp = jc + 1
    idiag = jc * jcp / 2
    jw = idiag + jcp
    temp = sqrt(rfac(jw - 1)**2 + rfac(jw)**2)
    cval = rfac(jw) / temp
    sval = rfac(jw - 1) / temp
    rfac(jw - 1) = sval * rfac(idiag)
    rfac(jw) = cval * rfac(idiag)
    rfac(idiag) = temp
    if (jcp < nact) then
        do j = jcp + 1, nact
            temp = sval * rfac(jw + jc) + cval * rfac(jw + jcp)
            rfac(jw + jcp) = cval * rfac(jw + jc) - sval * rfac(jw + jcp)
            rfac(jw + jc) = temp
            jw = jw + j
        end do
    end if
    jdiag = idiag - jc
    do i = 1, n
        if (i < jc) then
            temp = rfac(idiag + i)
            rfac(idiag + i) = rfac(jdiag + i)
            rfac(jdiag + i) = temp
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
nact = nact - 1
if (iflag == 1) then
    goto 50
elseif (iflag == 2) then
    goto 60
elseif (iflag == 3) then
    goto 280
end if
end subroutine getact


end module getact_mod
