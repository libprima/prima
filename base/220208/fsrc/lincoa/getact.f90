subroutine getact(n, m, amat, nact, iact, qfac, rfac, snorm, resnew, resact, g, dw, vlam, w)
use, non_intrinsic :: linalg_mod, only : inprod
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
!      DIMENSION AMAT(N,*),B(*),IACT(*),QFAC(N,*),RFAC(*),
dimension amat(n, m), iact(m), qfac(n, n), rfac(n * (n + 1) / 2), resnew(m), resact(m), g(n), dw(n), vlam(n), w(n)
!
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
!     Occasionally this projected direction is zero, and then the final value
!       of W(1) is set to zero. Otherwise, the direction itself is returned
!       in DW, and W(1) is set to the square of the length of the direction.
!
!     Set some constants and a temporary VLAM.
!
one = 1.0D0
tiny = 1.0D-60
zero = 0.0D0
tdel = 0.2D0 * snorm
ddsav = inprod(g, g) + inprod(g, g)
vlam = 0.0D0
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
if (nact == 0) then
    qfac = zero
    do i = 1, n
        qfac(i, i) = one
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
70 temp = zero
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
if (temp >= zero) goto 800
vlam(ic) = temp / rfac(idiag)
ic = ic - 1
if (ic > 0) goto 70
!
!     Set the new search direction D. Terminate if the 2-norm of D is zero
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
100 if (nact == n) goto 290
do j = nact + 1, n
    w(j) = zero
    do i = 1, n
        w(j) = w(j) + qfac(i, j) * g(i)
    end do
end do
dd = zero
do i = 1, n
    dw(i) = zero
    do j = nact + 1, n
        dw(i) = dw(i) - w(j) * qfac(i, j)
    end do
    dd = dd + dw(i)**2
end do
if (dd >= ddsav) goto 290
if (dd == zero) goto 300
ddsav = dd
dnorm = dsqrt(dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
l = 0
if (m > 0) then
    test = dnorm / snorm
    violmx = zero
    do j = 1, m
        if (resnew(j) > zero .and. resnew(j) <= tdel) then
            sum = zero
            do i = 1, n
                sum = sum + amat(i, j) * dw(i)
            end do
            if (sum > test * resnew(j)) then
                if (sum > violmx) then
                    l = j
                    violmx = sum
                end if
            end if
        end if
    end do
    ctol = zero
    temp = 0.01D0 * dnorm
    if (violmx > zero .and. violmx < temp) then
        if (nact > 0) then
            do k = 1, nact
                j = iact(k)
                sum = zero
                do i = 1, n
                    sum = sum + dw(i) * amat(i, j)
                end do
                ctol = dmax1(ctol, dabs(sum))
            end do
        end if
    end if
end if
w(1) = one
if (l == 0) goto 300
if (violmx <= 10.0D0 * ctol) goto 300
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
nactp = nact + 1
idiag = (nactp * nactp - nactp) / 2
rdiag = zero
do j = n, 1, -1
    sprod = zero
    do i = 1, n
        sprod = sprod + qfac(i, j) * amat(i, l)
    end do
    if (j <= nact) then
        rfac(idiag + j) = sprod
    else
        if (dabs(rdiag) <= 1.0D-20 * dabs(sprod)) then
            rdiag = sprod
        else
            temp = dsqrt(sprod * sprod + rdiag * rdiag)
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
if (rdiag < zero) then
    do i = 1, n
        qfac(i, nactp) = -qfac(i, nactp)
    end do
end if
rfac(idiag + nactp) = dabs(rdiag)
nact = nactp
iact(nact) = l
resact(nact) = resnew(l)
vlam(nact) = zero
resnew(l) = zero
!
!     Set the components of the vector VMU in W.
!
220 w(nact) = one / rfac((nact * nact + nact) / 2)**2
if (nact > 1) then
    do i = nact - 1, 1, -1
        idiag = (i * i + i) / 2
        jw = idiag + i
        sum = zero
        do j = i + 1, nact
            sum = sum - rfac(jw) * w(j)
            jw = jw + j
        end do
        w(i) = sum / rfac(idiag)
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
if (ic > 0) vlam(ic) = zero
violmx = dmax1(violmx - vmult, zero)
if (ic == 0) violmx = zero
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
iflag = 3
ic = nact
!!!! If NACT=0, then IC = 0, and hence IACT(IC) is undefined, which leads to memory error when
!RESNEW(IACT(IC)) is accessed.
270 if (vlam(ic) < zero) goto 280
resnew(iact(ic)) = dmax1(resact(ic), tiny)
goto 800
280 ic = ic - 1
if (ic > 0) goto 270
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
if (violmx > zero) goto 220
if (nact < n) goto 100
290 dd = zero
300 w(1) = dd
return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
800 resnew(iact(ic)) = dmax1(resact(ic), tiny)
jc = ic
810 if (jc < nact) then
    jcp = jc + 1
    idiag = jc * jcp / 2
    jw = idiag + jcp
    temp = dsqrt(rfac(jw - 1)**2 + rfac(jw)**2)
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
goto(50, 60, 280), iflag
end
