!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of getact.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 11-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==getact.f90  processed by SPAG 7.50RE at 17:53 on 31 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used
!      SUBROUTINE GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine GETACT(N, M, Amat, Nact, Iact, Qfac, Rfac, Snorm,&
     & Resnew, Resact, G, Dw, Vlam, W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--GETACT12
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer, intent(IN) :: N
            integer, intent(IN) :: M
            real*8, intent(IN), dimension(N, *) :: Amat
            integer, intent(INOUT) :: Nact
            integer, intent(INOUT), dimension(*) :: Iact
            real*8, intent(INOUT), dimension(N, *) :: Qfac
            real*8, intent(INOUT), dimension(*) :: Rfac
            real*8, intent(IN) :: Snorm
            real*8, intent(INOUT), dimension(*) :: Resnew
            real*8, intent(INOUT), dimension(*) :: Resact
            real*8, intent(IN), dimension(*) :: G
            real*8, intent(INOUT), dimension(*) :: Dw
            real*8, intent(INOUT), dimension(*) :: Vlam
            real*8, intent(INOUT), dimension(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: cosv, ctol, cval, dd, ddsav, dnorm, one, rdiag, si&
     &nv, sprod, sum, sval, tdel, temp, test, tiny, violmx, vmult, zero
            integer :: i, ic, idiag, iflag, j, jc, jcp, jdiag, jw, k, l,&
     & nactp
!*++
!*++ End of declarations rewritten by SPAG
!*++
!      DIMENSION AMAT(N,*),B(*),IACT(*),QFAC(N,*),RFAC(*),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            tdel = 0.2D0 * Snorm
            ddsav = zero
            do i = 1, N
                ddsav = ddsav + G(i)**2
                Vlam(i) = zero
            end do
            ddsav = ddsav + ddsav
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
            if (Nact == 0) then
                do i = 1, N
                    do j = 1, N
                        Qfac(i, j) = zero
                    end do
                    Qfac(i, i) = one
                end do
                goto 400
            end if
!
!     Remove any constraints from the initial active set whose residuals
!       exceed TDEL.
!
            iflag = 1
            ic = Nact
100   if (Resact(ic) > tdel) goto 900
200   ic = ic - 1
            if (ic > 0) goto 100
!
!     Remove any constraints from the initial active set whose Lagrange
!       multipliers are nonnegative, and set the surviving multipliers.
!
            iflag = 2
300   if (Nact /= 0) then
                ic = Nact
                do
                    temp = zero
                    do i = 1, N
                        temp = temp + Qfac(i, ic) * G(i)
                    end do
                    idiag = (ic * ic + ic) / 2
                    if (ic < Nact) then
                        jw = idiag + ic
                        do j = ic + 1, Nact
                            temp = temp - Rfac(jw) * Vlam(j)
                            jw = jw + j
                        end do
                    end if
                    if (temp >= zero) goto 900
                    Vlam(ic) = temp / Rfac(idiag)
                    ic = ic - 1
                    if (ic <= 0) exit
                end do
            end if
!
!     Set the new search direction D. Terminate if the 2-norm of D is zero
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
400   if (Nact == N) then
                dd = zero
                goto 800
            else
                do j = Nact + 1, N
                    W(j) = zero
                    do i = 1, N
                        W(j) = W(j) + Qfac(i, j) * G(i)
                    end do
                end do
                dd = zero
                do i = 1, N
                    Dw(i) = zero
                    do j = Nact + 1, N
                        Dw(i) = Dw(i) - W(j) * Qfac(i, j)
                    end do
                    dd = dd + Dw(i)**2
                end do
                if (dd >= ddsav) then
                    dd = zero
                    goto 800
                else
                    if (dd == zero) goto 800
                    ddsav = dd
                    dnorm = DSQRT(dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
                    l = 0
                    if (M > 0) then
                        test = dnorm / Snorm
                        violmx = zero
                        do j = 1, M
                            if (Resnew(j) > zero .and. Resnew(j) <= tdel&
     &) then
                                sum = zero
                                do i = 1, N
                                    sum = sum + Amat(i, j) * Dw(i)
                                end do
                                if (sum > test * Resnew(j)) then
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
                            if (Nact > 0) then
                                do k = 1, Nact
                                    j = Iact(k)
                                    sum = zero
                                    do i = 1, N
                                        sum = sum + Dw(i) * Amat(i, j)
                                    end do
                                    ctol = DMAX1(ctol, DABS(sum))
                                end do
                            end if
                        end if
                    end if
                    W(1) = one
                    if (l == 0) goto 800
                    if (violmx <= 10.0D0 * ctol) goto 800
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
                    nactp = Nact + 1
                    idiag = (nactp * nactp - nactp) / 2
                    rdiag = zero
                    do j = N, 1, -1
                        sprod = zero
                        do i = 1, N
                            sprod = sprod + Qfac(i, j) * Amat(i, l)
                        end do
                        if (j <= Nact) then
                            Rfac(idiag + j) = sprod
                        elseif (DABS(rdiag) <= 1.0D-20 * DABS(sprod)) th&
     &en
                            rdiag = sprod
                        else
                            temp = DSQRT(sprod * sprod + rdiag * rdiag)
                            cosv = sprod / temp
                            sinv = rdiag / temp
                            rdiag = temp
                            do i = 1, N
                                temp = cosv * Qfac(i, j) + sinv * Qfac(i&
     &, j + 1)
                                Qfac(i, j + 1) = -sinv * Qfac(i, j) + co&
     &sv * Qfac(i, j + 1)
                                Qfac(i, j) = temp
                            end do
                        end if
                    end do
                    if (rdiag < zero) then
                        do i = 1, N
                            Qfac(i, nactp) = -Qfac(i, nactp)
                        end do
                    end if
                    Rfac(idiag + nactp) = DABS(rdiag)
                    Nact = nactp
                    Iact(Nact) = l
                    Resact(Nact) = Resnew(l)
                    Vlam(Nact) = zero
                    Resnew(l) = zero
                end if
            end if
!
!     Set the components of the vector VMU in W.
!
500   W(Nact) = one / Rfac((Nact * Nact + Nact) / 2)**2
            if (Nact > 1) then
                do i = Nact - 1, 1, -1
                    idiag = (i * i + i) / 2
                    jw = idiag + i
                    sum = zero
                    do j = i + 1, Nact
                        sum = sum - Rfac(jw) * W(j)
                        jw = jw + j
                    end do
                    W(i) = sum / Rfac(idiag)
                end do
            end if
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
            vmult = violmx
            ic = 0
            j = 1
            do
                if (j < Nact) then
                    if (Vlam(j) >= vmult * W(j)) then
                        ic = j
                        vmult = Vlam(j) / W(j)
                    end if
                    j = j + 1
                    cycle
                end if
                do j = 1, Nact
                    Vlam(j) = Vlam(j) - vmult * W(j)
                end do
                if (ic > 0) Vlam(ic) = zero
                violmx = DMAX1(violmx - vmult, zero)
                if (ic == 0) violmx = zero
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
                iflag = 3
                ic = Nact
                exit
            end do
!!!! If NACT=0, then IC = 0, and hence IACT(IC) is undefined, which leads to memory error when
!RESNEW(IACT(IC)) is accessed.
600   if (Vlam(ic) >= zero) then
                Resnew(Iact(ic)) = DMAX1(Resact(ic), tiny)
                goto 900
            end if
700   ic = ic - 1
            if (ic > 0) goto 600
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
            if (violmx > zero) goto 500
            if (Nact < N) goto 400
            dd = zero
800   W(1) = dd
            return
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
900   Resnew(Iact(ic)) = DMAX1(Resact(ic), tiny)
            jc = ic
            do
                if (jc < Nact) then
                    jcp = jc + 1
                    idiag = jc * jcp / 2
                    jw = idiag + jcp
                    temp = DSQRT(Rfac(jw - 1)**2 + Rfac(jw)**2)
                    cval = Rfac(jw) / temp
                    sval = Rfac(jw - 1) / temp
                    Rfac(jw - 1) = sval * Rfac(idiag)
                    Rfac(jw) = cval * Rfac(idiag)
                    Rfac(idiag) = temp
                    if (jcp < Nact) then
                        do j = jcp + 1, Nact
                            temp = sval * Rfac(jw + jc) + cval * Rfac(jw&
     & + jcp)
                            Rfac(jw + jcp) = cval * Rfac(jw + jc) - sval&
     & * Rfac(jw + jcp)
                            Rfac(jw + jc) = temp
                            jw = jw + j
                        end do
                    end if
                    jdiag = idiag - jc
                    do i = 1, N
                        if (i < jc) then
                            temp = Rfac(idiag + i)
                            Rfac(idiag + i) = Rfac(jdiag + i)
                            Rfac(jdiag + i) = temp
                        end if
                        temp = sval * Qfac(i, jc) + cval * Qfac(i, jcp)
                        Qfac(i, jcp) = cval * Qfac(i, jc) - sval * Qfac(&
     &i, jcp)
                        Qfac(i, jc) = temp
                    end do
                    Iact(jc) = Iact(jcp)
                    Resact(jc) = Resact(jcp)
                    Vlam(jc) = Vlam(jcp)
                    jc = jcp
                    cycle
                end if
                Nact = Nact - 1
                if (iflag == 1) goto 200
                if (iflag == 2) goto 300
                if (iflag == 3) goto 700
                exit
            end do
            end subroutine GETACT