!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of qmstep.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 07-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==qmstep.f90  processed by SPAG 7.50RE at 17:53 on 31 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used
!      SUBROUTINE QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine QMSTEP(N, Npt, M, Amat, Xpt, Xopt, Nact, Iact, Re&
     &scon, Qfac, Kopt, Knew, Del, Step, Gl, Pqw, Rstat, W, Ifeas)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--QMSTEP12
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer, intent(IN) :: N
            integer, intent(IN) :: Npt
            integer, intent(IN) :: M
            real*8, intent(IN), dimension(N, *) :: Amat
            real*8, intent(IN), dimension(Npt, *) :: Xpt
            real*8, intent(IN), dimension(*) :: Xopt
            integer, intent(IN) :: Nact
            integer, intent(IN), dimension(*) :: Iact
            real*8, intent(IN), dimension(*) :: Rescon
            real*8, intent(IN), dimension(N, *) :: Qfac
            integer, intent(IN) :: Kopt
            integer, intent(IN) :: Knew
            real*8, intent(IN) :: Del
            real*8, intent(INOUT), dimension(*) :: Step
            real*8, intent(INOUT), dimension(*) :: Gl
            real*8, intent(IN), dimension(*) :: Pqw
            real*8, intent(INOUT), dimension(*) :: Rstat
            real*8, intent(INOUT), dimension(*) :: W
            integer, intent(INOUT) :: Ifeas
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: bigv, ctol, gg, ghg, half, one, resmax, sp, ss, st&
     &p, stpsav, sum, temp, tenth, test, vbig, vgrad, vlag, vnew, ww, ze&
     &ro
            integer :: i, j, jsav, k, ksav
!*++
!*++ End of declarations rewritten by SPAG
!*++
!      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),XOPT(*),IACT(*),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
!       same as the terms with these names in SUBROUTINE LINCOB.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DEL is the current restriction on the length of STEP, which is never
!       greater than the current trust region radius DELTA.
!     STEP will be set to the required step from XOPT to the new point.
!     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
!       is the KNEW-th Lagrange function. It is used also for some other
!       gradients of LFUNC.
!     PQW provides the second derivative parameters of LFUNC.
!     RSTAT and W are used for working space. Their lengths must be at least
!       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
!       J-th constraint is irrelevant, active, or both inactive and relevant,
!       respectively.
!     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
!
!     STEP is chosen to provide a relatively large value of the modulus of
!       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
!       calculated too, within the trust region, that does not alter the
!       residuals of the active constraints. The projected step is preferred
!       if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the
!       original one, but the greatest violation of a linear constraint must
!       be at least 0.2*DEL, in order to keep the interpolation points apart.
!       The remedy when the maximum constraint violation is too small is to
!       restore the original step, which is perturbed if necessary so that
!       its maximum constraint violation becomes 0.2*DEL.
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            tenth = 0.1D0
            zero = 0.0D0
            test = 0.2D0 * Del
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
            do k = 1, Npt
                temp = zero
                do j = 1, N
                    temp = temp + Xpt(k, j) * Xopt(j)
                end do
                temp = Pqw(k) * temp
                do i = 1, N
                    Gl(i) = Gl(i) + temp * Xpt(k, i)
                end do
            end do
            if (M > 0) then
                do j = 1, M
                    Rstat(j) = one
                    if (DABS(Rescon(j)) >= Del) Rstat(j) = -one
                end do
                do k = 1, Nact
                    Rstat(Iact(k)) = zero
                end do
            end if
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: IFLAG is never used
!      IFLAG=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            vbig = zero
            do k = 1, Npt
                if (k == Kopt) cycle
                ss = zero
                sp = zero
                do i = 1, N
                    temp = Xpt(k, i) - Xopt(i)
                    ss = ss + temp * temp
                    sp = sp + Gl(i) * temp
                end do
                stp = -Del / DSQRT(ss)
                if (k == Knew) then
                    if (sp * (sp - one) < zero) stp = -stp
                    vlag = DABS(stp * sp) + stp * stp * DABS(sp - one)
                else
                    vlag = DABS(stp * (one - stp) * sp)
                end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: With the original code, if either VLAG or VBIG is
! NaN, KSAV will not get a value. This may cause Segmentation Fault
! because XPT(KSAV, :) will later be accessed.
!      IF (VLAG .GT. VBIG) THEN
                if (.not. (vlag <= vbig)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ksav = k
                    stpsav = stp
                    vbig = vlag
                end if
            end do
!
!     Set STEP to the move that gives the greatest modulus calculated above.
!       This move may be replaced by a steepest ascent step from XOPT.
!
            gg = zero
            do i = 1, N
                gg = gg + Gl(i)**2
                Step(i) = stpsav * (Xpt(ksav, i) - Xopt(i))
            end do
            vgrad = Del * DSQRT(gg)
            if (vgrad > tenth * vbig) then
!
!     Make the replacement if it provides a larger value of VBIG.
!
                ghg = zero
                do k = 1, Npt
                    temp = zero
                    do j = 1, N
                        temp = temp + Xpt(k, j) * Gl(j)
                    end do
                    ghg = ghg + Pqw(k) * temp * temp
                end do
                vnew = vgrad + DABS(half * Del * Del * ghg / gg)
                if (vnew > vbig) then
                    vbig = vnew
                    stp = Del / DSQRT(gg)
                    if (ghg < zero) stp = -stp
                    do i = 1, N
                        Step(i) = stp * Gl(i)
                    end do
                end if
                if (Nact /= 0 .and. Nact /= N) then
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
                    do k = Nact + 1, N
                        W(k) = zero
                        do i = 1, N
                            W(k) = W(k) + Gl(i) * Qfac(i, k)
                        end do
                    end do
                    gg = zero
                    do i = 1, N
                        Gl(i) = zero
                        do k = Nact + 1, N
                            Gl(i) = Gl(i) + Qfac(i, k) * W(k)
                        end do
                        gg = gg + Gl(i)**2
                    end do
                    vgrad = Del * DSQRT(gg)
                    if (vgrad > tenth * vbig) then
                        ghg = zero
                        do k = 1, Npt
                            temp = zero
                            do j = 1, N
                                temp = temp + Xpt(k, j) * Gl(j)
                            end do
                            ghg = ghg + Pqw(k) * temp * temp
                        end do
                        vnew = vgrad + DABS(half * Del * Del * ghg / gg)
!
!     Set W to the possible move along the projected gradient.
!
                        stp = Del / DSQRT(gg)
                        if (ghg < zero) stp = -stp
                        ww = zero
                        do i = 1, N
                            W(i) = stp * Gl(i)
                            ww = ww + W(i)**2
                        end do
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
                        if (vnew / vbig >= 0.2D0) then
                            Ifeas = 1
                            bigv = zero
                            j = 0
                            do
                                j = j + 1
                                if (j <= M) then
                                    if (Rstat(j) == one) then
                                        temp = -Rescon(j)
                                        do i = 1, N
                                            temp = temp + W(i) * Amat(i,&
     & j)
                                        end do
                                        bigv = DMAX1(bigv, temp)
                                    end if
                                    if (bigv < test) cycle
                                    Ifeas = 0
                                end if
                                ctol = zero
                                temp = 0.01D0 * DSQRT(ww)
                                if (bigv > zero .and. bigv < temp) then
                                    do k = 1, Nact
                                        j = Iact(k)
                                        sum = zero
                                        do i = 1, N
                                            sum = sum + W(i) * Amat(i, j&
     &)
                                        end do
                                        ctol = DMAX1(ctol, DABS(sum))
                                    end do
                                end if
                                if (bigv <= 10.0D0 * ctol .or. bigv >= t&
     &est) then
                                    do i = 1, N
                                        Step(i) = W(i)
                                    end do
                                    goto 99999
                                end if
                                exit
                            end do
                        end if
                    end if
                end if
            end if
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
            Ifeas = 1
            bigv = zero
            resmax = zero
            j = 0
            do
                j = j + 1
                if (j <= M) then
                    if (Rstat(j) < zero) cycle
                    temp = -Rescon(j)
                    do i = 1, N
                        temp = temp + Step(i) * Amat(i, j)
                    end do
                    resmax = DMAX1(resmax, temp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .LT. TEST) THEN
                    if (.not. (temp >= test)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if (temp > bigv) then
                            bigv = temp
                            jsav = j
                            Ifeas = -1
                        end if
                        cycle
                    end if
                    Ifeas = 0
                end if
                if (Ifeas == -1) then
                    do i = 1, N
                        Step(i) = Step(i) + (test - bigv) * Amat(i, jsav&
     &)
                    end do
                    Ifeas = 0
                end if
                exit
            end do
!
!     Return the calculated STEP and the value of IFEAS.
!
      99999 end subroutine QMSTEP