!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trsbox.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 26-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==trsbox.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
            subroutine TRSBOX(N, Npt, Xpt, Xopt, Gopt, Hq, Pq, Sl, Su, D&
     &elta, Xnew, D, Gnew, Xbdi, S, Hs, Hred, Dsq, Crvmin)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--TRSBOX8
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer, intent(IN) :: N
            integer, intent(IN) :: Npt
            real*8, intent(IN), dimension(Npt, *) :: Xpt
            real*8, intent(IN), dimension(*) :: Xopt
            real*8, intent(IN), dimension(*) :: Gopt
            real*8, intent(IN), dimension(*) :: Hq
            real*8, intent(IN), dimension(*) :: Pq
            real*8, intent(IN), dimension(*) :: Sl
            real*8, intent(IN), dimension(*) :: Su
            real*8, intent(IN) :: Delta
            real*8, intent(INOUT), dimension(*) :: Xnew
            real*8, intent(INOUT), dimension(*) :: D
            real*8, intent(INOUT), dimension(*) :: Gnew
            real*8, intent(INOUT), dimension(*) :: Xbdi
            real*8, intent(INOUT), dimension(*) :: S
            real*8, intent(INOUT), dimension(*) :: Hs
            real*8, intent(INOUT), dimension(*) :: Hred
            real*8, intent(INOUT) :: Dsq
            real*8, intent(INOUT) :: Crvmin
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: angbd, angt, beta, blen, cth, delsq, dhd, dhs, dre&
     &dg, dredsq, ds, ggsav, gredsq, half, one, onemin, qred, rdnext, rd&
     &prev, redmax, rednew, redsav, resid, sdec, shs, sredg, ssq, stepsq&
     &, sth, stplen, temp, tempa, tempb, xsav, xsum, zero
            integer :: i, iact, ih, isav, itcsav, iterc, itermax, iu, j,&
     & k, nact
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
!       meanings as the corresponding arguments of BOBYQB.
!     DELTA is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance DELTA of
!       XOPT subject to the bounds on the variables.
!     XNEW will be set to a new vector of variables that is approximately
!       the one that minimizes the quadratic model within the trust region
!       subject to the SL and SU constraints on the variables. It satisfies
!       as equations the bounds that become active during the calculation.
!     D is the calculated trial step from XOPT, generated iteratively from an
!       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
!     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
!       when D is updated.
!     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
!       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
!       I-th variable has become fixed at a bound, the bound being SL(I) or
!       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
!       information is accumulated during the construction of XNEW.
!     The arrays S, HS and HRED are also used for working space. They hold the
!       current search direction, and the changes in the gradient of Q along S
!       and the reduced D, respectively, where the reduced D is the same as D,
!       except that the components of the fixed variables are zero.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-1.0D0 is set, however, if all of these searches are
!       constrained.
!
!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each one being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            onemin = -1.0D0
            zero = 0.0D0
!
!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at one of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.
!
            iterc = 0
            nact = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: SQSTP is never used
!      SQSTP=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i = 1, N
                Xbdi(i) = zero
                if (Xopt(i) <= Sl(i)) then
                    if (Gopt(i) >= zero) Xbdi(i) = onemin
                elseif (Xopt(i) >= Su(i)) then
                    if (Gopt(i) <= zero) Xbdi(i) = one
                end if
                if (Xbdi(i) /= zero) nact = nact + 1
                D(i) = zero
                Gnew(i) = Gopt(i)
            end do
            delsq = Delta * Delta
            qred = zero
            Crvmin = onemin
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
      100 beta = zero
      200 stepsq = zero
            do i = 1, N
                if (Xbdi(i) /= zero) then
                    S(i) = zero
                elseif (beta == zero) then
                    S(i) = -Gnew(i)
                else
                    S(i) = beta * S(i) - Gnew(i)
                end if
                stepsq = stepsq + S(i)**2
            end do
            if (stepsq == zero) goto 500
            if (beta == zero) then
                gredsq = stepsq
                itermax = iterc + N - nact
            end if
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
            if (gredsq * delsq > 1.0D-4 * qred * qred) goto 600
            goto 500
!
!     Prepare for the alternative iteration by calculating some scalars and
!     by multiplying the reduced D by the second derivative matrix of Q.
!
      300 if (nact >= N - 1) goto 500
            dredsq = zero
            dredg = zero
            gredsq = zero
            do i = 1, N
                if (Xbdi(i) == zero) then
                    dredsq = dredsq + D(i)**2
                    dredg = dredg + D(i) * Gnew(i)
                    gredsq = gredsq + Gnew(i)**2
                    S(i) = D(i)
                else
                    S(i) = zero
                end if
            end do
            itcsav = iterc
            goto 600
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
      400 iterc = iterc + 1
            temp = gredsq * dredsq - dredg * dredg
            if (temp > 1.0D-4 * qred * qred) then
                temp = DSQRT(temp)
                do i = 1, N
                    if (Xbdi(i) == zero) then
                        S(i) = (dredg * D(i) - dredsq * Gnew(i)) / temp
                    else
                        S(i) = zero
                    end if
                end do
                sredg = -temp
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
                angbd = one
                iact = 0
                do i = 1, N
                    if (Xbdi(i) == zero) then
                        tempa = Xopt(i) + D(i) - Sl(i)
                        tempb = Su(i) - Xopt(i) - D(i)
                        if (tempa <= zero) then
                            nact = nact + 1
                            Xbdi(i) = onemin
                            goto 300
                        elseif (tempb <= zero) then
                            nact = nact + 1
                            Xbdi(i) = one
                            goto 300
                        end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: RATIO is never used
!          RATIO=ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ssq = D(i)**2 + S(i)**2
                        temp = ssq - (Xopt(i) - Sl(i))**2
                        if (temp > zero) then
                            temp = DSQRT(temp) - S(i)
                            if (angbd * temp > tempa) then
                                angbd = tempa / temp
                                iact = i
                                xsav = onemin
                            end if
                        end if
                        temp = ssq - (Su(i) - Xopt(i))**2
                        if (temp > zero) then
                            temp = DSQRT(temp) + S(i)
                            if (angbd * temp > tempb) then
                                angbd = tempb / temp
                                iact = i
                                xsav = one
                            end if
                        end if
                    end if
                end do
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
                goto 600
            end if
      500 Dsq = zero
            do i = 1, N
                Xnew(i) = DMAX1(DMIN1(Xopt(i) + D(i), Su(i)), Sl(i))
                if (Xbdi(i) == onemin) Xnew(i) = Sl(i)
                if (Xbdi(i) == one) Xnew(i) = Su(i)
                D(i) = Xnew(i) - Xopt(i)
                Dsq = Dsq + D(i)**2
            end do
            return

!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
      600 ih = 0
            do j = 1, N
                Hs(j) = zero
                do i = 1, j
                    ih = ih + 1
                    if (i < j) Hs(j) = Hs(j) + Hq(ih) * S(i)
                    Hs(i) = Hs(i) + Hq(ih) * S(j)
                end do
            end do
            do k = 1, Npt
                if (Pq(k) /= zero) then
                    temp = zero
                    do j = 1, N
                        temp = temp + Xpt(k, j) * S(j)
                    end do
                    temp = temp * Pq(k)
                    do i = 1, N
                        Hs(i) = Hs(i) + temp * Xpt(k, i)
                    end do
                end if
            end do
            if (Crvmin /= zero) then
                resid = delsq
                ds = zero
                shs = zero
                do i = 1, N
                    if (Xbdi(i) == zero) then
                        resid = resid - D(i)**2
                        ds = ds + S(i) * D(i)
                        shs = shs + S(i) * Hs(i)
                    end if
                end do
                if (resid > zero) then
                    temp = DSQRT(stepsq * resid + ds * ds)
                    if (ds < zero) then
                        blen = (temp - ds) / stepsq
                    else
                        blen = resid / (temp + ds)
                    end if
                    stplen = blen
                    if (shs > zero) stplen = DMIN1(blen, gredsq / shs)

!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
                    iact = 0
                    do i = 1, N
                        if (S(i) /= zero) then
                            xsum = Xopt(i) + D(i)
                            if (S(i) > zero) then
                                temp = (Su(i) - xsum) / S(i)
                            else
                                temp = (Sl(i) - xsum) / S(i)
                            end if
                            if (temp < stplen) then
                                stplen = temp
                                iact = i
                            end if
                        end if
                    end do
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
                    sdec = zero
                    if (stplen > zero) then
                        iterc = iterc + 1
                        temp = shs / stepsq
                        if (iact == 0 .and. temp > zero) then
                            Crvmin = DMIN1(Crvmin, temp)
                            if (Crvmin == onemin) Crvmin = temp
                        end if
                        ggsav = gredsq
                        gredsq = zero
                        do i = 1, N
                            Gnew(i) = Gnew(i) + stplen * Hs(i)
                            if (Xbdi(i) == zero) gredsq = gredsq + Gnew(&
     &i)**2
                            D(i) = D(i) + stplen * S(i)
                        end do
                        sdec = DMAX1(stplen * (ggsav - half * stplen * s&
     &hs), zero)
                        qred = qred + sdec
                    end if
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
                    if (iact > 0) then
                        nact = nact + 1
                        Xbdi(iact) = one
                        if (S(iact) < zero) Xbdi(iact) = onemin
                        delsq = delsq - D(iact)**2
                        if (delsq > zero) goto 100
                        goto 650
                    end if
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
                    if (stplen < blen) then
                        if (iterc == itermax) goto 500
                        if (sdec <= 0.01D0 * qred) goto 500
                        beta = gredsq / ggsav
                        goto 200
                    end if
                end if
      650 Crvmin = zero
                goto 300
            elseif (iterc > itcsav) then
                shs = zero
                dhs = zero
                dhd = zero
                do i = 1, N
                    if (Xbdi(i) == zero) then
                        shs = shs + S(i) * Hs(i)
                        dhs = dhs + D(i) * Hs(i)
                        dhd = dhd + D(i) * Hred(i)
                    end if
                end do
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
                redmax = zero
                isav = 0
                redsav = zero
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IU=17.0D0*ANGBD+3.1D0
                iu = int(17.0D0 * angbd + 3.1D0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do i = 1, iu
                    angt = angbd * DFLOAT(i) / DFLOAT(iu)
                    sth = (angt + angt) / (one + angt * angt)
                    temp = shs + angt * (angt * dhd - dhs - dhs)
                    rednew = sth * (angt * dredg - sredg - half * sth * &
     &temp)
                    if (rednew > redmax) then
                        redmax = rednew
                        isav = i
                        rdprev = redsav
                    elseif (i == isav + 1) then
                        rdnext = rednew
                    end if
                    redsav = rednew
                end do
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
                if (isav == 0) goto 500
                if (isav < iu) then
                    temp = (rdnext - rdprev) / (redmax + redmax - rdprev&
     & - rdnext)
                    angt = angbd * (DFLOAT(isav) + half * temp) / DFLOAT&
     &(iu)
                end if
                cth = (one - angt * angt) / (one + angt * angt)
                sth = (angt + angt) / (one + angt * angt)
                temp = shs + angt * (angt * dhd - dhs - dhs)
                sdec = sth * (angt * dredg - sredg - half * sth * temp)
                if (sdec <= zero) goto 500
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
                dredg = zero
                gredsq = zero
                do i = 1, N
                    Gnew(i) = Gnew(i) + (cth - one) * Hred(i) + sth * Hs&
     &(i)
                    if (Xbdi(i) == zero) then
                        D(i) = cth * D(i) + sth * S(i)
                        dredg = dredg + D(i) * Gnew(i)
                        gredsq = gredsq + Gnew(i)**2
                    end if
                    Hred(i) = cth * Hred(i) + sth * Hs(i)
                end do
                qred = qred + sdec
                if (iact > 0 .and. isav == iu) then
                    nact = nact + 1
                    Xbdi(iact) = xsav
                    goto 300
                end if
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
                if (sdec <= 0.01D0 * qred) goto 500
                goto 400
            else
                do i = 1, N
                    Hred(i) = Hs(i)
                end do
                goto 400
            end if
            end subroutine TRSBOX