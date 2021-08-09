!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of altmov.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 10-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==altmov.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
            subroutine ALTMOV(N, Npt, Xpt, Xopt, Bmat, Zmat, Ndim, Sl, S&
     &u, Kopt, Knew, Adelt, Xnew, Xalt, Alpha, Cauchy, Glag, Hcol, W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--ALTMOV8
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer, intent(IN) :: N
            integer, intent(IN) :: Npt
            real*8, intent(IN), dimension(Npt, *) :: Xpt
            real*8, intent(IN), dimension(*) :: Xopt
            real*8, intent(IN), dimension(Ndim, *) :: Bmat
            real*8, intent(IN), dimension(Npt, *) :: Zmat
            integer, intent(IN) :: Ndim
            real*8, intent(IN), dimension(*) :: Sl
            real*8, intent(IN), dimension(*) :: Su
            integer, intent(IN) :: Kopt
            integer, intent(IN) :: Knew
            real*8, intent(IN) :: Adelt
            real*8, intent(OUT), dimension(*) :: Xnew
            real*8, intent(INOUT), dimension(*) :: Xalt
            real*8, intent(INOUT) :: Alpha
            real*8, intent(INOUT) :: Cauchy
            real*8, intent(INOUT), dimension(*) :: Glag
            real*8, intent(INOUT), dimension(*) :: Hcol
            real*8, intent(INOUT), dimension(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: bigstp, const, csave, curv, dderiv, diff, distsq, &
     &ggfree, gw, ha, half, one, predsq, presav, scale, slbd, step, stps&
     &av, subd, sumin, temp, tempa, tempb, tempd, vlag, wfixsq, wsqsav, &
     &zero
            integer :: i, ibdsav, iflag, ilbd, isbd, iubd, j, k, ksav
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!
            half = 0.5D0
            one = 1.0D0
            zero = 0.0D0
            const = one + DSQRT(2.0D0)
            do k = 1, Npt
                Hcol(k) = zero
            end do
            do j = 1, Npt - N - 1
                temp = Zmat(Knew, j)
                do k = 1, Npt
                    Hcol(k) = Hcol(k) + temp * Zmat(k, j)
                end do
            end do
            Alpha = Hcol(Knew)
            ha = half * Alpha
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
            do i = 1, N
                Glag(i) = Bmat(Knew, i)
            end do
            do k = 1, Npt
                temp = zero
                do j = 1, N
                    temp = temp + Xpt(k, j) * Xopt(j)
                end do
                temp = Hcol(k) * temp
                do i = 1, N
                    Glag(i) = Glag(i) + temp * Xpt(k, i)
                end do
            end do
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
            presav = zero
            do k = 1, Npt
                if (k == Kopt) cycle
                dderiv = zero
                distsq = zero
                do i = 1, N
                    temp = Xpt(k, i) - Xopt(i)
                    dderiv = dderiv + Glag(i) * temp
                    distsq = distsq + temp * temp
                end do
                subd = Adelt / DSQRT(distsq)
                slbd = -subd
                ilbd = 0
                iubd = 0
                sumin = DMIN1(one, subd)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
                do i = 1, N
                    temp = Xpt(k, i) - Xopt(i)
                    if (temp > zero) then
                        if (slbd * temp < Sl(i) - Xopt(i)) then
                            slbd = (Sl(i) - Xopt(i)) / temp
                            ilbd = -i
                        end if
                        if (subd * temp > Su(i) - Xopt(i)) then
                            subd = DMAX1(sumin, (Su(i) - Xopt(i)) / temp&
     &)
                            iubd = i
                        end if
                    elseif (temp < zero) then
                        if (slbd * temp > Su(i) - Xopt(i)) then
                            slbd = (Su(i) - Xopt(i)) / temp
                            ilbd = i
                        end if
                        if (subd * temp < Sl(i) - Xopt(i)) then
                            subd = DMAX1(sumin, (Sl(i) - Xopt(i)) / temp&
     &)
                            iubd = -i
                        end if
                    end if
                end do
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
                if (k == Knew) then
                    diff = dderiv - one
                    step = slbd
                    vlag = slbd * (dderiv - slbd * diff)
                    isbd = ilbd
                    temp = subd * (dderiv - subd * diff)
                    if (DABS(temp) > DABS(vlag)) then
                        step = subd
                        vlag = temp
                        isbd = iubd
                    end if
                    tempd = half * dderiv
                    tempa = tempd - diff * slbd
                    tempb = tempd - diff * subd
                    if (tempa * tempb < zero) then
                        temp = tempd * tempd / diff
                        if (DABS(temp) > DABS(vlag)) then
                            step = tempd / diff
                            vlag = temp
                            isbd = 0
                        end if
                    end if
!
!     Search along each of the other lines through XOPT and another point.
!
                else
                    step = slbd
                    vlag = slbd * (one - slbd)
                    isbd = ilbd
                    temp = subd * (one - subd)
                    if (DABS(temp) > DABS(vlag)) then
                        step = subd
                        vlag = temp
                        isbd = iubd
                    end if
                    if (subd > half) then
                        if (DABS(vlag) < 0.25D0) then
                            step = half
                            vlag = 0.25D0
                            isbd = 0
                        end if
                    end if
                    vlag = vlag * dderiv
                end if
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
                temp = step * (one - step) * distsq
                predsq = vlag * vlag * (vlag * vlag + ha * temp * temp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: With the original code, if either PREDSQ or PRESAV
! is NaN, KSAV/STPSAV/IBDSAV will not get a value. This may cause
! Segmentation Fault.
!      IF (PREDSQ .GT. PRESAV) THEN
                if (predsq > presav) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    presav = predsq
                    ksav = k
                    stpsav = step
                    ibdsav = isbd
                end if
            end do
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
            do i = 1, N
                temp = Xopt(i) + stpsav * (Xpt(ksav, i) - Xopt(i))
                Xnew(i) = DMAX1(Sl(i), DMIN1(Su(i), temp))
            end do
            if (ibdsav < 0) Xnew(-ibdsav) = Sl(-ibdsav)
            if (ibdsav > 0) Xnew(ibdsav) = Su(ibdsav)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
            bigstp = Adelt + Adelt
            iflag = 0
100   wfixsq = zero
            ggfree = zero
            do i = 1, N
                W(i) = zero
                tempa = DMIN1(Xopt(i) - Sl(i), Glag(i))
                tempb = DMAX1(Xopt(i) - Su(i), Glag(i))
                if (tempa > zero .or. tempb < zero) then
                    W(i) = bigstp
                    ggfree = ggfree + Glag(i)**2
                end if
            end do
            if (ggfree == zero) then
                Cauchy = zero
                goto 99999
            end if
            do
!
!     Investigate whether more components of W can be fixed.
!
                temp = Adelt * Adelt - wfixsq
                if (temp > zero) then
                    wsqsav = wfixsq
                    step = DSQRT(temp / ggfree)
                    ggfree = zero
                    do i = 1, N
                        if (W(i) == bigstp) then
                            temp = Xopt(i) - step * Glag(i)
                            if (temp <= Sl(i)) then
                                W(i) = Sl(i) - Xopt(i)
                                wfixsq = wfixsq + W(i)**2
                            elseif (temp >= Su(i)) then
                                W(i) = Su(i) - Xopt(i)
                                wfixsq = wfixsq + W(i)**2
                            else
                                ggfree = ggfree + Glag(i)**2
                            end if
                        end if
                    end do
                    if (wfixsq > wsqsav .and. ggfree > zero) cycle
                end if
                exit
            end do
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
            gw = zero
            do i = 1, N
                if (W(i) == bigstp) then
                    W(i) = -step * Glag(i)
                    Xalt(i) = DMAX1(Sl(i), DMIN1(Su(i), Xopt(i) + W(i)))
                elseif (W(i) == zero) then
                    Xalt(i) = Xopt(i)
                elseif (Glag(i) > zero) then
                    Xalt(i) = Sl(i)
                else
                    Xalt(i) = Su(i)
                end if
                gw = gw + Glag(i) * W(i)
            end do
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
            curv = zero
            do k = 1, Npt
                temp = zero
                do j = 1, N
                    temp = temp + Xpt(k, j) * W(j)
                end do
                curv = curv + Hcol(k) * temp * temp
            end do
            if (iflag == 1) curv = -curv
            if (curv > -gw .and. curv < -const * gw) then
                scale = -gw / curv
                do i = 1, N
                    temp = Xopt(i) + scale * W(i)
                    Xalt(i) = DMAX1(Sl(i), DMIN1(Su(i), temp))
                end do
                Cauchy = (half * gw * scale)**2
            else
                Cauchy = (gw + half * curv)**2
            end if
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
            if (iflag == 0) then
                do i = 1, N
                    Glag(i) = -Glag(i)
                    W(N + i) = Xalt(i)
                end do
                csave = Cauchy
                iflag = 1
                goto 100
            end if
            if (csave > Cauchy) then
                do i = 1, N
                    Xalt(i) = W(N + i)
                end do
                Cauchy = csave
            end if
99999 end subroutine ALTMOV