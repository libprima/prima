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
! Generated using the interform.m script by Zaikun ZHANG (www.zhangzk.net)
! on 11-Feb-2022.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine qmstep(n, npt, m, amat, xpt, xopt, nact, iact, rescon, &
     &qfac, kopt, knew, del, step, gl, pqw, rstat, w, ifeas)
!      IMPLICIT REAL*8 (A-H,O-Z)
      implicit real(kind(0.0D0)) (a - h, o - z)
      implicit integer(i - n)
!      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),XOPT(*),IACT(*),
      dimension amat(n, *), xpt(npt, *), xopt(*), iact(*), rescon(*), qf&
     &ac(n, *), step(*), gl(*), pqw(*), rstat(*), w(*)
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
      test = 0.2D0 * del
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
      do k = 1, npt
          temp = zero
          do j = 1, n
              temp = temp + xpt(k, j) * xopt(j)
          end do
          temp = pqw(k) * temp
          do i = 1, n
              gl(i) = gl(i) + temp * xpt(k, i)
          end do
      end do
      if (m > 0) then
          do j = 1, m
              rstat(j) = one
              if (dabs(rescon(j)) >= del) rstat(j) = -one
          end do
          do k = 1, nact
              rstat(iact(k)) = zero
          end do
      end if
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-15: IFLAG is never used
!      IFLAG=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      vbig = zero
      do k = 1, npt
          if (k == kopt) cycle
          ss = zero
          sp = zero
          do i = 1, n
              temp = xpt(k, i) - xopt(i)
              ss = ss + temp * temp
              sp = sp + gl(i) * temp
          end do
          stp = -del / dsqrt(ss)
          if (k == knew) then
              if (sp * (sp - one) < zero) stp = -stp
              vlag = dabs(stp * sp) + stp * stp * dabs(sp - one)
          else
              vlag = dabs(stp * (one - stp) * sp)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      do i = 1, n
          gg = gg + gl(i)**2
          step(i) = stpsav * (xpt(ksav, i) - xopt(i))
      end do
      vgrad = del * dsqrt(gg)
      if (vgrad <= tenth * vbig) goto 220
!
!     Make the replacement if it provides a larger value of VBIG.
!
      ghg = zero
      do k = 1, npt
          temp = zero
          do j = 1, n
              temp = temp + xpt(k, j) * gl(j)
          end do
          ghg = ghg + pqw(k) * temp * temp
      end do
      vnew = vgrad + dabs(half * del * del * ghg / gg)
      if (vnew > vbig) then
          vbig = vnew
          stp = del / dsqrt(gg)
          if (ghg < zero) stp = -stp
          do i = 1, n
              step(i) = stp * gl(i)
          end do
      end if
      if (nact == 0 .or. nact == n) goto 220
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
      do k = nact + 1, n
          w(k) = zero
          do i = 1, n
              w(k) = w(k) + gl(i) * qfac(i, k)
          end do
      end do
      gg = zero
      do i = 1, n
          gl(i) = zero
          do k = nact + 1, n
              gl(i) = gl(i) + qfac(i, k) * w(k)
          end do
          gg = gg + gl(i)**2
      end do
      vgrad = del * dsqrt(gg)
      if (vgrad <= tenth * vbig) goto 220
      ghg = zero
      do k = 1, npt
          temp = zero
          do j = 1, n
              temp = temp + xpt(k, j) * gl(j)
          end do
          ghg = ghg + pqw(k) * temp * temp
      end do
      vnew = vgrad + dabs(half * del * del * ghg / gg)
!
!     Set W to the possible move along the projected gradient.
!
      stp = del / dsqrt(gg)
      if (ghg < zero) stp = -stp
      ww = zero
      do i = 1, n
          w(i) = stp * gl(i)
          ww = ww + w(i)**2
      end do
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
      if (vnew / vbig >= 0.2D0) then
          ifeas = 1
          bigv = zero
          j = 0
170   j = j + 1
          if (j <= m) then
              if (rstat(j) == one) then
                  temp = -rescon(j)
                  do i = 1, n
                      temp = temp + w(i) * amat(i, j)
                  end do
                  bigv = dmax1(bigv, temp)
              end if
              if (bigv < test) goto 170
              ifeas = 0
          end if
          ctol = zero
          temp = 0.01D0 * dsqrt(ww)
          if (bigv > zero .and. bigv < temp) then
              do k = 1, nact
                  j = iact(k)
                  sum = zero
                  do i = 1, n
                      sum = sum + w(i) * amat(i, j)
                  end do
                  ctol = dmax1(ctol, dabs(sum))
              end do
          end if
          if (bigv <= 10.0D0 * ctol .or. bigv >= test) then
              do i = 1, n
                  step(i) = w(i)
              end do
              goto 260
          end if
      end if
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
220   ifeas = 1
      bigv = zero
      resmax = zero
      j = 0
230   j = j + 1
      if (j <= m) then
          if (rstat(j) < zero) goto 230
          temp = -rescon(j)
          do i = 1, n
              temp = temp + step(i) * amat(i, j)
          end do
          resmax = dmax1(resmax, temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IF (TEMP .LT. TEST) THEN
          if (.not. (temp >= test)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (temp <= bigv) goto 230
              bigv = temp
              jsav = j
              ifeas = -1
              goto 230
          end if
          ifeas = 0
      end if
      if (ifeas == -1) then
          do i = 1, n
              step(i) = step(i) + (test - bigv) * amat(i, jsav)
          end do
          ifeas = 0
      end if
!
!     Return the calculated STEP and the value of IFEAS.
!
260   return
      end