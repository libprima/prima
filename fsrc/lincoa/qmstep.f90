!*==qmstep.f90  processed by SPAG 7.50RE at 23:56 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used
!      SUBROUTINE QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE QMSTEP(N,Npt,M,Amat,Xpt,Xopt,Nact,Iact,Rescon,Qfac,    &
     &                  Kopt,Knew,Del,Step,Gl,Pqw,Rstat,W,Ifeas)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
      IMPLICIT NONE
!*--QMSTEP12
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Npt
      INTEGER , INTENT(IN) :: M
      REAL*8 , INTENT(IN) , DIMENSION(N,*) :: Amat
      REAL*8 , INTENT(IN) , DIMENSION(Npt,*) :: Xpt
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xopt
      INTEGER , INTENT(IN) :: Nact
      INTEGER , INTENT(IN) , DIMENSION(*) :: Iact
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Rescon
      REAL*8 , INTENT(IN) , DIMENSION(N,*) :: Qfac
      INTEGER , INTENT(IN) :: Kopt
      INTEGER , INTENT(IN) :: Knew
      REAL*8 , INTENT(IN) :: Del
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Step
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Gl
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Pqw
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Rstat
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER , INTENT(INOUT) :: Ifeas
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      REAL*8 :: bigv , ctol , gg , ghg , half , one , resmax , sp , ss ,  &
     &        stp , stpsav , sum , temp , tenth , test , vbig , vgrad , &
     &        vlag , vnew , ww , zero
      INTEGER :: i , j , jsav , k , ksav
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
      test = 0.2D0*Del
!
!     Replace GL by the gradient of LFUNC at the trust region centre, and
!       set the elements of RSTAT.
!
      DO k = 1 , Npt
         temp = zero
         DO j = 1 , N
            temp = temp + Xpt(k,j)*Xopt(j)
         ENDDO
         temp = Pqw(k)*temp
         DO i = 1 , N
            Gl(i) = Gl(i) + temp*Xpt(k,i)
         ENDDO
      ENDDO
      IF ( M>0 ) THEN
         DO j = 1 , M
            Rstat(j) = one
            IF ( DABS(Rescon(j))>=Del ) Rstat(j) = -one
         ENDDO
         DO k = 1 , Nact
            Rstat(Iact(k)) = zero
         ENDDO
      ENDIF
!
!     Find the greatest modulus of LFUNC on a line through XOPT and
!       another interpolation point within the trust region.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: IFLAG is never used
!      IFLAG=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      vbig = zero
      DO k = 1 , Npt
         IF ( k==Kopt ) CYCLE
         ss = zero
         sp = zero
         DO i = 1 , N
            temp = Xpt(k,i) - Xopt(i)
            ss = ss + temp*temp
            sp = sp + Gl(i)*temp
         ENDDO
         stp = -Del/DSQRT(ss)
         IF ( k==Knew ) THEN
            IF ( sp*(sp-one)<zero ) stp = -stp
            vlag = DABS(stp*sp) + stp*stp*DABS(sp-one)
         ELSE
            vlag = DABS(stp*(one-stp)*sp)
         ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: With the original code, if either VLAG or VBIG is
! NaN, KSAV will not get a value. This may cause Segmentation Fault
! because XPT(KSAV, :) will later be accessed.
!      IF (VLAG .GT. VBIG) THEN
         IF ( vlag>vbig ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ksav = k
            stpsav = stp
            vbig = vlag
         ENDIF
      ENDDO
!
!     Set STEP to the move that gives the greatest modulus calculated above.
!       This move may be replaced by a steepest ascent step from XOPT.
!
      gg = zero
      DO i = 1 , N
         gg = gg + Gl(i)**2
         Step(i) = stpsav*(Xpt(ksav,i)-Xopt(i))
      ENDDO
      vgrad = Del*DSQRT(gg)
      IF ( vgrad>tenth*vbig ) THEN
!
!     Make the replacement if it provides a larger value of VBIG.
!
         ghg = zero
         DO k = 1 , Npt
            temp = zero
            DO j = 1 , N
               temp = temp + Xpt(k,j)*Gl(j)
            ENDDO
            ghg = ghg + Pqw(k)*temp*temp
         ENDDO
         vnew = vgrad + DABS(half*Del*Del*ghg/gg)
         IF ( vnew>vbig ) THEN
            vbig = vnew
            stp = Del/DSQRT(gg)
            IF ( ghg<zero ) stp = -stp
            DO i = 1 , N
               Step(i) = stp*Gl(i)
            ENDDO
         ENDIF
         IF ( Nact/=0 .AND. Nact/=N ) THEN
!
!     Overwrite GL by its projection. Then set VNEW to the greatest
!       value of |LFUNC| on the projected gradient from XOPT subject to
!       the trust region bound. If VNEW is sufficiently large, then STEP
!       may be changed to a move along the projected gradient.
!
            DO k = Nact + 1 , N
               W(k) = zero
               DO i = 1 , N
                  W(k) = W(k) + Gl(i)*Qfac(i,k)
               ENDDO
            ENDDO
            gg = zero
            DO i = 1 , N
               Gl(i) = zero
               DO k = Nact + 1 , N
                  Gl(i) = Gl(i) + Qfac(i,k)*W(k)
               ENDDO
               gg = gg + Gl(i)**2
            ENDDO
            vgrad = Del*DSQRT(gg)
            IF ( vgrad>tenth*vbig ) THEN
               ghg = zero
               DO k = 1 , Npt
                  temp = zero
                  DO j = 1 , N
                     temp = temp + Xpt(k,j)*Gl(j)
                  ENDDO
                  ghg = ghg + Pqw(k)*temp*temp
               ENDDO
               vnew = vgrad + DABS(half*Del*Del*ghg/gg)
!
!     Set W to the possible move along the projected gradient.
!
               stp = Del/DSQRT(gg)
               IF ( ghg<zero ) stp = -stp
               ww = zero
               DO i = 1 , N
                  W(i) = stp*Gl(i)
                  ww = ww + W(i)**2
               ENDDO
!
!     Set STEP to W if W gives a sufficiently large value of the modulus
!       of the Lagrange function, and if W either preserves feasibility
!       or gives a constraint violation of at least 0.2*DEL. The purpose
!       of CTOL below is to provide a check on feasibility that includes
!       a tolerance for contributions from computer rounding errors.
!
               IF ( vnew/vbig>=0.2D0 ) THEN
                  Ifeas = 1
                  bigv = zero
                  j = 0
                  DO
                     j = j + 1
                     IF ( j<=M ) THEN
                        IF ( Rstat(j)==one ) THEN
                           temp = -Rescon(j)
                           DO i = 1 , N
                              temp = temp + W(i)*Amat(i,j)
                           ENDDO
                           bigv = DMAX1(bigv,temp)
                        ENDIF
                        IF ( bigv<test ) CYCLE
                        Ifeas = 0
                     ENDIF
                     ctol = zero
                     temp = 0.01D0*DSQRT(ww)
                     IF ( bigv>zero .AND. bigv<temp ) THEN
                        DO k = 1 , Nact
                           j = Iact(k)
                           sum = zero
                           DO i = 1 , N
                              sum = sum + W(i)*Amat(i,j)
                           ENDDO
                           ctol = DMAX1(ctol,DABS(sum))
                        ENDDO
                     ENDIF
                     IF ( bigv<=10.0D0*ctol .OR. bigv>=test ) THEN
                        DO i = 1 , N
                           Step(i) = W(i)
                        ENDDO
                        GOTO 99999
                     ENDIF
                     EXIT
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!
!     Calculate the greatest constraint violation at XOPT+STEP with STEP at
!       its original value. Modify STEP if this violation is unacceptable.
!
      Ifeas = 1
      bigv = zero
      resmax = zero
      j = 0
      DO
         j = j + 1
         IF ( j<=M ) THEN
            IF ( Rstat(j)<zero ) CYCLE
            temp = -Rescon(j)
            DO i = 1 , N
               temp = temp + Step(i)*Amat(i,j)
            ENDDO
            resmax = DMAX1(resmax,temp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .LT. TEST) THEN
            IF ( temp<test ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF ( temp>bigv ) THEN
                  bigv = temp
                  jsav = j
                  Ifeas = -1
               ENDIF
               CYCLE
            ENDIF
            Ifeas = 0
         ENDIF
         IF ( Ifeas==-1 ) THEN
            DO i = 1 , N
               Step(i) = Step(i) + (test-bigv)*Amat(i,jsav)
            ENDDO
            Ifeas = 0
         ENDIF
         EXIT
      ENDDO
!
!     Return the calculated STEP and the value of IFEAS.
!
99999 END SUBROUTINE QMSTEP
