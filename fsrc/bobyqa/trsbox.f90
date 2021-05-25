!*==trsbox.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
      SUBROUTINE TRSBOX(N,Npt,Xpt,Xopt,Gopt,Hq,Pq,Sl,Su,Delta,Xnew,D,   &
     &                  Gnew,Xbdi,S,Hs,Hred,Dsq,Crvmin)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
      USE F77KINDS                        
      IMPLICIT NONE
!*--TRSBOX8
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Npt
      REAL*8 , INTENT(IN) , DIMENSION(Npt,*) :: Xpt
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xopt
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Gopt
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Hq
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Pq
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Sl
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Su
      REAL*8 , INTENT(IN) :: Delta
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xnew
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: D
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Gnew
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xbdi
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: S
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Hs
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Hred
      REAL*8 , INTENT(INOUT) :: Dsq
      REAL*8 , INTENT(INOUT) :: Crvmin
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      REAL*8 :: angbd , angt , beta , blen , cth , delsq , dhd , dhs ,    &
     &        dredg , dredsq , ds , ggsav , gredsq , half , one ,       &
     &        onemin , qred , rdnext , rdprev , redmax , rednew ,       &
     &        redsav , resid , sdec , shs , sredg , ssq , stepsq , sth ,&
     &        stplen , temp , tempa , tempb , xsav , xsum , zero
      INTEGER :: i , iact , ih , isav , itcsav , iterc , itermax , iu , &
     &           j , k , nact
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
      DO i = 1 , N
         Xbdi(i) = zero
         IF ( Xopt(i)<=Sl(i) ) THEN
            IF ( Gopt(i)>=zero ) Xbdi(i) = onemin
         ELSEIF ( Xopt(i)>=Su(i) ) THEN
            IF ( Gopt(i)<=zero ) Xbdi(i) = one
         ENDIF
         IF ( Xbdi(i)/=zero ) nact = nact + 1
         D(i) = zero
         Gnew(i) = Gopt(i)
      ENDDO
      delsq = Delta*Delta
      qred = zero
      Crvmin = onemin
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
 100  beta = zero
 200  stepsq = zero
      DO i = 1 , N
         IF ( Xbdi(i)/=zero ) THEN
            S(i) = zero
         ELSEIF ( beta==zero ) THEN
            S(i) = -Gnew(i)
         ELSE
            S(i) = beta*S(i) - Gnew(i)
         ENDIF
         stepsq = stepsq + S(i)**2
      ENDDO
      IF ( stepsq==zero ) GOTO 500
      IF ( beta==zero ) THEN
         gredsq = stepsq
         itermax = iterc + N - nact
      ENDIF
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
      IF ( gredsq*delsq>1.0D-4*qred*qred ) GOTO 600
      GOTO 500
!
!     Prepare for the alternative iteration by calculating some scalars and
!     by multiplying the reduced D by the second derivative matrix of Q.
!
 300  IF ( nact>=N-1 ) GOTO 500
      dredsq = zero
      dredg = zero
      gredsq = zero
      DO i = 1 , N
         IF ( Xbdi(i)==zero ) THEN
            dredsq = dredsq + D(i)**2
            dredg = dredg + D(i)*Gnew(i)
            gredsq = gredsq + Gnew(i)**2
            S(i) = D(i)
         ELSE
            S(i) = zero
         ENDIF
      ENDDO
      itcsav = iterc
      GOTO 600
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
 400  iterc = iterc + 1
      temp = gredsq*dredsq - dredg*dredg
      IF ( temp>1.0D-4*qred*qred ) THEN
         temp = DSQRT(temp)
         DO i = 1 , N
            IF ( Xbdi(i)==zero ) THEN
               S(i) = (dredg*D(i)-dredsq*Gnew(i))/temp
            ELSE
               S(i) = zero
            ENDIF
         ENDDO
         sredg = -temp
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
         angbd = one
         iact = 0
         DO i = 1 , N
            IF ( Xbdi(i)==zero ) THEN
               tempa = Xopt(i) + D(i) - Sl(i)
               tempb = Su(i) - Xopt(i) - D(i)
               IF ( tempa<=zero ) THEN
                  nact = nact + 1
                  Xbdi(i) = onemin
                  GOTO 300
               ELSEIF ( tempb<=zero ) THEN
                  nact = nact + 1
                  Xbdi(i) = one
                  GOTO 300
               ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: RATIO is never used
!          RATIO=ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ssq = D(i)**2 + S(i)**2
               temp = ssq - (Xopt(i)-Sl(i))**2
               IF ( temp>zero ) THEN
                  temp = DSQRT(temp) - S(i)
                  IF ( angbd*temp>tempa ) THEN
                     angbd = tempa/temp
                     iact = i
                     xsav = onemin
                  ENDIF
               ENDIF
               temp = ssq - (Su(i)-Xopt(i))**2
               IF ( temp>zero ) THEN
                  temp = DSQRT(temp) + S(i)
                  IF ( angbd*temp>tempb ) THEN
                     angbd = tempb/temp
                     iact = i
                     xsav = one
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
         GOTO 600
      ENDIF
 500  Dsq = zero
      DO i = 1 , N
         Xnew(i) = DMAX1(DMIN1(Xopt(i)+D(i),Su(i)),Sl(i))
         IF ( Xbdi(i)==onemin ) Xnew(i) = Sl(i)
         IF ( Xbdi(i)==one ) Xnew(i) = Su(i)
         D(i) = Xnew(i) - Xopt(i)
         Dsq = Dsq + D(i)**2
      ENDDO
      RETURN
 
!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
 600  ih = 0
      DO j = 1 , N
         Hs(j) = zero
         DO i = 1 , j
            ih = ih + 1
            IF ( i<j ) Hs(j) = Hs(j) + Hq(ih)*S(i)
            Hs(i) = Hs(i) + Hq(ih)*S(j)
         ENDDO
      ENDDO
      DO k = 1 , Npt
         IF ( Pq(k)/=zero ) THEN
            temp = zero
            DO j = 1 , N
               temp = temp + Xpt(k,j)*S(j)
            ENDDO
            temp = temp*Pq(k)
            DO i = 1 , N
               Hs(i) = Hs(i) + temp*Xpt(k,i)
            ENDDO
         ENDIF
      ENDDO
      IF ( Crvmin/=zero ) THEN
         resid = delsq
         ds = zero
         shs = zero
         DO i = 1 , N
            IF ( Xbdi(i)==zero ) THEN
               resid = resid - D(i)**2
               ds = ds + S(i)*D(i)
               shs = shs + S(i)*Hs(i)
            ENDIF
         ENDDO
         IF ( resid>zero ) THEN
            temp = DSQRT(stepsq*resid+ds*ds)
            IF ( ds<zero ) THEN
               blen = (temp-ds)/stepsq
            ELSE
               blen = resid/(temp+ds)
            ENDIF
            stplen = blen
            IF ( shs>zero ) stplen = DMIN1(blen,gredsq/shs)
 
!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
            iact = 0
            DO i = 1 , N
               IF ( S(i)/=zero ) THEN
                  xsum = Xopt(i) + D(i)
                  IF ( S(i)>zero ) THEN
                     temp = (Su(i)-xsum)/S(i)
                  ELSE
                     temp = (Sl(i)-xsum)/S(i)
                  ENDIF
                  IF ( temp<stplen ) THEN
                     stplen = temp
                     iact = i
                  ENDIF
               ENDIF
            ENDDO
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
            sdec = zero
            IF ( stplen>zero ) THEN
               iterc = iterc + 1
               temp = shs/stepsq
               IF ( iact==0 .AND. temp>zero ) THEN
                  Crvmin = DMIN1(Crvmin,temp)
                  IF ( Crvmin==onemin ) Crvmin = temp
               ENDIF
               ggsav = gredsq
               gredsq = zero
               DO i = 1 , N
                  Gnew(i) = Gnew(i) + stplen*Hs(i)
                  IF ( Xbdi(i)==zero ) gredsq = gredsq + Gnew(i)**2
                  D(i) = D(i) + stplen*S(i)
               ENDDO
               sdec = DMAX1(stplen*(ggsav-half*stplen*shs),zero)
               qred = qred + sdec
            ENDIF
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
            IF ( iact>0 ) THEN
               nact = nact + 1
               Xbdi(iact) = one
               IF ( S(iact)<zero ) Xbdi(iact) = onemin
               delsq = delsq - D(iact)**2
               IF ( delsq>zero ) GOTO 100
               GOTO 650
            ENDIF
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
            IF ( stplen<blen ) THEN
               IF ( iterc==itermax ) GOTO 500
               IF ( sdec<=0.01D0*qred ) GOTO 500
               beta = gredsq/ggsav
               GOTO 200
            ENDIF
         ENDIF
 650     Crvmin = zero
         GOTO 300
      ELSEIF ( iterc>itcsav ) THEN
         shs = zero
         dhs = zero
         dhd = zero
         DO i = 1 , N
            IF ( Xbdi(i)==zero ) THEN
               shs = shs + S(i)*Hs(i)
               dhs = dhs + D(i)*Hs(i)
               dhd = dhd + D(i)*Hred(i)
            ENDIF
         ENDDO
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
         iu = INT(17.0D0*angbd+3.1D0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = 1 , iu
            angt = angbd*DFLOAT(i)/DFLOAT(iu)
            sth = (angt+angt)/(one+angt*angt)
            temp = shs + angt*(angt*dhd-dhs-dhs)
            rednew = sth*(angt*dredg-sredg-half*sth*temp)
            IF ( rednew>redmax ) THEN
               redmax = rednew
               isav = i
               rdprev = redsav
            ELSEIF ( i==isav+1 ) THEN
               rdnext = rednew
            ENDIF
            redsav = rednew
         ENDDO
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
         IF ( isav==0 ) GOTO 500
         IF ( isav<iu ) THEN
            temp = (rdnext-rdprev)/(redmax+redmax-rdprev-rdnext)
            angt = angbd*(DFLOAT(isav)+half*temp)/DFLOAT(iu)
         ENDIF
         cth = (one-angt*angt)/(one+angt*angt)
         sth = (angt+angt)/(one+angt*angt)
         temp = shs + angt*(angt*dhd-dhs-dhs)
         sdec = sth*(angt*dredg-sredg-half*sth*temp)
         IF ( sdec<=zero ) GOTO 500
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
         dredg = zero
         gredsq = zero
         DO i = 1 , N
            Gnew(i) = Gnew(i) + (cth-one)*Hred(i) + sth*Hs(i)
            IF ( Xbdi(i)==zero ) THEN
               D(i) = cth*D(i) + sth*S(i)
               dredg = dredg + D(i)*Gnew(i)
               gredsq = gredsq + Gnew(i)**2
            ENDIF
            Hred(i) = cth*Hred(i) + sth*Hs(i)
         ENDDO
         qred = qred + sdec
         IF ( iact>0 .AND. isav==iu ) THEN
            nact = nact + 1
            Xbdi(iact) = xsav
            GOTO 300
         ENDIF
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
         IF ( sdec<=0.01D0*qred ) GOTO 500
         GOTO 400
      ELSE
         DO i = 1 , N
            Hred(i) = Hs(i)
         ENDDO
         GOTO 400
      ENDIF
      END SUBROUTINE TRSBOX
