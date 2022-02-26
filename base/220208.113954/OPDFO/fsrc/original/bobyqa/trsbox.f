      SUBROUTINE TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,
     1  XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),
     1  XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
C
C     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
C       meanings as the corresponding arguments of BOBYQB.
C     DELTA is the trust region radius for the present calculation, which
C       seeks a small value of the quadratic model within distance DELTA of
C       XOPT subject to the bounds on the variables.
C     XNEW will be set to a new vector of variables that is approximately
C       the one that minimizes the quadratic model within the trust region
C       subject to the SL and SU constraints on the variables. It satisfies
C       as equations the bounds that become active during the calculation.
C     D is the calculated trial step from XOPT, generated iteratively from an
C       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
C     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
C       when D is updated.
C     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
C       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
C       I-th variable has become fixed at a bound, the bound being SL(I) or
C       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
C       information is accumulated during the construction of XNEW.
C     The arrays S, HS and HRED are also used for working space. They hold the
C       current search direction, and the changes in the gradient of Q along S
C       and the reduced D, respectively, where the reduced D is the same as D,
C       except that the components of the fixed variables are zero.
C     DSQ will be set to the square of the length of XNEW-XOPT.
C     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
C       it is set to the least curvature of H that occurs in the conjugate
C       gradient searches that are not restricted by any constraints. The
C       value CRVMIN=-1.0D0 is set, however, if all of these searches are
C       constrained.
C
C     A version of the truncated conjugate gradient is applied. If a line
C     search is restricted by a constraint, then the procedure is restarted,
C     the values of the variables that are at their bounds being fixed. If
C     the trust region boundary is reached, then further changes may be made
C     to D, each one being in the two dimensional space that is spanned
C     by the current D and the gradient of Q at XOPT+D, staying on the trust
C     region boundary. Termination occurs when the reduction in Q seems to
C     be close to the greatest reduction that can be achieved.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ONEMIN=-1.0D0
      ZERO=0.0D0
C
C     The sign of GOPT(I) gives the sign of the change to the I-th variable
C     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
C     or not to fix the I-th variable at one of its bounds initially, with
C     NACT being set to the number of fixed variables. D and GNEW are also
C     set for the first iteration. DELSQ is the upper bound on the sum of
C     squares of the free variables. QRED is the reduction in Q so far.
C
      ITERC=0
      NACT=0
      SQSTP=ZERO
      DO 10 I=1,N
      XBDI(I)=ZERO
      IF (XOPT(I) .LE. SL(I)) THEN
          IF (GOPT(I) .GE. ZERO) XBDI(I)=ONEMIN
      ELSE IF (XOPT(I) .GE. SU(I)) THEN
          IF (GOPT(I) .LE. ZERO) XBDI(I)=ONE
      END IF
      IF (XBDI(I) .NE. ZERO) NACT=NACT+1
      D(I)=ZERO
   10 GNEW(I)=GOPT(I)
      DELSQ=DELTA*DELTA
      QRED=ZERO
      CRVMIN=ONEMIN
C
C     Set the next search direction of the conjugate gradient method. It is
C     the steepest descent direction initially and when the iterations are
C     restarted because a variable has just been fixed by a bound, and of
C     course the components of the fixed variables are zero. ITERMAX is an
C     upper bound on the indices of the conjugate gradient iterations.
C
   20 BETA=ZERO
   30 STEPSQ=ZERO
      DO 40 I=1,N
      IF (XBDI(I) .NE. ZERO) THEN
          S(I)=ZERO
      ELSE IF (BETA .EQ. ZERO) THEN
          S(I)=-GNEW(I)
      ELSE
          S(I)=BETA*S(I)-GNEW(I)
      END IF
   40 STEPSQ=STEPSQ+S(I)**2
      IF (STEPSQ .EQ. ZERO) GOTO 190
      IF (BETA .EQ. ZERO) THEN
          GREDSQ=STEPSQ
          ITERMAX=ITERC+N-NACT
      END IF
      IF (GREDSQ*DELSQ .LE. 1.0D-4*QRED*QRED) GO TO 190
C
C     Multiply the search direction by the second derivative matrix of Q and
C     calculate some scalars for the choice of steplength. Then set BLEN to
C     the length of the the step to the trust region boundary and STPLEN to
C     the steplength, ignoring the simple bounds.
C
      GOTO 210
   50 RESID=DELSQ
      DS=ZERO
      SHS=ZERO
      DO 60 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          RESID=RESID-D(I)**2
          DS=DS+S(I)*D(I)
          SHS=SHS+S(I)*HS(I)
      END IF
   60 CONTINUE
      IF (RESID .LE. ZERO) GOTO 90
      TEMP=DSQRT(STEPSQ*RESID+DS*DS)
      IF (DS .LT. ZERO) THEN
          BLEN=(TEMP-DS)/STEPSQ
      ELSE
          BLEN=RESID/(TEMP+DS)
      END IF
      STPLEN=BLEN
      IF (SHS .GT. ZERO) THEN
          STPLEN=DMIN1(BLEN,GREDSQ/SHS)
      END IF
      
C
C     Reduce STPLEN if necessary in order to preserve the simple bounds,
C     letting IACT be the index of the new constrained variable.
C
      IACT=0
      DO 70 I=1,N
      IF (S(I) .NE. ZERO) THEN
          XSUM=XOPT(I)+D(I)
          IF (S(I) .GT. ZERO) THEN
              TEMP=(SU(I)-XSUM)/S(I)
          ELSE
              TEMP=(SL(I)-XSUM)/S(I)
          END IF
          IF (TEMP .LT. STPLEN) THEN
              STPLEN=TEMP
              IACT=I
          END IF
      END IF
   70 CONTINUE
C
C     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
C
      SDEC=ZERO
      IF (STPLEN .GT. ZERO) THEN
          ITERC=ITERC+1
          TEMP=SHS/STEPSQ
          IF (IACT .EQ. 0 .AND. TEMP .GT. ZERO) THEN
              CRVMIN=DMIN1(CRVMIN,TEMP)
              IF (CRVMIN .EQ. ONEMIN) CRVMIN=TEMP
          END IF 
          GGSAV=GREDSQ
          GREDSQ=ZERO
          DO 80 I=1,N
          GNEW(I)=GNEW(I)+STPLEN*HS(I)
          IF (XBDI(I) .EQ. ZERO) GREDSQ=GREDSQ+GNEW(I)**2
   80     D(I)=D(I)+STPLEN*S(I)
          SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
          QRED=QRED+SDEC
      END IF
C
C     Restart the conjugate gradient method if it has hit a new bound.
C
      IF (IACT .GT. 0) THEN
          NACT=NACT+1
          XBDI(IACT)=ONE
          IF (S(IACT) .LT. ZERO) XBDI(IACT)=ONEMIN
          DELSQ=DELSQ-D(IACT)**2
          IF (DELSQ .LE. ZERO) GOTO 90
          GOTO 20
      END IF
C
C     If STPLEN is less than BLEN, then either apply another conjugate
C     gradient iteration or RETURN.
C
      IF (STPLEN .LT. BLEN) THEN
          IF (ITERC .EQ. ITERMAX) GOTO 190
          IF (SDEC .LE. 0.01D0*QRED) GOTO 190
          BETA=GREDSQ/GGSAV
          GOTO 30
      END IF
   90 CRVMIN=ZERO
C
C     Prepare for the alternative iteration by calculating some scalars and
C     by multiplying the reduced D by the second derivative matrix of Q.
C
  100 IF (NACT .GE. N-1) GOTO 190
      DREDSQ=ZERO
      DREDG=ZERO
      GREDSQ=ZERO
      DO 110 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          DREDSQ=DREDSQ+D(I)**2
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
          S(I)=D(I)
      ELSE
          S(I)=ZERO
      END IF
  110 CONTINUE
      ITCSAV=ITERC
      GOTO 210
C
C     Let the search direction S be a linear combination of the reduced D
C     and the reduced G that is orthogonal to the reduced D.
C
  120 ITERC=ITERC+1
      TEMP=GREDSQ*DREDSQ-DREDG*DREDG
      IF (TEMP .LE. 1.0D-4*QRED*QRED) GOTO 190
      TEMP=DSQRT(TEMP)
      DO 130 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
      ELSE
          S(I)=ZERO
      END IF
  130 CONTINUE
      SREDG=-TEMP
C
C     By considering the simple bounds on the variables, calculate an upper
C     bound on the tangent of half the angle of the alternative iteration,
C     namely ANGBD, except that, if already a free variable has reached a
C     bound, there is a branch back to label 100 after fixing that variable.
C
      ANGBD=ONE
      IACT=0
      DO 140 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          TEMPA=XOPT(I)+D(I)-SL(I)
          TEMPB=SU(I)-XOPT(I)-D(I)
          IF (TEMPA .LE. ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONEMIN
              GOTO 100
          ELSE IF (TEMPB .LE. ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONE
              GOTO 100
          END IF
          RATIO=ONE
          SSQ=D(I)**2+S(I)**2
          TEMP=SSQ-(XOPT(I)-SL(I))**2
          IF (TEMP .GT. ZERO) THEN
              TEMP=DSQRT(TEMP)-S(I)
              IF (ANGBD*TEMP .GT. TEMPA) THEN
                  ANGBD=TEMPA/TEMP
                  IACT=I
                  XSAV=ONEMIN
              END IF
          END IF
          TEMP=SSQ-(SU(I)-XOPT(I))**2
          IF (TEMP .GT. ZERO) THEN
              TEMP=DSQRT(TEMP)+S(I)
              IF (ANGBD*TEMP .GT. TEMPB) THEN
                  ANGBD=TEMPB/TEMP
                  IACT=I
                  XSAV=ONE
              END IF
          END IF
      END IF
  140 CONTINUE
C
C     Calculate HHD and some curvatures for the alternative iteration.
C
      GOTO 210
  150 SHS=ZERO
      DHS=ZERO
      DHD=ZERO
      DO 160 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          SHS=SHS+S(I)*HS(I)
          DHS=DHS+D(I)*HS(I)
          DHD=DHD+D(I)*HRED(I)
      END IF
  160 CONTINUE
C
C     Seek the greatest reduction in Q for a range of equally spaced values
C     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
C     the alternative iteration.
C
      REDMAX=ZERO
      ISAV=0
      REDSAV=ZERO
      IU=17.0D0*ANGBD+3.1D0
      DO 170 I=1,IU
      ANGT=ANGBD*DFLOAT(I)/DFLOAT(IU)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (REDNEW .GT. REDMAX) THEN
          REDMAX=REDNEW
          ISAV=I
          RDPREV=REDSAV
      ELSE IF (I .EQ. ISAV+1) THEN
          RDNEXT=REDNEW
      END IF
  170 REDSAV=REDNEW
C
C     Return if the reduction is zero. Otherwise, set the sine and cosine
C     of the angle of the alternative iteration, and calculate SDEC.
C
      IF (ISAV .EQ. 0) GOTO 190
      IF (ISAV .LT. IU) THEN
          TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
          ANGT=ANGBD*(DFLOAT(ISAV)+HALF*TEMP)/DFLOAT(IU)
      END IF
      CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (SDEC .LE. ZERO) GOTO 190
C
C     Update GNEW, D and HRED. If the angle of the alternative iteration
C     is restricted by a bound on a free variable, that variable is fixed
C     at the bound.
C
      DREDG=ZERO
      GREDSQ=ZERO
      DO 180 I=1,N
      GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
      IF (XBDI(I) .EQ. ZERO) THEN
          D(I)=CTH*D(I)+STH*S(I)
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
      END IF
  180 HRED(I)=CTH*HRED(I)+STH*HS(I)
      QRED=QRED+SDEC
      IF (IACT .GT. 0 .AND. ISAV .EQ. IU) THEN
          NACT=NACT+1
          XBDI(IACT)=XSAV
          GOTO 100
      END IF
C
C     If SDEC is sufficiently small, then RETURN after setting XNEW to
C     XOPT+D, giving careful attention to the bounds.
C
      IF (SDEC .GT. 0.01D0*QRED) GOTO 120
  190 DSQ=ZERO
      DO 200 I=1,N
      XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
      IF (XBDI(I) .EQ. ONEMIN) XNEW(I)=SL(I)
      IF (XBDI(I) .EQ. ONE) XNEW(I)=SU(I)
      D(I)=XNEW(I)-XOPT(I)
  200 DSQ=DSQ+D(I)**2
      RETURN
 
C     The following instructions multiply the current S-vector by the second
C     derivative matrix of the quadratic model, putting the product in HS.
C     They are reached from three different parts of the software above and
C     they can be regarded as an external subroutine.
C
  210 IH=0
      DO 220 J=1,N
      HS(J)=ZERO
      DO 220 I=1,J
      IH=IH+1
      IF (I .LT. J) HS(J)=HS(J)+HQ(IH)*S(I)
  220 HS(I)=HS(I)+HQ(IH)*S(J)
      DO 250 K=1,NPT
      IF (PQ(K) .NE. ZERO) THEN
          TEMP=ZERO
          DO 230 J=1,N
  230     TEMP=TEMP+XPT(K,J)*S(J)
          TEMP=TEMP*PQ(K)
          DO 240 I=1,N
  240     HS(I)=HS(I)+TEMP*XPT(K,I)
      END IF
  250 CONTINUE
      IF (CRVMIN .NE. ZERO) GOTO 50
      IF (ITERC .GT. ITCSAV) GOTO 150
      DO 260 I=1,N
  260 HRED(I)=HS(I)
      GOTO 120
      END
