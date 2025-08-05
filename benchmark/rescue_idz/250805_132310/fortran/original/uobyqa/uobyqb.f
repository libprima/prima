      SUBROUTINE UOBYQB (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,NPT,XBASE,
     1  XOPT,XNEW,XPT,PQ,PL,H,G,D,VLAG,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),PQ(*),
     1  PL(NPT,*),H(N,*),G(*),D(*),VLAG(*),W(*)
C
C     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
C       the corresponding arguments in SUBROUTINE UOBYQA.
C     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
C     XBASE will contain a shift of origin that reduces the contributions from
C       rounding errors to values of the model and Lagrange functions.
C     XOPT will be set to the displacement from XBASE of the vector of
C       variables that provides the least calculated F so far.
C     XNEW will be set to the displacement from XBASE of the vector of
C       variables for the current calculation of F.
C     XPT will contain the interpolation point coordinates relative to XBASE.
C     PQ will contain the parameters of the quadratic model.
C     PL will contain the parameters of the Lagrange functions.
C     H will provide the second derivatives that TRSTEP and LAGMAX require.
C     G will provide the first derivatives that TRSTEP and LAGMAX require.
C     D is reserved for trial steps from XOPT, except that it will contain
C       diagonal second derivatives during the initialization procedure.
C     VLAG will contain the values of the Lagrange functions at a new point X.
C     The array W will be used for working space. Its length must be at least
C     max [ 6*N, ( N**2 + 3*N + 2 ) / 2 ].
C
C     Set some constants.
C
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      HALF=0.5D0
      TOL=0.01D0
      NNP=N+N+1
      NPTM=NPT-1
      NFTEST=MAX0(MAXFUN,1)
C
C     Initialization. NF is the number of function calculations so far.
C
      RHO=RHOBEG
      RHOSQ=RHO*RHO
      NF=0
      DO 10 I=1,N
      XBASE(I)=X(I)
      DO 10 K=1,NPT
   10 XPT(K,I)=ZERO
      DO 20 K=1,NPT
      DO 20 J=1,NPTM
   20 PL(K,J)=ZERO
C
C     The branch to label 120 obtains a new value of the objective function
C     and then there is a branch back to label 50, because the new function
C     value is needed to form the initial quadratic model. The least function
C     value so far and its index are noted below.
C
   30 DO 40 I=1,N
   40 X(I)=XBASE(I)+XPT(NF+1,I)
      GOTO 120
   50 IF (NF .EQ. 1) THEN
          FOPT=F
          KOPT=NF
          FBASE=F
          J=0
          JSWITCH=-1
          IH=N
      ELSE
          IF (F .LT. FOPT) THEN
              FOPT=F
              KOPT=NF
          END IF
      END IF
C
C     Form the gradient and diagonal second derivatives of the initial
C     quadratic model and Lagrange functions.
C
      IF (NF .LE. NNP) THEN
          JSWITCH=-JSWITCH
          IF (JSWITCH .GT. 0) THEN
              IF (J .GE. 1) THEN
                  IH=IH+J
                  IF (W(J) .LT. ZERO) THEN
                      D(J)=(FSAVE+F-TWO*FBASE)/RHOSQ
                      PQ(J)=(FSAVE-F)/(TWO*RHO)
                      PL(1,IH)=-TWO/RHOSQ
                      PL(NF-1,J)=HALF/RHO
                      PL(NF-1,IH)=ONE/RHOSQ
                  ELSE
                      PQ(J)=(4.0D0*FSAVE-3.0D0*FBASE-F)/(TWO*RHO)
                      D(J)=(FBASE+F-TWO*FSAVE)/RHOSQ
                      PL(1,J)=-1.5D0/RHO
                      PL(1,IH)=ONE/RHOSQ
                      PL(NF-1,J)=TWO/RHO
                      PL(NF-1,IH)=-TWO/RHOSQ
                  END IF
                  PQ(IH)=D(J)
                  PL(NF,J)=-HALF/RHO
                  PL(NF,IH)=ONE/RHOSQ
              END IF
C
C     Pick the shift from XBASE to the next initial interpolation point
C     that provides diagonal second derivatives.
C
              IF (J .LT. N) THEN
                  J=J+1
                  XPT(NF+1,J)=RHO
              END IF
          ELSE
              FSAVE=F
              IF (F .LT. FBASE) THEN
                  W(J)=RHO
                  XPT(NF+1,J)=TWO*RHO
              ELSE
                  W(J)=-RHO
                  XPT(NF+1,J)=-RHO
              END IF
          END IF
          IF (NF .LT. NNP) GOTO 30
C
C     Form the off-diagonal second derivatives of the initial quadratic model.
C
          IH=N
          IP=1
          IQ=2
      END IF
      IH=IH+1
      IF (NF .GT. NNP) THEN
          TEMP=ONE/(W(IP)*W(IQ))
          TEMPA=F-FBASE-W(IP)*PQ(IP)-W(IQ)*PQ(IQ)
          PQ(IH)=(TEMPA-HALF*RHOSQ*(D(IP)+D(IQ)))*TEMP
          PL(1,IH)=TEMP
          IW=IP+IP
          IF (W(IP) .LT. ZERO) IW=IW+1
          PL(IW,IH)=-TEMP
          IW=IQ+IQ
          IF (W(IQ) .LT. ZERO) IW=IW+1
          PL(IW,IH)=-TEMP
          PL(NF,IH)=TEMP
C
C     Pick the shift from XBASE to the next initial interpolation point
C     that provides off-diagonal second derivatives.
C
          IP=IP+1
      END IF
      IF (IP .EQ. IQ) THEN
          IH=IH+1
          IP=1
          IQ=IQ+1
      END IF
      IF (NF .LT. NPT) THEN
          XPT(NF+1,IP)=W(IP)
          XPT(NF+1,IQ)=W(IQ)
          GOTO 30
      END IF
C
C     Set parameters to begin the iterations for the current RHO.
C
      SIXTHM=ZERO
      DELTA=RHO
   60 TWORSQ=(TWO*RHO)**2
      RHOSQ=RHO*RHO
C
C     Form the gradient of the quadratic model at the trust region centre.
C
   70 KNEW=0
      IH=N
      DO 80 J=1,N
      XOPT(J)=XPT(KOPT,J)
      G(J)=PQ(J)
      DO 80 I=1,J
      IH=IH+1
      G(I)=G(I)+PQ(IH)*XOPT(J)
      IF (I .LT. J) G(J)=G(J)+PQ(IH)*XOPT(I)
   80 H(I,J)=PQ(IH)
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve conditioning,
C     and also calculate a lower bound on the Hessian term of the model Q.
C
      CALL TRSTEP (N,G,H,DELTA,TOL,D,W(1),W(N+1),W(2*N+1),W(3*N+1),
     1  W(4*N+1),W(5*N+1),EVALUE)
      TEMP=ZERO
      DO 90 I=1,N
   90 TEMP=TEMP+D(I)**2
      DNORM=DMIN1(DELTA,DSQRT(TEMP))
      ERRTOL=-ONE
      IF (DNORM .LT. HALF*RHO) THEN
          KNEW=-1
          ERRTOL=HALF*EVALUE*RHO*RHO
          IF (NF .LE. NPT+9) ERRTOL=ZERO
          GOTO 290
      END IF
C
C     Calculate the next value of the objective function.
C
  100 DO 110 I=1,N
      XNEW(I)=XOPT(I)+D(I)
  110 X(I)=XBASE(I)+XNEW(I)
  120 IF (NF .GE. NFTEST) THEN
          IF (IPRINT .GT. 0) PRINT 130
  130     FORMAT (/4X,'Return from UOBYQA because CALFUN has been',
     1      ' called MAXFUN times')
          GOTO 420
      END IF
      NF=NF+1
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 140, NF,F,(X(I),I=1,N)
  140      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .LE. NPT) GOTO 50
      IF (KNEW .EQ. -1) GOTO 420
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and find the values of the Lagrange functions at the new point.
C
      VQUAD=ZERO
      IH=N
      DO 150 J=1,N
      W(J)=D(J)
      VQUAD=VQUAD+W(J)*PQ(J)
      DO 150 I=1,J
      IH=IH+1
      W(IH)=D(I)*XNEW(J)+D(J)*XOPT(I)
      IF (I .EQ. J) W(IH)=HALF*W(IH)
  150 VQUAD=VQUAD+W(IH)*PQ(IH)
      DO 170 K=1,NPT
      TEMP=ZERO
      DO 160 J=1,NPTM
  160 TEMP=TEMP+W(J)*PL(K,J)
  170 VLAG(K)=TEMP
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     Update SIXTHM, which is a lower bound on one sixth of the greatest
C     third derivative of F.
C
      DIFF=F-FOPT-VQUAD
      SUM=ZERO
      DO 190 K=1,NPT
      TEMP=ZERO
      DO 180 I=1,N
  180 TEMP=TEMP+(XPT(K,I)-XNEW(I))**2
      TEMP=DSQRT(TEMP)
  190 SUM=SUM+DABS(TEMP*TEMP*TEMP*VLAG(K))
      SIXTHM=DMAX1(SIXTHM,DABS(DIFF)/SUM)
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. Then branch if D is not a trust region step.
C
      FSAVE=FOPT
      IF (F .LT. FOPT) THEN
          FOPT=F
          DO 200 I=1,N
  200     XOPT(I)=XNEW(I)
      END IF
      KSAVE=KNEW
      IF (KNEW .GT. 0) GOTO 240
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD .GE. ZERO) THEN
          IF (IPRINT .GT. 0) PRINT 210
  210     FORMAT (/4X,'Return from UOBYQA because a trust',
     1      ' region step has failed to reduce Q')
          GOTO 420
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO .LE. 0.1D0) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(DELTA,1.25D0*DNORM,DNORM+RHO)
      END IF
      IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      KTEMP=0
      DETRAT=ZERO
      IF (F .GE. FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO 230 K=1,NPT
      SUM=ZERO
      DO 220 I=1,N
  220 SUM=SUM+(XPT(K,I)-XOPT(I))**2
      TEMP=DABS(VLAG(K))
      IF (SUM .GT. RHOSQ) TEMP=TEMP*(SUM/RHOSQ)**1.5D0
      IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
          DETRAT=TEMP
          DDKNEW=SUM
          KNEW=K
      END IF
  230 CONTINUE
      IF (KNEW .EQ. 0) GOTO 290
C
C     Replace the interpolation point that has index KNEW by the point XNEW,
C     and also update the Lagrange functions and the quadratic model.
C
  240 DO 250 I=1,N
  250 XPT(KNEW,I)=XNEW(I)
      TEMP=ONE/VLAG(KNEW)
      DO 260 J=1,NPTM
      PL(KNEW,J)=TEMP*PL(KNEW,J)
  260 PQ(J)=PQ(J)+DIFF*PL(KNEW,J)
      DO 280 K=1,NPT
      IF (K .NE. KNEW) THEN
          TEMP=VLAG(K)
          DO 270 J=1,NPTM
  270     PL(K,J)=PL(K,J)-TEMP*PL(KNEW,J)
      END IF
  280 CONTINUE
C
C     Update KOPT if F is the least calculated value of the objective
C     function. Then branch for another trust region calculation. The
C     case KSAVE>0 indicates that a model step has just been taken.
C
      IF (F .LT. FSAVE) THEN
          KOPT=KNEW
          GOTO 70
      END IF
      IF (KSAVE .GT. 0) GOTO 70
      IF (DNORM .GT. TWO*RHO) GOTO 70
      IF (DDKNEW .GT. TWORSQ) GOTO 70
C
C     Alternatively, find out if the interpolation points are close
C     enough to the best point so far.
C
  290 DO 300 K=1,NPT
      W(K)=ZERO
      DO 300 I=1,N
  300 W(K)=W(K)+(XPT(K,I)-XOPT(I))**2
  310 KNEW=-1
      DISTEST=TWORSQ
      DO 320 K=1,NPT
      IF (W(K) .GT. DISTEST) THEN
          KNEW=K
          DISTEST=W(K)
      END IF
  320 CONTINUE
C
C     If a point is sufficiently far away, then set the gradient and Hessian
C     of its Lagrange function at the centre of the trust region, and find
C     half the sum of squares of components of the Hessian.
C
      IF (KNEW .GT. 0) THEN
          IH=N
          SUMH=ZERO
          DO 340 J=1,N
          G(J)=PL(KNEW,J)
          DO 330 I=1,J
          IH=IH+1
          TEMP=PL(KNEW,IH)
          G(J)=G(J)+TEMP*XOPT(I)
          IF (I .LT. J) THEN
              G(I)=G(I)+TEMP*XOPT(J)
              SUMH=SUMH+TEMP*TEMP
          END IF
  330     H(I,J)=TEMP
  340     SUMH=SUMH+HALF*TEMP*TEMP
C
C     If ERRTOL is positive, test whether to replace the interpolation point
C     with index KNEW, using a bound on the maximum modulus of its Lagrange
C     function in the trust region.
C
          IF (ERRTOL .GT. ZERO) THEN
              W(KNEW)=ZERO
              SUMG=ZERO
              DO 350 I=1,N
  350         SUMG=SUMG+G(I)**2
              ESTIM=RHO*(DSQRT(SUMG)+RHO*DSQRT(HALF*SUMH))
              WMULT=SIXTHM*DISTEST**1.5D0
              IF (WMULT*ESTIM .LE. ERRTOL) GOTO 310
          END IF
C
C     If the KNEW-th point may be replaced, then pick a D that gives a large
C     value of the modulus of its Lagrange function within the trust region.
C     Here the vector XNEW is used as temporary working space.
C
          CALL LAGMAX (N,G,H,RHO,D,XNEW,VMAX)
          IF (ERRTOL .GT. ZERO) THEN
              IF (WMULT*VMAX .LE. ERRTOL) GOTO 310
          END IF
          GOTO 100
      END IF
      IF (DNORM .GT. RHO) GOTO 70
C
C     Prepare to reduce RHO by shifting XBASE to the best point so far,
C     and make the corresponding changes to the gradients of the Lagrange
C     functions and the quadratic model.
C
      IF (RHO .GT. RHOEND) THEN
          IH=N
          DO 380 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
          DO 360 K=1,NPT
  360     XPT(K,J)=XPT(K,J)-XOPT(J)
          DO 380 I=1,J
          IH=IH+1
          PQ(I)=PQ(I)+PQ(IH)*XOPT(J)
          IF (I .LT. J) THEN
              PQ(J)=PQ(J)+PQ(IH)*XOPT(I)
              DO 370 K=1,NPT
  370         PL(K,J)=PL(K,J)+PL(K,IH)*XOPT(I)
          END IF
          DO 380 K=1,NPT
  380     PL(K,I)=PL(K,I)+PL(K,IH)*XOPT(J)
C
C     Pick the next values of RHO and DELTA.
C
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=0.1D0*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 390
  390         FORMAT (5X)
              PRINT 400, RHO,NF
  400         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 410, FOPT,(XBASE(I),I=1,N)
  410         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 60
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C     it is too short to have been tried before.
C
      IF (ERRTOL .GE. ZERO) GOTO 100
  420 IF (FOPT .LE. F) THEN
          DO 430 I=1,N
  430     X(I)=XBASE(I)+XOPT(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
          PRINT 440, NF
  440     FORMAT (/4X,'At the return from UOBYQA',5X,
     1      'Number of function values =',I6)
          PRINT 410, F,(X(I),I=1,N)
      END IF
      RETURN
      END
