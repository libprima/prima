      SUBROUTINE UOBYQB (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,NPT,XBASE,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     1  XOPT,XNEW,XPT,PQ,PL,H,G,D,VLAG,W)
     1  XOPT,XNEW,XPT,PQ,PL,H,G,D,VLAG,W,F,INFO,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      DO I=1,N
          XBASE(I)=X(I)
          DO K=1,NPT
              XPT(K,I)=ZERO
          END DO
      END DO
      DO K=1,NPT
          DO J=1,NPTM
              PL(K,J)=ZERO
          END DO
      END DO
C
C     The branch to label 120 obtains a new value of the objective function
C     and then there is a branch back to label 50, because the new function
C     value is needed to form the initial quadratic model. The least function
C     value so far and its index are noted below.
C
   30 DO I=1,N
          X(I)=XBASE(I)+XPT(NF+1,I)
      END DO
      GOTO 120
   50 IF (NF == 1) THEN
          FOPT=F
          KOPT=NF
          FBASE=F
          J=0
          JSWITCH=-1
          IH=N
      ELSE
          IF (F < FOPT) THEN
              FOPT=F
              KOPT=NF
          END IF
      END IF
C
C     Form the gradient and diagonal second derivatives of the initial
C     quadratic model and Lagrange functions.
C
      IF (NF <= NNP) THEN
          JSWITCH=-JSWITCH
          IF (JSWITCH > 0) THEN
              IF (J >= 1) THEN
                  IH=IH+J
                  IF (W(J) < ZERO) THEN
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
              IF (J < N) THEN
                  J=J+1
                  XPT(NF+1,J)=RHO
              END IF
          ELSE
              FSAVE=F
              IF (F < FBASE) THEN
                  W(J)=RHO
                  XPT(NF+1,J)=TWO*RHO
              ELSE
                  W(J)=-RHO
                  XPT(NF+1,J)=-RHO
              END IF
          END IF
          IF (NF < NNP) GOTO 30
C
C     Form the off-diagonal second derivatives of the initial quadratic model.
C
          IH=N
          IP=1
          IQ=2
      END IF
      IH=IH+1
      IF (NF > NNP) THEN
          TEMP=ONE/(W(IP)*W(IQ))
          TEMPA=F-FBASE-W(IP)*PQ(IP)-W(IQ)*PQ(IQ)
          PQ(IH)=(TEMPA-HALF*RHOSQ*(D(IP)+D(IQ)))*TEMP
          PL(1,IH)=TEMP
          IW=IP+IP
          IF (W(IP) < ZERO) IW=IW+1
          PL(IW,IH)=-TEMP
          IW=IQ+IQ
          IF (W(IQ) < ZERO) IW=IW+1
          PL(IW,IH)=-TEMP
          PL(NF,IH)=TEMP
C
C     Pick the shift from XBASE to the next initial interpolation point
C     that provides off-diagonal second derivatives.
C
          IP=IP+1
      END IF
      IF (IP == IQ) THEN
          IH=IH+1
          IP=1
          IQ=IQ+1
      END IF
      IF (NF < NPT) THEN
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
      DO J=1,N
          XOPT(J)=XPT(KOPT,J)
          G(J)=PQ(J)
          DO I=1,J
              IH=IH+1
              G(I)=G(I)+PQ(IH)*XOPT(J)
              IF (I < J) G(J)=G(J)+PQ(IH)*XOPT(I)
              H(I,J)=PQ(IH)
          END DO
      END DO
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve conditioning,
C     and also calculate a lower bound on the Hessian term of the model Q.
C
      CALL TRSTEP (N,G,H,DELTA,TOL,D,W(1),W(N+1),W(2*N+1),W(3*N+1),
     1  W(4*N+1),W(5*N+1),EVALUE)
      TEMP=ZERO
      DO I=1,N
          TEMP=TEMP+D(I)**2
      END DO
      DNORM=DMIN1(DELTA,DSQRT(TEMP))
      ERRTOL=-ONE
      IF (DNORM < HALF*RHO) THEN
          KNEW=-1
          ERRTOL=HALF*EVALUE*RHO*RHO
          IF (NF <= NPT+9) ERRTOL=ZERO
          GOTO 290
      END IF
C
C     Calculate the next value of the objective function.
C
  100 DO I=1,N
          XNEW(I)=XOPT(I)+D(I)
          X(I)=XBASE(I)+XNEW(I)
      END DO
  120 IF (NF >= NFTEST) THEN
          IF (IPRINT > 0) PRINT 130
  130     FORMAT (/4X,'Return from UOBYQA because CALFUN has been',
     1      ' called MAXFUN times')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
          INFO=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
          GOTO 420
      END IF
      NF=NF+1
      CALL CALFUN (N,X,F)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C     By Zaikun (commented on 02-06-2019; implemented in 2016):
C     Exit if F .LE. FTARGET.
      IF (F <= FTARGET) THEN
          INFO=1
          GOTO 436
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      IF (IPRINT == 3) THEN
          PRINT 140, NF,F,(X(I),I=1,N)
  140      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF <= NPT) GOTO 50
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      IF (KNEW .EQ. -1) GOTO 420
      IF (KNEW == -1) THEN
          INFO=0
          GOTO 420
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and find the values of the Lagrange functions at the new point.
C
      VQUAD=ZERO
      IH=N
      DO J=1,N
          W(J)=D(J)
          VQUAD=VQUAD+W(J)*PQ(J)
          DO I=1,J
              IH=IH+1
              W(IH)=D(I)*XNEW(J)+D(J)*XOPT(I)
              IF (I == J) W(IH)=HALF*W(IH)
              VQUAD=VQUAD+W(IH)*PQ(IH)
          END DO
      END DO
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,NPTM
              TEMP=TEMP+W(J)*PL(K,J)
          END DO
          VLAG(K)=TEMP
      END DO
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     Update SIXTHM, which is a lower bound on one sixth of the greatest
C     third derivative of F.
C
      DIFF=F-FOPT-VQUAD
      SUM=ZERO
      DO K=1,NPT
          TEMP=ZERO
          DO I=1,N
              TEMP=TEMP+(XPT(K,I)-XNEW(I))**2
          END DO
          TEMP=DSQRT(TEMP)
          SUM=SUM+DABS(TEMP*TEMP*TEMP*VLAG(K))
      END DO
      SIXTHM=DMAX1(SIXTHM,DABS(DIFF)/SUM)
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. Then branch if D is not a trust region step.
C
      FSAVE=FOPT
      IF (F < FOPT) THEN
          FOPT=F
          DO I=1,N
              XOPT(I)=XNEW(I)
          END DO
      END IF
      KSAVE=KNEW
      IF (KNEW > 0) GOTO 240
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD >= ZERO) THEN
          IF (IPRINT > 0) PRINT 210
  210     FORMAT (/4X,'Return from UOBYQA because a trust',
     1      ' region step has failed to reduce Q')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 420
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO <= 0.1D0) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(DELTA,1.25D0*DNORM,DNORM+RHO)
      END IF
      IF (DELTA <= 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      KTEMP=0
      DETRAT=ZERO
      IF (F >= FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO K=1,NPT
          SUM=ZERO
          DO I=1,N
              SUM=SUM+(XPT(K,I)-XOPT(I))**2
          END DO
          TEMP=DABS(VLAG(K))
          IF (SUM > RHOSQ) TEMP=TEMP*(SUM/RHOSQ)**1.5D0
          IF (TEMP > DETRAT .AND. K /= KTEMP) THEN
              DETRAT=TEMP
              DDKNEW=SUM
              KNEW=K
          END IF
      END DO
      IF (KNEW == 0) GOTO 290
C
C     Replace the interpolation point that has index KNEW by the point XNEW,
C     and also update the Lagrange functions and the quadratic model.
C
  240 DO I=1,N
          XPT(KNEW,I)=XNEW(I)
      END DO
      TEMP=ONE/VLAG(KNEW)
      DO J=1,NPTM
          PL(KNEW,J)=TEMP*PL(KNEW,J)
          PQ(J)=PQ(J)+DIFF*PL(KNEW,J)
      END DO
      DO K=1,NPT
          IF (K /= KNEW) THEN
              TEMP=VLAG(K)
              DO J=1,NPTM
                  PL(K,J)=PL(K,J)-TEMP*PL(KNEW,J)
              END DO
          END IF
      END DO
C
C     Update KOPT if F is the least calculated value of the objective
C     function. Then branch for another trust region calculation. The
C     case KSAVE>0 indicates that a model step has just been taken.
C
      IF (F < FSAVE) THEN
          KOPT=KNEW
          GOTO 70
      END IF
      IF (KSAVE > 0) GOTO 70
      IF (DNORM > TWO*RHO) GOTO 70
      IF (DDKNEW > TWORSQ) GOTO 70
C
C     Alternatively, find out if the interpolation points are close
C     enough to the best point so far.
C
  290 DO K=1,NPT
          W(K)=ZERO
          DO I=1,N
              W(K)=W(K)+(XPT(K,I)-XOPT(I))**2
          END DO
      END DO
  310 KNEW=-1
      DISTEST=TWORSQ
      DO K=1,NPT
          IF (W(K) > DISTEST) THEN
              KNEW=K
              DISTEST=W(K)
          END IF
      END DO
C
C     If a point is sufficiently far away, then set the gradient and Hessian
C     of its Lagrange function at the centre of the trust region, and find
C     half the sum of squares of components of the Hessian.
C
      IF (KNEW > 0) THEN
          IH=N
          SUMH=ZERO
          DO J=1,N
              G(J)=PL(KNEW,J)
              DO I=1,J
                  IH=IH+1
                  TEMP=PL(KNEW,IH)
                  G(J)=G(J)+TEMP*XOPT(I)
                  IF (I < J) THEN
                      G(I)=G(I)+TEMP*XOPT(J)
                      SUMH=SUMH+TEMP*TEMP
                  END IF
                  H(I,J)=TEMP
              END DO
              SUMH=SUMH+HALF*TEMP*TEMP
          END DO
C
C     If ERRTOL is positive, test whether to replace the interpolation point
C     with index KNEW, using a bound on the maximum modulus of its Lagrange
C     function in the trust region.
C
          IF (ERRTOL > ZERO) THEN
              W(KNEW)=ZERO
              SUMG=ZERO
              DO I=1,N
                  SUMG=SUMG+G(I)**2
              END DO
              ESTIM=RHO*(DSQRT(SUMG)+RHO*DSQRT(HALF*SUMH))
              WMULT=SIXTHM*DISTEST**1.5D0
              IF (WMULT*ESTIM <= ERRTOL) GOTO 310
          END IF
C
C     If the KNEW-th point may be replaced, then pick a D that gives a large
C     value of the modulus of its Lagrange function within the trust region.
C     Here the vector XNEW is used as temporary working space.
C
          CALL LAGMAX (N,G,H,RHO,D,XNEW,VMAX)
          IF (ERRTOL > ZERO) THEN
              IF (WMULT*VMAX <= ERRTOL) GOTO 310
          END IF
          GOTO 100
      END IF
      IF (DNORM > RHO) GOTO 70
C
C     Prepare to reduce RHO by shifting XBASE to the best point so far,
C     and make the corresponding changes to the gradients of the Lagrange
C     functions and the quadratic model.
C
      IF (RHO > RHOEND) THEN
          IH=N
          DO J=1,N
              XBASE(J)=XBASE(J)+XOPT(J)
              DO K=1,NPT
                  XPT(K,J)=XPT(K,J)-XOPT(J)
              END DO
              DO I=1,J
                  IH=IH+1
                  PQ(I)=PQ(I)+PQ(IH)*XOPT(J)
                  IF (I < J) THEN
                      PQ(J)=PQ(J)+PQ(IH)*XOPT(I)
                      DO K=1,NPT
                          PL(K,J)=PL(K,J)+PL(K,IH)*XOPT(I)
                      END DO
                  END IF
                  DO K=1,NPT
                      PL(K,I)=PL(K,I)+PL(K,IH)*XOPT(J)
                  END DO
              END DO
          END DO
C
C     Pick the next values of RHO and DELTA.
C
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO <= 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO <= 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=0.1D0*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT >= 2) THEN
              IF (IPRINT >= 3) PRINT 390
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INFO=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ERRTOL >= ZERO) GOTO 100
  420 IF (FOPT <= F) THEN
          DO I=1,N
              X(I)=XBASE(I)+XOPT(I)
          END DO
          F=FOPT
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (IPRINT .GE. 1) THEN
  436 IF (IPRINT >= 1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          PRINT 440, NF
  440     FORMAT (/4X,'At the return from UOBYQA',5X,
     1      'Number of function values =',I6)
          PRINT 410, F,(X(I),I=1,N)
      END IF
      RETURN
      END
