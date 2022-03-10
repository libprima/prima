      SUBROUTINE NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     1  XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W)
     1  XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W,F,INFO,
     2  FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),FVAL(*),
     1  GQ(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),VLAG(*),W(*)
C
C     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
C       to the corresponding arguments in SUBROUTINE NEWUOA.
C     XBASE will hold a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XOPT will be set to the displacement from XBASE of the vector of
C       variables that provides the least calculated F so far.
C     XNEW will be set to the displacement from XBASE of the vector of
C       variables for the current calculation of F.
C     XPT will contain the interpolation point coordinates relative to XBASE.
C     FVAL will hold the values of F at the interpolation points.
C     GQ will hold the gradient of the quadratic model at XBASE.
C     HQ will hold the explicit second derivatives of the quadratic model.
C     PQ will contain the parameters of the implicit second derivatives of
C       the quadratic model.
C     BMAT will hold the last N columns of H.
C     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
C       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
C       the elements of DZ are plus or minus one, as specified by IDZ.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     D is reserved for trial steps from XOPT.
C     VLAG will contain the values of the Lagrange functions at a new point X.
C       They are part of a product that requires VLAG to be of length NDIM.
C     The array W will be used for working space. Its length must be at least
C       10*NDIM = 10*(NPT+N).
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NH=(N*NP)/2
      NPTM=NPT-NP
      NFTEST=MAX0(MAXFUN,1)
C
C     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
C
      DO J=1,N
          XBASE(J)=X(J)
          DO K=1,NPT
              XPT(K,J)=ZERO
          END DO
          DO I=1,NDIM
              BMAT(I,J)=ZERO
          END DO
      END DO
      DO IH=1,NH
          HQ(IH)=ZERO
      END DO
      DO K=1,NPT
          PQ(K)=ZERO
          DO J=1,NPTM
              ZMAT(K,J)=ZERO
          END DO
      END DO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF,.).
C
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      NF=0
   50 NFM=NF
      NFMM=NF-N
      NF=NF+1
      IF (NFM <= 2*N) THEN
          IF (NFM >= 1 .AND. NFM <= N) THEN
              XPT(NF,NFM)=RHOBEG
          ELSE IF (NFM > N) THEN
              XPT(NF,NFMM)=-RHOBEG
          END IF
      ELSE
          ITEMP=(NFMM-1)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT > N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XIPT=RHOBEG
          IF (FVAL(IPT+NP) < FVAL(IPT+1)) XIPT=-XIPT
          XJPT=RHOBEG
          IF (FVAL(JPT+NP) < FVAL(JPT+1)) XJPT=-XJPT
          XPT(NF,IPT)=XIPT
          XPT(NF,JPT)=XJPT
      END IF
C
C     Calculate the next value of F, label 70 being reached immediately
C     after this calculation. The least function value so far and its index
C     are required.
C
      DO J=1,N
          X(J)=XPT(NF,J)+XBASE(J)
      END DO
      GOTO 310
   70 FVAL(NF)=F
      IF (NF == 1) THEN
          FBEG=F
          FOPT=F
          KOPT=1
      ELSE IF (F < FOPT) THEN
          FOPT=F
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in
C     the cases when NF is at most 2*N+1.
C
      IF (NFM <= 2*N) THEN
          IF (NFM >= 1 .AND. NFM <= N) THEN
              GQ(NFM)=(F-FBEG)/RHOBEG
              IF (NPT < NF+N) THEN
                  BMAT(1,NFM)=-ONE/RHOBEG
                  BMAT(NF,NFM)=ONE/RHOBEG
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NFM > N) THEN
              BMAT(NF-N,NFMM)=HALF/RHOBEG
              BMAT(NF,NFMM)=-HALF/RHOBEG
              ZMAT(1,NFMM)=-RECIQ-RECIQ
              ZMAT(NF-N,NFMM)=RECIQ
              ZMAT(NF,NFMM)=RECIQ
              IH=(NFMM*(NFMM+1))/2
              TEMP=(FBEG-F)/RHOBEG
              HQ(IH)=(GQ(NFMM)-TEMP)/RHOBEG
              GQ(NFMM)=HALF*(GQ(NFMM)+TEMP)
          END IF
C
C     Set the off-diagonal second derivatives of the Lagrange functions and
C     the initial quadratic model.
C
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          IF (XIPT < ZERO) IPT=IPT+N
          IF (XJPT < ZERO) JPT=JPT+N
          ZMAT(1,NFMM)=RECIP
          ZMAT(NF,NFMM)=RECIP
          ZMAT(IPT+1,NFMM)=-RECIP
          ZMAT(JPT+1,NFMM)=-RECIP
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/(XIPT*XJPT)
      END IF
      IF (NF < NPT) GOTO 50
C
C     Begin the iterative procedure, because the initial model is complete.
C
      RHO=RHOBEG
      DELTA=RHO
      IDZ=1
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      XOPTSQ=ZERO
      DO I=1,N
          XOPT(I)=XPT(KOPT,I)
          XOPTSQ=XOPTSQ+XOPT(I)**2
      END DO
   90 NFSAV=NF
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve the model.
C
  100 KNEW=0
      CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN)
      DSQ=ZERO
      DO I=1,N
          DSQ=DSQ+D(I)**2
      END DO
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM < HALF*RHO) THEN
          KNEW=-1
          DELTA=TENTH*DELTA
          RATIO=-1.0D0
          IF (DELTA <= 1.5D0*RHO) DELTA=RHO
          IF (NF <= NFSAV+2) GOTO 460
          TEMP=0.125D0*CRVMIN*RHO*RHO
          IF (TEMP <= DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
          GOTO 490
      END IF
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C     to BMAT that do not depend on ZMAT.
C
  120 IF (DSQ <= 1.0D-3*XOPTSQ) THEN
          TEMPQ=0.25D0*XOPTSQ
          DO K=1,NPT
              SUM=ZERO
              DO I=1,N
                  SUM=SUM+XPT(K,I)*XOPT(I)
              END DO
              TEMP=PQ(K)*SUM
              SUM=SUM-HALF*XOPTSQ
              W(NPT+K)=SUM
              DO I=1,N
                  GQ(I)=GQ(I)+TEMP*XPT(K,I)
                  XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
                  VLAG(I)=BMAT(K,I)
                  W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
                  IP=NPT+I
                  DO J=1,I
                      BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
                  END DO
              END DO
          END DO
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO K=1,NPTM
              SUMZ=ZERO
              DO I=1,NPT
                  SUMZ=SUMZ+ZMAT(I,K)
                  W(I)=W(NPT+I)*ZMAT(I,K)
              END DO
              DO J=1,N
                  SUM=TEMPQ*SUMZ*XOPT(J)
                  DO I=1,NPT
                      SUM=SUM+W(I)*XPT(I,J)
                  END DO
                  VLAG(J)=SUM
                  IF (K < IDZ) SUM=-SUM
                  DO I=1,NPT
                      BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
                  END DO
              END DO
              DO I=1,N
                  IP=I+NPT
                  TEMP=VLAG(I)
                  IF (K < IDZ) TEMP=-TEMP
                  DO J=1,I
                      BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
                  END DO
              END DO
          END DO
C
C     The following instructions complete the shift of XBASE, including
C     the changes to the parameters of the quadratic model.
C
          IH=0
          DO J=1,N
              W(J)=ZERO
              DO K=1,NPT
                  W(J)=W(J)+PQ(K)*XPT(K,J)
                  XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
              END DO
              DO I=1,J
                  IH=IH+1
                  IF (I < J) GQ(J)=GQ(J)+HQ(IH)*XOPT(I)
                  GQ(I)=GQ(I)+HQ(IH)*XOPT(J)
                  HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
                  BMAT(NPT+I,J)=BMAT(NPT+J,I)
              END DO
          END DO
          DO J=1,N
              XBASE(J)=XBASE(J)+XOPT(J)
              XOPT(J)=ZERO
          END DO
          XOPTSQ=ZERO
      END IF
C
C     Pick the model step if KNEW is positive. A different choice of D
C     may be made later, if the choice of D by BIGLAG causes substantial
C     cancellation in DENOM.
C
      IF (KNEW > 0) THEN
          CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,DSTEP,
     1      D,ALPHA,VLAG,VLAG(NPT+1),W,W(NP),W(NP+N))
      END IF
C
C     Calculate VLAG and BETA for the current choice of D. The first NPT
C     components of W_check will be held in W.
C
      DO K=1,NPT
          SUMA=ZERO
          SUMB=ZERO
          SUM=ZERO
          DO J=1,N
              SUMA=SUMA+XPT(K,J)*D(J)
              SUMB=SUMB+XPT(K,J)*XOPT(J)
              SUM=SUM+BMAT(K,J)*D(J)
          END DO
          W(K)=SUMA*(HALF*SUMA+SUMB)
          VLAG(K)=SUM
      END DO
      BETA=ZERO
      DO K=1,NPTM
          SUM=ZERO
          DO I=1,NPT
              SUM=SUM+ZMAT(I,K)*W(I)
          END DO
          IF (K < IDZ) THEN
              BETA=BETA+SUM*SUM
              SUM=-SUM
          ELSE
              BETA=BETA-SUM*SUM
          END IF
          DO I=1,NPT
              VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
          END DO
      END DO
      BSUM=ZERO
      DX=ZERO
      DO J=1,N
          SUM=ZERO
          DO I=1,NPT
              SUM=SUM+W(I)*BMAT(I,J)
          END DO
          BSUM=BSUM+SUM*D(J)
          JP=NPT+J
          DO K=1,N
              SUM=SUM+BMAT(JP,K)*D(K)
          END DO
          VLAG(JP)=SUM
          BSUM=BSUM+SUM*D(J)
          DX=DX+D(J)*XOPT(J)
      END DO
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is positive and if the cancellation in DENOM is unacceptable,
C     then BIGDEN calculates an alternative model step, XNEW being used for
C     working space.
C
      IF (KNEW > 0) THEN
          TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
          IF (DABS(TEMP) <= 0.8D0) THEN
              CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1          KNEW,D,W,VLAG,BETA,XNEW,W(NDIM+1),W(6*NDIM+1))
          END IF
      END IF
C
C     Calculate the next value of the objective function.
C
  290 DO I=1,N
          XNEW(I)=XOPT(I)+D(I)
          X(I)=XBASE(I)+XNEW(I)
      END DO
      NF=NF+1
  310 IF (NF > NFTEST) THEN
          NF=NF-1
          IF (IPRINT > 0) PRINT 320
  320     FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
     1      ' called MAXFUN times.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 530
      END IF
      CALL CALFUN (N,X,F)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Zaikun (commented on 02-06-2019; implemented in 2016):
C     Exit if F .LE. FTARGET.
      IF (F <= FTARGET) THEN
          INFO=1
          GOTO 546
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IPRINT == 3) THEN
          PRINT 330, NF,F,(X(I),I=1,N)
  330      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF <= NPT) GOTO 70
      IF (KNEW == -1) GOTO 530
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and set DIFF to the error of this prediction.
C
      VQUAD=ZERO
      IH=0
      DO J=1,N
          VQUAD=VQUAD+D(J)*GQ(J)
          DO I=1,J
              IH=IH+1
              TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
              IF (I == J) TEMP=HALF*TEMP
              VQUAD=VQUAD+TEMP*HQ(IH)
          END DO
      END DO
      DO K=1,NPT
          VQUAD=VQUAD+PQ(K)*W(K)
      END DO
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM > RHO) NFSAV=NF
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. The branch when KNEW is positive occurs if D is not
C     a trust region step.
C
      FSAVE=FOPT
      IF (F < FOPT) THEN
          FOPT=F
          XOPTSQ=ZERO
          DO I=1,N
              XOPT(I)=XNEW(I)
              XOPTSQ=XOPTSQ+XOPT(I)**2
          END DO
      END IF
      KSAVE=KNEW
      IF (KNEW > 0) GOTO 410
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD >= ZERO) THEN
          IF (IPRINT > 0) PRINT 370
  370     FORMAT (/4X,'Return from NEWUOA because a trust',
     1      ' region step has failed to reduce Q.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 530
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO <= TENTH) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO <= 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
      END IF
      IF (DELTA <= 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
      KTEMP=0
      DETRAT=ZERO
      IF (F >= FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO K=1,NPT
          HDIAG=ZERO
          DO J=1,NPTM
              TEMP=ONE
              IF (J < IDZ) TEMP=-ONE
              HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
          END DO
          TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
          DISTSQ=ZERO
          DO J=1,N
              DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
          END DO
          IF (DISTSQ > RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)**3
          IF (TEMP > DETRAT .AND. K /= KTEMP) THEN
              DETRAT=TEMP
              KNEW=K
          END IF
      END DO
      IF (KNEW == 0) GOTO 460
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C     can be moved. Begin the updating of the quadratic model, starting
C     with the explicit second derivative term.
C
  410 CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      FVAL(KNEW)=F
      IH=0
      DO I=1,N
          TEMP=PQ(KNEW)*XPT(KNEW,I)
          DO J=1,I
              IH=IH+1
              HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
          END DO
      END DO
      PQ(KNEW)=ZERO
C
C     Update the other second derivative parameters, and then the gradient
C     vector of the model. Also include the new interpolation point.
C
      DO J=1,NPTM
          TEMP=DIFF*ZMAT(KNEW,J)
          IF (J < IDZ) TEMP=-TEMP
          DO K=1,NPT
              PQ(K)=PQ(K)+TEMP*ZMAT(K,J)
          END DO
      END DO
      GQSQ=ZERO
      DO I=1,N
          GQ(I)=GQ(I)+DIFF*BMAT(KNEW,I)
          GQSQ=GQSQ+GQ(I)**2
          XPT(KNEW,I)=XNEW(I)
      END DO
C
C     If a trust region step makes a small change to the objective function,
C     then calculate the gradient of the least Frobenius norm interpolant at
C     XBASE, and store it in W, using VLAG for a vector of right hand sides.
C
      IF (KSAVE == 0 .AND. DELTA == RHO) THEN
          IF (DABS(RATIO) > 1.0D-2) THEN
              ITEST=0
          ELSE
              DO K=1,NPT
                  VLAG(K)=FVAL(K)-FVAL(KOPT)
              END DO
              GISQ=ZERO
              DO I=1,N
                  SUM=ZERO
                  DO K=1,NPT
                      SUM=SUM+BMAT(K,I)*VLAG(K)
                  END DO
                  GISQ=GISQ+SUM*SUM
                  W(I)=SUM
              END DO
C
C     Test whether to replace the new quadratic model by the least Frobenius
C     norm interpolant, making the replacement if the test is satisfied.
C
              ITEST=ITEST+1
              IF (GQSQ < 1.0D2*GISQ) ITEST=0
              IF (ITEST >= 3) THEN
                  DO I=1,N
                      GQ(I)=W(I)
                  END DO
                  DO IH=1,NH
                      HQ(IH)=ZERO
                  END DO
                  DO J=1,NPTM
                      W(J)=ZERO
                      DO K=1,NPT
                          W(J)=W(J)+VLAG(K)*ZMAT(K,J)
                      END DO
                      IF (J < IDZ) W(J)=-W(J)
                  END DO
                  DO K=1,NPT
                      PQ(K)=ZERO
                      DO J=1,NPTM
                          PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
                      END DO
                  END DO
                  ITEST=0
              END IF
          END IF
      END IF
      IF (F < FSAVE) KOPT=KNEW
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case KSAVE>0 occurs
C     when the new function value was calculated by a model step.
C
      IF (F <= FSAVE+TENTH*VQUAD) GOTO 100
      IF (KSAVE > 0) GOTO 100
C
C     Alternatively, find out if the interpolation points are close enough
C     to the best point so far.
C
      KNEW=0
  460 DISTSQ=4.0D0*DELTA*DELTA
      DO K=1,NPT
          SUM=ZERO
          DO J=1,N
              SUM=SUM+(XPT(K,J)-XOPT(J))**2
          END DO
          IF (SUM > DISTSQ) THEN
              KNEW=K
              DISTSQ=SUM
          END IF
      END DO
C
C     If KNEW is positive, then set DSTEP, and branch back for the next
C     iteration, which will generate a "model step".
C
      IF (KNEW > 0) THEN
          DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
          DSQ=DSTEP*DSTEP
          GOTO 120
      END IF
      IF (RATIO > ZERO) GOTO 100
      IF (DMAX1(DELTA,DNORM) > RHO) GOTO 100
C
C     The calculations with the current value of RHO are complete. Pick the
C     next values of RHO and DELTA.
C
  490 IF (RHO > RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO <= 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO <= 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT >= 2) THEN
              IF (IPRINT >= 3) PRINT 500
  500         FORMAT (5X)
              PRINT 510, RHO,NF
  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 520, FOPT,(XBASE(I)+XOPT(I),I=1,N)
  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 90
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ELSE
          INFO=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C     it is too short to have been tried before.
C
      IF (KNEW == -1) GOTO 290
  530 IF (FOPT <= F) THEN
          DO I=1,N
              X(I)=XBASE(I)+XOPT(I)
          END DO
          F=FOPT
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (IPRINT .GE. 1) THEN
  546 IF (IPRINT >= 1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
          PRINT 550, NF
  550     FORMAT (/4X,'At the return from NEWUOA',5X,
     1      'Number of function values =',I6)
          PRINT 520, F,(X(I),I=1,N)
      END IF
      RETURN
      END
