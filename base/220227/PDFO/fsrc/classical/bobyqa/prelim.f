      SUBROUTINE PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     1  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
     1  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,F,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),
     1  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
C
C     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
C       same as the corresponding arguments in SUBROUTINE BOBYQA.
C     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
C       are the same as the corresponding arguments in BOBYQB, the elements
C       of SL and SU being set in BOBYQA.
C     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
C       it is set by PRELIM to the gradient of the quadratic model at XBASE.
C       If XOPT is nonzero, BOBYQB will change it to its usual value later.
C     NF is maintaned as the number of calls of CALFUN so far.
C     KOPT will be such that the least calculated value of F so far is at
C       the point XPT(KOPT,.)+XBASE in the space of the variables.
C
C     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
C     BMAT and ZMAT for the first iteration, and it maintains the values of
C     NF and KOPT. The vector X is also changed by PRELIM.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      NP=N+1
C
C     Set XBASE to the initial vector of variables, and set the initial
C     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
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
      DO IH=1,(N*NP)/2
          HQ(IH)=ZERO
      END DO
      DO K=1,NPT
          PQ(K)=ZERO
          DO J=1,NPT-NP
              ZMAT(K,J)=ZERO
          END DO
      END DO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF+1,.).
C
      NF=0
   50 NFM=NF
      NFX=NF-N
      NF=NF+1
      IF (NFM <= 2*N) THEN
          IF (NFM >= 1 .AND. NFM <= N) THEN
              STEPA=RHOBEG
              IF (SU(NFM) == ZERO) STEPA=-STEPA
              XPT(NF,NFM)=STEPA
          ELSE IF (NFM > N) THEN
              STEPA=XPT(NF-N,NFX)
              STEPB=-RHOBEG
              IF (SL(NFX) == ZERO) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
              IF (SU(NFX) == ZERO) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
              XPT(NF,NFX)=STEPB
          END IF
      ELSE
          ITEMP=(NFM-NP)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT > N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XPT(NF,IPT)=XPT(IPT+1,IPT)
          XPT(NF,JPT)=XPT(JPT+1,JPT)
      END IF
C
C     Calculate the next value of F. The least function value so far and
C     its index are required.
C
      DO J=1,N
          X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
          IF (XPT(NF,J) == SL(J)) X(J)=XL(J)
          IF (XPT(NF,J) == SU(J)) X(J)=XU(J)
      END DO
      CALL CALFUN (N,X,F)
      IF (IPRINT == 3) THEN
          PRINT 70, NF,F,(X(I),I=1,N)
   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(NF)=F
      IF (NF == 1) THEN
          FBEG=F
          KOPT=1
      ELSE IF (F < FVAL(KOPT)) THEN
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in the
C     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
C     of the NF-th and (NF-N)-th interpolation points may be switched, in
C     order that the function value at the first of them contributes to the
C     off-diagonal second derivative terms of the initial quadratic model.
C
      IF (NF <= 2*N+1) THEN
          IF (NF >= 2 .AND. NF <= N+1) THEN
              GOPT(NFM)=(F-FBEG)/STEPA
              IF (NPT < NF+N) THEN
                  BMAT(1,NFM)=-ONE/STEPA
                  BMAT(NF,NFM)=ONE/STEPA
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NF >= N+2) THEN
              IH=(NFX*(NFX+1))/2
              TEMP=(F-FBEG)/STEPB
              DIFF=STEPB-STEPA
              HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
              GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
              IF (STEPA*STEPB < ZERO) THEN
                  IF (F < FVAL(NF-N)) THEN
                      FVAL(NF)=FVAL(NF-N)
                      FVAL(NF-N)=F
                      IF (KOPT == NF) KOPT=NF-N
                      XPT(NF-N,NFX)=STEPB
                      XPT(NF,NFX)=STEPA
                  END IF
              END IF
              BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
              BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
              BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
              ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
              ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
              ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
          END IF
C
C     Set the off-diagonal second derivatives of the Lagrange functions and
C     the initial quadratic model.
C
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          ZMAT(1,NFX)=RECIP
          ZMAT(NF,NFX)=RECIP
          ZMAT(IPT+1,NFX)=-RECIP
          ZMAT(JPT+1,NFX)=-RECIP
          TEMP=XPT(NF,IPT)*XPT(NF,JPT)
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom (on 04-06-2019):
C     If the target value is reached, stop the algorithm.
      IF (F <= FTARGET) GOTO 80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (NF < NPT .AND. NF < MAXFUN) GOTO 50
   80 RETURN
      END
