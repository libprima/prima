      SUBROUTINE PRELIM (N,NPT,M,AMAT,B,X,RHOBEG,IPRINT,XBASE,
     1  XPT,FVAL,XSAV,XOPT,GOPT,KOPT,HQ,PQ,BMAT,ZMAT,IDZ,NDIM,
     2  SP,RESCON,STEP,PQW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),X(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XSAV(*),XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
     2  SP(*),RESCON(*),STEP(*),PQW(*),W(*)
C
C     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
C       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
C       same as the corresponding arguments in SUBROUTINE LINCOB.
C     KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
C       region centre.
C     IDZ is going to be set to one, so that every element of Diag(DZ) is
C       one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
C       factorization of the leading NPT by NPT submatrix of H.
C     STEP, PQW and W are used for working space, the arrays STEP and PQW
C       being taken from LINCOB. The length of W must be at least N+NPT.
C
C     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
C       for the first iteration, an important feature being that, if any of
C       of the columns of XPT is an infeasible point, then the largest of
C       the constraint violations there is at least 0.2*RHOBEG. It also sets
C       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      TEST=0.2D0*RHOBEG
      IDZ=1
      KBASE=1
C
C     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. 
C
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 K=1,NPT
      SP(K)=ZERO
      DO 30 J=1,NPT-N-1
   30 ZMAT(K,J)=ZERO
C
C     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
C       but they may be altered later to make a constraint violation
C       sufficiently large. The initial nonzero elements of BMAT and of
C       the first min[N,NPT-N-1] columns of ZMAT are set also.
C
      DO 40 J=1,N
      XPT(J+1,J)=RHOBEG
      IF (J .LT. NPT-N) THEN
          JP=N+J+1
          XPT(JP,J)=-RHOBEG
          BMAT(J+1,J)=HALF/RHOBEG
          BMAT(JP,J)=-HALF/RHOBEG
          ZMAT(1,J)=-RECIQ-RECIQ
          ZMAT(J+1,J)=RECIQ
          ZMAT(JP,J)=RECIQ
      ELSE
          BMAT(1,J)=-ONE/RHOBEG
          BMAT(J+1,J)=ONE/RHOBEG
          BMAT(NPT+J,J)=-HALF*RHOSQ
      END IF
   40 CONTINUE
C
C     Set the remaining initial nonzero elements of XPT and ZMAT when the
C       number of interpolation points exceeds 2*N+1.
C
      IF (NPT .GT. 2*N+1) THEN
          DO 50 K=N+1,NPT-N-1
          ITEMP=(K-1)/N
          IPT=K-ITEMP*N
          JPT=IPT+ITEMP
          IF (JPT .GT. N) JPT=JPT-N
          XPT(N+K+1,IPT)=RHOBEG
          XPT(N+K+1,JPT)=RHOBEG
          ZMAT(1,K)=RECIP
          ZMAT(IPT+1,K)=-RECIP
          ZMAT(JPT+1,K)=-RECIP
   50     ZMAT(N+K+1,K)=RECIP
      END IF
C
C     Update the constraint right hand sides to allow for the shift XBASE.
C
      IF (M .GT. 0) THEN
          DO 70 J=1,M
          TEMP=ZERO
          DO 60 I=1,N
   60     TEMP=TEMP+AMAT(I,J)*XBASE(I)
   70     B(J)=B(J)-TEMP
      END IF
C
C     Go through the initial points, shifting every infeasible point if
C       necessary so that its constraint violation is at least 0.2*RHOBEG.
C
      DO 150 NF=1,NPT
      FEAS=ONE
      BIGV=ZERO
      J=0
   80 J=J+1
      IF (J .LE. M .AND. NF .GE. 2) THEN
          RESID=-B(J)
          DO 90 I=1,N
   90     RESID=RESID+XPT(NF,I)*AMAT(I,J)
          IF (RESID .LE. BIGV) GOTO 80
          BIGV=RESID
          JSAV=J
          IF (RESID .LE. TEST) THEN
              FEAS=-ONE
              GOTO 80
          END IF
          FEAS=ZERO
      END IF
      IF (FEAS .LT. ZERO) THEN
          DO 100 I=1,N
  100     STEP(I)=XPT(NF,I)+(TEST-BIGV)*AMAT(I,JSAV)
          DO 110 K=1,NPT
          SP(NPT+K)=ZERO
          DO 110 J=1,N
  110     SP(NPT+K)=SP(NPT+K)+XPT(K,J)*STEP(J)
          CALL UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,
     1      KBASE,NF,PQW,W)
          DO 120 I=1,N
  120     XPT(NF,I)=STEP(I)
      END IF
C
C     Calculate the objective function at the current interpolation point,
C       and set KOPT to the index of the first trust region centre.
C
      DO 130 J=1,N
  130 X(J)=XBASE(J)+XPT(NF,J)
      F=FEAS
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 140, NF,F,(X(I),I=1,N)
  140     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .EQ. 1) THEN
          KOPT=1
      ELSE IF (F .LT. FVAL(KOPT) .AND. FEAS .GT. ZERO) THEN
          KOPT=NF
      END IF
  150 FVAL(NF)=F
C
C     Set PQ for the first quadratic model.
C
      DO 160 J=1,NPTM
      W(J)=ZERO
      DO 160 K=1,NPT
  160 W(J)=W(J)+ZMAT(K,J)*FVAL(K)
      DO 170 K=1,NPT
      PQ(K)=ZERO
      DO 170 J=1,NPTM
  170 PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
C
C     Set XOPT, SP, GOPT and HQ for the first quadratic model.
C
      DO 180 J=1,N
      XOPT(J)=XPT(KOPT,J)
      XSAV(J)=XBASE(J)+XOPT(J)
  180 GOPT(J)=ZERO
      DO 200 K=1,NPT
      SP(K)=ZERO
      DO 190 J=1,N
  190 SP(K)=SP(K)+XPT(K,J)*XOPT(J)
      TEMP=PQ(K)*SP(K)
      DO 200 J=1,N
  200 GOPT(J)=GOPT(J)+FVAL(K)*BMAT(K,J)+TEMP*XPT(K,J)
      DO 210 I=1,(N*N+N)/2
  210 HQ(I)=ZERO
C
C     Set the initial elements of RESCON.
C
      DO 230 J=1,M
      TEMP=B(J)
      DO 220 I=1,N
  220 TEMP=TEMP-XOPT(I)*AMAT(I,J)
      TEMP=DMAX1(TEMP,ZERO)
      IF (TEMP .GE. RHOBEG) TEMP=-TEMP
  230 RESCON(J)=TEMP  
      RETURN
      END
