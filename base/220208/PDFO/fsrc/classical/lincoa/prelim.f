      SUBROUTINE PRELIM (N,NPT,M,AMAT,B,X,RHOBEG,IPRINT,XBASE,
     1  XPT,FVAL,XSAV,XOPT,GOPT,KOPT,HQ,PQ,BMAT,ZMAT,IDZ,NDIM,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  SP,RESCON,STEP,PQW,W)
     2  SP,RESCON,STEP,PQW,W,F,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      DO J=1,N
          XBASE(J)=X(J)
          DO K=1,NPT
              XPT(K,J)=ZERO
          END DO
          DO I=1,NDIM
              BMAT(I,J)=ZERO
          END DO
      END DO
      DO K=1,NPT
          SP(K)=ZERO
          DO J=1,NPT-N-1
              ZMAT(K,J)=ZERO
          END DO
      END DO
C
C     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
C       but they may be altered later to make a constraint violation
C       sufficiently large. The initial nonzero elements of BMAT and of
C       the first min[N,NPT-N-1] columns of ZMAT are set also.
C
      DO J=1,N
          XPT(J+1,J)=RHOBEG
          IF (J < NPT-N) THEN
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
      END DO
C
C     Set the remaining initial nonzero elements of XPT and ZMAT when the
C       number of interpolation points exceeds 2*N+1.
C
      IF (NPT > 2*N+1) THEN
          DO K=N+1,NPT-N-1
              ITEMP=(K-1)/N
              IPT=K-ITEMP*N
              JPT=IPT+ITEMP
              IF (JPT > N) JPT=JPT-N
              XPT(N+K+1,IPT)=RHOBEG
              XPT(N+K+1,JPT)=RHOBEG
              ZMAT(1,K)=RECIP
              ZMAT(IPT+1,K)=-RECIP
              ZMAT(JPT+1,K)=-RECIP
              ZMAT(N+K+1,K)=RECIP
          END DO
      END IF
C
C     Update the constraint right hand sides to allow for the shift XBASE.
C
      IF (M > 0) THEN
          DO J=1,M
              TEMP=ZERO
              DO I=1,N
                  TEMP=TEMP+AMAT(I,J)*XBASE(I)
              END DO
              B(J)=B(J)-TEMP
          END DO
      END IF
C
C     Go through the initial points, shifting every infeasible point if
C       necessary so that its constraint violation is at least 0.2*RHOBEG.
C
      DO NF=1,NPT
          FEAS=ONE
          BIGV=ZERO
          J=0
   80     J=J+1
          IF (J <= M .AND. NF >= 2) THEN
              RESID=-B(J)
              DO I=1,N
                  RESID=RESID+XPT(NF,I)*AMAT(I,J)
              END DO
              IF (RESID <= BIGV) GOTO 80
              BIGV=RESID
              JSAV=J
              IF (RESID <= TEST) THEN
                  FEAS=-ONE
                  GOTO 80
              END IF
              FEAS=ZERO
          END IF
          IF (FEAS < ZERO) THEN
              DO I=1,N
                  STEP(I)=XPT(NF,I)+(TEST-BIGV)*AMAT(I,JSAV)
              END DO
              DO K=1,NPT
                  SP(NPT+K)=ZERO
                  DO J=1,N
                      SP(NPT+K)=SP(NPT+K)+XPT(K,J)*STEP(J)
                  END DO
              END DO
              CALL UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,
     1          KBASE,NF,PQW,W)
              DO I=1,N
                  XPT(NF,I)=STEP(I)
              END DO
          END IF
C
C     Calculate the objective function at the current interpolation point,
C       and set KOPT to the index of the first trust region centre.
C
          DO J=1,N
              X(J)=XBASE(J)+XPT(NF,J)
          END DO
          F=FEAS
          CALL CALFUN (N,X,F)
          IF (IPRINT == 3) THEN
              PRINT 140, NF,F,(X(I),I=1,N)
  140         FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1          '    The corresponding X is:'/(2X,5D15.6))
          END IF
          IF (NF == 1) THEN
              KOPT=1
          ELSE IF (F < FVAL(KOPT) .AND. FEAS > ZERO) THEN
              KOPT=NF
          END IF
          FVAL(NF)=F
      END DO
C
C     Set PQ for the first quadratic model.
C
      DO J=1,NPTM
          W(J)=ZERO
          DO K=1,NPT
              W(J)=W(J)+ZMAT(K,J)*FVAL(K)
          END DO
      END DO
      DO K=1,NPT
          PQ(K)=ZERO
          DO J=1,NPTM
              PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
          END DO
      END DO
C
C     Set XOPT, SP, GOPT and HQ for the first quadratic model.
C
      DO J=1,N
          XOPT(J)=XPT(KOPT,J)
          XSAV(J)=XBASE(J)+XOPT(J)
          GOPT(J)=ZERO
      END DO
      DO K=1,NPT
          SP(K)=ZERO
          DO J=1,N
              SP(K)=SP(K)+XPT(K,J)*XOPT(J)
          END DO
          TEMP=PQ(K)*SP(K)
          DO J=1,N
              GOPT(J)=GOPT(J)+FVAL(K)*BMAT(K,J)+TEMP*XPT(K,J)
          END DO
      END DO
      DO I=1,(N*N+N)/2
          HQ(I)=ZERO
      END DO
C
C     Set the initial elements of RESCON.
C
      DO J=1,M
          TEMP=B(J)
          DO I=1,N
              TEMP=TEMP-XOPT(I)*AMAT(I,J)
          END DO
          TEMP=DMAX1(TEMP,ZERO)
          IF (TEMP >= RHOBEG) TEMP=-TEMP
          RESCON(J)=TEMP  
      END DO
      RETURN
      END
