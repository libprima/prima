      SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,
     1  DELTA,D,ALPHA,HCOL,GC,GD,S,W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  HCOL(*),GC(*),GD(*),S(*),W(*)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DELTA is the current trust region bound.
C     D will be set to the step from XOPT to the new point.
C     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
C     HCOL, GC, GD, S and W will be used for working space.
C
C     The step D is calculated in a way that attempts to maximize the modulus
C     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
C     the KNEW-th Lagrange function.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      DELSQ=DELTA*DELTA
      NPTM=NPT-N-1
C
C     Set the first NPT components of HCOL to the leading elements of the
C     KNEW-th column of H.
C
      ITERC=0
      DO K=1,NPT
          HCOL(K)=ZERO
      END DO
      DO J=1,NPTM
          TEMP=ZMAT(KNEW,J)
          IF (J < IDZ) TEMP=-TEMP
          DO K=1,NPT
              HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
          END DO
      END DO
      ALPHA=HCOL(KNEW)
C
C     Set the unscaled initial direction D. Form the gradient of LFUNC at
C     XOPT, and multiply D by the second derivative matrix of LFUNC.
C
      DD=ZERO
      DO I=1,N
          D(I)=XPT(KNEW,I)-XOPT(I)
          GC(I)=BMAT(KNEW,I)
          GD(I)=ZERO
          DD=DD+D(I)**2
      END DO
      DO K=1,NPT
          TEMP=ZERO
          SUM=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*XOPT(J)
              SUM=SUM+XPT(K,J)*D(J)
          END DO
          TEMP=HCOL(K)*TEMP
          SUM=HCOL(K)*SUM
          DO I=1,N
              GC(I)=GC(I)+TEMP*XPT(K,I)
              GD(I)=GD(I)+SUM*XPT(K,I)
          END DO
      END DO
C
C     Scale D and GD, with a sign change if required. Set S to another
C     vector in the initial two dimensional subspace.
C
      GG=ZERO
      SP=ZERO
      DHD=ZERO
      DO I=1,N
          GG=GG+GC(I)**2
          SP=SP+D(I)*GC(I)
          DHD=DHD+D(I)*GD(I)
      END DO
      SCALE=DELTA/DSQRT(DD)
      IF (SP*DHD < ZERO) SCALE=-SCALE
      TEMP=ZERO
      IF (SP*SP > 0.99D0*DD*GG) TEMP=ONE
      TAU=SCALE*(DABS(SP)+HALF*SCALE*DABS(DHD))
      IF (GG*DELSQ < 0.01D0*TAU*TAU) TEMP=ONE
      DO I=1,N
          D(I)=SCALE*D(I)
          GD(I)=SCALE*GD(I)
          S(I)=GC(I)+TEMP*GD(I)
      END DO
C
C     Begin the iteration by overwriting S with a vector that has the
C     required length and direction, except that termination occurs if
C     the given D and S are nearly parallel.
C
   80 ITERC=ITERC+1
      DD=ZERO
      SP=ZERO
      SS=ZERO
      DO I=1,N
          DD=DD+D(I)**2
          SP=SP+D(I)*S(I)
          SS=SS+S(I)**2
      END DO
      TEMP=DD*SS-SP*SP
      IF (TEMP <= 1.0D-8*DD*SS) GOTO 160
      DENOM=DSQRT(TEMP)
      DO I=1,N
          S(I)=(DD*S(I)-SP*D(I))/DENOM
          W(I)=ZERO
      END DO
C
C     Calculate the coefficients of the objective function on the circle,
C     beginning with the multiplication of S by the second derivative matrix.
C
      DO K=1,NPT
          SUM=ZERO
          DO J=1,N
              SUM=SUM+XPT(K,J)*S(J)
          END DO
          SUM=HCOL(K)*SUM
          DO I=1,N
              W(I)=W(I)+SUM*XPT(K,I)
          END DO
      END DO
      CF1=ZERO
      CF2=ZERO
      CF3=ZERO
      CF4=ZERO
      CF5=ZERO
      DO I=1,N
          CF1=CF1+S(I)*W(I)
          CF2=CF2+D(I)*GC(I)
          CF3=CF3+S(I)*GC(I)
          CF4=CF4+D(I)*GD(I)
          CF5=CF5+S(I)*GD(I)
      END DO
      CF1=HALF*CF1
      CF4=HALF*CF4-CF1
C
C     Seek the value of the angle that maximizes the modulus of TAU.
C
      TAUBEG=CF1+CF2+CF4
      TAUMAX=TAUBEG
      TAUOLD=TAUBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO I=1,IU
          ANGLE=DFLOAT(I)*TEMP
          CTH=DCOS(ANGLE)
          STH=DSIN(ANGLE)
          TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
          IF (DABS(TAU) > DABS(TAUMAX)) THEN
              TAUMAX=TAU
              ISAVE=I
              TEMPA=TAUOLD
          ELSE IF (I == ISAVE+1) THEN
              TEMPB=TAU
          END IF
          TAUOLD=TAU
      END DO
      IF (ISAVE == 0) TEMPA=TAU
      IF (ISAVE == IU) TEMPB=TAUBEG
      STEP=ZERO
      IF (TEMPA /= TEMPB) THEN
          TEMPA=TEMPA-TAUMAX
          TEMPB=TEMPB-TAUMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new D and GD. Then test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      DO I=1,N
          D(I)=CTH*D(I)+STH*S(I)
          GD(I)=CTH*GD(I)+STH*W(I)
          S(I)=GC(I)+GD(I)
      END DO
      IF (DABS(TAU) <= 1.1D0*DABS(TAUBEG)) GOTO 160
      IF (ITERC < N) GOTO 80
  160 RETURN
      END
