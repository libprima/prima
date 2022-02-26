      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,
     1  KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
C
C     The arrays BMAT and ZMAT are updated, as required by the new position
C     of the interpolation point that has the index KNEW. The vector VLAG has
C     N+NPT components, set on entry to the first NPT and last N components
C     of the product Hw in equation (4.11) of the Powell (2006) paper on
C     NEWUOA. Further, BETA is set on entry to the value of the parameter
C     with that name, and DENOM is set to the denominator of the updating
C     formula. Elements of ZMAT may be treated as zero if their moduli are
C     at most ZTEST. The first NDIM elements of W are used for working space.
C
C     Set some constants.
C
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      ZTEST=ZERO
      DO 10 K=1,NPT
      DO 10 J=1,NPTM
   10 ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
      ZTEST=1.0D-20*ZTEST
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      DO 30 J=2,NPTM
      IF (DABS(ZMAT(KNEW,J)) .GT. ZTEST) THEN
          TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,1)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 20 I=1,NPT
          TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
   20     ZMAT(I,1)=TEMP
      END IF
      ZMAT(KNEW,J)=ZERO
   30 CONTINUE
C
C     Put the first NPT components of the KNEW-th column of HLAG into W,
C     and calculate the parameters of the updating formula.
C
      DO 40 I=1,NPT
      W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
   40 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE
C
C     Complete the updating of ZMAT.
C
      TEMP=DSQRT(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      DO 50 I=1,NPT
   50 ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
C
C     Finally, update the matrix BMAT.
C
      DO 60 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 60 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   60 CONTINUE
      RETURN
      END
