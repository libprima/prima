      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,
     1  KNEW,W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      DO K=1,NPT
          DO J=1,NPTM
              ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
          END DO
      END DO
      ZTEST=1.0D-20*ZTEST
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C Zaikun 2019-08-15: JL is never used
C      JL=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO J=2,NPTM
          IF (DABS(ZMAT(KNEW,J)) > ZTEST) THEN
              TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
              TEMPA=ZMAT(KNEW,1)/TEMP
              TEMPB=ZMAT(KNEW,J)/TEMP
              DO I=1,NPT
                  TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
                  ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
                  ZMAT(I,1)=TEMP
              END DO
          END IF
          ZMAT(KNEW,J)=ZERO
      END DO
C
C     Put the first NPT components of the KNEW-th column of HLAG into W,
C     and calculate the parameters of the updating formula.
C
      DO I=1,NPT
          W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
      END DO
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE
C
C     Complete the updating of ZMAT.
C
      TEMP=DSQRT(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      DO I=1,NPT
          ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
      END DO
C
C     Finally, update the matrix BMAT.
C
      DO J=1,N
          JP=NPT+J
          W(JP)=BMAT(KNEW,J)
          TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
          TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
          DO I=1,JP
              BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
              IF (I > NPT) BMAT(JP,I-NPT)=BMAT(I,J)
          END DO
      END DO
      RETURN
      END
