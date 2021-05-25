      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
C
C     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
C     interpolation point that has index KNEW. On entry, VLAG contains the
C     components of the vector Theta*Wcheck+e_b of the updating formula
C     (6.11), and BETA holds the value of the parameter that has this name.
C     The vector W is used for working space.
C
C     Set some constants.
C
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      DO J=2,NPTM
          IF (J == IDZ) THEN
              JL=IDZ
          ELSE IF (ZMAT(KNEW,J) /= ZERO) THEN
              TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
              TEMPA=ZMAT(KNEW,JL)/TEMP
              TEMPB=ZMAT(KNEW,J)/TEMP
              DO I=1,NPT
                  TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J)
                  ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL)
                  ZMAT(I,JL)=TEMP
              END DO
              ZMAT(KNEW,J)=ZERO
          END IF
      END DO
C
C     Put the first NPT components of the KNEW-th column of HLAG into W,
C     and calculate the parameters of the updating formula.
C
      TEMPA=ZMAT(KNEW,1)
      IF (IDZ >= 2) TEMPA=-TEMPA
      IF (JL > 1) TEMPB=ZMAT(KNEW,JL)
      DO I=1,NPT
          W(I)=TEMPA*ZMAT(I,1)
          IF (JL > 1) W(I)=W(I)+TEMPB*ZMAT(I,JL)
      END DO
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      TAUSQ=TAU*TAU
      DENOM=ALPHA*BETA+TAUSQ
      VLAG(KNEW)=VLAG(KNEW)-ONE
C
C     Complete the updating of ZMAT when there is only one nonzero element
C     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
C     then the first column of ZMAT will be exchanged with another one later.
C
      IFLAG=0
      IF (JL == 1) THEN
          TEMP=DSQRT(DABS(DENOM))
          TEMPB=TEMPA/TEMP
          TEMPA=TAU/TEMP
          DO I=1,NPT
              ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
          END DO
          IF (IDZ == 1 .AND. TEMP < ZERO) IDZ=2
          IF (IDZ >= 2 .AND. TEMP >= ZERO) IFLAG=1
      ELSE
C
C     Complete the updating of ZMAT in the alternative case.
C
          JA=1
          IF (BETA >= ZERO) JA=JL
          JB=JL+1-JA
          TEMP=ZMAT(KNEW,JB)/DENOM
          TEMPA=TEMP*BETA
          TEMPB=TEMP*TAU
          TEMP=ZMAT(KNEW,JA)
          SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ)
          SCALB=SCALA*DSQRT(DABS(DENOM))
          DO I=1,NPT
              ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
              ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
          END DO
          IF (DENOM <= ZERO) THEN
              IF (BETA < ZERO) IDZ=IDZ+1
              IF (BETA >= ZERO) IFLAG=1
          END IF
      END IF
C
C     IDZ is reduced in the following case, and usually the first column
C     of ZMAT is exchanged with a later one.
C
      IF (IFLAG == 1) THEN
          IDZ=IDZ-1
          DO I=1,NPT
              TEMP=ZMAT(I,1)
              ZMAT(I,1)=ZMAT(I,IDZ)
              ZMAT(I,IDZ)=TEMP
          END DO
      END IF
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
