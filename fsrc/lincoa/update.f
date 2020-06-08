      SUBROUTINE UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,
     1  KOPT,KNEW,VLAG,W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),SP(*),STEP(*),
     2  VLAG(*),W(*)
C
C     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
C       identical to the corresponding arguments in SUBROUTINE LINCOB.
C     KOPT is such that XPT(KOPT,.) is the current trust region centre.
C     KNEW on exit is usually positive, and then it is the index of an
C       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
C       It is set on entry either to its final value or to 0. In the latter
C       case, the final value of KNEW is chosen to maximize the denominator
C       of the matrix updating formula times a weighting factor.
C     VLAG and W are used for working space, the first NPT+N elements of
C       both of these vectors being required.
C
C     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
C       the ones that are suitable after the shift of the KNEW-th point to
C       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero
C       occurs if the calculation fails due to a zero denominator in the
C       updating formula, which should never happen.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
C
C     Calculate VLAG and BETA for the current choice of STEP. The first NPT
C       elements of VLAG are set to the values of the Lagrange functions at
C       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
C       in W, where W_check is defined in a paper on the updating method.
C
      DO K=1,NPT
          W(K)=SP(NPT+K)*(HALF*SP(NPT+K)+SP(K))
          SUM=ZERO
          DO J=1,N
              SUM=SUM+BMAT(K,J)*STEP(J)
          END DO
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
      SSQ=ZERO
      DO J=1,N
          SUM=ZERO
          DO I=1,NPT
              SUM=SUM+W(I)*BMAT(I,J)
          END DO
          BSUM=BSUM+SUM*STEP(J)
          JP=NPT+J
          DO K=1,N
              SUM=SUM+BMAT(JP,K)*STEP(K)
          END DO
          VLAG(JP)=SUM
          BSUM=BSUM+SUM*STEP(J)
          DX=DX+STEP(J)*XPT(KOPT,J)
          SSQ=SSQ+STEP(J)**2
      END DO
      BETA=DX*DX+SSQ*(SP(KOPT)+DX+DX+HALF*SSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is zero initially, then pick the index of the interpolation
C       point to be deleted, by maximizing the absolute value of the
C       denominator of the updating formula times a weighting factor.
C       
C
      IF (KNEW == 0) THEN
          DENMAX=ZERO
          DO K=1,NPT
              HDIAG=ZERO
              DO J=1,NPTM
                  TEMP=ONE
                  IF (J < IDZ) TEMP=-ONE
                  HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
              END DO
              DENABS=DABS(BETA*HDIAG+VLAG(K)**2)
              DISTSQ=ZERO
              DO J=1,N
                  DISTSQ=DISTSQ+(XPT(K,J)-XPT(KOPT,J))**2
              END DO
              TEMP=DENABS*DISTSQ*DISTSQ
              IF (TEMP > DENMAX) THEN
                  DENMAX=TEMP
                  KNEW=K
              END IF
          END DO
      END IF
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      IF (NPTM >= 2) THEN
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
      END IF
C
C     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
C       into W, and calculate the parameters of the updating formula.
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
      IF (DENOM == ZERO) THEN
          KNEW=0
          GOTO 180
      END IF
      SQRTDN=DSQRT(DABS(DENOM))
C
C     Complete the updating of ZMAT when there is only one nonzero element
C       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
C       the value of IDZ is going to be reduced.
C
      IFLAG=0
      IF (JL == 1) THEN
          TEMPA=TAU/SQRTDN
          TEMPB=ZMAT(KNEW,1)/SQRTDN
          DO I=1,NPT
              ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
          END DO
          IF (DENOM < ZERO) THEN
              IF (IDZ == 1) THEN
                  IDZ=2
              ELSE
                  IFLAG=1
              END IF
          END IF
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
          SCALB=SCALA*SQRTDN
          DO I=1,NPT
              ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
              ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
          END DO
          IF (DENOM <= ZERO) THEN
              IF (BETA < ZERO) THEN
                  IDZ=IDZ+1
              ELSE
                  IFLAG=1
              END IF
          END IF
      END IF
C
C     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
C       ZMAT^T factorization gains another positive element. Then exchange
C       the first and IDZ-th columns of ZMAT.
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
  180 RETURN
      END
