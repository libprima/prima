      SUBROUTINE TRSTLP (N,M,A,B,RHO,DX,IFULL,IACT,Z,ZDOTA,VMULTC,
     1  SDIRN,DXNEW,VMULTD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION A(N,*),B(*),DX(*),IACT(*),Z(N,*),ZDOTA(*),
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*)
     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*),DSAV(N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     This subroutine calculates an N-component vector DX by applying the
C     following two stages. In the first stage, DX is set to the shortest
C     vector that minimizes the greatest violation of the constraints
C       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
C     subject to the Euclidean length of DX being at most RHO. If its length is
C     strictly less than RHO, then we use the resultant freedom in DX to
C     minimize the objective function
C              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
C     subject to no increase in any greatest constraint violation. This
C     notation allows the gradient of the objective function to be regarded as
C     the gradient of a constraint. Therefore the two stages are distinguished
C     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
C     degeneracy may prevent DX from attaining the target length RHO. Then the
C     value IFULL=0 would be set, but usually IFULL=1 on return.
C
C     In general NACT is the number of constraints in the active set and
C     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
C     contains a permutation of the remaining constraint indices. Further, Z is
C     an orthogonal matrix whose first NACT columns can be regarded as the
C     result of Gram-Schmidt applied to the active constraint gradients. For
C     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
C     column of Z with the gradient of the J-th active constraint. DX is the
C     current vector of variables and here the residuals of the active
C     constraints should be zero. Further, the active constraints have
C     nonnegative Lagrange multipliers that are held at the beginning of
C     VMULTC. The remainder of this vector holds the residuals of the inactive
C     constraints at DX, the ordering of the components of VMULTC being in
C     agreement with the permutation of the indices of the constraints that is
C     in IACT. All these residuals are nonnegative, which is achieved by the
C     shift RESMAX that makes the least residual zero.
C
C     Initialize Z and some other variables. The value of RESMAX will be
C     appropriate to DX=0, while ICON will be the index of a most violated
C     constraint if RESMAX is positive. Usually during the first stage the
C     vector SDIRN gives a search direction that reduces all the active
C     constraint violations by one simultaneously.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 26-06-2019: See the code below line number 80
      ITERC=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IFULL=1
      MCON=M
      NACT=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RESMAX=0.0
      RESMAX=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          DO J=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   10 Z(I,J)=0.0
C      Z(I,I)=1.0
C   20 DX(I)=0.0
              Z(I,J)=0.0D0
          END DO
          Z(I,I)=1.0D0
          DX(I)=0.0D0
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (M >= 1) THEN
          DO K=1,M
              IF (B(K) > RESMAX) THEN
                  RESMAX=B(K)
                  ICON=K
              END IF
          END DO
          DO K=1,M
              IACT(K)=K
              VMULTC(K)=RESMAX-B(K)
          END DO
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (RESMAX .EQ. 0.0) GOTO 480
      IF (RESMAX == 0.0D0) GOTO 480
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   50 SDIRN(I)=0.0
          SDIRN(I)=0.0D0
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 26-06-2019: See the code below line number 80
      DO I = 1, N
         DSAV(I) = DX(I)
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     End the current stage of the calculation if 3 consecutive iterations
C     have either failed to reduce the best calculated value of the objective
C     function or to increase the number of active constraints since the best
C     value was calculated. This strategy prevents cycling, but there is a
C     remote possibility that it will cause premature termination.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   60 OPTOLD=0.0
   60 OPTOLD=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICOUNT=0
   70 IF (MCON == M) THEN
          OPTNEW=RESMAX
      ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          OPTNEW=0.0
          OPTNEW=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              OPTNEW=OPTNEW-DX(I)*A(I,MCON)
          END DO
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 26-06-2019
C The original code can still encounter infinite cycling,
C which did happen when testing the following CUTEst problems:
C DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
C Indeed, in all these cases, Inf and NaN appear in D due to extremely
C large values in A (up to 10^219). 
C To avoid wasting energy, we do the following.
      SUMD = 0.0D0 
      DO I = 1, N
          SUMD = SUMD + DABS(DX(I))
      END DO
      IF (SUMD >= 1.0D100 .OR. SUMD /= SUMD) THEN
          DO I = 1, N
             DX(I) = DSAV(I)
          END DO
          GOTO 500 
      ELSE
          DO I = 1, N
             DSAV(I) = DX(I)
          END DO
      END IF
      ITERC = ITERC + 1
      IF (ITERC > MIN(10000, 100*N)) THEN
          GOTO 500
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ICOUNT == 0 .OR. OPTNEW < OPTOLD) THEN
          OPTOLD=OPTNEW
          NACTX=NACT
          ICOUNT=3
      ELSE IF (NACT > NACTX) THEN
          NACTX=NACT
          ICOUNT=3
      ELSE
          ICOUNT=ICOUNT-1
          IF (ICOUNT == 0) GOTO 490
      END IF
C
C     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
C     the active set. Apply Givens rotations so that the last N-NACT-1 columns
C     of Z are orthogonal to the gradient of the new constraint, a scalar
C     product being set to zero if its nonzero value could be due to computer
C     rounding errors. The array DXNEW is used for working space.
C
      IF (ICON <= NACT) GOTO 260
      KK=IACT(ICON)
      DO I=1,N
          DXNEW(I)=A(I,KK)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TOT=0.0
      TOT=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      K=N
  100 IF (K > NACT) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          SP=0.0
C          SPABS=0.0
          SP=0.0D0
          SPABS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              TEMP=Z(I,K)*DXNEW(I)
              SP=SP+TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  110     SPABS=SPABS+ABS(TEMP)
C          ACCA=SPABS+0.1*ABS(SP)
C          ACCB=SPABS+0.2*ABS(SP)
C          IF (SPABS .GE. ACCA .OR. ACCA .GE. ACCB) SP=0.0
C          IF (TOT .EQ. 0.0) THEN
              SPABS=SPABS+DABS(TEMP)
          END DO
          ACCA=SPABS+0.1D0*DABS(SP)
          ACCB=SPABS+0.2D0*DABS(SP)
          IF (SPABS >= ACCA .OR. ACCA >= ACCB) SP=0.0D0
          IF (TOT == 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TOT=SP
          ELSE
              KP=K+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              TEMP=SQRT(SP*SP+TOT*TOT)
              TEMP=DSQRT(SP*SP+TOT*TOT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ALPHA=SP/TEMP
              BETA=TOT/TEMP
              TOT=TEMP
              DO I=1,N
                  TEMP=ALPHA*Z(I,K)+BETA*Z(I,KP)
                  Z(I,KP)=ALPHA*Z(I,KP)-BETA*Z(I,K)
                  Z(I,K)=TEMP
              END DO
          END IF
          K=K-1
          GOTO 100
      END IF
C
C     Add the new constraint if this can be done without a deletion from the
C     active set.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (TOT .NE. 0.0) THEN
      IF (TOT /= 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          NACT=NACT+1
          ZDOTA(NACT)=TOT
          VMULTC(ICON)=VMULTC(NACT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          VMULTC(NACT)=0.0
          VMULTC(NACT)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 210
      END IF
C
C     The next instruction is reached if a deletion has to be made from the
C     active set in order to make room for the new active constraint, because
C     the new constraint gradient is a linear combination of the gradients of
C     the old active constraints. Set the elements of VMULTD to the multipliers
C     of the linear combination. Further, set IOUT to the index of the
C     constraint to be deleted, but branch if no suitable index can be found.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RATIO=-1.0
      RATIO=-1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      K=NACT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  130 ZDOTV=0.0
C      ZDVABS=0.0
  130 ZDOTV=0.0D0
      ZDVABS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=Z(I,K)*DXNEW(I)
          ZDOTV=ZDOTV+TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  140 ZDVABS=ZDVABS+ABS(TEMP)
C      ACCA=ZDVABS+0.1*ABS(ZDOTV)
C      ACCB=ZDVABS+0.2*ABS(ZDOTV)
          ZDVABS=ZDVABS+DABS(TEMP)
      END DO
      ACCA=ZDVABS+0.1D0*DABS(ZDOTV)
      ACCB=ZDVABS+0.2D0*DABS(ZDOTV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ZDVABS < ACCA .AND. ACCA < ACCB) THEN
          TEMP=ZDOTV/ZDOTA(K)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TEMP .GT. 0.0 .AND. IACT(K) .LE. M) THEN
          IF (TEMP > 0.0D0 .AND. IACT(K) <= M) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TEMPA=VMULTC(K)/TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              IF (RATIO .LT. 0.0 .OR. TEMPA .LT. RATIO) THEN
              IF (RATIO < 0.0D0 .OR. TEMPA < RATIO) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  RATIO=TEMPA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-15: IOUT is never used
C                  IOUT=K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              END IF
           END IF
          IF (K >= 2) THEN
              KW=IACT(K)
              DO I=1,N
                  DXNEW(I)=DXNEW(I)-TEMP*A(I,KW)
              END DO
          END IF
          VMULTD(K)=TEMP
      ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          VMULTD(K)=0.0
          VMULTD(K)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
      K=K-1
      IF (K > 0) GOTO 130
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (RATIO .LT. 0.0) GOTO 490
      IF (RATIO < 0.0D0) GOTO 490
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Revise the Lagrange multipliers and reorder the active constraints so
C     that the one to be replaced is at the end of the list. Also calculate the
C     new value of ZDOTA(NACT) and branch if it is not acceptable.
C
      DO K=1,NACT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  160 VMULTC(K)=AMAX1(0.0,VMULTC(K)-RATIO*VMULTD(K))
          VMULTC(K)=DMAX1(0.0D0,VMULTC(K)-RATIO*VMULTD(K))
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ICON < NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  170     KP=K+1
          KW=IACT(KP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          SP=0.0
          SP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SP=SP+Z(I,K)*A(I,KW)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          TEMP=DSQRT(SP*SP+ZDOTA(KP)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ALPHA=ZDOTA(KP)/TEMP
          BETA=SP/TEMP
          ZDOTA(KP)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO I=1,N
              TEMP=ALPHA*Z(I,KP)+BETA*Z(I,K)
              Z(I,KP)=ALPHA*Z(I,K)-BETA*Z(I,KP)
              Z(I,K)=TEMP
          END DO
          IACT(K)=KW
          VMULTC(K)=VMULTC(KP)
          K=KP
          IF (K < NACT) GOTO 170
          IACT(K)=ISAVE
          VMULTC(K)=VSAVE
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=TEMP+Z(I,NACT)*A(I,KK)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (TEMP .EQ. 0.0) GOTO 490
      IF (TEMP == 0.0D0) GOTO 490
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ZDOTA(NACT)=TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      VMULTC(ICON)=0.0
      VMULTC(ICON)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VMULTC(NACT)=RATIO
C
C     Update IACT and ensure that the objective function continues to be
C     treated as the last active constraint when MCON>M.
C
  210 IACT(ICON)=IACT(NACT)
      IACT(NACT)=KK
      IF (MCON > M .AND. KK /= MCON) THEN
          K=NACT-1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          SP=0.0
          SP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SP=SP+Z(I,K)*A(I,KK)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
          TEMP=DSQRT(SP*SP+ZDOTA(NACT)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ALPHA=ZDOTA(NACT)/TEMP
          BETA=SP/TEMP
          ZDOTA(NACT)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO I=1,N
              TEMP=ALPHA*Z(I,NACT)+BETA*Z(I,K)
              Z(I,NACT)=ALPHA*Z(I,K)-BETA*Z(I,NACT)
              Z(I,K)=TEMP
          END DO
          IACT(NACT)=IACT(K)
          IACT(K)=KK
          TEMP=VMULTC(K)
          VMULTC(K)=VMULTC(NACT)
          VMULTC(NACT)=TEMP
      END IF
C
C     If stage one is in progress, then set SDIRN to the direction of the next
C     change to the current vector of variables.
C
      IF (MCON > M) GOTO 320
      KK=IACT(NACT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=TEMP+SDIRN(I)*A(I,KK)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=TEMP-1.0
      TEMP=TEMP-1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TEMP=TEMP/ZDOTA(NACT)
      DO I=1,N
          SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT)
      END DO
      GOTO 340
C
C     Delete the constraint that has the index IACT(ICON) from the active set.
C
  260 IF (ICON < NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  270     KP=K+1
          KK=IACT(KP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          SP=0.0
          SP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SP=SP+Z(I,K)*A(I,KK)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          TEMP=DSQRT(SP*SP+ZDOTA(KP)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ALPHA=ZDOTA(KP)/TEMP
          BETA=SP/TEMP
          ZDOTA(KP)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO I=1,N
              TEMP=ALPHA*Z(I,KP)+BETA*Z(I,K)
              Z(I,KP)=ALPHA*Z(I,K)-BETA*Z(I,KP)
              Z(I,K)=TEMP
          END DO
          IACT(K)=KK
          VMULTC(K)=VMULTC(KP)
          K=KP
          IF (K < NACT) GOTO 270
          IACT(K)=ISAVE
          VMULTC(K)=VSAVE
      END IF
      NACT=NACT-1
C
C     If stage one is in progress, then set SDIRN to the direction of the next
C     change to the current vector of variables.
C
      IF (MCON > M) GOTO 320
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=TEMP+SDIRN(I)*Z(I,NACT+1)
      END DO
      DO I=1,N
          SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT+1)
      END DO
      GO TO 340
C
C     Pick the next search direction of stage two.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  320 TEMP=1.0/ZDOTA(NACT)
  320 TEMP=1.0D0/ZDOTA(NACT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          SDIRN(I)=TEMP*Z(I,NACT)
      END DO
C
C     Calculate the step to the boundary of the trust region or take the step
C     that reduces RESMAX to zero. The two statements below that include the
C     factor 1.0E-6 prevent some harmless underflows that occurred in a test
C     calculation. Further, we skip the step if it could be zero within a
C     reasonable tolerance for computer rounding errors.
C
  340 DD=RHO*RHO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SD=0.0
C      SS=0.0
      SD=0.0D0
      SS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
          IF (DABS(DX(I)) >= 1.0D-6*RHO) DD=DD-DX(I)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SD=SD+DX(I)*SDIRN(I)
          SS=SS+SDIRN(I)**2
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (DD .LE. 0.0) GOTO 490
C      TEMP=SQRT(SS*DD)
C      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
      IF (DD <= 0.0D0) GOTO 490
      TEMP=DSQRT(SS*DD)
      IF (DABS(SD) >= 1.0D-6*TEMP) TEMP=DSQRT(SS*DD+SD*SD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      STPFUL=DD/(TEMP+SD)
      STEP=STPFUL
      IF (MCON == M) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          ACCA=STEP+0.1*RESMAX
C          ACCB=STEP+0.2*RESMAX
          ACCA=STEP+0.1D0*RESMAX
          ACCB=STEP+0.2D0*RESMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (STEP >= ACCA .OR. ACCA >= ACCB) GOTO 480
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          STEP=AMIN1(STEP,RESMAX)
          STEP=DMIN1(STEP,RESMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
C
C     Set DXNEW to the new variables if STEP is the steplength, and reduce
C     RESMAX to the corresponding maximum residual if stage one is being done.
C     Because DXNEW will be changed during the calculation of some Lagrange
C     multipliers, it will be restored to the following value later.
C
      DO I=1,N
          DXNEW(I)=DX(I)+STEP*SDIRN(I)
      END DO
      IF (MCON == M) THEN
          RESOLD=RESMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          RESMAX=0.0
          RESMAX=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO K=1,NACT
              KK=IACT(K)
              TEMP=B(KK)
              DO I=1,N
                  TEMP=TEMP-A(I,KK)*DXNEW(I)
              END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          RESMAX=AMAX1(RESMAX,TEMP)
              RESMAX=DMAX1(RESMAX,TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END DO
      END IF
C
C     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
C     device is included to force VMULTD(K)=0.0 if deviations from this value
C     can be attributed to computer rounding errors. First calculate the new
C     Lagrange multipliers.
C
      K=NACT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  390 ZDOTW=0.0
C      ZDWABS=0.0
  390 ZDOTW=0.0D0
      ZDWABS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=Z(I,K)*DXNEW(I)
          ZDOTW=ZDOTW+TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  400 ZDWABS=ZDWABS+ABS(TEMP)
C      ACCA=ZDWABS+0.1*ABS(ZDOTW)
C      ACCB=ZDWABS+0.2*ABS(ZDOTW)
C      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0
          ZDWABS=ZDWABS+DABS(TEMP)
      END DO
      ACCA=ZDWABS+0.1D0*DABS(ZDOTW)
      ACCB=ZDWABS+0.2D0*DABS(ZDOTW)
      IF (ZDWABS >= ACCA .OR. ACCA >= ACCB) ZDOTW=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VMULTD(K)=ZDOTW/ZDOTA(K)
      IF (K >= 2) THEN
          KK=IACT(K)
          DO I=1,N
              DXNEW(I)=DXNEW(I)-VMULTD(K)*A(I,KK)
          END DO
          K=K-1
          GOTO 390
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (MCON .GT. M) VMULTD(NACT)=AMAX1(0.0,VMULTD(NACT))
      IF (MCON > M) VMULTD(NACT)=DMAX1(0.0D0,VMULTD(NACT))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Complete VMULTC by finding the new constraint residuals.
C
      DO I=1,N
          DXNEW(I)=DX(I)+STEP*SDIRN(I)
      END DO
      IF (MCON > NACT) THEN
          KL=NACT+1
          DO K=KL,MCON
              KK=IACT(K)
              SUM=RESMAX-B(KK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          SUMABS=RESMAX+ABS(B(KK))
              SUMABS=RESMAX+DABS(B(KK))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO I=1,N
                  TEMP=A(I,KK)*DXNEW(I)
                  SUM=SUM+TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  430     SUMABS=SUMABS+ABS(TEMP)
C          ACCA=SUMABS+0.1*ABS(SUM)
C          ACCB=SUMABS+0.2*ABS(SUM)
C          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
                  SUMABS=SUMABS+DABS(TEMP)
              END DO
              ACCA=SUMABS+0.1D0*DABS(SUM)
              ACCB=SUMABS+0.2D0*DABS(SUM)
              IF (SUMABS >= ACCA .OR. ACCA >= ACCB) SUM=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              VMULTD(K)=SUM
          END DO
      END IF
C
C     Calculate the fraction of the step from DX to DXNEW that will be taken.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RATIO=1.0
      RATIO=1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICON=0
      DO K=1,MCON
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (VMULTD(K) .LT. 0.0) THEN
          IF (VMULTD(K) < 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TEMP=VMULTC(K)/(VMULTC(K)-VMULTD(K))
              IF (TEMP < RATIO) THEN
                  RATIO=TEMP
                  ICON=K
              END IF
          END IF
      END DO
C
C     Update DX, VMULTC and RESMAX.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=1.0-RATIO
      TEMP=1.0D0-RATIO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          DX(I)=TEMP*DX(I)+RATIO*DXNEW(I)
      END DO
      DO K=1,MCON
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
          VMULTC(K)=DMAX1(0.0D0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (MCON == M) RESMAX=RESOLD+RATIO*(RESMAX-RESOLD)
C
C     If the full step is not acceptable then begin another iteration.
C     Otherwise switch to stage two or end the calculation.
C
      IF (ICON > 0) GOTO 70
      IF (STEP == STPFUL) GOTO 500
  480 MCON=M+1
      ICON=MCON
      IACT(MCON)=MCON
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      VMULTC(MCON)=0.0
      VMULTC(MCON)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      GOTO 60
C
C     We employ any freedom that may be available to reduce the objective
C     function before returning a DX whose length is less than RHO.
C
  490 IF (MCON == M) GOTO 480
      IFULL=0
  500 RETURN
      END
