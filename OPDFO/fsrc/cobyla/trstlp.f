      SUBROUTINE TRSTLP (N,M,A,B,RHO,DX,IFULL,IACT,Z,ZDOTA,VMULTC,
     1  SDIRN,DXNEW,VMULTD)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!-----------------------!!!!!!
      USE DIRTY_TEMPORARY_MOD4POWELL_MOD!
      !!!!!!-----------------------!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION A(N,*),B(*),DX(*),IACT(*),Z(N,*),ZDOTA(*),
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*)
     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*),dold(N)
      integer :: stage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This subroutine calculates an N-component vector DX by applying the
!     following two stages. In the first stage, DX is set to the shortest
!     vector that minimizes the greatest violation of the constraints
!       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
!     subject to the Euclidean length of DX being at most RHO. If its length is
!     strictly less than RHO, then we use the resultant freedom in DX to
!     minimize the objective function
!              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
!     subject to no increase in any greatest constraint violation. This
!     notation allows the gradient of the objective function to be regarded as
!     the gradient of a constraint. Therefore the two stages are distinguished
!     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
!     degeneracy may prevent DX from attaining the target length RHO. Then the
!     value IFULL=0 would be set, but usually IFULL=1 on return.
!
!     In general NACT is the number of constraints in the active set and
!     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
!     contains a permutation of the remaining constraint indices. Further, Z is
!     an orthogonal matrix whose first NACT columns can be regarded as the
!     result of Gram-Schmidt applied to the active constraint gradients. For
!     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
!     column of Z with the gradient of the J-th active constraint. DX is the
!     current vector of variables and here the residuals of the active
!     constraints should be zero. Further, the active constraints have
!     nonnegative Lagrange multipliers that are held at the beginning of
!     VMULTC. The remainder of this vector holds the residuals of the inactive
!     constraints at DX, the ordering of the components of VMULTC being in
!     agreement with the permutation of the indices of the constraints that is
!     in IACT. All these residuals are nonnegative, which is achieved by the
!     shift RESMAX that makes the least residual zero.
!
!     Initialize Z and some other variables. The value of RESMAX will be
!     appropriate to DX=0, while ICON will be the index of a most violated
!     constraint if RESMAX is positive. Usually during the first stage the
!     vector SDIRN gives a search direction that reduces all the active
!     constraint violations by one simultaneously.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      stage = 1
      ITERC=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IFULL=1
      MCON=M
      NACT=0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESMAX=0.0
      RESMAX=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          DO J=1,N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   10 Z(I,J)=0.0
!      Z(I,I)=1.0
!   20 DX(I)=0.0
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


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RESMAX .EQ. 0.0) GOTO 480
!      IF (RESMAX == 0.0D0) GOTO 480
      IF (RESMAX <= 0.0D0) GOTO 480
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   50 SDIRN(I)=0.0
          SDIRN(I)=0.0D0
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      !DO I = 1, N
      !   DSAV(I) = DX(I)
      !END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     End the current stage of the calculation if 3 consecutive iterations
!     have either failed to reduce the best calculated value of the objective
!     function or to increase the number of active constraints since the best
!     value was calculated. This strategy prevents cycling, but there is a
!     remote possibility that it will cause premature termination.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   60 OPTOLD=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICOUNT=0

   70 IF (MCON == M) THEN
!      IF (MCON == M) THEN
          OPTNEW=RESMAX
      ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          OPTNEW=0.0
          OPTNEW=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              OPTNEW=OPTNEW-DX(I)*A(I,MCON)
          END DO
      END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019
! The original code can still encounter infinite cycling,
! which did happen when testing the following CUTEst problems:
! DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in A (up to 10^219).
! To avoid wasting energy, we do the following.
      !SUMD = 0.0D0
      !DO I = 1, N
      !    SUMD = SUMD + DABS(DX(I))
      !END DO
      !IF (SUMD > HUGE(SUMD) .OR. SUMD /= SUMD) THEN
      !    DO I = 1, N
      !       DX(I) = DSAV(I)
      !    END DO
      !    GOTO 490
      !ELSE
      !    DO I = 1, N
      !       DSAV(I) = DX(I)
      !    END DO
      !END IF
      ITERC = ITERC + 1
      IF (ITERC > Min(50000, 100*max(M, N))) THEN
          GOTO 490
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
!
!     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
!     the active set. Apply Givens rotations so that the last N-NACT-1 columns
!     of Z are orthogonal to the gradient of the new constraint, a scalar
!     product being set to zero if its nonzero value could be due to computer
!     rounding errors. The array DXNEW is used for working space.
!
      IF (ICON <= NACT) GOTO 260
      KK=IACT(ICON)
      DO I=1,N
          DXNEW(I)=A(I,KK)
      END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TOT=0.0
      TOT=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      K=N
  100 IF (K > NACT) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
!          SPABS=0.0
          SP=0.0D0
          SPABS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              TEMP=Z(I,K)*DXNEW(I)
              SP=SP+TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  110     SPABS=SPABS+ABS(TEMP)
!          ACCA=SPABS+0.1*ABS(SP)
!          ACCB=SPABS+0.2*ABS(SP)
!          IF (SPABS .GE. ACCA .OR. ACCA .GE. ACCB) SP=0.0
!          IF (TOT .EQ. 0.0) THEN
              SPABS=SPABS+DABS(TEMP)
          END DO
          ACCA=SPABS+0.1D0*DABS(SP)
          ACCB=SPABS+2.0D0*0.1D0*DABS(SP)
          IF (SPABS >= ACCA .OR. ACCA >= ACCB) SP=0.0D0
          !IF (TOT == 0.0D0) THEN
          IF (TOT == 0.0D0 .or. TOT /= TOT) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TOT=SP
          ELSE
              KP=K+1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=SQRT(SP*SP+TOT*TOT)
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
!
!     Add the new constraint if this can be done without a deletion from the
!     active set.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TOT .NE. 0.0) THEN
      !IF (TOT /= 0.0D0) THEN
      IF (ABS(TOT) > 0.0D0) THEN  ! When TOT = NaN, the conditions work differently
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          NACT=NACT+1
          ZDOTA(NACT)=TOT
          VMULTC(ICON)=VMULTC(NACT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTC(NACT)=0.0
          VMULTC(NACT)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 210
      END IF

!!!!!!>>>>>>>>>>>> Zaikun 20211112
      !if (nact > 0) then
      !    zdasav = zdota(nact)
      !    TEMP=0.0D0
      !    DO I=1,N
      !        TEMP=TEMP+Z(I,NACT)*A(I,KK)
      !    END DO
      !    if (isminor(temp,inprod(abs(Z(1:N,nact)),abs(A(1:N,kk)))))then
      !        temp = 0.0d0
      !    end if
      ! IF (.not. ABS(temp) > 0.0D0) THEN
      !    GOTO 490
      !  end if
      !end if

      if (nact == 0) GOTO 490
!!!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!     The next instruction is reached if a deletion has to be made from the
!     active set in order to make room for the new active constraint, because
!     the new constraint gradient is a linear combination of the gradients of
!     the old active constraints. Set the elements of VMULTD to the multipliers
!     of the linear combination. Further, set IOUT to the index of the
!     constraint to be deleted, but branch if no suitable index can be found.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=-1.0
      RATIO=-1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      K=NACT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  130 ZDOTV=0.0
!      ZDVABS=0.0
  130 ZDOTV=0.0D0
      ZDVABS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=Z(I,K)*DXNEW(I)
          ZDOTV=ZDOTV+TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  140 ZDVABS=ZDVABS+ABS(TEMP)
!      ACCA=ZDVABS+0.1*ABS(ZDOTV)
!      ACCB=ZDVABS+0.2*ABS(ZDOTV)
          ZDVABS=ZDVABS+DABS(TEMP)
      END DO
      ACCA=ZDVABS+0.1D0*DABS(ZDOTV)
      ACCB=ZDVABS+2.0D0*0.1D0*DABS(ZDOTV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (ZDVABS < ACCA .AND. ACCA < ACCB) THEN
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          TEMP=ZDOTV/ZDOTA(K)
          !!!temp=zdotv/inprod(A(:, iact(k)), z(:, k))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .GT. 0.0 .AND. IACT(K) .LE. M) THEN
          IF (TEMP > 0.0D0 .AND. IACT(K) <= M) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TEMPA=VMULTC(K)/TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              IF (RATIO .LT. 0.0 .OR. TEMPA .LT. RATIO) THEN
              IF (RATIO < 0.0D0 .OR. TEMPA < RATIO) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  RATIO=TEMPA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: IOUT is never used
!                  IOUT=K
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTD(K)=0.0
          VMULTD(K)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
      K=K-1
      IF (K > 0) GOTO 130

      !write(17,*) 'vmultd1', iterc, vmultd(1:nact)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RATIO .LT. 0.0) GOTO 490
      IF (RATIO < 0.0D0) GOTO 490
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Revise the Lagrange multipliers and reorder the active constraints so
!     that the one to be replaced is at the end of the list. Also calculate the
!     new value of ZDOTA(NACT) and branch if it is not acceptable.
!
      DO K=1,NACT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  160 VMULTC(K)=AMAX1(0.0,VMULTC(K)-RATIO*VMULTD(K))
          VMULTC(K)=DMAX1(0.0D0,VMULTC(K)-RATIO*VMULTD(K))
      END DO

!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> This part cannot be reached!
      IF (ICON < NACT) THEN
          stop "icon<nact" !!!!
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  170     KP=K+1
          KW=IACT(KP)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          SP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SP=SP+Z(I,K)*A(I,KW)
          END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
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
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< This part cannot be reached
!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=TEMP+Z(I,NACT)*A(I,KK)
      END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TEMP .EQ. 0.0) GOTO 490
!      IF (TEMP == 0.0D0) GOTO 490
!      IF (.not. ABS(TEMP) > 0.0D0) GOTO 490
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Zaikun 20211112
      if (isminor(temp,inprod(abs(Z(1:N,nact)),abs(A(1:N,kk)))))then
          temp = 0.0d0
      end if
      IF (.not. ABS(TEMP) > 0.0D0) GOTO 490
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ZDOTA(NACT)=TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(ICON)=0.0
      VMULTC(ICON)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VMULTC(NACT)=RATIO
!
!     Update IACT and ensure that the objective function continues to be
!     treated as the last active constraint when MCON>M.
!
  210 IACT(ICON)=IACT(NACT)
      IACT(NACT)=KK
      IF (MCON > M .AND. KK /= MCON) THEN
          K=NACT-1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          SP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SP=SP+Z(I,K)*A(I,KK)
          END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
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
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
      IF (MCON > M) GOTO 320
      KK=IACT(NACT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=TEMP+SDIRN(I)*A(I,KK)
      END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=TEMP-1.0
      TEMP=TEMP-1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TEMP=TEMP/ZDOTA(NACT)
      DO I=1,N
          SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT)
      END DO
      GOTO 340
!
!     Delete the constraint that has the index IACT(ICON) from the active set.
!
  260 IF (ICON < NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  270     KP=K+1
          KK=IACT(KP)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          SP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SP=SP+Z(I,K)*A(I,KK)
          END DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
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
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
      IF (MCON > M) GOTO 320
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=TEMP+SDIRN(I)*Z(I,NACT+1)
      END DO
      DO I=1,N
          SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT+1)
      END DO
      GO TO 340
!
!     Pick the next search direction of stage two.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  320 TEMP=1.0/ZDOTA(NACT)
  320 TEMP=1.0D0/ZDOTA(NACT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          SDIRN(I)=TEMP*Z(I,NACT)
      END DO
      !write(17,*) iterc, 'vmultc', vmultc(1:mcon)
!
!     Calculate the step to the boundary of the trust region or take the step
!     that reduces RESMAX to zero. The two statements below that include the
!     factor 1.0E-6 prevent some harmless underflows that occurred in a test
!     calculation. Further, we skip the step if it could be zero within a
!     reasonable tolerance for computer rounding errors.
!
!  340 DD=RHO*RHO
  340 DD=0.0D0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SD=0.0
!      SS=0.0
      SD=0.0D0
      SS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
          !IF (DABS(DX(I)) >= 1.0D-6*RHO) DD=DD+DX(I)**2
          IF (DABS(DX(I)) >= EPS*RHO) DD=DD+DX(I)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SD=SD+DX(I)*SDIRN(I)
          SS=SS+SDIRN(I)**2
      END DO
      DD = RHO**2-DD
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (DD .LE. 0.0) GOTO 490
!      TEMP=SQRT(SS*DD)
!      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
      IF (DD <= 0.0D0) GOTO 490
      TEMP=DSQRT(SS*DD)
      !IF (DABS(SD) >= 1.0D-6*TEMP) TEMP=DSQRT(SS*DD+SD*SD)
!      IF (DABS(SD) >= EPS*TEMP) TEMP=DSQRT(SS*DD+SD*SD)
      IF (DABS(SD) >= EPS*TEMP) TEMP=DSQRT(SS*DD+SD**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      STPFUL=DD/(TEMP+SD)
      STEP=STPFUL
      IF (MCON == M) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          ACCA=STEP+0.1*RESMAX
!          ACCB=STEP+0.2*RESMAX
          ACCA=STEP+0.1D0*RESMAX
          ACCB=STEP+2.0D0*0.1D0*RESMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (STEP >= ACCA .OR. ACCA >= ACCB) THEN
              GOTO 480
          END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          STEP=AMIN1(STEP,RESMAX)
          STEP=DMIN1(STEP,RESMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
!
!     Set DXNEW to the new variables if STEP is the steplength, and reduce
!     RESMAX to the corresponding maximum residual if stage one is being done.
!     Because DXNEW will be changed during the calculation of some Lagrange
!     multipliers, it will be restored to the following value later.
!
      DO I=1,N
          DXNEW(I)=DX(I)+STEP*SDIRN(I)
      END DO
      IF (MCON == M) THEN
          RESOLD=RESMAX
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=0.0
          RESMAX=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO K=1,NACT
              KK=IACT(K)
              !TEMP=B(KK)
              TEMP=0.0D0
              DO I=1,N
                  !TEMP=TEMP-A(I,KK)*DXNEW(I)
                  TEMP=TEMP+A(I,KK)*DXNEW(I)
              END DO
              TEMP = B(KK) - TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=AMAX1(RESMAX,TEMP)
              RESMAX=DMAX1(RESMAX,TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END DO

          !resmax = maxval([b(1:m)-matprod(dxnew(1:n),A(1:n,1:m)),0.0D0])

      END IF

      !write(17,*), mcon-m+1, iterc, 'oovd', vmultd(1:mcon)
!     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
!     device is included to force VMULTD(K)=0.0 if deviations from this value
!     can be attributed to computer rounding errors. First calculate the new
!     Lagrange multipliers.

      K=NACT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  390 ZDOTW=0.0
!      ZDWABS=0.0
  390 ZDOTW=0.0D0
      ZDWABS=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          TEMP=Z(I,K)*DXNEW(I)
          ZDOTW=ZDOTW+TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  400 ZDWABS=ZDWABS+ABS(TEMP)
!      ACCA=ZDWABS+0.1*ABS(ZDOTW)
!      ACCB=ZDWABS+0.2*ABS(ZDOTW)
!      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0
          ZDWABS=ZDWABS+DABS(TEMP)
      END DO
      ACCA=ZDWABS+0.1D0*DABS(ZDOTW)
      ACCB=ZDWABS+2.0D0*0.1D0*DABS(ZDOTW)
!      IF (ZDWABS >= ACCA .OR. ACCA >= ACCB) ZDOTW=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      VMULTD(K)=ZDOTW/ZDOTA(K)
      IF (ZDWABS >= ACCA .OR. ACCA >= ACCB) then
        VMULTD(K)=0.0D0
      else
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        VMULTD(K)=ZDOTW/ZDOTA(K)
        !!!vmultd(k)=zdotw/inprod(A(:, iact(k)), Z(:, k))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
      IF (K >= 2) THEN
          KK=IACT(K)
          DO I=1,N
              DXNEW(I)=DXNEW(I)-VMULTD(K)*A(I,KK)
          END DO
          K=K-1
          GOTO 390
      END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (MCON .GT. M) VMULTD(NACT)=AMAX1(0.0,VMULTD(NACT))
      IF (MCON > M) VMULTD(NACT)=DMAX1(0.0D0,VMULTD(NACT))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !write(17,*) 'vmultd21', iterc, vmultd(1:nact)
!
!     Complete VMULTC by finding the new constraint residuals.
!
      DO I=1,N
          DXNEW(I)=DX(I)+STEP*SDIRN(I)
      END DO
      !write(17,*) iterc, 'ovd', vmultd(1:MCON)
      !write(17,*) iterc, nact, mcon
      !write(17,*) 'ad', A(1:n, 1:mcon)
      !write(17,*) 'ad', dxnew(1:n)
      !write(17,*) 'b', b(1:mcon)
      IF (MCON > NACT) THEN
          KL=NACT+1
          DO K=KL,MCON
              KK=IACT(K)
              !SUM=RESMAX-B(KK)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SUMABS=RESMAX+ABS(B(KK))
              !SUMABS=RESMAX+DABS(B(KK))
              SUM = 0.0D0
              SUMABS = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO I=1,N
                  TEMP=A(I,KK)*DXNEW(I)
                  SUM=SUM+TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  430     SUMABS=SUMABS+ABS(TEMP)
!          ACCA=SUMABS+0.1*ABS(SUM)
!          ACCB=SUMABS+0.2*ABS(SUM)
!          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
                  SUMABS=SUMABS+DABS(TEMP)
              END DO
              SUM = SUM - B(KK) + RESMAX
              SUMABS = SUMABS + ABS(B(KK)) + RESMAX
              ACCA=SUMABS+0.1D0*DABS(SUM)
              ACCB=SUMABS+2.0D0*0.1D0*DABS(SUM)
              IF (SUMABS >= ACCA .OR. ACCA >= ACCB) SUM=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              VMULTD(K)=SUM
          END DO
      END IF
      !write(17,*) 'a', iterc, a(1:n, 1:mcon)
      !write(17,*) 'b', iterc, b(1:mcon)
      !write(17,*) 'c', iterc, resmax
      !write(17,*) 'iact', iterc, iact(1:mcon)
      !write(17,*) 'vmultd22', iterc, vmultd(nact+1:mcon)
!
!     Calculate the fraction of the step from DX to DXNEW that will be taken.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=1.0
      RATIO=1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICON=0
      !write(17,*) iterc, 'vc', vmultc(1:MCON)
!      write(17,*) iterc, 'vd', vmultd(nact+1:MCON)
      !write(17,*) iterc, 'vd',
      !1 matprod(dxnew(1:n), A(1:n, iact(nact+1:mcon)))
      !1 - b(iact(nact+1:mcon)) + resmax
      DO K=1,MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (VMULTD(K) .LT. 0.0) THEN
          IF (VMULTD(K) < 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TEMP=VMULTC(K)/(VMULTC(K)-VMULTD(K))
              IF (TEMP < RATIO) THEN
                  RATIO=TEMP
                  ICON=K
              END IF
          END IF
      END DO
!
!     Update DX, VMULTC and RESMAX.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0-RATIO
      TEMP=1.0D0-RATIO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dold = dx(1:n)
      !write(17,*) 'iter', iterc, RATIO, DX(1:N), DXNEW(1:N)
      DO I=1,N
          DX(I)=TEMP*DX(I)+RATIO*DXNEW(I)
      END DO
      !write(17,*) 'iterd', iterc, DX(1:N)

      ! Exit in case of Inf/NaN in D.
      sumd  = 0.0D0
      do i = 1, n
          sumd = sumd + abs(dx(i))
      end do
      if (sumd /= sumd .or. sumd > huge(0.0d0)) then
        dx(1:n) = dold(1:n)
        goto 490
      end if


      DO K=1,MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
          VMULTC(K)=DMAX1(0.0D0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !IF (MCON == M) RESMAX=RESOLD+RATIO*(RESMAX-RESOLD)
      IF (MCON == M) RESMAX=TEMP*RESOLD+RATIO*RESMAX
      !if (mcon == m)  then
      !    resmax = maxval([b(1:m) - matprod(dx(1:n), A(:, 1:m)), 0.0D0])
      !end if
!
!     If the full step is not acceptable then begin another iteration.
!     Otherwise switch to stage two or end the calculation.


      IF (ICON > 0) GOTO 70
      !IF (STEP == STPFUL) THEN
      if (mcon > m .or. dot_product(dx(1:n),dx(1:n))>=rho**2) goto 500
480   MCON=M+1
      ICON=MCON
      IACT(MCON)=MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(MCON)=0.0
      VMULTC(MCON)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !write(17,*) 'stage', stage, 'd', dx(1:n)
      stage = 2  !!! Stage 2 starts here.
      ITERC = 0
      resmax = maxval([b(1:m) - matprod(dx(1:n), A(:, 1:m)), 0.0D0])  ! Do not use MATMUL! Do not use ZERO unless it is defined!
      zdota(1:nact) = [(inprod(z(:, i), A(:, iact(i))), i=1, nact)]
!      write(17,*) 'stage2', vmultc(1:m)
!      iact(1:mcon), nact, A(1:n, 1:mcon), b(1:mcon), rho,
      GOTO 60
!
!     We employ any freedom that may be available to reduce the objective
!     function before returning a DX whose length is less than RHO.
!
  490 IF (MCON == M) then
      GOTO 480
      end if
      IFULL=0
  500 continue
      !write(17,*) 'stage', stage, 'd', dx(1:n)
      RETURN
      END
