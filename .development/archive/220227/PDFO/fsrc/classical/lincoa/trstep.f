      SUBROUTINE TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
     1  QFAC,RFAC,SNORM,STEP,G,RESNEW,RESACT,D,DW,W)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),HQ(*),PQ(*),IACT(*),
     1  RESCON(*),QFAC(N,*),RFAC(*),STEP(*),G(*),RESNEW(*),RESACT(*),
     2  D(*),DW(*),W(*)
C
C     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
C       are the same as the terms with these names in LINCOB. If RESCON(J)
C       is negative, then |RESCON(J)| must be no less than the trust region
C       radius, so that the J-th constraint can be ignored.
C     SNORM is set to the trust region radius DELTA initially. On the
C       return, however, it is the length of the calculated STEP, which is
C       set to zero if the constraints do not allow a long enough step.
C     STEP is the total calculated step so far from the trust region centre,
C       its final value being given by the sequence of CG iterations, which
C       terminate if the trust region boundary is reached.
C     G must be set on entry to the gradient of the quadratic model at the
C       trust region centre. It is used as working space, however, and is
C       always the gradient of the model at the current STEP, except that
C       on return the value of G(1) is set to ONE instead of to ZERO if
C       and only if GETACT is called more than once.
C     RESNEW, RESACT, D, DW and W are used for working space. A negative
C       value of RESNEW(J) indicates that the J-th constraint does not
C       restrict the CG steps of the current trust region calculation, a
C       zero value of RESNEW(J) indicates that the J-th constraint is active,
C       and otherwise RESNEW(J) is set to the greater of TINY and the actual
C       residual of the J-th constraint for the current STEP. RESACT holds
C       the residuals of the active constraints, which may be positive.
C       D is the search direction of each line search. DW is either another
C       search direction or the change in gradient along D. The length of W
C       must be at least MAX[M,2*N].
C
C     Set some numbers for the conjugate gradient iterations.
C
      HALF=0.5D0
      ONE=1.0D0
      TINY=1.0D-60
      ZERO=0.0D0
      CTEST=0.01D0
      SNSQ=SNORM*SNORM
C
C     Set the initial elements of RESNEW, RESACT and STEP.
C
      IF (M > 0) THEN
          DO J=1,M
              RESNEW(J)=RESCON(J)
              IF (RESCON(J) >= SNORM) THEN
                  RESNEW(J)=-ONE
              ELSE IF (RESCON(J) >= ZERO) THEN
                  RESNEW(J)=DMAX1(RESNEW(J),TINY)
              END IF
          END DO
          IF (NACT > 0) THEN
              DO K=1,NACT
                  RESACT(K)=RESCON(IACT(K))
                  RESNEW(IACT(K))=ZERO
              END DO
          END IF
      END IF
      DO I=1,N
          STEP(I)=ZERO
      END DO
      SS=ZERO
      REDUCT=ZERO
      NCALL=0
C
C     GETACT picks the active set for the current STEP. It also sets DW to
C       the vector closest to -G that is orthogonal to the normals of the
C       active constraints. DW is scaled to have length 0.2*SNORM, as then
C       a move of DW from STEP is allowed by the linear constraints.
C
   40 NCALL=NCALL+1
      CALL GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,RESNEW,
     1  RESACT,G,DW,W,W(N+1))
      IF (W(N+1) == ZERO) GOTO 320
      SCALE=0.2D0*SNORM/DSQRT(W(N+1))
      DO I=1,N
          DW(I)=SCALE*DW(I)
      END DO
C
C     If the modulus of the residual of an active constraint is substantial,
C       then set D to the shortest move from STEP to the boundaries of the
C       active constraints.
C
      RESMAX=ZERO
      IF (NACT > 0) THEN
          DO K=1,NACT
              RESMAX=DMAX1(RESMAX,RESACT(K))
          END DO
      END IF
      GAMMA=ZERO
      IF (RESMAX > 1.0D-4*SNORM) THEN
          IR=0
          DO K=1,NACT
              TEMP=RESACT(K)
              IF (K >= 2) THEN
                  DO I=1,K-1
                      IR=IR+1
                      TEMP=TEMP-RFAC(IR)*W(I)
                  END DO
              END IF
              IR=IR+1
              W(K)=TEMP/RFAC(IR)
          END DO
          DO I=1,N
              D(I)=ZERO
              DO K=1,NACT
                  D(I)=D(I)+W(K)*QFAC(I,K)
              END DO
          END DO
C
C     The vector D that has just been calculated is also the shortest move
C       from STEP+DW to the boundaries of the active constraints. Set GAMMA
C       to the greatest steplength of this move that satisfies the trust
C       region bound.
C
          RHS=SNSQ
          DS=ZERO
          DD=ZERO
          DO I=1,N
              SUM=STEP(I)+DW(I)
              RHS=RHS-SUM*SUM
              DS=DS+D(I)*SUM
              DD=DD+D(I)**2
          END DO
          IF (RHS > ZERO) THEN
              TEMP=DSQRT(DS*DS+DD*RHS)
              IF (DS <= ZERO) THEN
                  GAMMA=(TEMP-DS)/DD
              ELSE
                  GAMMA=RHS/(TEMP+DS)
              END IF
          END IF
C
C     Reduce the steplength GAMMA if necessary so that the move along D
C       also satisfies the linear constraints.
C
          J=0
  110     IF (GAMMA > ZERO) THEN
              J=J+1
              IF (RESNEW(J) > ZERO) THEN
                  AD=ZERO
                  ADW=ZERO
                  DO I=1,N
                      AD=AD+AMAT(I,J)*D(I)
                      ADW=ADW+AMAT(I,J)*DW(I)
                  END DO
                  IF (AD > ZERO) THEN
                      TEMP=DMAX1((RESNEW(J)-ADW)/AD,ZERO)
                      GAMMA=DMIN1(GAMMA,TEMP)
                  END IF
              END IF
              IF (J < M) GOTO 110
          END IF
          GAMMA=DMIN1(GAMMA,ONE)
      END IF
C
C     Set the next direction for seeking a reduction in the model function
C       subject to the trust region bound and the linear constraints.
C
      IF (GAMMA <= ZERO) THEN
          DO I=1,N
              D(I)=DW(I)
          END DO
          ICOUNT=NACT
      ELSE
          DO I=1,N
              D(I)=DW(I)+GAMMA*D(I)
          END DO
          ICOUNT=NACT-1
      END IF
      ALPBD=ONE
C
C     Set ALPHA to the steplength from STEP along D to the trust region
C       boundary. Return if the first derivative term of this step is
C       sufficiently small or if no further progress is possible.
C
  150 ICOUNT=ICOUNT+1
      RHS=SNSQ-SS
      IF (RHS <= ZERO) GOTO 320
      DG=ZERO
      DS=ZERO
      DD=ZERO
      DO I=1,N
          DG=DG+D(I)*G(I)
          DS=DS+D(I)*STEP(I)
          DD=DD+D(I)**2
      END DO
      IF (DG >= ZERO) GOTO 320
      TEMP=DSQRT(RHS*DD+DS*DS)
      IF (DS <= ZERO) THEN
          ALPHA=(TEMP-DS)/DD
      ELSE
          ALPHA=RHS/(TEMP+DS)
      END IF
      IF (-ALPHA*DG <= CTEST*REDUCT) GOTO 320
C
C     Set DW to the change in gradient along D.
C
      IH=0
      DO J=1,N
          DW(J)=ZERO
          DO I=1,J
              IH=IH+1
              IF (I < J) DW(J)=DW(J)+HQ(IH)*D(I)
              DW(I)=DW(I)+HQ(IH)*D(J)
          END DO
      END DO
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*D(J)
          END DO
          TEMP=PQ(K)*TEMP
          DO I=1,N
              DW(I)=DW(I)+TEMP*XPT(K,I)
          END DO
      END DO
C
C     Set DGD to the curvature of the model along D. Then reduce ALPHA if
C       necessary to the value that minimizes the model.
C
      DGD=ZERO
      DO I=1,N
          DGD=DGD+D(I)*DW(I)
      END DO
      ALPHT=ALPHA
      IF (DG+ALPHA*DGD > ZERO) THEN
          ALPHA=-DG/DGD
      END IF
C
C     Make a further reduction in ALPHA if necessary to preserve feasibility,
C       and put some scalar products of D with constraint gradients in W.
C
      ALPHM=ALPHA
      JSAV=0
      IF (M > 0) THEN
          DO J=1,M
              AD=ZERO
              IF (RESNEW(J) > ZERO) THEN
                  DO I=1,N
                      AD=AD+AMAT(I,J)*D(I)
                  END DO
                  IF (ALPHA*AD > RESNEW(J)) THEN
                      ALPHA=RESNEW(J)/AD
                      JSAV=J
                  END IF
              END IF
              W(J)=AD
          END DO
      END IF
      ALPHA=DMAX1(ALPHA,ALPBD)
      ALPHA=DMIN1(ALPHA,ALPHM)
      IF (ICOUNT == NACT) ALPHA=DMIN1(ALPHA,ONE)
C
C     Update STEP, G, RESNEW, RESACT and REDUCT.
C
      SS=ZERO
      DO I=1,N
          STEP(I)=STEP(I)+ALPHA*D(I)
          SS=SS+STEP(I)**2
          G(I)=G(I)+ALPHA*DW(I)
      END DO
      IF (M > 0) THEN
          DO J=1,M
              IF (RESNEW(J) > ZERO) THEN
                  RESNEW(J)=DMAX1(RESNEW(J)-ALPHA*W(J),TINY)
              END IF
          END DO
      END IF
      IF (ICOUNT == NACT .AND. NACT > 0) THEN
          DO K=1,NACT
              RESACT(K)=(ONE-GAMMA)*RESACT(K)
          END DO
      END IF
      REDUCT=REDUCT-ALPHA*(DG+HALF*ALPHA*DGD)
C
C     Test for termination. Branch to label 40 if there is a new active
C       constraint and if the distance from STEP to the trust region
C       boundary is at least 0.2*SNORM.
C
      IF (ALPHA == ALPHT) GOTO 320
      TEMP=-ALPHM*(DG+HALF*ALPHM*DGD)
      IF (TEMP <= CTEST*REDUCT) GOTO 320
      IF (JSAV > 0) THEN
          IF (SS <= 0.64D0*SNSQ) GOTO 40
          GOTO 320
      END IF
      IF (ICOUNT == N) GOTO 320
C
C     Calculate the next search direction, which is conjugate to the
C       previous one except in the case ICOUNT=NACT.
C
      IF (NACT > 0) THEN
          DO J=NACT+1,N
              W(J)=ZERO
              DO I=1,N
                  W(J)=W(J)+G(I)*QFAC(I,J)
              END DO
          END DO
          DO I=1,N
              TEMP=ZERO
              DO J=NACT+1,N
                  TEMP=TEMP+QFAC(I,J)*W(J)
              END DO
              W(N+I)=TEMP
          END DO
      ELSE
          DO I=1,N
              W(N+I)=G(I)
          END DO
      END IF
      IF (ICOUNT == NACT) THEN
          BETA=ZERO
      ELSE
          WGD=ZERO
          DO I=1,N
              WGD=WGD+W(N+I)*DW(I)
          END DO
          BETA=WGD/DGD
      END IF
      DO I=1,N
          D(I)=-W(N+I)+BETA*D(I)
      END DO
      ALPBD=ZERO
      GOTO 150
C
C     Return from the subroutine.
C
  320 SNORM=ZERO
      IF (REDUCT > ZERO) SNORM=DSQRT(SS)
      G(1)=ZERO
      IF (NCALL > 1) G(1)=ONE
      RETURN
      END
