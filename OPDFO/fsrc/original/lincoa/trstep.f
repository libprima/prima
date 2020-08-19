      SUBROUTINE TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
     1  QFAC,RFAC,SNORM,STEP,G,RESNEW,RESACT,D,DW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
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
      IF (M .GT. 0) THEN
          DO 10 J=1,M
          RESNEW(J)=RESCON(J)
          IF (RESCON(J) .GE. SNORM) THEN
              RESNEW(J)=-ONE
          ELSE IF (RESCON(J) .GE. ZERO) THEN
              RESNEW(J)=DMAX1(RESNEW(J),TINY)
          END IF
   10     CONTINUE
          IF (NACT .GT. 0) THEN
              DO 20 K=1,NACT
              RESACT(K)=RESCON(IACT(K))
   20         RESNEW(IACT(K))=ZERO
          END IF
      END IF
      DO 30 I=1,N
   30 STEP(I)=ZERO
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
      IF (W(N+1) .EQ. ZERO) GOTO 320
      SCALE=0.2D0*SNORM/DSQRT(W(N+1))
      DO 50 I=1,N
   50 DW(I)=SCALE*DW(I)
C
C     If the modulus of the residual of an active constraint is substantial,
C       then set D to the shortest move from STEP to the boundaries of the
C       active constraints.
C
      RESMAX=ZERO
      IF (NACT .GT. 0) THEN
          DO 60 K=1,NACT
   60     RESMAX=DMAX1(RESMAX,RESACT(K))
      END IF
      GAMMA=ZERO
      IF (RESMAX .GT. 1.0D-4*SNORM) THEN
          IR=0
          DO 80 K=1,NACT
          TEMP=RESACT(K)
          IF (K .GE. 2) THEN
              DO 70 I=1,K-1
              IR=IR+1
   70         TEMP=TEMP-RFAC(IR)*W(I)
          END IF
          IR=IR+1
   80     W(K)=TEMP/RFAC(IR)
          DO 90 I=1,N
          D(I)=ZERO
          DO 90 K=1,NACT
   90     D(I)=D(I)+W(K)*QFAC(I,K)
C
C     The vector D that has just been calculated is also the shortest move
C       from STEP+DW to the boundaries of the active constraints. Set GAMMA
C       to the greatest steplength of this move that satisfies the trust
C       region bound.
C
          RHS=SNSQ
          DS=ZERO
          DD=ZERO
          DO 100 I=1,N
          SUM=STEP(I)+DW(I)
          RHS=RHS-SUM*SUM
          DS=DS+D(I)*SUM
  100     DD=DD+D(I)**2
          IF (RHS .GT. ZERO) THEN
              TEMP=DSQRT(DS*DS+DD*RHS)
              IF (DS .LE. ZERO) THEN
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
  110     IF (GAMMA .GT. ZERO) THEN
              J=J+1
              IF (RESNEW(J) .GT. ZERO) THEN
                  AD=ZERO
                  ADW=ZERO
                  DO 120 I=1,N
                  AD=AD+AMAT(I,J)*D(I)
  120             ADW=ADW+AMAT(I,J)*DW(I)
                  IF (AD .GT. ZERO) THEN
                      TEMP=DMAX1((RESNEW(J)-ADW)/AD,ZERO)
                      GAMMA=DMIN1(GAMMA,TEMP)
                  END IF
              END IF
              IF (J .LT. M) GOTO 110
          END IF
          GAMMA=DMIN1(GAMMA,ONE)
      END IF
C
C     Set the next direction for seeking a reduction in the model function
C       subject to the trust region bound and the linear constraints.
C
      IF (GAMMA .LE. ZERO) THEN
          DO 130 I=1,N
  130     D(I)=DW(I)
          ICOUNT=NACT
      ELSE
          DO 140 I=1,N
  140     D(I)=DW(I)+GAMMA*D(I)
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
      IF (RHS .LE. ZERO) GOTO 320
      DG=ZERO
      DS=ZERO
      DD=ZERO
      DO 160 I=1,N
      DG=DG+D(I)*G(I)
      DS=DS+D(I)*STEP(I)
  160 DD=DD+D(I)**2
      IF (DG .GE. ZERO) GOTO 320
      TEMP=DSQRT(RHS*DD+DS*DS)
      IF (DS .LE. ZERO) THEN
          ALPHA=(TEMP-DS)/DD
      ELSE
          ALPHA=RHS/(TEMP+DS)
      END IF
      IF (-ALPHA*DG .LE. CTEST*REDUCT) GOTO 320
C
C     Set DW to the change in gradient along D.
C
      IH=0
      DO 170 J=1,N
      DW(J)=ZERO
      DO 170 I=1,J
      IH=IH+1
      IF (I .LT. J) DW(J)=DW(J)+HQ(IH)*D(I)
  170 DW(I)=DW(I)+HQ(IH)*D(J)
      DO 190 K=1,NPT
      TEMP=ZERO
      DO 180 J=1,N
  180 TEMP=TEMP+XPT(K,J)*D(J)
      TEMP=PQ(K)*TEMP
      DO 190 I=1,N
  190 DW(I)=DW(I)+TEMP*XPT(K,I)
C
C     Set DGD to the curvature of the model along D. Then reduce ALPHA if
C       necessary to the value that minimizes the model.
C
      DGD=ZERO
      DO 200 I=1,N
  200 DGD=DGD+D(I)*DW(I)
      ALPHT=ALPHA
      IF (DG+ALPHA*DGD .GT. ZERO) THEN
          ALPHA=-DG/DGD
      END IF
C
C     Make a further reduction in ALPHA if necessary to preserve feasibility,
C       and put some scalar products of D with constraint gradients in W.
C
      ALPHM=ALPHA
      JSAV=0
      IF (M .GT. 0) THEN
          DO 220 J=1,M
          AD=ZERO
          IF (RESNEW(J) .GT. ZERO) THEN
              DO 210 I=1,N
  210         AD=AD+AMAT(I,J)*D(I)
              IF (ALPHA*AD .GT. RESNEW(J)) THEN
                  ALPHA=RESNEW(J)/AD
                  JSAV=J
              END IF
          END IF
  220     W(J)=AD
      END IF
      ALPHA=DMAX1(ALPHA,ALPBD)
      ALPHA=DMIN1(ALPHA,ALPHM)
      IF (ICOUNT .EQ. NACT) ALPHA=DMIN1(ALPHA,ONE)
C
C     Update STEP, G, RESNEW, RESACT and REDUCT.
C
      SS=ZERO
      DO 230 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      SS=SS+STEP(I)**2
  230 G(I)=G(I)+ALPHA*DW(I)
      IF (M .GT. 0) THEN
          DO 240 J=1,M
          IF (RESNEW(J) .GT. ZERO) THEN
              RESNEW(J)=DMAX1(RESNEW(J)-ALPHA*W(J),TINY)
          END IF
  240     CONTINUE
      END IF
      IF (ICOUNT .EQ. NACT .AND. NACT .GT. 0) THEN
          DO 250 K=1,NACT
  250     RESACT(K)=(ONE-GAMMA)*RESACT(K)
      END IF
      REDUCT=REDUCT-ALPHA*(DG+HALF*ALPHA*DGD)
C
C     Test for termination. Branch to label 40 if there is a new active
C       constraint and if the distance from STEP to the trust region
C       boundary is at least 0.2*SNORM.
C
      IF (ALPHA .EQ. ALPHT) GOTO 320
      TEMP=-ALPHM*(DG+HALF*ALPHM*DGD)
      IF (TEMP .LE. CTEST*REDUCT) GOTO 320
      IF (JSAV .GT. 0) THEN
          IF (SS .LE. 0.64D0*SNSQ) GOTO 40
          GOTO 320
      END IF
      IF (ICOUNT .EQ. N) GOTO 320
C
C     Calculate the next search direction, which is conjugate to the
C       previous one except in the case ICOUNT=NACT.
C
      IF (NACT .GT. 0) THEN
          DO 260 J=NACT+1,N
          W(J)=ZERO
          DO 260 I=1,N
  260     W(J)=W(J)+G(I)*QFAC(I,J)
          DO 280 I=1,N
          TEMP=ZERO
          DO 270 J=NACT+1,N
  270     TEMP=TEMP+QFAC(I,J)*W(J)
  280     W(N+I)=TEMP
      ELSE
          DO 290 I=1,N
  290     W(N+I)=G(I)
      END IF
      IF (ICOUNT .EQ. NACT) THEN
          BETA=ZERO
      ELSE
          WGD=ZERO
          DO 300 I=1,N
  300     WGD=WGD+W(N+I)*DW(I)
          BETA=WGD/DGD
      END IF
      DO 310 I=1,N
  310 D(I)=-W(N+I)+BETA*D(I)
      ALPBD=ZERO
      GOTO 150
C
C     Return from the subroutine.
C
  320 SNORM=ZERO
      IF (REDUCT .GT. ZERO) SNORM=DSQRT(SS)
      G(1)=ZERO
      IF (NCALL .GT. 1) G(1)=ONE
      RETURN
      END
