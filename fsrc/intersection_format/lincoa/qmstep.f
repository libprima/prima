CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: B is never used 
C      SUBROUTINE QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,
      SUBROUTINE QMSTEP (N,NPT,M,AMAT,XPT,XOPT,NACT,IACT,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     1  RESCON,QFAC,KOPT,KNEW,DEL,STEP,GL,PQW,RSTAT,W,IFEAS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),XOPT(*),IACT(*),
      DIMENSION AMAT(N,*),XPT(NPT,*),XOPT(*),IACT(*),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     1  RESCON(*),QFAC(N,*),STEP(*),GL(*),PQW(*),RSTAT(*),W(*)
C
C     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
C       same as the terms with these names in SUBROUTINE LINCOB.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DEL is the current restriction on the length of STEP, which is never
C       greater than the current trust region radius DELTA.
C     STEP will be set to the required step from XOPT to the new point.
C     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
C       is the KNEW-th Lagrange function. It is used also for some other
C       gradients of LFUNC.
C     PQW provides the second derivative parameters of LFUNC.
C     RSTAT and W are used for working space. Their lengths must be at least
C       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
C       J-th constraint is irrelevant, active, or both inactive and relevant,
C       respectively.
C     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
C
C     STEP is chosen to provide a relatively large value of the modulus of
C       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
C       calculated too, within the trust region, that does not alter the
C       residuals of the active constraints. The projected step is preferred
C       if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the
C       original one, but the greatest violation of a linear constraint must
C       be at least 0.2*DEL, in order to keep the interpolation points apart.
C       The remedy when the maximum constraint violation is too small is to
C       restore the original step, which is perturbed if necessary so that
C       its maximum constraint violation becomes 0.2*DEL.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      TEST=0.2D0*DEL
C
C     Replace GL by the gradient of LFUNC at the trust region centre, and
C       set the elements of RSTAT.
C
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*XOPT(J)
          END DO
          TEMP=PQW(K)*TEMP
          DO I=1,N
              GL(I)=GL(I)+TEMP*XPT(K,I)
          END DO
      END DO
      IF (M > 0) THEN
          DO J=1,M
              RSTAT(J)=ONE
              IF (DABS(RESCON(J)) >= DEL) RSTAT(J)=-ONE
          END DO
          DO K=1,NACT
              RSTAT(IACT(K))=ZERO
          END DO
      END IF
C
C     Find the greatest modulus of LFUNC on a line through XOPT and
C       another interpolation point within the trust region.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-15: IFLAG is never used 
C      IFLAG=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VBIG=ZERO
      DO K=1,NPT
          IF (K == KOPT) CYCLE 
          SS=ZERO
          SP=ZERO
          DO I=1,N
              TEMP=XPT(K,I)-XOPT(I)
              SS=SS+TEMP*TEMP
              SP=SP+GL(I)*TEMP
          END DO
          STP=-DEL/DSQRT(SS)
          IF (K == KNEW) THEN
              IF (SP*(SP-ONE) < ZERO) STP=-STP
              VLAG=DABS(STP*SP)+STP*STP*DABS(SP-ONE)
          ELSE
              VLAG=DABS(STP*(ONE-STP)*SP)
          END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: With the original code, if either VLAG or VBIG is 
C NaN, KSAV will not get a value. This may cause Segmentation Fault
C because XPT(KSAV, :) will later be accessed. 
C      IF (VLAG .GT. VBIG) THEN
          IF (.NOT. (VLAG <= VBIG)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              KSAV=K
              STPSAV=STP
              VBIG=VLAG
          END IF
      END DO
C
C     Set STEP to the move that gives the greatest modulus calculated above.
C       This move may be replaced by a steepest ascent step from XOPT.
C
      GG=ZERO
      DO I=1,N
          GG=GG+GL(I)**2
          STEP(I)=STPSAV*(XPT(KSAV,I)-XOPT(I))
      END DO
      VGRAD=DEL*DSQRT(GG)
      IF (VGRAD <= TENTH*VBIG) GOTO 220
C
C     Make the replacement if it provides a larger value of VBIG.
C
      GHG=ZERO
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*GL(J)
          END DO
          GHG=GHG+PQW(K)*TEMP*TEMP
      END DO
      VNEW=VGRAD+DABS(HALF*DEL*DEL*GHG/GG)
      IF (VNEW > VBIG) THEN
          VBIG=VNEW
          STP=DEL/DSQRT(GG)
          IF (GHG < ZERO) STP=-STP
          DO I=1,N
              STEP(I)=STP*GL(I)
          END DO
      END IF
      IF (NACT == 0 .OR. NACT == N) GOTO 220
C
C     Overwrite GL by its projection. Then set VNEW to the greatest
C       value of |LFUNC| on the projected gradient from XOPT subject to
C       the trust region bound. If VNEW is sufficiently large, then STEP
C       may be changed to a move along the projected gradient.
C
      DO K=NACT+1,N
          W(K)=ZERO
          DO I=1,N
              W(K)=W(K)+GL(I)*QFAC(I,K)
          END DO
      END DO
      GG=ZERO
      DO I=1,N
          GL(I)=ZERO
          DO K=NACT+1,N
              GL(I)=GL(I)+QFAC(I,K)*W(K)
          END DO
          GG=GG+GL(I)**2
      END DO
      VGRAD=DEL*DSQRT(GG)
      IF (VGRAD <= TENTH*VBIG) GOTO 220
      GHG=ZERO
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*GL(J)
          END DO
          GHG=GHG+PQW(K)*TEMP*TEMP
      END DO
      VNEW=VGRAD+DABS(HALF*DEL*DEL*GHG/GG)
C
C     Set W to the possible move along the projected gradient.
C
      STP=DEL/DSQRT(GG)
      IF (GHG < ZERO) STP=-STP
      WW=ZERO
      DO I=1,N
          W(I)=STP*GL(I)
          WW=WW+W(I)**2
      END DO
C
C     Set STEP to W if W gives a sufficiently large value of the modulus
C       of the Lagrange function, and if W either preserves feasibility
C       or gives a constraint violation of at least 0.2*DEL. The purpose
C       of CTOL below is to provide a check on feasibility that includes
C       a tolerance for contributions from computer rounding errors.
C
      IF (VNEW/VBIG >= 0.2D0) THEN
          IFEAS=1
          BIGV=ZERO
          J=0
  170     J=J+1
          IF (J <= M) THEN
              IF (RSTAT(J) == ONE) THEN
                  TEMP=-RESCON(J)
                  DO I=1,N
                      TEMP=TEMP+W(I)*AMAT(I,J)
                  END DO
                  BIGV=DMAX1(BIGV,TEMP)
              END IF
              IF (BIGV < TEST) GOTO 170
              IFEAS=0
          END IF
          CTOL=ZERO
          TEMP=0.01D0*DSQRT(WW)
          IF (BIGV > ZERO .AND. BIGV < TEMP) THEN
              DO K=1,NACT
                  J=IACT(K)
                  SUM=ZERO
                  DO I=1,N
                      SUM=SUM+W(I)*AMAT(I,J)
                  END DO
                  CTOL=DMAX1(CTOL,DABS(SUM))
              END DO
          END IF
          IF (BIGV <= 10.0D0*CTOL .OR. BIGV >= TEST) THEN
              DO I=1,N
                  STEP(I)=W(I)
              END DO
              GOTO 260
          END IF
      END IF
C
C     Calculate the greatest constraint violation at XOPT+STEP with STEP at
C       its original value. Modify STEP if this violation is unacceptable.
C
  220 IFEAS=1
      BIGV=ZERO
      RESMAX=ZERO
      J=0
  230 J=J+1
      IF (J <= M) THEN
          IF (RSTAT(J) < ZERO) GOTO 230
          TEMP=-RESCON(J)
          DO I=1,N
              TEMP=TEMP+STEP(I)*AMAT(I,J)
          END DO
          RESMAX=DMAX1(RESMAX,TEMP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TEMP .LT. TEST) THEN
          IF (.NOT. (TEMP >= TEST)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
              IF (TEMP <= BIGV) GOTO 230
              BIGV=TEMP
              JSAV=J
              IFEAS=-1
              GOTO 230
          END IF
          IFEAS=0
      END IF
      IF (IFEAS == -1) THEN
          DO I=1,N
              STEP(I)=STEP(I)+(TEST-BIGV)*AMAT(I,JSAV)
          END DO
          IFEAS=0
      END IF
C
C     Return the calculated STEP and the value of IFEAS.
C
  260 RETURN
      END
