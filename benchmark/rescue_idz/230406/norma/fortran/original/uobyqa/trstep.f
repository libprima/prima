      SUBROUTINE TRSTEP (N,G,H,DELTA,TOL,D,GG,TD,TN,W,PIV,Z,EVALUE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION G(*),H(N,*),D(*),GG(*),TD(*),TN(*),W(*),PIV(*),Z(*)
C
C     N is the number of variables of a quadratic objective function, Q say.
C     G is the gradient of Q at the origin.
C     H is the Hessian matrix of Q. Only the upper triangular and diagonal
C       parts need be set. The lower triangular part is used to store the
C       elements of a Householder similarity transformation.
C     DELTA is the trust region radius, and has to be positive.
C     TOL is the value of a tolerance from the open interval (0,1).
C     D will be set to the calculated vector of variables.
C     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
C     EVALUE will be set to the least eigenvalue of H if and only if D is a
C     Newton-Raphson step. Then EVALUE will be positive, but otherwise it
C     will be set to zero.
C
C     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| .LEQ. DELTA,
C     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
C     We take the view that any D is acceptable if it has the properties
C
C             ||D|| .LEQ. DELTA  and  ACTRED .LEQ. (1-TOL)*MAXRED.
C
C     The calculation of D is done by the method of Section 2 of the paper
C     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
C     after transforming H to tridiagonal form.
C
C     Initialization.
C
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      DELSQ=DELTA*DELTA
      EVALUE=ZERO
      NM=N-1
      DO 10 I=1,N
      D(I)=ZERO
      TD(I)=H(I,I)
      DO 10 J=1,I
   10 H(I,J)=H(J,I)
C
C     Apply Householder transformations to obtain a tridiagonal matrix that
C     is similar to H, and put the elements of the Householder vectors in
C     the lower triangular part of H. Further, TD and TN will contain the
C     diagonal and other nonzero elements of the tridiagonal matrix.
C
      DO 80 K=1,NM
      KP=K+1
      SUM=ZERO
      IF (KP .LT. N) THEN
          KPP=KP+1
          DO 20 I=KPP,N
   20     SUM=SUM+H(I,K)**2
      END IF
      IF (SUM .EQ. ZERO) THEN
          TN(K)=H(KP,K)
          H(KP,K)=ZERO
      ELSE
          TEMP=H(KP,K)
          TN(K)=DSIGN(DSQRT(SUM+TEMP*TEMP),TEMP)
          H(KP,K)=-SUM/(TEMP+TN(K))
          TEMP=DSQRT(TWO/(SUM+H(KP,K)**2))
          DO 30 I=KP,N
          W(I)=TEMP*H(I,K)
          H(I,K)=W(I)
   30     Z(I)=TD(I)*W(I)
          WZ=ZERO
          DO 50 J=KP,NM
          JP=J+1
          DO 40 I=JP,N
          Z(I)=Z(I)+H(I,J)*W(J)
   40     Z(J)=Z(J)+H(I,J)*W(I)
   50     WZ=WZ+W(J)*Z(J)
          WZ=WZ+W(N)*Z(N)
          DO 70 J=KP,N
          TD(J)=TD(J)+W(J)*(WZ*W(J)-TWO*Z(J))
          IF (J .LT. N) THEN
              JP=J+1
              DO 60 I=JP,N
   60         H(I,J)=H(I,J)-W(I)*Z(J)-W(J)*(Z(I)-WZ*W(I))
          END IF
   70     CONTINUE
      END IF
   80 CONTINUE
C
C     Form GG by applying the similarity transformation to G.
C
      GSQ=ZERO
      DO 90 I=1,N
      GG(I)=G(I)
   90 GSQ=GSQ+G(I)**2
      GNORM=DSQRT(GSQ)
      DO 110 K=1,NM
      KP=K+1
      SUM=ZERO
      DO 100 I=KP,N
  100 SUM=SUM+GG(I)*H(I,K)
      DO 110 I=KP,N
  110 GG(I)=GG(I)-SUM*H(I,K)
C
C     Begin the trust region calculation with a tridiagonal matrix by
C     calculating the norm of H. Then treat the case when H is zero.
C
      HNORM=DABS(TD(1))+DABS(TN(1))
      TDMIN=TD(1)
      TN(N)=ZERO
      DO 120 I=2,N
      TEMP=DABS(TN(I-1))+DABS(TD(I))+DABS(TN(I))
      HNORM=DMAX1(HNORM,TEMP)
  120 TDMIN=DMIN1(TDMIN,TD(I))
      IF (HNORM .EQ. ZERO) THEN
          IF (GNORM .EQ. ZERO) GOTO 400
          SCALE=DELTA/GNORM
          DO 130 I=1,N
  130     D(I)=-SCALE*GG(I)
          GOTO 370
      END IF
C
C     Set the initial values of PAR and its bounds.
C
      PARL=DMAX1(ZERO,-TDMIN,GNORM/DELTA-HNORM)
      PARLEST=PARL
      PAR=PARL
      PARU=ZERO
      PARUEST=ZERO
      POSDEF=ZERO
      ITERC=0
C
C     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
C
  140 ITERC=ITERC+1
      KSAV=0
      PIV(1)=TD(1)+PAR
      K=1
  150 IF (PIV(K) .GT. ZERO) THEN
          PIV(K+1)=TD(K+1)+PAR-TN(K)**2/PIV(K)
      ELSE
          IF (PIV(K) .LT. ZERO .OR. TN(K) .NE. ZERO) GOTO 160
          KSAV=K
          PIV(K+1)=TD(K+1)+PAR
      END IF
      K=K+1
      IF (K .LT. N) GOTO 150
      IF (PIV(K) .LT. ZERO) GOTO 160
      IF (PIV(K) .EQ. ZERO) KSAV=K
C
C     Branch if all the pivots are positive, allowing for the case when
C     G is zero.
C
      IF (KSAV .EQ. 0 .AND. GSQ .GT. ZERO) GOTO 230
      IF (GSQ .EQ. ZERO) THEN
          IF (PAR .EQ. ZERO) GOTO 370
          PARU=PAR
          PARUEST=PAR
          IF (KSAV .EQ. 0) GOTO 190
      END IF
      K=KSAV
C
C     Set D to a direction of nonpositive curvature of the given tridiagonal
C     matrix, and thus revise PARLEST.
C
  160 D(K)=ONE
      IF (DABS(TN(K)) .LE. DABS(PIV(K))) THEN
          DSQ=ONE
          DHD=PIV(K)
      ELSE
          TEMP=TD(K+1)+PAR
          IF (TEMP .LE. DABS(PIV(K))) THEN
              D(K+1)=DSIGN(ONE,-TN(K))
              DHD=PIV(K)+TEMP-TWO*DABS(TN(K))
          ELSE
              D(K+1)=-TN(K)/TEMP
              DHD=PIV(K)+TN(K)*D(K+1)
          END IF
          DSQ=ONE+D(K+1)**2
      END IF
  170 IF (K .GT. 1) THEN
          K=K-1
          IF (TN(K) .NE. ZERO) THEN
              D(K)=-TN(K)*D(K+1)/PIV(K)
              DSQ=DSQ+D(K)**2
              GOTO 170
          END IF
          DO 180 I=1,K
  180     D(I)=ZERO
      END IF
      PARL=PAR
      PARLEST=PAR-DHD/DSQ
C
C     Terminate with D set to a multiple of the current D if the following
C     test suggests that it suitable to do so.
C
  190 TEMP=PARUEST
      IF (GSQ .EQ. ZERO) TEMP=TEMP*(ONE-TOL)
      IF (PARUEST .GT. ZERO .AND. PARLEST .GE. TEMP) THEN
          DTG=ZERO
          DO 200 I=1,N
  200     DTG=DTG+D(I)*GG(I)
          SCALE=-DSIGN(DELTA/DSQRT(DSQ),DTG)
          DO 210 I=1,N
  210     D(I)=SCALE*D(I)
          GOTO 370
      END IF
C
C     Pick the value of PAR for the next iteration.
C
  220 IF (PARU .EQ. ZERO) THEN
          PAR=TWO*PARLEST+GNORM/DELTA
      ELSE
          PAR=0.5D0*(PARL+PARU)
          PAR=DMAX1(PAR,PARLEST)
      END IF
      IF (PARUEST .GT. ZERO) PAR=DMIN1(PAR,PARUEST)
      GOTO 140
C
C     Calculate D for the current PAR in the positive definite case.
C
  230 W(1)=-GG(1)/PIV(1)
      DO 240 I=2,N
  240 W(I)=(-GG(I)-TN(I-1)*W(I-1))/PIV(I)
      D(N)=W(N)
      DO 250 I=NM,1,-1
  250 D(I)=W(I)-TN(I)*D(I+1)/PIV(I)
C
C     Branch if a Newton-Raphson step is acceptable.
C
      DSQ=ZERO
      WSQ=ZERO
      DO 260 I=1,N
      DSQ=DSQ+D(I)**2
  260 WSQ=WSQ+PIV(I)*W(I)**2
      IF (PAR .EQ. ZERO .AND. DSQ .LE. DELSQ) GOTO 320
C
C     Make the usual test for acceptability of a full trust region step.
C
      DNORM=DSQRT(DSQ)
      PHI=ONE/DNORM-ONE/DELTA
      TEMP=TOL*(ONE+PAR*DSQ/WSQ)-DSQ*PHI*PHI
      IF (TEMP .GE. ZERO) THEN
          SCALE=DELTA/DNORM
          DO 270 I=1,N
  270     D(I)=SCALE*D(I)
          GOTO 370
      END IF
      IF (ITERC .GE. 2 .AND. PAR .LE. PARL) GOTO 370
      IF (PARU .GT. ZERO .AND. PAR .GE. PARU) GOTO 370
C
C     Complete the iteration when PHI is negative.
C
      IF (PHI .LT. ZERO) THEN
          PARLEST=PAR
          IF (POSDEF. EQ. ONE) THEN
              IF (PHI .LE. PHIL) GOTO 370
              SLOPE=(PHI-PHIL)/(PAR-PARL)
              PARLEST=PAR-PHI/SLOPE
          END IF
          SLOPE=ONE/GNORM
          IF (PARU .GT. ZERO) SLOPE=(PHIU-PHI)/(PARU-PAR)
          TEMP=PAR-PHI/SLOPE
          IF (PARUEST .GT. ZERO) TEMP=DMIN1(TEMP,PARUEST)
          PARUEST=TEMP
          POSDEF=ONE
          PARL=PAR
          PHIL=PHI
          GOTO 220
      END IF
C
C     If required, calculate Z for the alternative test for convergence.
C
      IF (POSDEF .EQ. ZERO) THEN
          W(1)=ONE/PIV(1)
          DO 280 I=2,N
          TEMP=-TN(I-1)*W(I-1)
  280     W(I)=(DSIGN(ONE,TEMP)+TEMP)/PIV(I)
          Z(N)=W(N)
          DO 290 I=NM,1,-1
  290     Z(I)=W(I)-TN(I)*Z(I+1)/PIV(I)
          WWSQ=ZERO
          ZSQ=ZERO
          DTZ=ZERO
          DO 300 I=1,N
          WWSQ=WWSQ+PIV(I)*W(I)**2
          ZSQ=ZSQ+Z(I)**2
  300     DTZ=DTZ+D(I)*Z(I)
C
C     Apply the alternative test for convergence.
C
          TEMPA=DABS(DELSQ-DSQ)
          TEMPB=DSQRT(DTZ*DTZ+TEMPA*ZSQ)
          GAM=TEMPA/(DSIGN(TEMPB,DTZ)+DTZ)
          TEMP=TOL*(WSQ+PAR*DELSQ)-GAM*GAM*WWSQ
          IF (TEMP .GE. ZERO) THEN
              DO 310 I=1,N
  310         D(I)=D(I)+GAM*Z(I)
              GOTO 370
          END IF
          PARLEST=DMAX1(PARLEST,PAR-WWSQ/ZSQ)
      END IF
C
C     Complete the iteration when PHI is positive.
C
      SLOPE=ONE/GNORM
      IF (PARU .GT. ZERO) THEN
          IF (PHI .GE. PHIU) GOTO 370
          SLOPE=(PHIU-PHI)/(PARU-PAR)
      END IF
      PARLEST=DMAX1(PARLEST,PAR-PHI/SLOPE)
      PARUEST=PAR
      IF (POSDEF .EQ. ONE) THEN
          SLOPE=(PHI-PHIL)/(PAR-PARL)
          PARUEST=PAR-PHI/SLOPE
      END IF
      PARU=PAR
      PHIU=PHI
      GOTO 220
C
C     Set EVALUE to the least eigenvalue of the second derivative matrix if
C     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
C
  320 SHFMIN=ZERO
      PIVOT=TD(1)
      SHFMAX=PIVOT
      DO 330 K=2,N
      PIVOT=TD(K)-TN(K-1)**2/PIVOT
  330 SHFMAX=DMIN1(SHFMAX,PIVOT)
C
C     Find EVALUE by a bisection method, but occasionally SHFMAX may be
C     adjusted by the rule of false position.
C
      KSAVE=0
  340 SHIFT=0.5D0*(SHFMIN+SHFMAX)
      K=1
      TEMP=TD(1)-SHIFT
  350 IF (TEMP .GT. ZERO) THEN
          PIV(K)=TEMP
          IF (K .LT. N) THEN
              TEMP=TD(K+1)-SHIFT-TN(K)**2/TEMP
              K=K+1
              GOTO 350
          END IF
          SHFMIN=SHIFT
      ELSE
          IF (K .LT. KSAVE) GOTO 360
          IF (K .EQ. KSAVE) THEN
              IF (PIVKSV .EQ. ZERO) GOTO 360
              IF (PIV(K)-TEMP .LT. TEMP-PIVKSV) THEN
                  PIVKSV=TEMP
                  SHFMAX=SHIFT
              ELSE
                  PIVKSV=ZERO
                  SHFMAX=(SHIFT*PIV(K)-SHFMIN*TEMP)/(PIV(K)-TEMP)
              END IF
          ELSE
              KSAVE=K
              PIVKSV=TEMP
              SHFMAX=SHIFT
          END IF
      END IF
      IF (SHFMIN .LE. 0.99D0*SHFMAX) GOTO 340
  360 EVALUE=SHFMIN
C
C     Apply the inverse Householder transformations to D.
C
  370 NM=N-1
      DO 390 K=NM,1,-1
      KP=K+1
      SUM=ZERO
      DO 380 I=KP,N
  380 SUM=SUM+D(I)*H(I,K)
      DO 390 I=KP,N
  390 D(I)=D(I)-SUM*H(I,K)
C
C     Return from the subroutine.
C
  400 RETURN
      END
