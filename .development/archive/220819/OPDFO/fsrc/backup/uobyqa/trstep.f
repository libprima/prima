      SUBROUTINE TRSTEP (N,G,H,DELTA,TOL,D,GG,TD,TN,W,PIV,Z,EVALUE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      DIMENSION G(*),H(N,*),D(*),GG(*),TD(*),TN(*),W(*),PIV(*),Z(*) 
      DIMENSION G(*),H(N,*),D(*),GG(*),TD(*),TN(*),W(*),PIV(*),Z(*), 
     +  DSAV(N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      DO I=1,N
          D(I)=ZERO
          TD(I)=H(I,I)
          DO J=1,I
              H(I,J)=H(J,I)
          END DO
      END DO
C
C     Apply Householder transformations to obtain a tridiagonal matrix that
C     is similar to H, and put the elements of the Householder vectors in
C     the lower triangular part of H. Further, TD and TN will contain the
C     diagonal and other nonzero elements of the tridiagonal matrix.
C
      DO K=1,NM
          KP=K+1
          SUM=ZERO
          IF (KP < N) THEN
              KPP=KP+1
              DO I=KPP,N
                  SUM=SUM+H(I,K)**2
              END DO
          END IF
          IF (SUM == ZERO) THEN
              TN(K)=H(KP,K)
              H(KP,K)=ZERO
          ELSE
              TEMP=H(KP,K)
              TN(K)=DSIGN(DSQRT(SUM+TEMP*TEMP),TEMP)
              H(KP,K)=-SUM/(TEMP+TN(K))
              TEMP=DSQRT(TWO/(SUM+H(KP,K)**2))
              DO I=KP,N
                  W(I)=TEMP*H(I,K)
                  H(I,K)=W(I)
                  Z(I)=TD(I)*W(I)
              END DO
              WZ=ZERO
              DO J=KP,NM
                  JP=J+1
                  DO I=JP,N
                      Z(I)=Z(I)+H(I,J)*W(J)
                      Z(J)=Z(J)+H(I,J)*W(I)
                  END DO
                  WZ=WZ+W(J)*Z(J)
              END DO
              WZ=WZ+W(N)*Z(N)
              DO J=KP,N
                  TD(J)=TD(J)+W(J)*(WZ*W(J)-TWO*Z(J))
                  IF (J < N) THEN
                      JP=J+1
                      DO I=JP,N
                          H(I,J)=H(I,J)-W(I)*Z(J)-W(J)*(Z(I)-WZ*W(I))
                      END DO
                  END IF
              END DO
          END IF
      END DO
C
C     Form GG by applying the similarity transformation to G.
C
      GSQ=ZERO
      DO I=1,N
          GG(I)=G(I)
          GSQ=GSQ+G(I)**2
      END DO
      GNORM=DSQRT(GSQ)
      DO K=1,NM
          KP=K+1
          SUM=ZERO
          DO I=KP,N
              SUM=SUM+GG(I)*H(I,K)
          END DO
          DO I=KP,N
              GG(I)=GG(I)-SUM*H(I,K)
          END DO
      END DO
C
C     Begin the trust region calculation with a tridiagonal matrix by
C     calculating the norm of H. Then treat the case when H is zero.
C
      HNORM=DABS(TD(1))+DABS(TN(1))
      TDMIN=TD(1)
      TN(N)=ZERO
      DO I=2,N
          TEMP=DABS(TN(I-1))+DABS(TD(I))+DABS(TN(I))
          HNORM=DMAX1(HNORM,TEMP)
          TDMIN=DMIN1(TDMIN,TD(I))
      END DO
      IF (HNORM == ZERO) THEN
          IF (GNORM == ZERO) GOTO 400
          SCALE=DELTA/GNORM
          DO I=1,N
              D(I)=-SCALE*GG(I)
          END DO
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 26-06-2019: See the lines below line number 140
      DO I = 1, N
         DSAV(I) = D(I)
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
C
  140 ITERC=ITERC+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 26-06-2019
C The original code can encounter infinite cycling, which did happen
C when testing the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS. 
C Indeed, in all these cases, Inf and NaN appear in D due to extremely
C large values in H (up to 10^219). 
C To avoid wasting energy, we do the following 
      SUMD = ZERO
      DO I = 1, N
          SUMD = SUMD + DABS(D(I))
      END DO
      IF (SUMD >= 1.0D100 .OR. SUMD /= SUMD) THEN
          DO I = 1, N
             D(I) = DSAV(I)
          END DO
          GOTO 370
      ELSE
          DO I = 1, N
             DSAV(I) = D(I)
          END DO
      END IF
      IF (ITERC > MIN(10000, 100*N)) THEN
          GOTO 370
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      KSAV=0
      PIV(1)=TD(1)+PAR
      K=1
  150 IF (PIV(K) > ZERO) THEN
          PIV(K+1)=TD(K+1)+PAR-TN(K)**2/PIV(K)
      ELSE
          IF (PIV(K) < ZERO .OR. TN(K) /= ZERO) GOTO 160
          KSAV=K
          PIV(K+1)=TD(K+1)+PAR
      END IF
      K=K+1
      IF (K < N) GOTO 150
      IF (PIV(K) < ZERO) GOTO 160
      IF (PIV(K) == ZERO) KSAV=K
C
C     Branch if all the pivots are positive, allowing for the case when
C     G is zero.
C
      IF (KSAV == 0 .AND. GSQ > ZERO) GOTO 230
      IF (GSQ == ZERO) THEN
          IF (PAR == ZERO) GOTO 370
          PARU=PAR
          PARUEST=PAR
          IF (KSAV == 0) GOTO 190
      END IF
      K=KSAV
C
C     Set D to a direction of nonpositive curvature of the given tridiagonal
C     matrix, and thus revise PARLEST.
C
  160 D(K)=ONE
      IF (DABS(TN(K)) <= DABS(PIV(K))) THEN
          DSQ=ONE
          DHD=PIV(K)
      ELSE
          TEMP=TD(K+1)+PAR
          IF (TEMP <= DABS(PIV(K))) THEN
              D(K+1)=DSIGN(ONE,-TN(K))
              DHD=PIV(K)+TEMP-TWO*DABS(TN(K))
          ELSE
              D(K+1)=-TN(K)/TEMP
              DHD=PIV(K)+TN(K)*D(K+1)
          END IF
          DSQ=ONE+D(K+1)**2
      END IF
  170 IF (K > 1) THEN
          K=K-1
          IF (TN(K) /= ZERO) THEN
              D(K)=-TN(K)*D(K+1)/PIV(K)
              DSQ=DSQ+D(K)**2
              GOTO 170
          END IF
          DO I=1,K
              D(I)=ZERO
          END DO
      END IF
      PARL=PAR
      PARLEST=PAR-DHD/DSQ
C
C     Terminate with D set to a multiple of the current D if the following
C     test suggests that it suitable to do so.
C
  190 TEMP=PARUEST
      IF (GSQ == ZERO) TEMP=TEMP*(ONE-TOL)
      IF (PARUEST > ZERO .AND. PARLEST >= TEMP) THEN
          DTG=ZERO
          DO I=1,N
              DTG=DTG+D(I)*GG(I)
          END DO
          SCALE=-DSIGN(DELTA/DSQRT(DSQ),DTG)
          DO I=1,N
              D(I)=SCALE*D(I)
          END DO
          GOTO 370
      END IF
C
C     Pick the value of PAR for the next iteration.
C
  220 IF (PARU == ZERO) THEN
          PAR=TWO*PARLEST+GNORM/DELTA
      ELSE
          PAR=0.5D0*(PARL+PARU)
          PAR=DMAX1(PAR,PARLEST)
      END IF
      IF (PARUEST > ZERO) PAR=DMIN1(PAR,PARUEST)
      GOTO 140
C
C     Calculate D for the current PAR in the positive definite case.
C
  230 W(1)=-GG(1)/PIV(1)
      DO I=2,N
          W(I)=(-GG(I)-TN(I-1)*W(I-1))/PIV(I)
      END DO
      D(N)=W(N)
      DO I=NM,1,-1
          D(I)=W(I)-TN(I)*D(I+1)/PIV(I)
      END DO
C
C     Branch if a Newton-Raphson step is acceptable.
C
      DSQ=ZERO
      WSQ=ZERO
      DO I=1,N
          DSQ=DSQ+D(I)**2
          WSQ=WSQ+PIV(I)*W(I)**2
      END DO
      IF (PAR == ZERO .AND. DSQ <= DELSQ) GOTO 320
C
C     Make the usual test for acceptability of a full trust region step.
C
      DNORM=DSQRT(DSQ)
      PHI=ONE/DNORM-ONE/DELTA
      TEMP=TOL*(ONE+PAR*DSQ/WSQ)-DSQ*PHI*PHI
      IF (TEMP >= ZERO) THEN
          SCALE=DELTA/DNORM
          DO I=1,N
              D(I)=SCALE*D(I)
          END DO
          GOTO 370
      END IF
      IF (ITERC >= 2 .AND. PAR <= PARL) GOTO 370
      IF (PARU > ZERO .AND. PAR >= PARU) GOTO 370
C
C     Complete the iteration when PHI is negative.
C
      IF (PHI < ZERO) THEN
          PARLEST=PAR
          IF (POSDEF. EQ. ONE) THEN
              IF (PHI <= PHIL) GOTO 370
              SLOPE=(PHI-PHIL)/(PAR-PARL)
              PARLEST=PAR-PHI/SLOPE
          END IF
          SLOPE=ONE/GNORM
          IF (PARU > ZERO) SLOPE=(PHIU-PHI)/(PARU-PAR)
          TEMP=PAR-PHI/SLOPE
          IF (PARUEST > ZERO) TEMP=DMIN1(TEMP,PARUEST)
          PARUEST=TEMP
          POSDEF=ONE
          PARL=PAR
          PHIL=PHI
          GOTO 220
      END IF
C
C     If required, calculate Z for the alternative test for convergence.
C
      IF (POSDEF == ZERO) THEN
          W(1)=ONE/PIV(1)
          DO I=2,N
              TEMP=-TN(I-1)*W(I-1)
              W(I)=(DSIGN(ONE,TEMP)+TEMP)/PIV(I)
          END DO
          Z(N)=W(N)
          DO I=NM,1,-1
              Z(I)=W(I)-TN(I)*Z(I+1)/PIV(I)
          END DO
          WWSQ=ZERO
          ZSQ=ZERO
          DTZ=ZERO
          DO I=1,N
              WWSQ=WWSQ+PIV(I)*W(I)**2
              ZSQ=ZSQ+Z(I)**2
              DTZ=DTZ+D(I)*Z(I)
          END DO
C
C     Apply the alternative test for convergence.
C
          TEMPA=DABS(DELSQ-DSQ)
          TEMPB=DSQRT(DTZ*DTZ+TEMPA*ZSQ)
          GAM=TEMPA/(DSIGN(TEMPB,DTZ)+DTZ)
          TEMP=TOL*(WSQ+PAR*DELSQ)-GAM*GAM*WWSQ
          IF (TEMP >= ZERO) THEN
              DO I=1,N
                  D(I)=D(I)+GAM*Z(I)
              END DO
              GOTO 370
          END IF
          PARLEST=DMAX1(PARLEST,PAR-WWSQ/ZSQ)
      END IF
C
C     Complete the iteration when PHI is positive.
C
      SLOPE=ONE/GNORM
      IF (PARU > ZERO) THEN
          IF (PHI >= PHIU) GOTO 370
          SLOPE=(PHIU-PHI)/(PARU-PAR)
      END IF
      PARLEST=DMAX1(PARLEST,PAR-PHI/SLOPE)
      PARUEST=PAR
      IF (POSDEF == ONE) THEN
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
      DO K=2,N
          PIVOT=TD(K)-TN(K-1)**2/PIVOT
          SHFMAX=DMIN1(SHFMAX,PIVOT)
      END DO
C
C     Find EVALUE by a bisection method, but occasionally SHFMAX may be
C     adjusted by the rule of false position.
C
      KSAVE=0
  340 SHIFT=0.5D0*(SHFMIN+SHFMAX)
      K=1
      TEMP=TD(1)-SHIFT
  350 IF (TEMP > ZERO) THEN
          PIV(K)=TEMP
          IF (K < N) THEN
              TEMP=TD(K+1)-SHIFT-TN(K)**2/TEMP
              K=K+1
              GOTO 350
          END IF
          SHFMIN=SHIFT
      ELSE
          IF (K < KSAVE) GOTO 360
          IF (K == KSAVE) THEN
              IF (PIVKSV == ZERO) GOTO 360
              IF (PIV(K)-TEMP < TEMP-PIVKSV) THEN
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
      IF (SHFMIN <= 0.99D0*SHFMAX) GOTO 340
  360 EVALUE=SHFMIN
C
C     Apply the inverse Householder transformations to D.
C
  370 NM=N-1
      DO K=NM,1,-1
          KP=K+1
          SUM=ZERO
          DO I=KP,N
              SUM=SUM+D(I)*H(I,K)
          END DO
          DO I=KP,N
              D(I)=D(I)-SUM*H(I,K)
          END DO
      END DO
C
C     Return from the subroutine.
C
  400 RETURN
      END
