      SUBROUTINE TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
     1  D,G,HD,HS,CRVMIN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*),
     1  D(*),G(*),HD(*),HS(*)
C
C     N is the number of variables of a quadratic objective function, Q say.
C     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
C       in order to define the current quadratic model Q.
C     DELTA is the trust region radius, and has to be positive.
C     STEP will be set to the calculated trial step.
C     The arrays D, G, HD and HS will be used for working space.
C     CRVMIN will be set to the least curvature of H along the conjugate
C       directions that occur, except that it is set to zero if STEP goes
C       all the way to the trust region boundary.
C
C     The calculation of STEP begins with the truncated conjugate gradient
C     method. If the boundary of the trust region is reached, then further
C     changes to STEP may be made, each one being in the 2D space spanned
C     by the current STEP and the corresponding gradient of Q. Thus STEP
C     should provide a substantial reduction to Q within the trust region.
C
C     Initialization, which includes setting HD to H times XOPT.
C
      HALF=0.5D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(1.0D0)
      DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      DO I=1,N
          D(I)=XOPT(I)
      END DO
      GOTO 170
C
C     Prepare for the first line search.
C
   20 QRED=ZERO
      DD=ZERO
      DO I=1,N
          STEP(I)=ZERO
          HS(I)=ZERO
          G(I)=GQ(I)+HD(I)
          D(I)=-G(I)
          DD=DD+D(I)**2
      END DO
      CRVMIN=ZERO
      IF (DD == ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
C
C     Calculate the step to the trust region boundary and the product HD.
C
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
      GOTO 170
   50 DHD=ZERO
      DO J=1,N
          DHD=DHD+D(J)*HD(J)
      END DO
C
C     Update CRVMIN and set the step-length ALPHA.
C
      ALPHA=BSTEP
      IF (DHD > ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC == 1) CRVMIN=TEMP
          CRVMIN=DMIN1(CRVMIN,TEMP)
          ALPHA=DMIN1(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
C
C     Update STEP and HS.
C
      GGSAV=GG
      GG=ZERO
      DO I=1,N
          STEP(I)=STEP(I)+ALPHA*D(I)
          HS(I)=HS(I)+ALPHA*HD(I)
          GG=GG+(G(I)+HS(I))**2
      END DO
C
C     Begin another conjugate direction iteration if required.
C
      IF (ALPHA < BSTEP) THEN
          IF (QADD <= 0.01D0*QRED) GOTO 160
          IF (GG <= 1.0D-4*GGBEG) GOTO 160
          IF (ITERC == ITERMAX) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO I=1,N
              D(I)=TEMP*D(I)-G(I)-HS(I)
              DD=DD+D(I)**2
              DS=DS+D(I)*STEP(I)
              SS=SS+STEP(I)**2
          END DO
          IF (DS <= ZERO) GOTO 160
          IF (SS < DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
C
C     Test whether an alternative iteration is required.
C
   90 IF (GG <= 1.0D-4*GGBEG) GOTO 160
      SG=ZERO
      SHS=ZERO
      DO I=1,N
          SG=SG+STEP(I)*G(I)
          SHS=SHS+STEP(I)*HS(I)
      END DO
      SGK=SG+SHS
      ANGTEST=SGK/DSQRT(GG*DELSQ)
      IF (ANGTEST <= -0.99D0) GOTO 160
C
C     Begin the alternative iteration by calculating D and HD and some
C     scalar products.
C
      ITERC=ITERC+1
      TEMP=DSQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO I=1,N
          D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      END DO
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO I=1,N
          DG=DG+D(I)*G(I)
          DHD=DHD+HD(I)*D(I)
          DHS=DHS+HD(I)*STEP(I)
      END DO
C
C     Seek the value of the angle that minimizes Q.
C
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO I=1,IU
          ANGLE=DFLOAT(I)*TEMP
          CTH=DCOS(ANGLE)
          STH=DSIN(ANGLE)
          QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
          IF (QNEW < QMIN) THEN
              QMIN=QNEW
              ISAVE=I
              TEMPA=QSAV
          ELSE IF (I == ISAVE+1) THEN
              TEMPB=QNEW
          END IF
          QSAV=QNEW
      END DO
      IF (ISAVE == ZERO) TEMPA=QNEW
      IF (ISAVE == IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA /= TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+ANGLE)
C
C     Calculate the new STEP and HS. Then test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO I=1,N
          STEP(I)=CTH*STEP(I)+STH*D(I)
          HS(I)=CTH*HS(I)+STH*HD(I)
          GG=GG+(G(I)+HS(I))**2
      END DO
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
      IF (ITERC < ITERMAX .AND. RATIO > 0.01D0) GOTO 90
  160 RETURN
C
C     The following instructions act as a subroutine for setting the vector
C     HD to the vector D multiplied by the second derivative matrix of Q.
C     They are called from three different places, which are distinguished
C     by the value of ITERC.
C
  170 DO I=1,N
          HD(I)=ZERO
      END DO
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*D(J)
          END DO
          TEMP=TEMP*PQ(K)
          DO I=1,N
              HD(I)=HD(I)+TEMP*XPT(K,I)
          END DO
      END DO
      IH=0
      DO J=1,N
          DO I=1,J
              IH=IH+1
              IF (I < J) HD(J)=HD(J)+HQ(IH)*D(I)
              HD(I)=HD(I)+HQ(IH)*D(J)
          END DO
      END DO
      IF (ITERC == 0) GOTO 20
      IF (ITERC <= ITERSW) GOTO 50
      GOTO 120
      END
