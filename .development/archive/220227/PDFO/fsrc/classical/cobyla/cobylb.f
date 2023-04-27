CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
C     1  CON,SIM,SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT)
      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,CON,SIM,
     1  SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT,F,INFO,FTARGET,RESMAX)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION X(*),CON(*),SIM(N,*),SIMI(N,*),DATMAT(MPP,*),
     1  A(N,*),VSIG(*),VETA(*),SIGBAR(*),DX(*),W(*),IACT(*),CONSAV(M)
C
C     Set the initial values of some parameters. The last column of SIM holds
C     the optimal vertex of the current simplex, and the preceding N columns
C     hold the displacements from the optimal vertex to the other vertices.
C     Further, SIMI holds the inverse of the matrix that is contained in the
C     first N columns of SIM.
C
      IPTEM=MIN0(N,5)
      IPTEMP=IPTEM+1
      NP=N+1
      MP=M+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ALPHA=0.25
C      BETA=2.1
C      GAMMA=0.5
C      DELTA=1.1
      ALPHA=0.25D0
      BETA=2.1D0
      GAMMA=0.5D0
      DELTA=1.1D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RHO=RHOBEG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      PARMU=0.0
      PARMU=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (IPRINT >= 2) PRINT 10, RHO
   10 FORMAT (/3X,'The initial value of RHO is',1PE13.6,2X,
     1  'and PARMU is set to zero.')
      NFVALS=0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=1.0/RHO
      TEMP=1.0D0/RHO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          SIM(I,NP)=X(I)
          DO J=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SIM(I,J)=0.0
C   20 SIMI(I,J)=0.0
              SIM(I,J)=0.0D0
              SIMI(I,J)=0.0D0
          END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SIM(I,I)=RHO
          SIMI(I,I)=TEMP
      END DO
      JDROP=NP
      IBRNCH=0
C
C     Make the next call of the user-supplied subroutine CALCFC. These
C     instructions are also used for calling CALCFC during the iterations of
C     the algorithm.
C
   40 IF (NFVALS >= MAXFUN .AND. NFVALS > 0) THEN
          IF (IPRINT >= 1) PRINT 50
   50     FORMAT (/3X,'Return from subroutine COBYLA because the ',
     1      'MAXFUN limit has been reached.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 600
      END IF
      NFVALS=NFVALS+1
      CALL CALCFC (N,M,X,F,CON)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RESMAX=0.0
      RESMAX=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (M > 0) THEN
          DO K=1,M
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   60     RESMAX=AMAX1(RESMAX,-CON(K))
              RESMAX=DMAX1(RESMAX,-CON(K))
          END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
      IF (NFVALS == IPRINT-1 .OR. IPRINT == 3) THEN
          PRINT 70, NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
   70     FORMAT (/3X,'NFVALS =',I5,3X,'F =',1PE13.6,4X,'MAXCV =',
     1      1PE13.6/3X,'X =',1PE13.6,1P4E15.6)
          IF (IPTEM < N) PRINT 80, (X(I),I=IPTEMP,N)
   80     FORMAT (1PE19.6,1P4E15.6)
      END IF
      CON(MP)=F
      CON(MPP)=RESMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom 2020/05/15:
C     The values in CON are revised when the linear models are
C     calculated. We store, therefore, its values into CONSAV to access
C     after the execution the true constraint function evaluations.
      DO K=1,M
          CONSAV(K)=CON(K)
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IBRNCH == 1) GOTO 440
C
C     Set the recently calculated function values in a column of DATMAT. This
C     array has a column for each vertex of the current simplex, the entries of
C     each column being the values of the constraint functions (if any)
C     followed by the objective function and the greatest constraint violation
C     at the vertex.
C
      DO K=1,MPP
          DATMAT(K,JDROP)=CON(K)
      END DO
      IF (NFVALS > NP) GOTO 130
C
C     Exchange the new vertex of the initial simplex with the optimal vertex if
C     necessary. Then, if the initial simplex is not complete, pick its next
C     vertex and calculate the function values there.
C
      IF (JDROP <= N) THEN
          IF (DATMAT(MP,NP) <= F) THEN
              X(JDROP)=SIM(JDROP,NP)
          ELSE
              SIM(JDROP,NP)=X(JDROP)
              DO K=1,MPP
                  DATMAT(K,JDROP)=DATMAT(K,NP)
                  DATMAT(K,NP)=CON(K)
              END DO
              DO K=1,JDROP
                  SIM(JDROP,K)=-RHO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              TEMP=0.0
                  TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  DO I=K,JDROP
                      TEMP=TEMP-SIMI(I,K)
                  END DO
                  SIMI(JDROP,K)=TEMP
              END DO
          END IF
      END IF
      IF (NFVALS <= N) THEN
          JDROP=NFVALS
          X(JDROP)=X(JDROP)+RHO
          GOTO 40
      END IF
  130 IBRNCH=1
C
C     Identify the optimal vertex of the current simplex.
C
  140 PHIMIN=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
      NBEST=NP
      DO J=1,N
          TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
          IF (TEMP < PHIMIN) THEN
              NBEST=J
              PHIMIN=TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ELSE IF (TEMP .EQ. PHIMIN .AND. PARMU .EQ. 0.0) THEN
          ELSE IF (TEMP == PHIMIN .AND. PARMU == 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (DATMAT(MPP,J) < DATMAT(MPP,NBEST)) NBEST=J
          END IF
      END DO
C
C     Switch the best vertex into pole position if it is not there already,
C     and also update SIM, SIMI and DATMAT.
C
      IF (NBEST <= N) THEN
          DO I=1,MPP
              TEMP=DATMAT(I,NP)
              DATMAT(I,NP)=DATMAT(I,NBEST)
              DATMAT(I,NBEST)=TEMP
          END DO
          DO I=1,N
              TEMP=SIM(I,NBEST)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          SIM(I,NBEST)=0.0
              SIM(I,NBEST)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              SIM(I,NP)=SIM(I,NP)+TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMPA=0.0
              TEMPA=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO K=1,N
                  SIM(I,K)=SIM(I,K)-TEMP
                  TEMPA=TEMPA-SIMI(K,I)
              END DO
              SIMI(NBEST,I)=TEMPA
          END DO
      END IF
C
C     Make an error return if SIGI is a poor approximation to the inverse of
C     the leading N by N submatrix of SIG.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ERROR=0.0
      ERROR=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          DO J=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
C      IF (I .EQ. J) TEMP=TEMP-1.0
              TEMP=0.0D0
              IF (I == J) TEMP=TEMP-1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO K=1,N
                  TEMP=TEMP+SIMI(I,K)*SIM(K,J)
              END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  200 ERROR=AMAX1(ERROR,ABS(TEMP))
C      IF (ERROR .GT. 0.1) THEN
              ERROR=DMAX1(ERROR,DABS(TEMP))
          END DO
      END DO
      IF (ERROR > 0.1D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (IPRINT >= 1) PRINT 210
  210     FORMAT (/3X,'Return from subroutine COBYLA because ',
     1      'rounding errors are becoming damaging.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO = 7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 600
      END IF
C
C     Calculate the coefficients of the linear approximations to the objective
C     and constraint functions, placing minus the objective function gradient
C     after the constraint gradients in the array A. The vector W is used for
C     working space.
C
      DO K=1,MP
          CON(K)=-DATMAT(K,NP)
          DO J=1,N
              W(J)=DATMAT(K,J)+CON(K)
          END DO
          DO I=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
              TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO J=1,N
                  TEMP=TEMP+W(J)*SIMI(J,I)
              END DO
              IF (K == MP) TEMP=-TEMP
              A(I,K)=TEMP
          END DO
      END DO
C
C     Calculate the values of sigma and eta, and set IFLAG=0 if the current
C     simplex is not acceptable.
C
      IFLAG=1
      PARSIG=ALPHA*RHO
      PARETA=BETA*RHO
      DO J=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      WSIG=0.0
C      WETA=0.0
          WSIG=0.0D0
          WETA=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              WSIG=WSIG+SIMI(J,I)**2
              WETA=WETA+SIM(I,J)**2
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      VSIG(J)=1.0/SQRT(WSIG)
C      VETA(J)=SQRT(WETA)
          VSIG(J)=1.0D0/DSQRT(WSIG)
          VETA(J)=DSQRT(WETA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (VSIG(J) < PARSIG .OR. VETA(J) > PARETA) IFLAG=0
      END DO
C
C     If a new vertex is needed to improve acceptability, then decide which
C     vertex to drop from the simplex.
C
      IF (IBRNCH == 1 .OR. IFLAG == 1) GOTO 370
      JDROP=0
      TEMP=PARETA
      DO J=1,N
          IF (VETA(J) > TEMP) THEN
              JDROP=J
              TEMP=VETA(J)
          END IF
      END DO
      IF (JDROP == 0) THEN
          DO J=1,N
              IF (VSIG(J) < TEMP) THEN
                  JDROP=J
                  TEMP=VSIG(J)
              END IF
          END DO
      END IF
C
C     Calculate the step to the new vertex and its sign.
C
      TEMP=GAMMA*RHO*VSIG(JDROP)
      DO I=1,N
          DX(I)=TEMP*SIMI(JDROP,I)
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      CVMAXP=0.0
C      CVMAXM=0.0
      CVMAXP=0.0D0
      CVMAXM=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO K=1,MP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SUM=0.0
          SUM=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SUM=SUM+A(I,K)*DX(I)
          END DO
          IF (K < MP) THEN
              TEMP=DATMAT(K,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CVMAXP=AMAX1(CVMAXP,-SUM-TEMP)
C          CVMAXM=AMAX1(CVMAXM,SUM-TEMP)
              CVMAXP=DMAX1(CVMAXP,-SUM-TEMP)
              CVMAXM=DMAX1(CVMAXM,SUM-TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          END IF
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      DXSIGN=1.0
C      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUM+SUM) DXSIGN=-1.0
      DXSIGN=1.0D0
      IF (PARMU*(CVMAXP-CVMAXM) > SUM+SUM) DXSIGN=-1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
C
C     Update the elements of SIM and SIMI, and set the next X.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          DX(I)=DXSIGN*DX(I)
          SIM(I,JDROP)=DX(I)
          TEMP=TEMP+SIMI(JDROP,I)*DX(I)
      END DO
      DO I=1,N
          SIMI(JDROP,I)=SIMI(JDROP,I)/TEMP
      END DO
      DO J=1,N
          IF (J /= JDROP) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=0.0
              TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO I=1,N
                  TEMP=TEMP+SIMI(J,I)*DX(I)
              END DO
              DO I=1,N
                  SIMI(J,I)=SIMI(J,I)-TEMP*SIMI(JDROP,I)
              END DO
          END IF
          X(J)=SIM(J,NP)+DX(J)
      END DO
      GOTO 40
C
C     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO.
C
  370 IZ=1
      IZDOTA=IZ+N*N
      IVMC=IZDOTA+N
      ISDIRN=IVMC+MP
      IDXNEW=ISDIRN+N
      IVMD=IDXNEW+N
      CALL TRSTLP (N,M,A,CON,RHO,DX,IFULL,IACT,W(IZ),W(IZDOTA),
     1  W(IVMC),W(ISDIRN),W(IDXNEW),W(IVMD))
      IF (IFULL == 0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=0.0
          TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              TEMP=TEMP+DX(I)**2
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TEMP .LT. 0.25*RHO*RHO) THEN
          IF (TEMP < 0.25D0*RHO*RHO) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IBRNCH=1
              GOTO 550
          END IF
      END IF
C
C     Predict the change to F and the new maximum constraint violation if the
C     variables are altered from x(0) to x(0)+DX.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RESNEW=0.0
C      CON(MP)=0.0
      RESNEW=0.0D0
      CON(MP)=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      DO K=1,MP
          SUM=CON(K)
          DO I=1,N
              SUM=SUM-A(I,K)*DX(I)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUM)
          IF (K < MP) RESNEW=DMAX1(RESNEW,SUM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      END DO
C
C     Increase PARMU if necessary and branch back if this change alters the
C     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
C     reductions in the merit function and the maximum constraint violation
C     respectively.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      BARMU=0.0
      BARMU=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PREREC=DATMAT(MPP,NP)-RESNEW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (PREREC .GT. 0.0) BARMU=SUM/PREREC
C      IF (PARMU .LT. 1.5*BARMU) THEN
C          PARMU=2.0*BARMU
      IF (PREREC > 0.0D0) BARMU=SUM/PREREC
      IF (PARMU < 1.5D0*BARMU) THEN
          PARMU=2.0D0*BARMU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (IPRINT >= 2) PRINT 410, PARMU
  410     FORMAT (/3X,'Increase in PARMU to',1PE13.6)
          PHI=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
          DO J=1,N
              TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
              IF (TEMP < PHI) GOTO 140
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TEMP .EQ. PHI .AND. PARMU .EQ. 0.0) THEN
              IF (TEMP == PHI .AND. PARMU == 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IF (DATMAT(MPP,J) < DATMAT(MPP,NP)) GOTO 140
              END IF
          END DO
      END IF
      PREREM=PARMU*PREREC-SUM
C
C     Calculate the constraint and objective functions at x(*). Then find the
C     actual reduction in the merit function.
C
      DO I=1,N
          X(I)=SIM(I,NP)+DX(I)
      END DO
      IBRNCH=1
      GOTO 40
  440 VMOLD=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
      VMNEW=F+PARMU*RESMAX
      TRURED=VMOLD-VMNEW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (PARMU .EQ. 0.0 .AND. F .EQ. DATMAT(MP,NP)) THEN
      IF (PARMU == 0.0D0 .AND. F == DATMAT(MP,NP)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          PREREM=PREREC
          TRURED=DATMAT(MPP,NP)-RESMAX
      END IF
C
C     Begin the operations that decide whether x(*) should replace one of the
C     vertices of the current simplex, the change being mandatory if TRURED is
C     positive. Firstly, JDROP is set to the index of the vertex that is to be
C     replaced.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RATIO=0.0
C      IF (TRURED .LE. 0.0) RATIO=1.0
      RATIO=0.0D0
      IF (TRURED <= 0.0D0) RATIO=1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      JDROP=0
      DO J=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
          TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              TEMP=TEMP+SIMI(J,I)*DX(I)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=ABS(TEMP)
          TEMP=DABS(TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (TEMP > RATIO) THEN
              JDROP=J
              RATIO=TEMP
          END IF
          SIGBAR(J)=TEMP*VSIG(J)
      END DO
C
C     Calculate the value of ell.
C
      EDGMAX=DELTA*RHO
      L=0
      DO J=1,N
          IF (SIGBAR(J) >= PARSIG .OR. SIGBAR(J) >= VSIG(J)) THEN
              TEMP=VETA(J)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TRURED .GT. 0.0) THEN
C              TEMP=0.0
              IF (TRURED > 0.0D0) THEN
                  TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  DO I=1,N
                      TEMP=TEMP+(DX(I)-SIM(I,J))**2
                  END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              TEMP=SQRT(TEMP)
                  TEMP=DSQRT(TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              END IF
              IF (TEMP > EDGMAX) THEN
                  L=J
                  EDGMAX=TEMP
              END IF
          END IF
      END DO
      IF (L > 0) JDROP=L
      IF (JDROP == 0) GOTO 550
C
C     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          SIM(I,JDROP)=DX(I)
          TEMP=TEMP+SIMI(JDROP,I)*DX(I)
      END DO
      DO I=1,N
          SIMI(JDROP,I)=SIMI(JDROP,I)/TEMP
      END DO
      DO J=1,N
          IF (J /= JDROP) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=0.0
              TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO I=1,N
                  TEMP=TEMP+SIMI(J,I)*DX(I)
              END DO
              DO I=1,N
                  SIMI(J,I)=SIMI(J,I)-TEMP*SIMI(JDROP,I)
              END DO
          END IF
      END DO
      DO K=1,MPP
          DATMAT(K,JDROP)=CON(K)
      END DO
C
C     Branch back for further iterations with the current RHO.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (TRURED .GT. 0.0 .AND. TRURED .GE. 0.1*PREREM) GOTO 140
      IF (TRURED > 0.0D0 .AND. TRURED >= 0.1D0*PREREM) GOTO 140
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  550 IF (IFLAG == 0) THEN
          IBRNCH=0
          GOTO 140
      END IF
C
C     Otherwise reduce RHO if it is not at its least value and reset PARMU.
C
      IF (RHO > RHOEND) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          RHO=0.5*RHO
C          IF (RHO .LE. 1.5*RHOEND) RHO=RHOEND
C          IF (PARMU .GT. 0.0) THEN
C              DENOM=0.0
          RHO=0.5D0*RHO
          IF (RHO <= 1.5D0*RHOEND) RHO=RHOEND
          IF (PARMU > 0.0D0) THEN
              DENOM=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO K=1,MP
                  CMIN=DATMAT(K,NP)
                  CMAX=CMIN
                  DO I=1,N
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              CMIN=AMIN1(CMIN,DATMAT(K,I))
C  560         CMAX=AMAX1(CMAX,DATMAT(K,I))
C              IF (K .LE. M .AND. CMIN .LT. 0.5*CMAX) THEN
C                  TEMP=AMAX1(CMAX,0.0)-CMIN
C                  IF (DENOM .LE. 0.0) THEN
                      CMIN=DMIN1(CMIN,DATMAT(K,I))
                      CMAX=DMAX1(CMAX,DATMAT(K,I))
                  END DO
                  IF (K <= M .AND. CMIN < 0.5D0*CMAX) THEN
                      TEMP=DMAX1(CMAX,0.0D0)-CMIN
                      IF (DENOM <= 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          DENOM=TEMP
                      ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                      DENOM=AMIN1(DENOM,TEMP)
                          DENOM=DMIN1(DENOM,TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      END IF
                  END IF
              END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              IF (DENOM .EQ. 0.0) THEN
C                  PARMU=0.0
              IF (DENOM == 0.0D0) THEN
                  PARMU=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ELSE IF (CMAX-CMIN < PARMU*DENOM) THEN
                  PARMU=(CMAX-CMIN)/DENOM
              END IF
          END IF
          IF (IPRINT >= 2) PRINT 580, RHO,PARMU
  580     FORMAT (/3X,'Reduction in RHO to',1PE13.6,'  and PARMU =',
     1      1PE13.6)
          IF (IPRINT == 2) THEN
              PRINT 70, NFVALS,DATMAT(MP,NP),DATMAT(MPP,NP),
     1          (SIM(I,NP),I=1,IPTEM)
              IF (IPTEM < N) PRINT 80, (X(I),I=IPTEMP,N)
          END IF
          GOTO 140
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ELSE
          INFO = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
C
C     Return the best calculated values of the variables.
C
      IF (IPRINT >= 1) PRINT 590
  590 FORMAT (/3X,'Normal return from subroutine COBYLA')
      IF (IFULL == 1) GOTO 620
  600 DO I=1,N
          X(I)=SIM(I,NP)
      END DO
      F=DATMAT(MP,NP)
      RESMAX=DATMAT(MPP,NP)
  620 IF (IPRINT >= 1) THEN
          PRINT 70, NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
          IF (IPTEM < N) PRINT 80, (X(I),I=IPTEMP,N)
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom 2020/05/15:
C     The true values of the constraint function evaluations, stored in
C     CONSAV are dumped into CON to be returned by COBYLA.
      DO K=1,M
          CON(K)=CONSAV(K)
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MAXFUN=NFVALS
      RETURN
      END
