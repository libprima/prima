      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
     1  CON,SIM,SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT)
      DIMENSION X(*),CON(*),SIM(N,*),SIMI(N,*),DATMAT(MPP,*),
     1  A(N,*),VSIG(*),VETA(*),SIGBAR(*),DX(*),W(*),IACT(*)
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
      ALPHA=0.25
      BETA=2.1
      GAMMA=0.5
      DELTA=1.1
      RHO=RHOBEG
      PARMU=0.0
      IF (IPRINT .GE. 2) PRINT 10, RHO
   10 FORMAT (/3X,'The initial value of RHO is',1PE13.6,2X,
     1  'and PARMU is set to zero.')
      NFVALS=0
      TEMP=1.0/RHO
      DO 30 I=1,N
      SIM(I,NP)=X(I)
      DO 20 J=1,N
      SIM(I,J)=0.0
   20 SIMI(I,J)=0.0
      SIM(I,I)=RHO
   30 SIMI(I,I)=TEMP
      JDROP=NP
      IBRNCH=0
C
C     Make the next call of the user-supplied subroutine CALCFC. These
C     instructions are also used for calling CALCFC during the iterations of
C     the algorithm.
C
   40 IF (NFVALS .GE. MAXFUN .AND. NFVALS .GT. 0) THEN
          IF (IPRINT .GE. 1) PRINT 50
   50     FORMAT (/3X,'Return from subroutine COBYLA because the ',
     1      'MAXFUN limit has been reached.')
          GOTO 600
      END IF
      NFVALS=NFVALS+1
      CALL CALCFC (N,M,X,F,CON)
      RESMAX=0.0
      IF (M .GT. 0) THEN
          DO 60 K=1,M
   60     RESMAX=AMAX1(RESMAX,-CON(K))
      END IF
      IF (NFVALS .EQ. IPRINT-1 .OR. IPRINT .EQ. 3) THEN
          PRINT 70, NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
   70     FORMAT (/3X,'NFVALS =',I5,3X,'F =',1PE13.6,4X,'MAXCV =',
     1      1PE13.6/3X,'X =',1PE13.6,1P4E15.6)
          IF (IPTEM .LT. N) PRINT 80, (X(I),I=IPTEMP,N)
   80     FORMAT (1PE19.6,1P4E15.6)
      END IF
      CON(MP)=F
      CON(MPP)=RESMAX
      IF (IBRNCH .EQ. 1) GOTO 440
C
C     Set the recently calculated function values in a column of DATMAT. This
C     array has a column for each vertex of the current simplex, the entries of
C     each column being the values of the constraint functions (if any)
C     followed by the objective function and the greatest constraint violation
C     at the vertex.
C
      DO 90 K=1,MPP
   90 DATMAT(K,JDROP)=CON(K)
      IF (NFVALS .GT. NP) GOTO 130
C
C     Exchange the new vertex of the initial simplex with the optimal vertex if
C     necessary. Then, if the initial simplex is not complete, pick its next
C     vertex and calculate the function values there.
C
      IF (JDROP .LE. N) THEN
          IF (DATMAT(MP,NP) .LE. F) THEN
              X(JDROP)=SIM(JDROP,NP)
          ELSE
              SIM(JDROP,NP)=X(JDROP)
              DO 100 K=1,MPP
              DATMAT(K,JDROP)=DATMAT(K,NP)
  100         DATMAT(K,NP)=CON(K)
              DO 120 K=1,JDROP
              SIM(JDROP,K)=-RHO
              TEMP=0.0
              DO 110 I=K,JDROP
  110         TEMP=TEMP-SIMI(I,K)
  120         SIMI(JDROP,K)=TEMP
          END IF
      END IF
      IF (NFVALS .LE. N) THEN
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
      DO 150 J=1,N
      TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
      IF (TEMP .LT. PHIMIN) THEN
          NBEST=J
          PHIMIN=TEMP
      ELSE IF (TEMP .EQ. PHIMIN .AND. PARMU .EQ. 0.0) THEN
          IF (DATMAT(MPP,J) .LT. DATMAT(MPP,NBEST)) NBEST=J
      END IF
  150 CONTINUE
C
C     Switch the best vertex into pole position if it is not there already,
C     and also update SIM, SIMI and DATMAT.
C
      IF (NBEST .LE. N) THEN
          DO 160 I=1,MPP
          TEMP=DATMAT(I,NP)
          DATMAT(I,NP)=DATMAT(I,NBEST)
  160     DATMAT(I,NBEST)=TEMP
          DO 180 I=1,N
          TEMP=SIM(I,NBEST)
          SIM(I,NBEST)=0.0
          SIM(I,NP)=SIM(I,NP)+TEMP
          TEMPA=0.0
          DO 170 K=1,N
          SIM(I,K)=SIM(I,K)-TEMP
  170     TEMPA=TEMPA-SIMI(K,I)
  180     SIMI(NBEST,I)=TEMPA
      END IF
C
C     Make an error return if SIGI is a poor approximation to the inverse of
C     the leading N by N submatrix of SIG.
C
      ERROR=0.0
      DO 200 I=1,N
      DO 200 J=1,N
      TEMP=0.0
      IF (I .EQ. J) TEMP=TEMP-1.0
      DO 190 K=1,N
  190 TEMP=TEMP+SIMI(I,K)*SIM(K,J)
  200 ERROR=AMAX1(ERROR,ABS(TEMP))
      IF (ERROR .GT. 0.1) THEN
          IF (IPRINT .GE. 1) PRINT 210
  210     FORMAT (/3X,'Return from subroutine COBYLA because ',
     1      'rounding errors are becoming damaging.')
          GOTO 600
      END IF
C
C     Calculate the coefficients of the linear approximations to the objective
C     and constraint functions, placing minus the objective function gradient
C     after the constraint gradients in the array A. The vector W is used for
C     working space.
C
      DO 240 K=1,MP
      CON(K)=-DATMAT(K,NP)
      DO 220 J=1,N
  220 W(J)=DATMAT(K,J)+CON(K)
      DO 240 I=1,N
      TEMP=0.0
      DO 230 J=1,N
  230 TEMP=TEMP+W(J)*SIMI(J,I)
      IF (K .EQ. MP) TEMP=-TEMP
  240 A(I,K)=TEMP
C
C     Calculate the values of sigma and eta, and set IFLAG=0 if the current
C     simplex is not acceptable.
C
      IFLAG=1
      PARSIG=ALPHA*RHO
      PARETA=BETA*RHO
      DO 260 J=1,N
      WSIG=0.0
      WETA=0.0
      DO 250 I=1,N
      WSIG=WSIG+SIMI(J,I)**2
  250 WETA=WETA+SIM(I,J)**2
      VSIG(J)=1.0/SQRT(WSIG)
      VETA(J)=SQRT(WETA)
      IF (VSIG(J) .LT. PARSIG .OR. VETA(J) .GT. PARETA) IFLAG=0
  260 CONTINUE
C
C     If a new vertex is needed to improve acceptability, then decide which
C     vertex to drop from the simplex.
C
      IF (IBRNCH .EQ. 1 .OR. IFLAG .EQ. 1) GOTO 370
      JDROP=0
      TEMP=PARETA
      DO 270 J=1,N
      IF (VETA(J) .GT. TEMP) THEN
          JDROP=J
          TEMP=VETA(J)
      END IF
  270 CONTINUE
      IF (JDROP .EQ. 0) THEN
          DO 280 J=1,N
          IF (VSIG(J) .LT. TEMP) THEN
              JDROP=J
              TEMP=VSIG(J)
          END IF
  280     CONTINUE
      END IF
C
C     Calculate the step to the new vertex and its sign.
C
      TEMP=GAMMA*RHO*VSIG(JDROP)
      DO 290 I=1,N
  290 DX(I)=TEMP*SIMI(JDROP,I)
      CVMAXP=0.0
      CVMAXM=0.0
      DO 310 K=1,MP
      SUM=0.0
      DO 300 I=1,N
  300 SUM=SUM+A(I,K)*DX(I)
      IF (K .LT. MP) THEN
          TEMP=DATMAT(K,NP)
          CVMAXP=AMAX1(CVMAXP,-SUM-TEMP)
          CVMAXM=AMAX1(CVMAXM,SUM-TEMP)
      END IF
  310 CONTINUE
      DXSIGN=1.0
      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUM+SUM) DXSIGN=-1.0
C
C     Update the elements of SIM and SIMI, and set the next X.
C
      TEMP=0.0
      DO 320 I=1,N
      DX(I)=DXSIGN*DX(I)
      SIM(I,JDROP)=DX(I)
  320 TEMP=TEMP+SIMI(JDROP,I)*DX(I)
      DO 330 I=1,N
  330 SIMI(JDROP,I)=SIMI(JDROP,I)/TEMP
      DO 360 J=1,N
      IF (J .NE. JDROP) THEN
          TEMP=0.0
          DO 340 I=1,N
  340     TEMP=TEMP+SIMI(J,I)*DX(I)
          DO 350 I=1,N
  350     SIMI(J,I)=SIMI(J,I)-TEMP*SIMI(JDROP,I)
      END IF
  360 X(J)=SIM(J,NP)+DX(J)
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
      IF (IFULL .EQ. 0) THEN
          TEMP=0.0
          DO 380 I=1,N
  380     TEMP=TEMP+DX(I)**2
          IF (TEMP .LT. 0.25*RHO*RHO) THEN
              IBRNCH=1
              GOTO 550
          END IF
      END IF
C
C     Predict the change to F and the new maximum constraint violation if the
C     variables are altered from x(0) to x(0)+DX.
C
      RESNEW=0.0
      CON(MP)=0.0
      DO 400 K=1,MP
      SUM=CON(K)
      DO 390 I=1,N
  390 SUM=SUM-A(I,K)*DX(I)
      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUM)
  400 CONTINUE
C
C     Increase PARMU if necessary and branch back if this change alters the
C     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
C     reductions in the merit function and the maximum constraint violation
C     respectively.
C
      BARMU=0.0
      PREREC=DATMAT(MPP,NP)-RESNEW
      IF (PREREC .GT. 0.0) BARMU=SUM/PREREC
      IF (PARMU .LT. 1.5*BARMU) THEN
          PARMU=2.0*BARMU
          IF (IPRINT .GE. 2) PRINT 410, PARMU
  410     FORMAT (/3X,'Increase in PARMU to',1PE13.6)
          PHI=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
          DO 420 J=1,N
          TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
          IF (TEMP .LT. PHI) GOTO 140
          IF (TEMP .EQ. PHI .AND. PARMU .EQ. 0.0) THEN
              IF (DATMAT(MPP,J) .LT. DATMAT(MPP,NP)) GOTO 140
          END IF
  420     CONTINUE
      END IF
      PREREM=PARMU*PREREC-SUM
C
C     Calculate the constraint and objective functions at x(*). Then find the
C     actual reduction in the merit function.
C
      DO 430 I=1,N
  430 X(I)=SIM(I,NP)+DX(I)
      IBRNCH=1
      GOTO 40
  440 VMOLD=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
      VMNEW=F+PARMU*RESMAX
      TRURED=VMOLD-VMNEW
      IF (PARMU .EQ. 0.0 .AND. F .EQ. DATMAT(MP,NP)) THEN
          PREREM=PREREC
          TRURED=DATMAT(MPP,NP)-RESMAX
      END IF
C
C     Begin the operations that decide whether x(*) should replace one of the
C     vertices of the current simplex, the change being mandatory if TRURED is
C     positive. Firstly, JDROP is set to the index of the vertex that is to be
C     replaced.
C
      RATIO=0.0
      IF (TRURED .LE. 0.0) RATIO=1.0
      JDROP=0
      DO 460 J=1,N
      TEMP=0.0
      DO 450 I=1,N
  450 TEMP=TEMP+SIMI(J,I)*DX(I)
      TEMP=ABS(TEMP)
      IF (TEMP .GT. RATIO) THEN
          JDROP=J
          RATIO=TEMP
      END IF
  460 SIGBAR(J)=TEMP*VSIG(J)
C
C     Calculate the value of ell.
C
      EDGMAX=DELTA*RHO
      L=0
      DO 480 J=1,N
      IF (SIGBAR(J) .GE. PARSIG .OR. SIGBAR(J) .GE. VSIG(J)) THEN
          TEMP=VETA(J)
          IF (TRURED .GT. 0.0) THEN
              TEMP=0.0
              DO 470 I=1,N
  470         TEMP=TEMP+(DX(I)-SIM(I,J))**2
              TEMP=SQRT(TEMP)
          END IF
          IF (TEMP .GT. EDGMAX) THEN
              L=J
              EDGMAX=TEMP
          END IF
      END IF
  480 CONTINUE
      IF (L .GT. 0) JDROP=L
      IF (JDROP .EQ. 0) GOTO 550
C
C     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
C
      TEMP=0.0
      DO 490 I=1,N
      SIM(I,JDROP)=DX(I)
  490 TEMP=TEMP+SIMI(JDROP,I)*DX(I)
      DO 500 I=1,N
  500 SIMI(JDROP,I)=SIMI(JDROP,I)/TEMP
      DO 530 J=1,N
      IF (J .NE. JDROP) THEN
          TEMP=0.0
          DO 510 I=1,N
  510     TEMP=TEMP+SIMI(J,I)*DX(I)
          DO 520 I=1,N
  520     SIMI(J,I)=SIMI(J,I)-TEMP*SIMI(JDROP,I)
      END IF
  530 CONTINUE
      DO 540 K=1,MPP
  540 DATMAT(K,JDROP)=CON(K)
C
C     Branch back for further iterations with the current RHO.
C
      IF (TRURED .GT. 0.0 .AND. TRURED .GE. 0.1*PREREM) GOTO 140
  550 IF (IFLAG .EQ. 0) THEN
          IBRNCH=0
          GOTO 140
      END IF
C
C     Otherwise reduce RHO if it is not at its least value and reset PARMU.
C
      IF (RHO .GT. RHOEND) THEN
          RHO=0.5*RHO
          IF (RHO .LE. 1.5*RHOEND) RHO=RHOEND
          IF (PARMU .GT. 0.0) THEN
              DENOM=0.0
              DO 570 K=1,MP
              CMIN=DATMAT(K,NP)
              CMAX=CMIN
              DO 560 I=1,N
              CMIN=AMIN1(CMIN,DATMAT(K,I))
  560         CMAX=AMAX1(CMAX,DATMAT(K,I))
              IF (K .LE. M .AND. CMIN .LT. 0.5*CMAX) THEN
                  TEMP=AMAX1(CMAX,0.0)-CMIN
                  IF (DENOM .LE. 0.0) THEN
                      DENOM=TEMP
                  ELSE
                      DENOM=AMIN1(DENOM,TEMP)
                  END IF
              END IF
  570         CONTINUE
              IF (DENOM .EQ. 0.0) THEN
                  PARMU=0.0
              ELSE IF (CMAX-CMIN .LT. PARMU*DENOM) THEN
                  PARMU=(CMAX-CMIN)/DENOM
              END IF
          END IF
          IF (IPRINT .GE. 2) PRINT 580, RHO,PARMU
  580     FORMAT (/3X,'Reduction in RHO to',1PE13.6,'  and PARMU =',
     1      1PE13.6)
          IF (IPRINT .EQ. 2) THEN
              PRINT 70, NFVALS,DATMAT(MP,NP),DATMAT(MPP,NP),
     1          (SIM(I,NP),I=1,IPTEM)
              IF (IPTEM .LT. N) PRINT 80, (X(I),I=IPTEMP,N)
          END IF
          GOTO 140
      END IF
C
C     Return the best calculated values of the variables.
C
      IF (IPRINT .GE. 1) PRINT 590
  590 FORMAT (/3X,'Normal return from subroutine COBYLA')
      IF (IFULL .EQ. 1) GOTO 620
  600 DO 610 I=1,N
  610 X(I)=SIM(I,NP)
      F=DATMAT(MP,NP)
      RESMAX=DATMAT(MPP,NP)
  620 IF (IPRINT .GE. 1) THEN
          PRINT 70, NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
          IF (IPTEM .LT. N) PRINT 80, (X(I),I=IPTEMP,N)
      END IF
      MAXFUN=NFVALS
      RETURN
      END
