!     IFLAG = 0: The current simplex is not acceptable
!     IFLAG = 1: The current simplex is acceptable
!     IBRNCH = 0: The next iteration is geometry step
!     IBRNCH = 1: The next iteration is trust-region step
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
C     1  CON,SIM,SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT)
      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,CON,SIM,
     1  SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT,F,INFO,FTARGET,RESMAX,
     1 NW)



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!-----------------------!!!!!!
        USE DIRTY_TEMPORARY_MOD4POWELL_MOD!
        !!!!!!-----------------------!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL BETTER
      PARAMETER (NSMAX = 2000)
C NSMAX is the maximal number of "dropped X" to save (see comments below
C line number 480)
      PARAMETER (CTOL = EPSILON(1.0D0))
C CTOL is the tolerance for consraint violation. A point X is considered
C to be feasible if its constraint violation (RESMAX) is less than CTOL.
C EPSILON(1.0D0) returns the machine epsilon corresponding to 1.0D0,
C which is expected to be about 2.0D-16.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DIMENSION X(*),CON(*),SIM(N,*),SIMI(N,*),DATMAT(MPP,*),
!     1  A(N,*),VSIG(*),VETA(*),SIGBAR(*),DX(*),W(*),IACT(*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION X(N),CON(MPP),SIM(N,N+1),SIMI(N,N),DATMAT(MPP,N+1),
     1 A(N,M+1),VSIG(N),VETA(N),SIGBAR(N),DX(N),W(NW),IACT(M+1),
     1 CONSAV(MPP),XSAV(N,NSMAX),DATSAV(MPP,NSMAX),XDROP(N),DATDROP(MPP)
     1 ,fhist(n+2+NSMAX), chist(n+2+NSMAX)
     1 , DATMAT_OLD(M+2, N+1), SIM_OLD(N, N+1), SIMI_OLD(N, N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Set the initial values of some parameters. The last column of SIM holds
C     the optimal vertex of the current simplex, and the preceding N columns
C     hold the displacements from the optimal vertex to the other vertices.
C     Further, SIMI holds the inverse of the matrix that is contained in the
C     first N columns of SIM.
C
      INFO = 2147483647
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

      DATMAT(1:M+2,1:N+1) = HUGENUM
      DATSAV(1:M+2,1:NSMAX) = HUGENUM
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NSAV = 0
      DATSAV = HUGENUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Make the next call of the user-supplied subroutine CALCFC. These
C     instructions are also used for calling CALCFC during the iterations of
C     the algorithm.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   40 IF (NFVALS .GE. MAXFUN .AND. NFVALS .GT. 0) THEN
C          IF (IPRINT .GE. 1) PRINT 50
C   50     FORMAT (/3X,'Return from subroutine COBYLA because the ',
C     1      'MAXFUN limit has been reached.')
C          GOTO 600
C      END IF
C      NFVALS=NFVALS+1
C      CALL CALCFC (N,M,X,F,CON)
C
C     By Zaikun (02-06-2019):
   40 DO I=1,N
          IF (X(I) /= X(I)) THEN
              F=X(I) ! Set F to NaN
              INFO = -1
              GOTO 600
          END IF
      END DO

      CALL CALCFC (N,M,X,F,CON)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (F /= F) F = HUGEFUN
      where (CON(1:M) /= CON(1:M)) CON(1:M) = HUGECON
      F = MIN(F, HUGEFUN)
      CON(1:M) = MAX(CON(1:M), -HUGECON)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NFVALS=NFVALS+1


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
C By Zaikun 20190819:
C CONSAV always containts the constraint value of the current x.
C CON, however, will be changed during the calculation (see the lines
C above line number 220).
      DO K = 1, MPP
          CONSAV(K) = CON(K)
      END DO
      CALL SAVEX (X(1:N), CON(1:MPP), XSAV(1:N, 1:NSMAX),
     1     DATSAV(1:MPP, 1:NSMAX), N, M, NSAV, NSMAX, CTOL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom/Zaikun (on 04-06-2019/10-06-2019):
C     CSUMM containts the sum of the absolute value of the constraints to
C     check whether it contains a NaN value.
      CSUMM=0.0D0
      DO K=1,M
          CSUMM=CSUMM+DABS(CON(K))
      END DO
      IF (CSUMM /= CSUMM) THEN
          RESMAX = CSUMM ! Set RESMAX to NaN
          INFO = -2
          GOTO 600
      END IF
C     If the objective function value or the constraints contain a NaN or an
C     infinite value, the algorithm stops.
      IF (F /= F .OR. F > HUGE(0.0D0)) THEN
          INFO = -2
          GOTO 600
      END IF
C     If the objective function achieves the target value at a feasible
C     point, then exit.
C      IF (F .LE. FTARGET .AND. RESMAX .LE. 0.0D0) THEN
      IF (F <= FTARGET .AND. RESMAX <= CTOL) THEN
C         The feasibility is guarantee because RESMAX .LE. CTOL
          INFO = 1
          GOTO 600
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Zaikun (on 06-06-2019)
C     The following code was placed before "CALL CALCFC (N,M,X,F,CON)".
C     This led to a bug, because F may not equal F(X) if the subroutine
C     exits due to NFVALS .GE. MAXFUN (X is updated but F is not evaluated
C     at X). Similar thing can be said about RESMAX.
      IF (NFVALS >= MAXFUN .AND. NFVALS > 0) THEN
          IF (IPRINT >= 1) PRINT 50
   50     FORMAT (/3X,'Return from subroutine COBYLA because the ',
     1      'MAXFUN limit has been reached.')
          INFO = 3
      END IF
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
      ! IF we do not go to 130 but continue to below, then NFVALS <= NP.
      ! Thus NFVALS may be NP = N+1 > N.
C
C     Exchange the new vertex of the initial simplex with the optimal vertex if
C     necessary. Then, if the initial simplex is not complete, pick its next
C     vertex and calculate the function values there.
C
      IF (JDROP <= N) THEN
          IF (DATMAT(MP,NP) <= F) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C When NFVALS<=N, it is necessary to update X(JDROP) because next X will
C be calculated based on the current one (see the code below line number
C 120). The purpose of this update is to make this X hold the variable that
C has the smallest function value up to now. The next X is defined as
C a perturbation of this optimal X, which is reasonable.
C However, this update leads to inconsistency between X and [F, RESMAX],
C meaning that X is not necessarily the X corresponding to [F, RESMAX].
C This can cause COBYLA return inconsistent [X, F, RESMAX].
C Fortunately, this is not a problem if NFVALS <= N. Because, if COBYLA
C returns with a NFVALS <= N, then X contains NaN or F = NaN or nearly
C Inf or the constraint contains NaN, all of which would lead to an
C immediate jump to line 600 without coming here. Therefore, if we
C restrict the update to only the case with NFVALS <= N, ther will be no
C inconsistency at return.
C With the original code, if COBYLA returns with NFVALS = NP (this can
C happen if the initial TR problem constantly produces too short steps
C so that RHO is reduced to RHOEND without producing any acceptable trial
C step; see the code below line number 380), then, when the code arrives
C at line number 600, X and [F, RESMAX] may be inconsistent. However,
C recall that the inconsistency occurs only if SIM(:, NP) has a lower
C function value than X (before updated). Therefore, as long as the code
C takes SIM(:, NP) instead of X, no harm would be done. It is the case
C likely the in the original code, because if COBYLA returns with
C NFVALS=NP, then PARMU=0, and hence SIM(:, NP) will be returned.
C
C              X(JDROP)=SIM(JDROP,NP)
               IF (NFVALS <= N)  X(JDROP)=SIM(JDROP,NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

! 120
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DATMAT_OLD = DATMAT(1:M+2,1:N+1)
      SIM_OLD  = SIM(1:N, 1:N+1)
      SIMI_OLD = SIMI(1:N, 1:N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2021-05-30
      IF (INFO == 3) GOTO 600
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      IF (.NOT. (ERROR <= 0.1D0)) THEN
          DATMAT(1:M+2,1:N+1) = DATMAT_OLD
          SIM(1:N, 1:N+1) = SIM_OLD
          SIMI(1:N, 1:N) = SIMI_OLD
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
! 220
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun 2021-07-19: With the following line, the geometry step
      ! will not be taken if the current simplex is acceptable but then
      ! becomes unacceptable due to the update of PARMU in the lines
      ! above the line number 410.
      IF (IFLAG == 1) IBRNCH = 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 20190822: If VETA or VSIG become NaN due to rounding errors,
C JDROP may end up being 0. If we continue, then a Segmentation Fault
C will happen because we will read SIM(:, JDROP) and VSIG(JDROP).
!      summ = sum(veta(1:n))
!      if (summ /= summ) jdrop = 0
!      summ = sum(vsig(1:n))
!      if (summ /= summ) jdrop = 0
      IF (JDROP == 0) THEN
          IF (IPRINT >= 1) PRINT 286
  286     FORMAT (/3X,'Return from subroutine COBYLA because ',
     1      'rounding errors are becoming damaging.')
          INFO = 7
          GOTO 600
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 20190820: See the comments below line number 480
      DO I = 1, N
          XDROP(I) = SIM(I, NP) + SIM(I, JDROP) ! JDROP<NP is guaranteed
      END DO
      DO K = 1, MPP
          DATDROP(K) = DATMAT(K, JDROP)
      END DO
!      CALL SAVEX (XDROP(1:N), DATDROP(1:MPP), XSAV(1:N, 1:NSMAX),
!     1     DATSAV(1:MPP, 1:NSMAX), N, M, NSAV, NSMAX, CTOL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
C      SUMM=0.0
          SUMM=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              SUMM=SUMM+A(I,K)*DX(I)
          END DO
          IF (K < MP) THEN
              TEMP=DATMAT(K,NP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          CVMAXP=AMAX1(CVMAXP,-SUMM-TEMP)
C          CVMAXM=AMAX1(CVMAXM,SUMM-TEMP)
              CVMAXP=DMAX1(CVMAXP,-SUMM-TEMP)
              CVMAXM=DMAX1(CVMAXM,SUMM-TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          END IF
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      DXSIGN=1.0
C      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUMM+SUMM) DXSIGN=-1.0
      DXSIGN=1.0D0
      IF (PARMU*(CVMAXP-CVMAXM) > 2.0D0*SUMM) DXSIGN=-1.0D0
      !DXSIGN=SIGN(1.0D0, 2.0D0*SUMM-PARMU*(CVMAXP-CVMAXM))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Update the elements of SIM and SIMI, and set the next X.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      TEMP=0.0
      TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          !DX(I)=SIGN(1.0D0, 2.0D0*SUMM-PARMU*(CVMAXP-CVMAXM))*DX(I)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
C models. In such a case, we terminate the code. Otherwise, the behavior
C of TRSTLP is not predictable, and Segmentation Fault or infinite
C cycling may happen. This is because any equality/inequality comparison
C involving NaN returns FALSE, which can lead to unintended behavior of
C the code, including uninitialized indices.
      DO J = 1, N
          DO I = 1, N
              IF (SIMI(I, J) /= SIMI(I, J)) THEN
                  IF (IPRINT >= 1) PRINT 376
  376             FORMAT (/3X,'Return from subroutine COBYLA because ',
     1            'rounding errors are becoming damaging.')
                  INFO = 7
                  GOTO 600
              END IF
          END DO
      END DO
      DO J = 1, MP
          DO I = 1, N
              IF (A(I, J) /= A(I, J)) THEN
                  INFO = -3
                  GOTO 600
              END IF
          END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL TRSTLP (N,M,A,CON,RHO,DX,IFULL,IACT,W(IZ),W(IZDOTA),
     1  W(IVMC),W(ISDIRN),W(IDXNEW),W(IVMD))

      !IF (IFULL == 0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          TEMP=0.0
          TEMP=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              TEMP=TEMP+DX(I)**2
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TEMP .LT. 0.25*RHO*RHO) THEN
!          IF (TEMP < 0.25D0*RHO*RHO) THEN
          IF (TEMP < 0.25D0*RHO**2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IBRNCH=1
              GOTO 550
          END IF
      !END IF
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
        SUMM = 0.0D0
        do I = 1, N
            SUMM = SUMM + A(I, K) * DX(I)
        end do
        SUMM = CON(K) - SUMM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUMM)
          IF (K < MP) RESNEW=DMAX1(RESNEW,SUMM)
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
C      IF (PREREC .GT. 0.0) BARMU=SUMM/PREREC
C      IF (PARMU .LT. 1.5*BARMU) THEN
C          PARMU=2.0*BARMU
      IF (PREREC > 0.0D0) BARMU=SUMM/PREREC
      IF (PARMU < 1.5D0*BARMU) THEN
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Zaikun 20211108
          !PARMU=2.0D0*BARMU
          PARMU=min(2.0D0*BARMU, HUGENUM)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (IPRINT >= 2) PRINT 410, PARMU
  410     FORMAT (/3X,'Increase in PARMU to',1PE13.6)
          PHI=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
          DO J=1,N
              TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
              IF (TEMP < PHI) THEN
                  GOTO 140
              END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (TEMP .EQ. PHI .AND. PARMU .EQ. 0.0) THEN
              IF (TEMP == PHI .AND. PARMU == 0.0D0) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IF (DATMAT(MPP,J) < DATMAT(MPP,NP)) THEN
                      GOTO 140
                  END IF
              END IF
          END DO
      END IF
      PREREM=PARMU*PREREC-SUMM
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun 20211008
      if (is_nan(trured)) then trured = -hugenum
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C     Begin the operations that decide whether x(*) should replace one of the
C     vertices of the current simplex, the change being mandatory if TRURED is
C     positive. Firstly, JDROP is set to the index of the vertex that is to be
C     replaced.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      RATIO=0.0
C      IF (TRURED .LE. 0.0) RATIO=1.0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Zaikun 20211108: The follow two versions differ when TRUERED=NaN
      !RATIO=0.0D0
      !IF (TRURED <= 0.0D0) RATIO=1.0D0
      RATIO = 1.0D0
      IF (TRURED > 0.0D0) RATIO=0.0D0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !! Zaikun 20210812
              !if (temp /= temp) then
              !    jdrop = 0
              !    l = 0
              !    exit
              !end if
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              IF (TEMP > EDGMAX) THEN
                  L=J
                  EDGMAX=TEMP
              END IF
          END IF
      END DO
      IF (L > 0) JDROP=L



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun 20211107
      ! With the following, JDROP should not be 0 unless VETA is all NaN,
      ! which should not happen if X0 does not contain NaN, trust-region/
      ! geometry step never contains NaN, and we exit once we arrive at
      ! an iterate containing Inf (due to overflow).
      ! Nevertheless, we still keep the code from 20190820 that handles
      ! the case JDROP = 0 for robustness.
      if (trured > 0 .and. jdrop < 1) then
      jdrop=int(maxloc(veta(1:n),mask=(.not.is_nan(veta(1:n))), dim=1))
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !summ = sum(abs(sigbar(1:n)))
      !if (summ /= summ) jdrop = 0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 20190820:
C When JDROP=0, the algorithm decides not to include the trust-region
C trial point X into the simplex, because X is not good enough according
C to the merit function PHI = F + PARMU*RESMAX. In this case, X will
C simply be discarded in the original code. However, this decision
C depends on the value of PARMU. When PARMU is updated later, the
C discarded X might turn out better, sometimes even better than
C SIM(:, NP), which is supposed to be the best point in the simplex.
C For this reason, we save the to-be-discarded X in XSAV and compare
C them with SIM(:, NP) right before exiting. If a vector in XSAV turns
C out better than SIM(:, NP), we replace SIM(:, NP) by this vector
C
C When JDROP > 0, SIM(:, JDROP) will be removed from the simplex
C according to PHI with the current PARMU. Similar to X, SIM(:, JDROP)
C may turn out better when PARMU is updated. Therefore, XSAV also takes
C SIM(:, JDROP) into account.
C
C We save at most NSMAX to-be-discarded X.
C
      IF (JDROP == 0) THEN
          DO I = 1, N
              XDROP(I) = X(I)
          END DO
          DO K = 1, MPP
              DATDROP(K) = CONSAV(K)
          END DO
      ELSE ! JDROP < NP is guaranteed
          DO I = 1, N
              XDROP(I) = SIM(I, NP) + SIM(I, JDROP)
          END DO
          DO K = 1, MPP
              DATDROP(K) = DATMAT(K, JDROP)
          END DO
      END IF
!      CALL SAVEX (XDROP(1:N), DATDROP(1:MPP), XSAV(1:N, 1:NSMAX),
!     1     DATSAV(1:MPP, 1:NSMAX), N, M, NSAV, NSMAX, CTOL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (JDROP == 0) then
          GOTO 550
      end if
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
      IF (NFVALS >= MAXFUN) GOTO 600
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
!              ELSE IF (CMAX-CMIN < PARMU*DENOM) THEN
!                  PARMU=(CMAX-CMIN)/DENOM
              ELSE
                  PARMU = min(PARMU, (CMAX-CMIN)/DENOM)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (IFULL .EQ. 1) GOTO 620
C  600 DO 610 I=1,N
C  610 X(I)=SIM(I,NP)
C      F=DATMAT(MP,NP)
C      RESMAX=DATMAT(MPP,NP)
C
C      Zaikun 01-06-2019:
C      Why go to 620 directly without setting X and F? This seems
C      INCORRECT, because it may lead to a return with X and F that
C      are not the best available.
C      The following code defines X as an "optimal" one among the
C      vectors:
C      DATSAV(:, 1:NSAV), (when NSAV >= 1),
C      SIM(:, NP), and SIM(:, 1:min(NP-1, NFVALS-2)) (when NFVALS>=2).
C      Here, X being "optimal" means
C      1. the constraint violation of X is at most RESREF
C      2. no other vector is better than X according to ISBETTER with
C      the current PARMU.
C
C      Note:
C      0. The last evaluated X and its function/constraint information
C      are saved in [X, CONSAV, F, RESMAX].
C      1. When NFVALS=1, SIM and DATMAT have not been initialized yet.
C      2. When 2<=NFVALS<=NP, the first evaluated X are saved in
C      SIM(:, NP), its function/constraint in DATMAT(:, NP), while the
C      other X are saved in SIM(:, NFVALS-1), its function/constraint
C      in DATMAT(:, NFVALS-1). However, when the code arrives at line 600,
C      [X, CON, F, RESMAX] may have not been saved into SIM(:, NFVALS-1)
C      and DATMAT(:, NFVALS-1) yet. That is why we check SIM up to
C      NFVALS-2 instead of NFVALS-1.
  600 DO K = 1, M
           CON(K) = CONSAV(K)
      END DO

!      CALL SAVEX (X(1:N), CONSAV(1:MPP), XSAV(1:N, 1:NSMAX),
!     1     DATSAV(1:MPP, 1:NSMAX), N, M, NSAV, NSMAX, CTOL)

      !fhist(1:n+2+nsav) = [datmat(mp, 1:np), DATSAV(mp, 1:nsav), f]
      !chist(1:n+2+nsav) = [datmat(mpp, 1:np),DATSAV(mpp,1:nsav),resmax]
      fhist(1:nsav) = DATSAV(mp, 1:nsav)
      chist(1:nsav) = DATSAV(mpp,1:nsav)
      if (any(fhist(1:nsav)< HUGEFUN .and. chist(1:nsav)
     1 < HUGECON)) then
          fref = HUGEFUN
          resref = HUGECON
      elseif (any(fhist(1:nsav) < HUGENUM .and. chist(1:nsav)
     1     < HUGECON)) then
          fref = hugenum
          resref = hugecon
      elseif (any(fhist(1:nsav) < HUGEFUN .and. chist(1:nsav)
     1     < HUGENUM)) then
          fref = hugefun
          resref = hugenum
      else
          fref = hugenum
          resref = hugenum
      end if


      if (any (fhist(1:nsav) < fref .and. chist(1:nsav)
     1 < resref)) then
      cmin = minval(chist(1:nsav),mask= (fhist(1:nsav) < FREF))
      PARMU = MAX(PARMU, 1.0D8)
      RESREF = 2.0D0*max(cmin - ctol, 0.0D0)
      IF (NFVALS >= 2) THEN ! See the comments above for why NFVALS>2
          resmax = huge(0.0D0)
          f = huge(0.0D0)
          !if (.not. (resmax-ctol <= resref .and. f < FREF)) then
          !    resmax = huge(0.0D0)
          !    f = huge(0.0D0)
          !end if
          !DO J = 1, NP
!C See the comments above for why to check these J
          !    IF (max(DATMAT(MPP, J)-ctol, 0.0D0) <= RESREF .and.
      !1            datmat(mp, J) < FREF) THEN
          !        IF (datmat(mp, j) + parmu*max(datmat(mpp, j)-ctol,
      !1             0.0D0) < f + parmu*max(resmax-ctol, 0.0D0)) THEN
          !            if (J <= N) then
          !                DO I = 1, N
          !                    X(I) = SIM(I, J) + SIM(I, NP)
          !                END DO
          !            else
          !                X(1:N) = SIM(1:N, NP)
          !            end if
          !            F = DATMAT(MP, J)
          !            RESMAX = DATMAT(MPP, J)
          !            DO K = 1, M
          !                CON(K) = DATMAT(K, J)
          !            END DO
          !        END IF
          !    END IF
          !END DO
          DO J = 1, NSAV
              IF (max(DATSAV(MPP, J)-ctol, 0.0D0) <= RESREF .and.
     1         datsav(mp, j) < FREF) THEN
                  IF (datsav(mp, j)+ parmu*max(datsav(mpp, j)-ctol,
     1              0.0D0) < f + parmu*max(resmax-ctol, 0.0D0)) THEN
                      DO I = 1, N
                          X(I) = XSAV(I, J)
                      END DO
                      F = DATSAV(MP, J)
                      RESMAX = DATSAV(MPP, J)
                      DO K = 1, M
                          CON(K) = DATSAV(K, J)
                      END DO
                  END IF
              ENDIF
          END DO


          PHIMIN =  f + parmu*max(resmax-ctol, 0.0D0)

          !DO J = 1, NP
!C See the comments above for why to check these J
          !    IF (max(DATMAT(MPP, J)-ctol, 0.0D0) <= resref .and.
      !1         datmat(mp, j) < fref .and.
      !1             datmat(mp, j) + parmu*max(datmat(mpp, j)-ctol,
      !1              0.0D0) <= phimin .and. datmat(mpp, j) < resmax) THEN
          !            if (J <= N) then
          !                DO I = 1, N
          !                    X(I) = SIM(I, J) + SIM(I, NP)
          !                END DO
          !            else
          !                X(1:N) = SIM(1:N, NP)
          !            end if
          !            F = DATMAT(MP, J)
          !            RESMAX = DATMAT(MPP, J)
          !            DO K = 1, M
          !                CON(K) = DATMAT(K, J)
          !            END DO
          !    END IF
          !END DO
          DO J = 1, NSAV
              IF (max(DATSAV(MPP, J)-ctol, 0.0D0) <= resref .and.
     1        datsav(mp, j) < fref .and.
     1             datsav(mp, j) + parmu*max(datsav(mpp, j)-ctol,
     1              0.0D0) <= phimin .and. datsav(mpp, j) < resmax) THEN
                      DO I = 1, N
                          X(I) = XSAV(I, J)
                      END DO
                      F = DATSAV(MP, J)
                      RESMAX = DATSAV(MPP, J)
                      DO K = 1, M
                          CON(K) = DATSAV(K, J)
                      END DO
              ENDIF
          END DO


          RESREF = resmax

          !DO J = 1, NP
!C See the comments above for why to check these J
          !    IF (DATMAT(MPP, J) <= RESREF .and. datmat(mp, j) < f) THEN
          !            if (J <= N) then
          !                DO I = 1, N
          !                    X(I) = SIM(I, J) + SIM(I, NP)
          !                END DO
          !            else
          !                X(1:N) = SIM(1:N, NP)
          !            end if
          !            F = DATMAT(MP, J)
          !            RESMAX = DATMAT(MPP, J)
          !            DO K = 1, M
          !                CON(K) = DATMAT(K, J)
          !            END DO
          !    END IF
          !END DO
          DO J = 1, NSAV
              IF (DATSAV(MPP, J) <= RESREF .and. datsav(mp, j) < f) THEN
                      DO I = 1, N
                          X(I) = XSAV(I, J)
                      END DO
                      F = DATSAV(MP, J)
                      RESMAX = DATSAV(MPP, J)
                      DO K = 1, M
                          CON(K) = DATSAV(K, J)
                      END DO
              ENDIF
          END DO

!
!
!          FREF = F
!          DO J = 1, NP
!C See the comments above for why to check these J
!              IF (DATMAT(MPP, J) < RESmax .and.
!     1             datmat(mp, j) + parmu*max(datmat(mpp, j)-ctol,
!     1             0.0D0) <= phimin .and. datmat(mp, j) <= fref) then
!                      if (J <= N) then
!                          DO I = 1, N
!                              X(I) = SIM(I, J) + SIM(I, NP)
!                          END DO
!                      else
!                          X(1:N) = SIM(1:N, NP)
!                      end if
!                      F = DATMAT(MP, J)
!                      RESMAX = DATMAT(MPP, J)
!                      DO K = 1, M
!                          CON(K) = DATMAT(K, J)
!                      END DO
!              END IF
!          END DO
!          DO J = 1, NSAV
!              IF (DATSAV(MPP, J) < RESmax .and.
!     1             datsav(mp, j) + parmu*max(datsav(mpp, j)-ctol,
!     1              0.0D0) <= phimin .and. datsav(mp, j) <= fref) THEN
!                      DO I = 1, N
!                          X(I) = XSAV(I, J)
!                      END DO
!                      F = DATSAV(MP, J)
!                      RESMAX = DATSAV(MPP, J)
!                      DO K = 1, M
!                          CON(K) = DATSAV(K, J)
!                      END DO
!              ENDIF
!          END DO



      END IF
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  620 IF (IPRINT >= 1) THEN
          PRINT 70, NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
          IF (IPTEM < N) PRINT 80, (X(I),I=IPTEMP,N)
      END IF
      MAXFUN=NFVALS

      close(17)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 20190820: See the comments below line number 480
      SUBROUTINE SAVEX (XDROP, DATDROP, XSAV, DATSAV, N, M, NSAV,
     +    NSMAX, CTOL)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!-----------------------!!!!!!
      USE DIRTY_TEMPORARY_MOD4POWELL_MOD!
      !!!!!!-----------------------!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless
C XDROP is dominated by a vector in XSAV(:, 1:NSAV). If XDROP dominates
C some vectors in XSAV(:, 1:NSAV), then these vectors will be removed.
C If XDROP does not dominate any of XSAV(:, 1:NSAV) but NSAV=NSMAX,
C then we remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NSAV).
C
C When COBYLA calls this subroutine, XDROP is a vector to be "dropped",
C and  DATDROP contains its function/constraint inforation (in particular,
C DATDROP(MP) = F(XDROP), and DATDROP(MPP) = RESMAX(X)). XSAV and DATSAV
C save at most NSMAX vectors "dropped" by COBYLB and their function/constraint
C information. Only XSAV(:, 1:NSAV) and DATSAV(:, 1:NSAV) contains such
C vectors, while XSAV(:, NSAV+1:NSMAX) and DATSAV(:, NSAV+1:NSMAX) are
C not initialized yet.
C
C Note: X dominates Y if and only if the function/constraint of X is
C better than the function/constraint of Y accoring to the ISBETTER
C subroutine with PARMU = -1.0D0. Indeed, PARMU can be any negative
C number. This is because, due to the implementation of ISBETTER,
C X dominates Y (i.e., X is better than Y with PARMU < 0)
C ==> X is better than Y with any PARMU >= 0,
C ==> X is better than Y regardless of PARMU.
C Fot this reason, it is sufficient to save all the "dropped" X that are not
C dominated by any vectors in XSAV (as we do in this subroutine),
C because the other X are always worse than some vector in XSAV.
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, M, NSMAX
      INTEGER, INTENT(INOUT) :: NSAV
      REAL(KIND(0.0D0)), INTENT(IN) :: XDROP(N), DATDROP(M+2), CTOL
      REAL(KIND(0.0D0)), INTENT(INOUT) :: XSAV(N, NSMAX)
      REAL(KIND(0.0D0)), INTENT(INOUT) :: DATSAV(M+2, NSMAX)
      REAL(KIND(0.0D0)) :: PARMU
      INTEGER :: MP, MPP, I, J, K, L, IREMOVE(NSMAX), NREMOVE
      LOGICAL :: BETTER

      IF (NSMAX <= 0) RETURN ! Do nothing if NSMAX=0

      MP = M + 1
      MPP = M + 2
      PARMU = -1.0D0 ! See the comments above for why PARMU = -1

C IREMOVE: indices of vectors to remove from XSAV
C NREMOVE: number of vectors to remove from XSAV
      DO J = 1, NSMAX
C It is not enough to initialize IREMOVE(1:NSAV), because NSAV may be
C incremented by 1 latter, and then IREMOVE(NSAV+1) will be accessed.
          IREMOVE(J) = -1
      END DO
      NREMOVE = 0
      DO I = 1, NSAV
C If XDROP is dominated by XSAV(:, I), then return immediately,
C because XDROP should not be inluded into XSAV.
          CALL ISBETTER (DATDROP(MP), DATDROP(MPP), DATSAV(MP, I),
     1         DATSAV(MPP, I), PARMU, CTOL, BETTER)
          IF (BETTER) RETURN
C If XDROP dominates XSAV(:, I), then increment NREMOVE by 1 and save
C I as IREMOVE(NREMOVE).
          CALL ISBETTER (DATSAV(MP, I), DATSAV(MPP, I), DATDROP(MP),
     1         DATDROP(MPP), PARMU, CTOL, BETTER)
          IF (BETTER) THEN
              NREMOVE = NREMOVE + 1
              IREMOVE(NREMOVE) = I
          END IF
      END DO

C The code did not return and NREMOVE=0 (no vector to remove from XSAV).
C If NSAV=NSMAX, then remove XSAV(:, 1); otherwise, increment NSAV by
C 1 and then "remove" XSAV(:, NSAV) (though there is no vector saved there)
      IF (NREMOVE == 0) THEN
          IF (NSAV == NSMAX) THEN
              IREMOVE(1) = 1
          ELSE
              NSAV = NSAV + 1
              IREMOVE(1) = NSAV
          END IF
          NREMOVE = 1
      END IF

C Remove from XSAV the vectors indexed by IREMOVE
      J = 1 ! Index of IREMOVE
      K = 1 ! Index of the new XSAV
      DO I = 1, NSAV ! Index of the old XSAV
          IF (I == IREMOVE(J)) THEN
              J = J + 1 ! Do nothing but incrementing J by 1
          ELSE ! Set XSAV(:, K) = XSAV(:, I)
              DO L = 1, N
                  XSAV(L, K) = XSAV(L, I)
              END DO
              DO L = 1, MPP
                  DATSAV(L, K) = DATSAV(L, I)
              END DO
              K = K + 1 ! Increment K by 1
          END IF
      END DO

C Set the number of vectors in the new XSAV
      NSAV = NSAV - NREMOVE + 1

C Save XDROP in XSAV(:, NSAV) (with NSAV updated as above)
      IF (NSAV >= 1 .AND. NSAV <= NSMAX) THEN
          ! This inequlity is not guaranteed if NSMAX=0, where NSAV will
          ! be 0 and hence a Segmentation Fault when accessing
          ! XSAV(:, NSAV). Although we return immediately if NSMAX=0,
          ! we still check this inequlity to be safe.
          DO L = 1, N
              XSAV(L, NSAV) = XDROP(L)
          END DO
          DO L = 1, MPP
              DATSAV(L, NSAV) = DATDROP(L)
          END DO
      END IF

      END SUBROUTINE SAVEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 20190820:
      SUBROUTINE ISBETTER (F0, CONV0, F, CONV, PARMU, CTOL, BETTER)
C This subroutine compares whether (F, R) is (strictly) better than
C (F0, R0) in the sense of decreasing the merit function PHI = F + PARMU*R.
C It takes care of the cases where some of these values are NaN or Inf.
C At return, BETTER=true iff (F, R) is better than (F0, R0).
C Note:
C 1. A = Inf if and only if A .GT. HUGENUM (defined below);
C 2. A = NaN if and only if A .NE. A;
C 3. If A = NaN, then any comparison (except .NE.) with another number B
C    (can be Inf or NaN) returns false.
C

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!-----------------------!!!!!!
        USE DIRTY_TEMPORARY_MOD4POWELL_MOD!
        !!!!!!-----------------------!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      real(RP), intent(IN) :: f0
      real(RP), intent(IN) :: conv0
      real(RP), intent(IN) :: f
      real(RP), intent(IN) :: conv
      real(RP), intent(IN) :: parmu
      real(RP), intent(IN) :: ctol
      logical :: better

      ! Local variables
      real(RP) :: cref
      logical :: f0infnan
      logical :: finfnan
      logical :: fle
      logical :: flt
      logical :: c0infnan
      logical :: cinfnan
      logical :: cle
      logical :: clt

      better = .false.
      better = better .or. (.not. is_nan(f) .and. .not. is_nan(conv)
     & .and. (is_nan(f0) .or. is_nan(conv0)))
      better = better .or. (f <= f0 .and. conv < conv0)
      better = better .or. (f < f0 .and. conv <= conv0)

      cref = 10.0D0*max(ctol, epsilon(ctol))
      better = better .or. (f < HUGENUM .and. conv <= ctol .and.
     1 (conv0 > cref .or. is_nan(conv0)))

      END SUBROUTINE
