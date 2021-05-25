      SUBROUTINE LINCOB (N,NPT,M,AMAT,B,X,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XSAV,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  STEP,SP,XNEW,IACT,RESCON,QFAC,RFAC,PQW,W)
     2  STEP,SP,XNEW,IACT,RESCON,QFAC,RFAC,PQW,W,F,INFO,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION AMAT(N,*),B(*),X(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XSAV(*),XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),
     2  ZMAT(NPT,*),STEP(*),SP(*),XNEW(*),IACT(*),RESCON(*),
     3  QFAC(N,*),RFAC(*),PQW(*),W(*)
C
C     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
C       identical to the corresponding arguments in SUBROUTINE LINCOA.
C     AMAT is a matrix whose columns are the constraint gradients, scaled
C       so that they have unit length.
C     B contains on entry the right hand sides of the constraints, scaled
C       as above, but later B is modified for variables relative to XBASE.
C     XBASE holds a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XPT contains the interpolation point coordinates relative to XBASE.
C     FVAL holds the values of F at the interpolation points.
C     XSAV holds the best feasible vector of variables so far, without any
C       shift of origin.
C     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
C       the feasible vector of variables that provides the least calculated
C       F so far, this vector being the current trust region centre.
C     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
C     HQ holds the explicit second derivatives of the quadratic model.
C     PQ contains the parameters of the implicit second derivatives of the
C       quadratic model.
C     BMAT holds the last N columns of the big inverse matrix H.
C     ZMAT holds the factorization of the leading NPT by NPT submatrix
C       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
C       where the elements of DZ are plus or minus one, as specified by IDZ.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     STEP is employed for trial steps from XOPT. It is also used for working
C       space when XBASE is shifted and in PRELIM.
C     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
C       followed by STEP^T XPT(K,.), K=1,2,...,NPT.
C     XNEW is the displacement from XBASE of the vector of variables for
C       the current calculation of F, except that SUBROUTINE TRSTEP uses it
C       for working space.
C     IACT is an integer array for the indices of the active constraints.
C     RESCON holds useful information about the constraint residuals. Every
C       nonnegative RESCON(J) is the residual of the J-th constraint at the
C       current trust region centre. Otherwise, if RESCON(J) is negative, the
C       J-th constraint holds as a strict inequality at the trust region
C       centre, its residual being at least |RESCON(J)|; further, the value
C       of |RESCON(J)| is at least the current trust region radius DELTA.
C     QFAC is the orthogonal part of the QR factorization of the matrix of
C       active constraint gradients, these gradients being ordered in
C       accordance with IACT. When NACT is less than N, columns are added
C       to QFAC to complete an N by N orthogonal matrix, which is important
C       for keeping calculated steps sufficiently close to the boundaries
C       of the active constraints.
C     RFAC is the upper triangular part of this QR factorization, beginning
C       with the first diagonal element, followed by the two elements in the
C       upper triangular part of the second column and so on.
C     PQW is used for working space, mainly for storing second derivative
C       coefficients of quadratic functions. Its length is NPT+N.
C     The array W is also used for working space. The required number of
C       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NH=(N*NP)/2
      NPTM=NPT-NP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ALMOST_INFINITY=HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 15-08-2019
C See the comments below line number 210
      IMPRV = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
C       ZMAT and SP for the first iteration. An important feature is that,
C       if the interpolation point XPT(K,.) is not feasible, where K is any
C       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
C       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
C       is set so that XPT(KOPT,.) is the initial trust region centre.
C
      CALL PRELIM (N,NPT,M,AMAT,B,X,RHOBEG,IPRINT,XBASE,XPT,FVAL,
     1  XSAV,XOPT,GOPT,KOPT,HQ,PQ,BMAT,ZMAT,IDZ,NDIM,SP,RESCON,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  STEP,PQW,W)
     2  STEP,PQW,W,F,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom (on 04-06-2019):
      IF (F /= F .OR. F > ALMOST_INFINITY) THEN
          FOPT=FVAL(KOPT)
          INFO=-2
          GOTO 600
      END IF
C     By Tom/Zaikun (on 04-06-2019/07-06-2019):
C     Note that we should NOT compare F and FTARGET, because X may not
C     be feasible at the exit of PRELIM.
      IF (FVAL(KOPT) <= FTARGET) THEN
          F=FVAL(KOPT)
          X(1:N)=XSAV(1:N)
          INFO=1
          GOTO 616
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Begin the iterative procedure.
C
      NF=NPT
      FOPT=FVAL(KOPT)
      RHO=RHOBEG
      DELTA=RHO
      IFEAS=0
      NACT=0
      ITEST=3
   10 KNEW=0
      NVALA=0
      NVALB=0
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C       to BMAT that do not depend on ZMAT.
C
   20 FSAVE=FOPT
      XOPTSQ=ZERO
      DO I=1,N
          XOPTSQ=XOPTSQ+XOPT(I)**2
      END DO
      IF (XOPTSQ >= 1.0D4*DELTA*DELTA) THEN
          QOPTSQ=0.25D0*XOPTSQ
          DO K=1,NPT
              SUM=ZERO
              DO I=1,N
                  SUM=SUM+XPT(K,I)*XOPT(I)
              END DO
              SUM=SUM-HALF*XOPTSQ
              W(NPT+K)=SUM
              SP(K)=ZERO
              DO I=1,N
                  XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
                  STEP(I)=BMAT(K,I)
                  W(I)=SUM*XPT(K,I)+QOPTSQ*XOPT(I)
                  IP=NPT+I
                  DO J=1,I
                      BMAT(IP,J)=BMAT(IP,J)+STEP(I)*W(J)+W(I)*STEP(J)
                  END DO
              END DO
          END DO
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO K=1,NPTM
              SUMZ=ZERO
              DO I=1,NPT
                  SUMZ=SUMZ+ZMAT(I,K)
                  W(I)=W(NPT+I)*ZMAT(I,K)
              END DO
              DO J=1,N
                  SUM=QOPTSQ*SUMZ*XOPT(J)
                  DO I=1,NPT
                      SUM=SUM+W(I)*XPT(I,J)
                  END DO
                  STEP(J)=SUM
                  IF (K < IDZ) SUM=-SUM
                  DO I=1,NPT
                      BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
                  END DO
              END DO
              DO I=1,N
                  IP=I+NPT
                  TEMP=STEP(I)
                  IF (K < IDZ) TEMP=-TEMP
                  DO J=1,I
                      BMAT(IP,J)=BMAT(IP,J)+TEMP*STEP(J)
                  END DO
              END DO
          END DO
C
C     Update the right hand sides of the constraints.
C
          IF (M > 0) THEN
              DO J=1,M
                  TEMP=ZERO
                  DO I=1,N
                      TEMP=TEMP+AMAT(I,J)*XOPT(I)
                  END DO
                  B(J)=B(J)-TEMP
              END DO
          END IF
C
C     The following instructions complete the shift of XBASE, including the
C       changes to the parameters of the quadratic model.
C
          IH=0
          DO J=1,N
              W(J)=ZERO
              DO K=1,NPT
                  W(J)=W(J)+PQ(K)*XPT(K,J)
                  XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
              END DO
              DO I=1,J
                  IH=IH+1
                  HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
                  BMAT(NPT+I,J)=BMAT(NPT+J,I)
              END DO
          END DO
          DO J=1,N
              XBASE(J)=XBASE(J)+XOPT(J)
              XOPT(J)=ZERO
              XPT(KOPT,J)=ZERO
          END DO
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 21-03-2020
C Exit if BMAT or ZMAT contians NaN
      DO J = 1,N
          DO I = 1,NDIM
              IF (BMAT(I,J) /= BMAT(I,J)) THEN
                  INFO = -3
                  GOTO 600
              END IF
          END DO
      END DO
      DO J = 1,NPTM
          DO I = 1,NPT
              IF (ZMAT(I,J) /= ZMAT(I,J)) THEN
                  INFO = -3
                  GOTO 600
              END IF
          END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C     In the case KNEW=0, generate the next trust region step by calling
C       TRSTEP, where SNORM is the current trust region radius initially.
C       The final value of SNORM is the length of the calculated step,
C       except that SNORM is zero on return if the projected gradient is
C       unsuitable for starting the conjugate gradient iterations.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
C models. In such a case, we terminate the code. Otherwise, the behavior
C of TRSTEP or QMSTEP is not predictable, and Segmentation Fault or
C infinite cycling may happen. This is because any equality/inequality
C comparison involving NaN returns FALSE, which can lead to unintended
C behavior of the code, including uninitialized indices, which can lead
C to segmentation faults. 
      DO J = 1,N
          IF (GOPT(J) /= GOPT(J)) THEN
              INFO = -3
              GOTO 600
          END IF
      END DO
      DO I = 1, NH
          IF (HQ(I) /= HQ(I)) THEN
              INFO = -3
              GOTO 600
          END IF
      END DO
      DO I = 1, NPT
          IF (PQ(I) /= PQ(I)) THEN
              INFO = -3
              GOTO 600
          END IF
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DELSAV=DELTA
      KSAVE=KNEW
      IF (KNEW == 0) THEN
          SNORM=DELTA
          DO I=1,N
              XNEW(I)=GOPT(I)
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 19-03-2020: B is never used in TRSTEP
C          CALL TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
          CALL TRSTEP (N,NPT,M,AMAT,XPT,HQ,PQ,NACT,IACT,RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     1      QFAC,RFAC,SNORM,STEP,XNEW,W,W(M+1),PQW,PQW(NP),W(M+NP))
C
C     A trust region step is applied whenever its length, namely SNORM, is at
C       least HALF*DELTA. It is also applied if its length is at least 0.1999
C       times DELTA and if a line search of TRSTEP has caused a change to the
C       active set. Otherwise there is a branch below to label 530 or 560.
C
          TEMP=HALF*DELTA
          IF (XNEW(1) >= HALF) TEMP=0.1999D0*DELTA
          IF (SNORM <= TEMP) THEN
              DELTA=HALF*DELTA
              IF (DELTA <= 1.4D0*RHO) DELTA=RHO
              NVALA=NVALA+1
              NVALB=NVALB+1
              TEMP=SNORM/RHO
              IF (DELSAV > RHO) TEMP=ONE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 24-07-2019
C              IF (TEMP .GE. HALF) NVALA=ZERO
C              IF (TEMP .GE. TENTH) NVALB=ZERO
              IF (TEMP >= HALF) NVALA=0
              IF (TEMP >= TENTH) NVALB=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (DELSAV > RHO) GOTO 530
              IF (NVALA < 5 .AND. NVALB < 3) GOTO 530
              IF (SNORM > ZERO) KSAVE=-1
              GOTO 560
          END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 24-07-2019
C          NVALA=ZERO
C          NVALB=ZERO
          NVALA=0
          NVALB=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Alternatively, KNEW is positive. Then the model step is calculated
C       within a trust region of radius DEL, after setting the gradient at
C       XBASE and the second derivative parameters of the KNEW-th Lagrange
C       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
C
      ELSE
          DEL=DMAX1(TENTH*DELTA,RHO)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: See the comments below line number 140
C          DO 160 I=1,N
C  160     W(I)=BMAT(KNEW,I)
          DO I = 1, N
              W(I)=BMAT(KNEW,I)
              IF (W(I) /= W(I)) THEN
                  INFO = -3
                  GOTO 600
              END IF
          END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO K=1,NPT
              PQW(K)=ZERO
          END DO
          DO J=1,NPTM
              TEMP=ZMAT(KNEW,J)
              IF (J < IDZ) TEMP=-TEMP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: See the comments below line number 140
C Note that the data in PQW is used in QMSTEP below
C          DO 180 K=1,NPT
C  180     PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
              DO K = 1, NPT
                  PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
                  IF (PQW(K) /= PQW(K)) THEN
                      INFO = -3
                      GOTO 600
                  END IF
              END DO
          END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: B is never used in QMSTEP
C          CALL QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,RESCON,
          CALL QMSTEP (N,NPT,M,AMAT,XPT,XOPT,NACT,IACT,RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     1      QFAC,KOPT,KNEW,DEL,STEP,W,PQW,W(NP),W(NP+M),IFEAS)
      END IF
C
C     Set VQUAD to the change to the quadratic model when the move STEP is
C       made from XOPT. If STEP is a trust region step, then VQUAD should be
C       negative. If it is nonnegative due to rounding errors in this case,
C       there is a branch to label 530 to try to improve the model.
C
      VQUAD=ZERO
      IH=0
      DO J=1,N
          VQUAD=VQUAD+STEP(J)*GOPT(J)
          DO I=1,J
              IH=IH+1
              TEMP=STEP(I)*STEP(J)
              IF (I == J) TEMP=HALF*TEMP
              VQUAD=VQUAD+TEMP*HQ(IH)
          END DO
      END DO
      DO K=1,NPT
          TEMP=ZERO
          DO J=1,N
              TEMP=TEMP+XPT(K,J)*STEP(J)
              SP(NPT+K)=TEMP
          END DO
          VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 15-08-2019
C Although very rarely, with the original code, an infinite loop can occur
C in the following scenario.
C Suppose that, at an certain iteration,
C KNEW = 0, SNORM > 0.5*DELTA > RHO, VQUAD >= 0, and
C \sum_{K=1}^NPT ||XPT(K,:)-XOPT(:)||^2 < DELTA^2
C (i.e., DELTA is large and SNORM is not small, yet VQUAD >= 0 due to
C rounding errors and XPT are not far from XOPT).
C Then the program will goto 530 and then goto 20, where XBASE may be
C shifted to the current best point, in the hope of reducing rounding
C errors and 'improve' the model. Afterwards, another trust region step
C is produced by the 'improved' model. Note that DELTA remains unchanged
C in this process. If the new trust region step turns out to satisfy
C SNORM > 0.5*DELTA and VQUAD >= 0 again (i.e., the 'improved' model
C still suffers from rounding errors), then the program will goto 530
C and then goto 20, where shifting will not happen because either XBASE
C was already shifted to the current best point in last step, or XBASE
C is close to the current best point. Consequently, the model will
C remain unchanged, and produce the same trust region step again. This
C leads to an infinite loop.
C The infinite loop did happen when the MATLAB interface was applied to
C min atan(x+100) s.t. x<=-99 (x0=-99, npt=3, rhobeg=1, rhoend=1e-6).
C The problem does not exist in NEWUOA or BOBYQA, where the program will
C exit immediately when VQUAD >= 0.
C To prevent such a loop, here we use IMPRV to record whether the path
C 530 --> 20 has already happened for last trust region step. IMPRV=1
C implies that last trust region step satisfies VQUAD >= 0 and followed
C 530 --> 20. With IMPRV=1, if VQUAD is again nonnegative for the new trust
C region step, we should not goto 530 but goto 560, where IMPRV will be
C set to 0 and DELTA will be reduced. Otherwise, an infinite loop would happen.
C      IF (KSAVE .EQ. 0 .AND. VQUAD .GE. ZERO) GOTO 530
      IF (KSAVE == 0 .AND. .NOT. (VQUAD < ZERO)) THEN
          IF (IMPRV == 1) THEN
              GOTO 560
          ELSE
              IMPRV = 1
              GOTO 530
          END IF
      ELSE
          IMPRV = 0
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     Calculate the next value of the objective function. The difference
C       between the actual new value of F and the value predicted by the
C       model is recorded in DIFF.
C
  220 NF=NF+1
      IF (NF > MAXFUN) THEN
          NF=NF-1
          IF (IPRINT > 0) PRINT 230
  230     FORMAT (/4X,'Return from LINCOA because CALFUN has been',
     1      ' called MAXFUN times.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 600
      END IF
      XDIFF=ZERO
      DO I=1,N
          XNEW(I)=XOPT(I)+STEP(I)
          X(I)=XBASE(I)+XNEW(I)
          XDIFF=XDIFF+(X(I)-XSAV(I))**2
      END DO
      XDIFF=DSQRT(XDIFF)
      IF (KSAVE == -1) XDIFF=RHO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (XDIFF .LE. TENTH*RHO .OR. XDIFF .GE. DELTA+DELTA) THEN
      IF (.NOT.(XDIFF > TENTH*RHO .AND. XDIFF <DELTA+DELTA)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IFEAS=0
          IF (IPRINT > 0) PRINT 250
  250     FORMAT (/4X,'Return from LINCOA because rounding errors',
     1      ' prevent reasonable changes to X.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 600
      END IF
      IF (KSAVE <= 0) IFEAS=1
      F=DFLOAT(IFEAS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO I=1,N
          IF (X(I) /= X(I)) THEN
              F = X(I) ! Set F to NaN
              IF (NF == 1) THEN
                  FOPT = F
                  XOPT(1:N) = ZERO
              END IF
              INFO=-1
              GOTO 600
          END IF
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL CALFUN (N,X,F)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom (on 04-06-2019):
      IF (F /= F .OR. F > ALMOST_INFINITY) THEN
          IF (NF == 1) THEN
              FOPT = F
              XOPT(1:N) = ZERO
          END IF
          INFO=-2
          GOTO 600
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (IPRINT == 3) THEN
          PRINT 260, NF,F,(X(I),I=1,N)
  260     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (KSAVE .EQ. -1) GOTO 600
      IF (KSAVE == -1) THEN
          INFO=0
          GOTO 600
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIFF=F-FOPT-VQUAD
C
C     If X is feasible, then set DFFALT to the difference between the new
C       value of F and the value predicted by the alternative model.
C
      IF (IFEAS == 1 .AND. ITEST < 3) THEN
          DO K=1,NPT
              PQW(K)=ZERO
              W(K)=FVAL(K)-FVAL(KOPT)
          END DO
          DO J=1,NPTM
              SUM=ZERO
              DO I=1,NPT
                  SUM=SUM+W(I)*ZMAT(I,J)
              END DO
              IF (J < IDZ) SUM=-SUM
              DO K=1,NPT
                  PQW(K)=PQW(K)+SUM*ZMAT(K,J)
              END DO
          END DO
          VQALT=ZERO
          DO K=1,NPT
              SUM=ZERO
              DO J=1,N
                  SUM=SUM+BMAT(K,J)*STEP(J)
              END DO
              VQALT=VQALT+SUM*W(K)
              VQALT=VQALT+PQW(K)*SP(NPT+K)*(HALF*SP(NPT+K)+SP(K))
          END DO
          DFFALT=F-FOPT-VQALT
      END IF
      IF (ITEST == 3) THEN
          DFFALT=DIFF
          ITEST=0
      END IF
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (KSAVE == 0) THEN
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO <= TENTH) THEN
              DELTA=HALF*DELTA
          ELSE IF (RATIO <= 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,SNORM)
          ELSE
              TEMP=DSQRT(2.0D0)*DELTA
              DELTA=DMAX1(HALF*DELTA,SNORM+SNORM)
              DELTA=DMIN1(DELTA,TEMP)
          END IF
          IF (DELTA <= 1.4D0*RHO) DELTA=RHO
      END IF
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C       can be moved. If STEP is a trust region step, then KNEW is zero at
C       present, but a positive value is picked by subroutine UPDATE.
C
      CALL UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,KOPT,
     1  KNEW,PQW,W)
      IF (KNEW == 0) THEN
          IF (IPRINT > 0) PRINT 320
  320     FORMAT (/4X,'Return from LINCOA because the denominator'
     1      ' of the updating formula is zero.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 600
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 19-03-2020
C Exit if BMAT or ZMAT contians NaN
      DO J = 1,N
          DO I = 1,NDIM
              IF (BMAT(I,J) /= BMAT(I,J)) THEN
                  INFO = -3
                  GOTO 600
              END IF
          END DO
      END DO
      DO J = 1,NPTM
          DO I = 1,NPT
              IF (ZMAT(I,J) /= ZMAT(I,J)) THEN
                  INFO = -3
                  GOTO 600
              END IF
          END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C     If ITEST is increased to 3, then the next quadratic model is the
C       one whose second derivative matrix is least subject to the new
C       interpolation conditions. Otherwise the new model is constructed
C       by the symmetric Broyden method in the usual way.
C
      IF (IFEAS == 1) THEN
          ITEST=ITEST+1
          IF (DABS(DFFALT) >= TENTH*DABS(DIFF)) ITEST=0
      END IF
C
C     Update the second derivatives of the model by the symmetric Broyden
C       method, using PQW for the second derivative parameters of the new
C       KNEW-th Lagrange function. The contribution from the old parameter
C       PQ(KNEW) is included in the second derivative matrix HQ. W is used
C       later for the gradient of the new KNEW-th Lagrange function.
C
      IF (ITEST < 3) THEN
          DO K=1,NPT
              PQW(K)=ZERO
          END DO
          DO J=1,NPTM
              TEMP=ZMAT(KNEW,J)
              IF (TEMP /= ZERO) THEN
                  IF (J < IDZ) TEMP=-TEMP
                  DO K=1,NPT
                      PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
                  END DO
              END IF
          END DO
          IH=0
          DO I=1,N
              W(I)=BMAT(KNEW,I)
              TEMP=PQ(KNEW)*XPT(KNEW,I)
              DO J=1,I
                  IH=IH+1
                  HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
              END DO
          END DO
          PQ(KNEW)=ZERO
          DO K=1,NPT
              PQ(K)=PQ(K)+DIFF*PQW(K)
          END DO
      END IF
C
C     Include the new interpolation point with the corresponding updates of
C       SP. Also make the changes of the symmetric Broyden method to GOPT at
C       the old XOPT if ITEST is less than 3.
C
      FVAL(KNEW)=F
      SP(KNEW)=SP(KOPT)+SP(NPT+KOPT)
      SSQ=ZERO
      DO I=1,N
          XPT(KNEW,I)=XNEW(I)
          SSQ=SSQ+STEP(I)**2
      END DO
      SP(NPT+KNEW)=SP(NPT+KOPT)+SSQ
      IF (ITEST < 3) THEN
          DO K=1,NPT
              TEMP=PQW(K)*SP(K)
              DO I=1,N
                  W(I)=W(I)+TEMP*XPT(K,I)
              END DO
          END DO
          DO I=1,N
              GOPT(I)=GOPT(I)+DIFF*W(I)
          END DO
      END IF
C
C     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
C       least calculated value so far with a feasible vector of variables.
C
      IF (F < FOPT .AND. IFEAS == 1) THEN
          FOPT=F
          DO J=1,N
              XSAV(J)=X(J)
              XOPT(J)=XNEW(J)
          END DO
          KOPT=KNEW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         By Tom (on 04-06-2019):
          IF (FOPT <= FTARGET) THEN
              INFO=1
              GOTO 616
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SNORM=DSQRT(SSQ)
          DO J=1,M
              IF (RESCON(J) >= DELTA+SNORM) THEN
                  RESCON(J)=SNORM-RESCON(J)
              ELSE
                  RESCON(J)=RESCON(J)+SNORM
                  IF (RESCON(J)+DELTA > ZERO) THEN
                      TEMP=B(J)
                      DO I=1,N
                          TEMP=TEMP-XOPT(I)*AMAT(I,J)
                      END DO
                      TEMP=DMAX1(TEMP,ZERO)
                      IF (TEMP >= DELTA) TEMP=-TEMP
                      RESCON(J)=TEMP
                  END IF
              END IF
          END DO
          DO K=1,NPT
              SP(K)=SP(K)+SP(NPT+K)
          END DO
C
C     Also revise GOPT when symmetric Broyden updating is applied.
C
          IF (ITEST < 3) THEN
              IH=0
              DO J=1,N
                  DO I=1,J
                      IH=IH+1
                      IF (I < J) GOPT(J)=GOPT(J)+HQ(IH)*STEP(I)
                      GOPT(I)=GOPT(I)+HQ(IH)*STEP(J)
                  END DO
              END DO
              DO K=1,NPT
                  TEMP=PQ(K)*SP(NPT+K)
                  DO I=1,N
                      GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
                  END DO
              END DO
          END IF
      END IF
C
C     Replace the current model by the least Frobenius norm interpolant if
C       this interpolant gives substantial reductions in the predictions
C       of values of F at feasible points.
C
      IF (ITEST == 3) THEN
          DO K=1,NPT
              PQ(K)=ZERO
              W(K)=FVAL(K)-FVAL(KOPT)
          END DO
          DO J=1,NPTM
              SUM=ZERO
              DO I=1,NPT
                  SUM=SUM+W(I)*ZMAT(I,J)
              END DO
              IF (J < IDZ) SUM=-SUM
              DO K=1,NPT
                  PQ(K)=PQ(K)+SUM*ZMAT(K,J)
              END DO
          END DO
          DO J=1,N
              GOPT(J)=ZERO
              DO I=1,NPT
                  GOPT(J)=GOPT(J)+W(I)*BMAT(I,J)
              END DO
          END DO
          DO K=1,NPT
              TEMP=PQ(K)*SP(K)
              DO I=1,N
                  GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
              END DO
          END DO
          DO IH=1,NH
              HQ(IH)=ZERO
          END DO
      END IF
C
C     If a trust region step has provided a sufficient decrease in F, then
C       branch for another trust region calculation. Every iteration that
C       takes a model step is followed by an attempt to take a trust region
C       step.
C
      KNEW=0
      IF (KSAVE > 0) GOTO 20
      IF (RATIO >= TENTH) GOTO 20
C
C     Alternatively, find out if the interpolation points are close enough
C       to the best point so far.
C
  530 DISTSQ=DMAX1(DELTA*DELTA,4.0D0*RHO*RHO)
      DO K=1,NPT
          SUM=ZERO
          DO J=1,N
              SUM=SUM+(XPT(K,J)-XOPT(J))**2
          END DO
          IF (SUM > DISTSQ) THEN
              KNEW=K
              DISTSQ=SUM
          END IF
      END DO
C
C     If KNEW is positive, then branch back for the next iteration, which
C       will generate a "model step". Otherwise, if the current iteration
C       has reduced F, or if DELTA was above its lower bound when the last
C       trust region step was calculated, then try a "trust region" step
C       instead.
C
      IF (KNEW > 0) GOTO 20
      KNEW=0
      IF (FOPT < FSAVE) GOTO 20
      IF (DELSAV > RHO) GOTO 20
C
C     The calculations with the current value of RHO are complete.
C       Pick the next value of RHO.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 15-08-2019
C See the comments below line number 210
C  560 IF (RHO .GT. RHOEND) THEN
  560 IMPRV = 0
      IF (RHO > RHOEND) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DELTA=HALF*RHO
          IF (RHO > 250.0D0*RHOEND) THEN
              RHO=TENTH*RHO
          ELSE IF (RHO <= 16.0D0*RHOEND) THEN
              RHO=RHOEND
          ELSE
              RHO=DSQRT(RHO*RHOEND)
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT >= 2) THEN
              IF (IPRINT >= 3) PRINT 570
  570         FORMAT (5X)
              PRINT 580, RHO,NF
  580         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 590, FOPT,(XBASE(I)+XOPT(I),I=1,N)
  590         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 10
      END IF
C
C     Return from the calculation, after branching to label 220 for another
C       Newton-Raphson step if it has not been tried before.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INFO=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (KSAVE == -1) GOTO 220
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  600 IF (FOPT .LE. F .OR. IFEAS .EQ. 0) THEN
  600 IF (FOPT <= F .OR. IFEAS == 0 .OR. F /= F) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              X(I)=XSAV(I)
          END DO
          F=FOPT
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (IPRINT .GE. 1) THEN
  616 IF (IPRINT >= 1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          PRINT 620, NF
  620     FORMAT (/4X,'At the return from LINCOA',5X,
     1      'Number of function values =',I6)
          PRINT 590, F,(X(I),I=1,N)
      END IF
      W(1)=F
      W(2)=DFLOAT(NF)+HALF
      RETURN
      END
