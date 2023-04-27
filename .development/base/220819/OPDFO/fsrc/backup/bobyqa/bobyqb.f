      SUBROUTINE BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  SL,SU,XNEW,XALT,D,VLAG,W)
     2  SL,SU,XNEW,XALT,D,VLAG,W,F,INFO,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
     2  SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
C
C     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
C       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
C     XBASE holds a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XPT is a two-dimensional array that holds the coordinates of the
C       interpolation points relative to XBASE.
C     FVAL holds the values of F at the interpolation points.
C     XOPT is set to the displacement from XBASE of the trust region centre.
C     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
C     HQ holds the explicit second derivatives of the quadratic model.
C     PQ contains the parameters of the implicit second derivatives of the
C       quadratic model.
C     BMAT holds the last N columns of H.
C     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
C       this factorization being ZMAT times ZMAT^T, which provides both the
C       correct rank and positive semi-definiteness.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
C       All the components of every XOPT are going to satisfy the bounds
C       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
C       XOPT is on a constraint boundary.
C     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
C       vector of variables for the next call of CALFUN. XNEW also satisfies
C       the SL and SU constraints in the way that has just been mentioned.
C     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
C       in order to increase the denominator in the updating of UPDATE.
C     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
C     VLAG contains the values of the Lagrange functions at a new point X.
C       They are part of a product that requires VLAG to be of length NDIM.
C     W is a one-dimensional array that is used for working space. Its length
C       must be at least 3*NDIM = 3*(NPT+N).
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TEN=10.0D0
      TENTH=0.1D0
      TWO=2.0D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NH=(N*NP)/2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      ALMOST_INFINITY=HUGE(0.0D0)/2.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
C
C     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
C     BMAT and ZMAT for the first iteration, with the corresponding values of
C     of NF and KOPT, which are the number of calls of CALFUN so far and the
C     index of the interpolation point at the trust region centre. Then the
C     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
C     less than NPT. GOPT will be updated if KOPT is different from KBASE.
C
      CALL PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     1  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
     1  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,F,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      XOPTSQ=ZERO
      DO I=1,N
          XOPT(I)=XPT(KOPT,I)
          XOPTSQ=XOPTSQ+XOPT(I)**2
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: FSAVE is not needed any more. See line number 720.
C      FSAVE=FVAL(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom/Zaikun (on 04-06-2019/07-06-2019):
      IF (F /= F .OR. F > ALMOST_INFINITY) THEN
          INFO=-2
          GOTO 720
      END IF
C     By Tom (on 04-06-2019):
C     If F reached the target function, PRELIM will stop and BOBYQB
C     should stop here.
      IF (F <= FTARGET) THEN
          INFO=1
          GOTO 736
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (NF < NPT) THEN
          IF (IPRINT > 0) PRINT 390
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 720
      END IF
      KBASE=1
C
C     Complete the settings that are required for the iterative procedure.
C
      RHO=RHOBEG
      DELTA=RHO
      NRESC=NF
      NTRITS=0
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      NFSAV=NF
C
C     Update GOPT if necessary before the first iteration and after each
C     call of RESCUE that makes a call of CALFUN.
C
   20 IF (KOPT /= KBASE) THEN
          IH=0
          DO J=1,N
              DO I=1,J
                  IH=IH+1
                  IF (I < J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
                  GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
              END DO
          END DO
          IF (NF > NPT) THEN
              DO K=1,NPT
                  TEMP=ZERO
                  DO J=1,N
                      TEMP=TEMP+XPT(K,J)*XOPT(J)
                  END DO
                  TEMP=PQ(K)*TEMP
                  DO I=1,N
                      GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
                  END DO
              END DO
          END IF
      END IF
C
C     Generate the next point in the trust region that provides a small value
C     of the quadratic model subject to the constraints on the variables.
C     The integer NTRITS is set to the number "trust region" iterations that
C     have occurred since the last "alternative" iteration. If the length
C     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
C     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
C models. In such a case, we terminate the code. Otherwise, the behavior
C of TRBOX, ALTMOV, or RESCUE is not predictable, and Segmentation Fault or
C infinite cycling may happen. This is because any equality/inequality
C comparison involving NaN returns FALSE, which can lead to unintended
C behavior of the code, including uninitialized indices.
C
C   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
   60 DO I = 1, N
          IF (GOPT(I) /= GOPT(I)) THEN
              INFO = -3
              GOTO 720
          END IF 
      END DO
      DO I = 1, NH
          IF (HQ(I) /= HQ(I)) THEN
              INFO = -3
              GOTO 720
          END IF 
      END DO 
      DO I = 1, NPT
          IF (PQ(I) /= PQ(I)) THEN
              INFO = -3
              GOTO 720
          END IF
      END DO
      CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     1  W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM < HALF*RHO) THEN
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          IF (NF <= NFSAV+2) GOTO 650
C
C     The following choice between labels 650 and 680 depends on whether or
C     not our work with the current RHO seems to be complete. Either RHO is
C     decreased or termination occurs if the errors in the quadratic model at
C     the last three interpolation points compare favourably with predictions
C     of likely improvements to the model within distance HALF*RHO of XOPT.
C
          ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125D0*RHO*RHO
          IF (CRVMIN > ZERO .AND. ERRBIG > FRHOSQ*CRVMIN)
     1       GOTO 650
          BDTOL=ERRBIG/RHO
          DO J=1,N
              BDTEST=BDTOL
              IF (XNEW(J) == SL(J)) BDTEST=W(J)
              IF (XNEW(J) == SU(J)) BDTEST=-W(J)
              IF (BDTEST < BDTOL) THEN
                  CURV=HQ((J+J*J)/2)
                  DO K=1,NPT
                      CURV=CURV+PQ(K)*XPT(K,J)**2
                  END DO
                  BDTEST=BDTEST+HALF*CURV*RHO
                  IF (BDTEST < BDTOL) GOTO 650
              END IF
          END DO
          GOTO 680
      END IF
      NTRITS=NTRITS+1
C
C     Severe cancellation is likely to occur if XOPT is too far from XBASE.
C     If the following test holds, then XBASE is shifted so that XOPT becomes
C     zero. The appropriate changes are made to BMAT and to the second
C     derivatives of the current model, beginning with the changes to BMAT
C     that do not depend on ZMAT. VLAG is used temporarily for working space.
C
   90 IF (DSQ <= 1.0D-3*XOPTSQ) THEN
          FRACSQ=0.25D0*XOPTSQ
          SUMPQ=ZERO
          DO K=1,NPT
              SUMPQ=SUMPQ+PQ(K)
              SUM=-HALF*XOPTSQ
              DO I=1,N
                  SUM=SUM+XPT(K,I)*XOPT(I)
              END DO
              W(NPT+K)=SUM
              TEMP=FRACSQ-HALF*SUM
              DO I=1,N
                  W(I)=BMAT(K,I)
                  VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
                  IP=NPT+I
                  DO J=1,I
                      BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
                  END DO
              END DO
          END DO
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO JJ=1,NPTM
              SUMZ=ZERO
              SUMW=ZERO
              DO K=1,NPT
                  SUMZ=SUMZ+ZMAT(K,JJ)
                  VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
                  SUMW=SUMW+VLAG(K)
              END DO
              DO J=1,N
                  SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
                  DO K=1,NPT
                      SUM=SUM+VLAG(K)*XPT(K,J)
                  END DO
                  W(J)=SUM
                  DO K=1,NPT
                      BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
                  END DO
              END DO
              DO I=1,N
                  IP=I+NPT
                  TEMP=W(I)
                  DO J=1,I
                      BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
                  END DO
              END DO
          END DO
C
C     The following instructions complete the shift, including the changes
C     to the second derivative parameters of the quadratic model.
C
          IH=0
          DO J=1,N
              W(J)=-HALF*SUMPQ*XOPT(J)
              DO K=1,NPT
                  W(J)=W(J)+PQ(K)*XPT(K,J)
                  XPT(K,J)=XPT(K,J)-XOPT(J)
              END DO
              DO I=1,J
                  IH=IH+1
                  HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
                  BMAT(NPT+I,J)=BMAT(NPT+J,I)
              END DO
          END DO
          DO I=1,N
              XBASE(I)=XBASE(I)+XOPT(I)
              XNEW(I)=XNEW(I)-XOPT(I)
              SL(I)=SL(I)-XOPT(I)
              SU(I)=SU(I)-XOPT(I)
              XOPT(I)=ZERO
          END DO
          XOPTSQ=ZERO
      END IF
      IF (NTRITS == 0) GOTO 210
      GOTO 230
C
C     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
C     more expensive than the previous shift, because new matrices BMAT and
C     ZMAT are generated from scratch, which may include the replacement of
C     interpolation points whose positions seem to be causing near linear
C     dependence in the interpolation conditions. Therefore RESCUE is called
C     only if rounding errors have reduced by at least a factor of two the
C     denominator of the formula for updating the H matrix. It provides a
C     useful safeguard, but is not invoked in most applications of BOBYQA.
C
  190 NFSAV=NF
      KBASE=KOPT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29
C See the comments above line number 60. 
      DO I = 1, N
          IF (GOPT(I) /= GOPT(I)) THEN
              INFO = -3
              GOTO 720
          END IF 
      END DO
      DO I = 1, NH
          IF (HQ(I) /= HQ(I)) THEN
              INFO = -3
              GOTO 720
          END IF 
      END DO 
      DO I = 1, NPT
          IF (PQ(I) /= PQ(I)) THEN
              INFO = -3
              GOTO 720
          END IF
      END DO
      DO J = 1, N
          DO I = 1, NDIM
              IF (BMAT(I,J) /= BMAT(I,J)) THEN
                  INFO = -3
                  GOTO 720
              END IF
          END DO
      END DO
      DO J = 1, NPTM
          DO I = 1, NPT
              IF (ZMAT(I,J) /= ZMAT(I,J)) THEN
                  INFO = -3
                  GOTO 720
              END IF
          END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,
     1  XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  VLAG,W,W(N+NP),W(NDIM+NP))
     2  VLAG,W,W(N+NP),W(NDIM+NP),F,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     XOPT is updated now in case the branch below to label 720 is taken.
C     Any updating of GOPT occurs after the branch below to label 20, which
C     leads to a trust region iteration as does the branch to label 60.
C
      XOPTSQ=ZERO
      IF (KOPT /= KBASE) THEN
          DO I=1,N
              XOPT(I)=XPT(KOPT,I)
              XOPTSQ=XOPTSQ+XOPT(I)**2
          END DO
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom/Zaikun (on 04-06-2019/07-06-2019):
      IF (F /= F .OR. F > ALMOST_INFINITY) THEN
          INFO=-2
          GOTO 720
      END IF
C     By Tom (on 04-06-2019):
C     If F reached the target function, RESCUE will stop and BOBYQB
C     should stop here.
      IF (F <= FTARGET) THEN
          INFO=1
          GOTO 736
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (NF < 0) THEN
          NF=MAXFUN
          IF (IPRINT > 0) PRINT 390
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 720
      END IF
      NRESC=NF
      IF (NFSAV < NF) THEN
          NFSAV=NF
          GOTO 20
      END IF
      IF (NTRITS > 0) GOTO 60
C
C     Pick two alternative vectors of variables, relative to XBASE, that
C     are suitable as new positions of the KNEW-th interpolation point.
C     Firstly, XNEW is set to the point on a line through XOPT and another
C     interpolation point that minimizes the predicted value of the next
C     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
C     and SU bounds. Secondly, XALT is set to the best feasible point on
C     a constrained version of the Cauchy step of the KNEW-th Lagrange
C     function, the corresponding value of the square of this function
C     being returned in CAUCHY. The choice between these alternatives is
C     going to be made when the denominator is calculated.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Zaikun 23-07-2019:
C  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
C     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
C
C  Although very rare, NaN can sometimes occur in BMAT or ZMAT. If it
C  happens, we terminate the code. See the comments above line number 60.
C  Indeed, if ALTMOV is called with such matrices, then altmov.f will 
C  encounter a memory error at lines 173--174. This is because the first 
C  value of PREDSQ in ALTOMOV (see line 159 of altmov.f) will be NaN, line 
C  164 will not be reached, and hence no value will be assigned to IBDSAV.
C  
C  Such an error was observed when BOBYQA was (mistakenly) tested on CUTEst 
C  problem CONCON. CONCON is a nonlinearly constrained problem with
C  bounds. By mistake, BOBYQA was called to solve this problem,
C  neglecting all the constraints other than bounds. With only the bound
C  constraints, the objective function turned to be unbounded from
C  below, which led to abnormal values in BMAT (indeed, BETA defined in
C  lines 366--389 took NaN/infinite values).
C
  210 DO J = 1,N
          DO I = 1,NDIM
              IF (BMAT(I,J) /= BMAT(I,J)) THEN
                  INFO = -3
                  GOTO 720
              END IF
          END DO
      END DO
      DO J = 1,NPTM
          DO I = 1,NPT
              IF (ZMAT(I,J) /= ZMAT(I,J)) THEN
                  INFO = -3
                  GOTO 720
              END IF
          END DO
      END DO
      CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I=1,N
          D(I)=XNEW(I)-XOPT(I)
      END DO
C
C     Calculate VLAG and BETA for the current choice of D. The scalar
C     product of D with XPT(K,.) is going to be held in W(NPT+K) for
C     use when VQUAD is calculated.
C
  230 DO K=1,NPT
          SUMA=ZERO
          SUMB=ZERO
          SUM=ZERO
          DO J=1,N
              SUMA=SUMA+XPT(K,J)*D(J)
              SUMB=SUMB+XPT(K,J)*XOPT(J)
              SUM=SUM+BMAT(K,J)*D(J)
          END DO
          W(K)=SUMA*(HALF*SUMA+SUMB)
          VLAG(K)=SUM
          W(NPT+K)=SUMA
      END DO
      BETA=ZERO
      DO JJ=1,NPTM
          SUM=ZERO
          DO K=1,NPT
              SUM=SUM+ZMAT(K,JJ)*W(K)
          END DO
          BETA=BETA-SUM*SUM
          DO K=1,NPT
              VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
          END DO
      END DO
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      DO J=1,N
          DSQ=DSQ+D(J)**2
          SUM=ZERO
          DO K=1,NPT
              SUM=SUM+W(K)*BMAT(K,J)
          END DO
          BSUM=BSUM+SUM*D(J)
          JP=NPT+J
          DO I=1,N
              SUM=SUM+BMAT(JP,I)*D(I)
          END DO
          VLAG(JP)=SUM
          BSUM=BSUM+SUM*D(J)
          DX=DX+D(J)*XOPT(J)
      END DO
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If NTRITS is zero, the denominator may be increased by replacing
C     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
C     rounding errors have damaged the chosen denominator.
C
      IF (NTRITS == 0) THEN
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          IF (DENOM < CAUCHY .AND. CAUCHY > ZERO) THEN
              DO I=1,N
                  XNEW(I)=XALT(I)
                  D(I)=XNEW(I)-XOPT(I)
              END DO
              CAUCHY=ZERO
              GO TO 230
          END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
          IF (.NOT. (DENOM > HALF*VLAG(KNEW)**2)) THEN
!111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (NF > NRESC) GOTO 190
              IF (IPRINT > 0) PRINT 320
  320         FORMAT (/5X,'Return from BOBYQA because of much',
     1          ' cancellation in a denominator.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              INFO=4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GOTO 720
          END IF
C
C     Alternatively, if NTRITS is positive, then set KNEW to the index of
C     the next interpolation point to be deleted to make room for a trust
C     region step. Again RESCUE may be called if rounding errors have damaged
C     the chosen denominator, which is the reason for attempting to select
C     KNEW before calculating the next value of the objective function.
C
      ELSE
          DELSQ=DELTA*DELTA
          SCADEN=ZERO
          BIGLSQ=ZERO
          KNEW=0
          DO K=1,NPT
              IF (K == KOPT) CYCLE 
              HDIAG=ZERO
              DO JJ=1,NPTM
                  HDIAG=HDIAG+ZMAT(K,JJ)**2
              END DO
              DEN=BETA*HDIAG+VLAG(K)**2
              DISTSQ=ZERO
              DO J=1,N
                  DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
              END DO
              TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
              IF (TEMP*DEN > SCADEN) THEN
                  SCADEN=TEMP*DEN
                  KNEW=K
                  DENOM=DEN
              END IF
              BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
          END DO
          IF (SCADEN <= HALF*BIGLSQ) THEN
              IF (NF > NRESC) GOTO 190
              IF (IPRINT > 0) PRINT 320
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              INFO=4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GOTO 720
          END IF
      END IF
C
C     Put the variables for the next calculation of the objective function
C       in XNEW, with any adjustments for the bounds.
C
C
C     Calculate the value of the objective function at XBASE+XNEW, unless
C       the limit on the number of calculations of F has been reached.
C
  360 DO I=1,N
          X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
          IF (XNEW(I) == SL(I)) X(I)=XL(I)
          IF (XNEW(I) == SU(I)) X(I)=XU(I)
      END DO
      IF (NF >= MAXFUN) THEN
          IF (IPRINT > 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been',
     1      ' called MAXFUN times.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 720
      END IF
      NF=NF+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO I=1,N
          IF (X(I) /= X(I)) THEN
              F=X(I) ! Set F to NaN
              IF (NF == 1) THEN
                  FOPT=F
                  XOPT(1:N)=ZERO
              END IF
              INFO=-1
              GOTO 720
          END IF
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL CALFUN (N,X,F)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom (on 04-06-2019):
      IF (F /= F .OR. F > ALMOST_INFINITY) THEN
          IF (NF == 1) THEN
              FOPT=F
              XOPT(1:N)=ZERO
          END IF
          INFO=-2
          GOTO 720
      END IF
C     By Tom (on 04-06-2019):
C     If F achieves the function value, the algorithm exits.
      IF (F <= FTARGET) THEN
          INFO=1
          GOTO 736
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IPRINT == 3) THEN
          PRINT 400, NF,F,(X(I),I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NTRITS == -1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: FSAVE is not needed any more. See line number 720.
C          FSAVE=F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GOTO 720
      END IF
C
C     Use the quadratic model to predict the change in F due to the step D,
C       and set DIFF to the error of this prediction.
C
      FOPT=FVAL(KOPT)
      VQUAD=ZERO
      IH=0
      DO J=1,N
          VQUAD=VQUAD+D(J)*GOPT(J)
          DO I=1,J
              IH=IH+1
              TEMP=D(I)*D(J)
              IF (I == J) TEMP=HALF*TEMP
              VQUAD=VQUAD+HQ(IH)*TEMP
          END DO
      END DO
      DO K=1,NPT
          VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
      END DO
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM > RHO) NFSAV=NF
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (NTRITS > 0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IF (VQUAD .GE. ZERO) THEN
          IF (.NOT. (VQUAD < ZERO)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (IPRINT > 0) PRINT 430
  430         FORMAT (/4X,'Return from BOBYQA because a trust',
     1          ' region step has failed to reduce Q.')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              INFO=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              GOTO 720
          END IF
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO <= TENTH) THEN
              DELTA=DMIN1(HALF*DELTA,DNORM)
          ELSE IF (RATIO <= 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,DNORM)
          ELSE
              DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
          END IF
          IF (DELTA <= 1.5D0*RHO) DELTA=RHO
C
C     Recalculate KNEW and DENOM if the new F is less than FOPT.
C
          IF (F < FOPT) THEN
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              DO K=1,NPT
                  HDIAG=ZERO
                  DO JJ=1,NPTM
                      HDIAG=HDIAG+ZMAT(K,JJ)**2
                  END DO
                  DEN=BETA*HDIAG+VLAG(K)**2
                  DISTSQ=ZERO
                  DO J=1,N
                      DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
                  END DO
                  TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
                  IF (TEMP*DEN > SCADEN) THEN
                      SCADEN=TEMP*DEN
                      KNEW=K
                      DENOM=DEN
                  END IF
                  BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
              END DO
              IF (SCADEN <= HALF*BIGLSQ) THEN
                  KNEW=KSAV
                  DENOM=DENSAV
              END IF
          END IF
      END IF
C
C     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
C     moved. Also update the second derivative terms of the model.
C
      CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      IH=0
      PQOLD=PQ(KNEW)
      PQ(KNEW)=ZERO
      DO I=1,N
          TEMP=PQOLD*XPT(KNEW,I)
          DO J=1,I
              IH=IH+1
              HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
          END DO
      END DO
      DO JJ=1,NPTM
          TEMP=DIFF*ZMAT(KNEW,JJ)
          DO K=1,NPT
              PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
          END DO
      END DO
C
C     Include the new interpolation point, and make the changes to GOPT at
C     the old XOPT that are caused by the updating of the quadratic model.
C
      FVAL(KNEW)=F
      DO I=1,N
          XPT(KNEW,I)=XNEW(I)
          W(I)=BMAT(KNEW,I)
      END DO
      DO K=1,NPT
          SUMA=ZERO
          DO JJ=1,NPTM
              SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
          END DO
          SUMB=ZERO
          DO J=1,N
              SUMB=SUMB+XPT(K,J)*XOPT(J)
          END DO
          TEMP=SUMA*SUMB
          DO I=1,N
              W(I)=W(I)+TEMP*XPT(K,I)
          END DO
      END DO
      DO I=1,N
          GOPT(I)=GOPT(I)+DIFF*W(I)
      END DO
C
C     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
C
      IF (F < FOPT) THEN
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          DO J=1,N
              XOPT(J)=XNEW(J)
              XOPTSQ=XOPTSQ+XOPT(J)**2
              DO I=1,J
                  IH=IH+1
                  IF (I < J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
                  GOPT(I)=GOPT(I)+HQ(IH)*D(J)
              END DO
          END DO
          DO K=1,NPT
              TEMP=ZERO
              DO J=1,N
                  TEMP=TEMP+XPT(K,J)*D(J)
              END DO
              TEMP=PQ(K)*TEMP
              DO I=1,N
                  GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
              END DO
          END DO
      END IF
C
C     Calculate the parameters of the least Frobenius norm interpolant to
C     the current data, the gradient of this interpolant at XOPT being put
C     into VLAG(NPT+I), I=1,2,...,N.
C
      IF (NTRITS > 0) THEN
          DO K=1,NPT
              VLAG(K)=FVAL(K)-FVAL(KOPT)
              W(K)=ZERO
          END DO
          DO J=1,NPTM
              SUM=ZERO
              DO K=1,NPT
                  SUM=SUM+ZMAT(K,J)*VLAG(K)
              END DO
              DO K=1,NPT
                  W(K)=W(K)+SUM*ZMAT(K,J)
              END DO
          END DO
          DO K=1,NPT
              SUM=ZERO
              DO J=1,N
                  SUM=SUM+XPT(K,J)*XOPT(J)
              END DO
              W(K+NPT)=W(K)
              W(K)=SUM*W(K)
          END DO
          GQSQ=ZERO
          GISQ=ZERO
          DO I=1,N
              SUM=ZERO
              DO K=1,NPT
                  SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
              END DO
              IF (XOPT(I) == SL(I)) THEN
                  GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
                  GISQ=GISQ+DMIN1(ZERO,SUM)**2
              ELSE IF (XOPT(I) == SU(I)) THEN
                  GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
                  GISQ=GISQ+DMAX1(ZERO,SUM)**2
              ELSE
                  GQSQ=GQSQ+GOPT(I)**2
                  GISQ=GISQ+SUM*SUM
              END IF
              VLAG(NPT+I)=SUM
          END DO
C
C     Test whether to replace the new quadratic model by the least Frobenius
C     norm interpolant, making the replacement if the test is satisfied.
C
          ITEST=ITEST+1
          IF (GQSQ < TEN*GISQ) ITEST=0
          IF (ITEST >= 3) THEN
              DO I=1,MAX0(NPT,NH)
                  IF (I <= N) GOPT(I)=VLAG(NPT+I)
                  IF (I <= NPT) PQ(I)=W(NPT+I)
                  IF (I <= NH) HQ(I)=ZERO
                  ITEST=0
              END DO
          END IF
      END IF
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case NTRITS=0 occurs
C     when the new interpolation point was reached by an alternative step.
C
      IF (NTRITS == 0) GOTO 60
      IF (F <= FOPT+TENTH*VQUAD) GOTO 60
C
C     Alternatively, find out if the interpolation points are close enough
C       to the best point so far.
C
      DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
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
C     If KNEW is positive, then ALTMOV finds alternative new positions for
C     the KNEW-th interpolation point within distance ADELT of XOPT. It is
C     reached via label 90. Otherwise, there is a branch to label 60 for
C     another trust region iteration, unless the calculations with the
C     current RHO are complete.
C
      IF (KNEW > 0) THEN
          DIST=DSQRT(DISTSQ)
          IF (NTRITS == -1) THEN
              DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
              IF (DELTA <= 1.5D0*RHO) DELTA=RHO
          END IF
          NTRITS=0
          ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          GOTO 90
      END IF
      IF (NTRITS == -1) GOTO 680
      IF (RATIO > ZERO) GOTO 60
      IF (DMAX1(DELTA,DNORM) > RHO) GOTO 60
C
C     The calculations with the current value of RHO are complete. Pick the
C       next values of RHO and DELTA.
C
  680 IF (RHO > RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO <= 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO <= 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT >= 2) THEN
              IF (IPRINT >= 3) PRINT 690
  690         FORMAT (5X)
              PRINT 700, RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
  710         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          NTRITS=0
          NFSAV=NF
          GOTO 60
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ELSE
          INFO=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C       it is too short to have been tried before.
C
      IF (NTRITS == -1) GOTO 360
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  720 IF (FVAL(KOPT) .LE. FSAVE) THEN
C  Why update X only when FVAL(KOPT) .LE. FSAVE? This seems INCORRECT, 
C  because it may lead to a return with F and X that are not the best
C  available. 
  720 IF (FVAL(KOPT) <= F .OR. F /= F) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DO I=1,N
              X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
              IF (XOPT(I) == SL(I)) X(I)=XL(I)
              IF (XOPT(I) == SU(I)) X(I)=XU(I)
          END DO
          F=FVAL(KOPT)
      END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IF (IPRINT .GE. 1) THEN
  736 IF (IPRINT >= 1) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X,
     1      'Number of function values =',I6)
          PRINT 710, F,(X(I),I=1,N)
      END IF
      RETURN
      END
