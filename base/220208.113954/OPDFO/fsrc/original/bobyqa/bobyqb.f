      SUBROUTINE BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
     2  SL,SU,XNEW,XALT,D,VLAG,W)
      IMPLICIT REAL*8 (A-H,O-Z)
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
C
C     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
C     BMAT and ZMAT for the first iteration, with the corresponding values of
C     of NF and KOPT, which are the number of calls of CALFUN so far and the
C     index of the interpolation point at the trust region centre. Then the
C     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
C     less than NPT. GOPT will be updated if KOPT is different from KBASE.
C
      CALL PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
      XOPTSQ=ZERO
      DO 10 I=1,N
      XOPT(I)=XPT(KOPT,I)
   10 XOPTSQ=XOPTSQ+XOPT(I)**2
      FSAVE=FVAL(1)
      IF (NF .LT. NPT) THEN
          IF (IPRINT .GT. 0) PRINT 390
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
   20 IF (KOPT .NE. KBASE) THEN
          IH=0
          DO 30 J=1,N
          DO 30 I=1,J
          IH=IH+1
          IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
   30     GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
          IF (NF .GT. NPT) THEN
              DO 50 K=1,NPT
              TEMP=ZERO
              DO 40 J=1,N
   40         TEMP=TEMP+XPT(K,J)*XOPT(J)
              TEMP=PQ(K)*TEMP
              DO 50 I=1,N
   50         GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
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
   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
     1  W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          IF (NF .LE. NFSAV+2) GOTO 650
C
C     The following choice between labels 650 and 680 depends on whether or
C     not our work with the current RHO seems to be complete. Either RHO is
C     decreased or termination occurs if the errors in the quadratic model at
C     the last three interpolation points compare favourably with predictions
C     of likely improvements to the model within distance HALF*RHO of XOPT.
C
          ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125D0*RHO*RHO
          IF (CRVMIN .GT. ZERO .AND. ERRBIG .GT. FRHOSQ*CRVMIN)
     1       GOTO 650
          BDTOL=ERRBIG/RHO
          DO 80 J=1,N
          BDTEST=BDTOL
          IF (XNEW(J) .EQ. SL(J)) BDTEST=W(J)
          IF (XNEW(J) .EQ. SU(J)) BDTEST=-W(J)
          IF (BDTEST .LT. BDTOL) THEN
              CURV=HQ((J+J*J)/2)
              DO 70 K=1,NPT
   70         CURV=CURV+PQ(K)*XPT(K,J)**2
              BDTEST=BDTEST+HALF*CURV*RHO
              IF (BDTEST .LT. BDTOL) GOTO 650
          END IF
   80     CONTINUE
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
   90 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          FRACSQ=0.25D0*XOPTSQ
          SUMPQ=ZERO
          DO 110 K=1,NPT
          SUMPQ=SUMPQ+PQ(K)
          SUM=-HALF*XOPTSQ
          DO 100 I=1,N
  100     SUM=SUM+XPT(K,I)*XOPT(I)
          W(NPT+K)=SUM
          TEMP=FRACSQ-HALF*SUM
          DO 110 I=1,N
          W(I)=BMAT(K,I)
          VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
          IP=NPT+I
          DO 110 J=1,I
  110     BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 150 JJ=1,NPTM
          SUMZ=ZERO
          SUMW=ZERO
          DO 120 K=1,NPT
          SUMZ=SUMZ+ZMAT(K,JJ)
          VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
  120     SUMW=SUMW+VLAG(K)
          DO 140 J=1,N
          SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
          DO 130 K=1,NPT
  130     SUM=SUM+VLAG(K)*XPT(K,J)
          W(J)=SUM
          DO 140 K=1,NPT
  140     BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
          DO 150 I=1,N
          IP=I+NPT
          TEMP=W(I)
          DO 150 J=1,I
  150     BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
C
C     The following instructions complete the shift, including the changes
C     to the second derivative parameters of the quadratic model.
C
          IH=0
          DO 170 J=1,N
          W(J)=-HALF*SUMPQ*XOPT(J)
          DO 160 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  160     XPT(K,J)=XPT(K,J)-XOPT(J)
          DO 170 I=1,J
          IH=IH+1
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  170     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 180 I=1,N
          XBASE(I)=XBASE(I)+XOPT(I)
          XNEW(I)=XNEW(I)-XOPT(I)
          SL(I)=SL(I)-XOPT(I)
          SU(I)=SU(I)-XOPT(I)
  180     XOPT(I)=ZERO
          XOPTSQ=ZERO
      END IF
      IF (NTRITS .EQ. 0) GOTO 210
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
      CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,
     1  XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,
     2  VLAG,W,W(N+NP),W(NDIM+NP))
C
C     XOPT is updated now in case the branch below to label 720 is taken.
C     Any updating of GOPT occurs after the branch below to label 20, which
C     leads to a trust region iteration as does the branch to label 60.
C
      XOPTSQ=ZERO
      IF (KOPT .NE. KBASE) THEN
          DO 200 I=1,N
          XOPT(I)=XPT(KOPT,I)
  200     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      IF (NF .LT. 0) THEN
          NF=MAXFUN
          IF (IPRINT .GT. 0) PRINT 390
          GOTO 720
      END IF
      NRESC=NF
      IF (NFSAV .LT. NF) THEN
          NFSAV=NF
          GOTO 20
      END IF
      IF (NTRITS .GT. 0) GOTO 60
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
  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
      DO 220 I=1,N
  220 D(I)=XNEW(I)-XOPT(I)
C
C     Calculate VLAG and BETA for the current choice of D. The scalar
C     product of D with XPT(K,.) is going to be held in W(NPT+K) for
C     use when VQUAD is calculated.
C
  230 DO 250 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 240 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  240 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
      VLAG(K)=SUM
  250 W(NPT+K)=SUMA
      BETA=ZERO
      DO 270 JJ=1,NPTM
      SUM=ZERO
      DO 260 K=1,NPT
  260 SUM=SUM+ZMAT(K,JJ)*W(K)
      BETA=BETA-SUM*SUM
      DO 270 K=1,NPT
  270 VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      DO 300 J=1,N
      DSQ=DSQ+D(J)**2
      SUM=ZERO
      DO 280 K=1,NPT
  280 SUM=SUM+W(K)*BMAT(K,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 290 I=1,N
  290 SUM=SUM+BMAT(JP,I)*D(I)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  300 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If NTRITS is zero, the denominator may be increased by replacing
C     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
C     rounding errors have damaged the chosen denominator.
C
      IF (NTRITS .EQ. 0) THEN
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          IF (DENOM .LT. CAUCHY .AND. CAUCHY .GT. ZERO) THEN
              DO 310 I=1,N
              XNEW(I)=XALT(I)
  310         D(I)=XNEW(I)-XOPT(I)
              CAUCHY=ZERO
              GO TO 230
          END IF
          IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
              IF (NF .GT. NRESC) GOTO 190
              IF (IPRINT .GT. 0) PRINT 320
  320         FORMAT (/5X,'Return from BOBYQA because of much',
     1          ' cancellation in a denominator.')
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
          DO 350 K=1,NPT
          IF (K .EQ. KOPT) GOTO 350
          HDIAG=ZERO
          DO 330 JJ=1,NPTM
  330     HDIAG=HDIAG+ZMAT(K,JJ)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          DISTSQ=ZERO
          DO 340 J=1,N
  340     DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
          TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
          IF (TEMP*DEN .GT. SCADEN) THEN
              SCADEN=TEMP*DEN
              KNEW=K
              DENOM=DEN
          END IF
          BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
  350     CONTINUE
          IF (SCADEN .LE. HALF*BIGLSQ) THEN
              IF (NF .GT. NRESC) GOTO 190
              IF (IPRINT .GT. 0) PRINT 320
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
  360 DO 380 I=1,N
      X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
      IF (XNEW(I) .EQ. SL(I)) X(I)=XL(I)
      IF (XNEW(I) .EQ. SU(I)) X(I)=XU(I)
  380 CONTINUE
      IF (NF .GE. MAXFUN) THEN
          IF (IPRINT .GT. 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 720
      END IF
      NF=NF+1
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 400, NF,F,(X(I),I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NTRITS .EQ. -1) THEN
          FSAVE=F
          GOTO 720
      END IF
C
C     Use the quadratic model to predict the change in F due to the step D,
C       and set DIFF to the error of this prediction.
C
      FOPT=FVAL(KOPT)
      VQUAD=ZERO
      IH=0
      DO 410 J=1,N
      VQUAD=VQUAD+D(J)*GOPT(J)
      DO 410 I=1,J
      IH=IH+1
      TEMP=D(I)*D(J)
      IF (I .EQ. J) TEMP=HALF*TEMP
  410 VQUAD=VQUAD+HQ(IH)*TEMP
      DO 420 K=1,NPT
  420 VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (NTRITS .GT. 0) THEN
          IF (VQUAD .GE. ZERO) THEN
              IF (IPRINT .GT. 0) PRINT 430
  430         FORMAT (/4X,'Return from BOBYQA because a trust',
     1          ' region step has failed to reduce Q.')
              GOTO 720
          END IF
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO .LE. TENTH) THEN
              DELTA=DMIN1(HALF*DELTA,DNORM)
          ELSE IF (RATIO .LE. 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,DNORM)
          ELSE
              DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
          END IF
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Recalculate KNEW and DENOM if the new F is less than FOPT.
C
          IF (F .LT. FOPT) THEN
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              DO 460 K=1,NPT
              HDIAG=ZERO
              DO 440 JJ=1,NPTM
  440         HDIAG=HDIAG+ZMAT(K,JJ)**2
              DEN=BETA*HDIAG+VLAG(K)**2
              DISTSQ=ZERO
              DO 450 J=1,N
  450         DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
              TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
              IF (TEMP*DEN .GT. SCADEN) THEN
                  SCADEN=TEMP*DEN
                  KNEW=K
                  DENOM=DEN
              END IF
  460         BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
              IF (SCADEN .LE. HALF*BIGLSQ) THEN
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
      DO 470 I=1,N
      TEMP=PQOLD*XPT(KNEW,I)
      DO 470 J=1,I
      IH=IH+1
  470 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      DO 480 JJ=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,JJ)
      DO 480 K=1,NPT
  480 PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
C
C     Include the new interpolation point, and make the changes to GOPT at
C     the old XOPT that are caused by the updating of the quadratic model.
C
      FVAL(KNEW)=F
      DO 490 I=1,N
      XPT(KNEW,I)=XNEW(I)
  490 W(I)=BMAT(KNEW,I)
      DO 520 K=1,NPT
      SUMA=ZERO
      DO 500 JJ=1,NPTM
  500 SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
      SUMB=ZERO
      DO 510 J=1,N
  510 SUMB=SUMB+XPT(K,J)*XOPT(J)
      TEMP=SUMA*SUMB
      DO 520 I=1,N
  520 W(I)=W(I)+TEMP*XPT(K,I)
      DO 530 I=1,N
  530 GOPT(I)=GOPT(I)+DIFF*W(I)
C
C     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
C
      IF (F .LT. FOPT) THEN
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          DO 540 J=1,N
          XOPT(J)=XNEW(J)
          XOPTSQ=XOPTSQ+XOPT(J)**2
          DO 540 I=1,J
          IH=IH+1
          IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
  540     GOPT(I)=GOPT(I)+HQ(IH)*D(J)
          DO 560 K=1,NPT
          TEMP=ZERO
          DO 550 J=1,N
  550     TEMP=TEMP+XPT(K,J)*D(J)
          TEMP=PQ(K)*TEMP
          DO 560 I=1,N
  560     GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
      END IF
C
C     Calculate the parameters of the least Frobenius norm interpolant to
C     the current data, the gradient of this interpolant at XOPT being put
C     into VLAG(NPT+I), I=1,2,...,N.
C
      IF (NTRITS .GT. 0) THEN
          DO 570 K=1,NPT
          VLAG(K)=FVAL(K)-FVAL(KOPT)
  570     W(K)=ZERO
          DO 590 J=1,NPTM
          SUM=ZERO
          DO 580 K=1,NPT
  580     SUM=SUM+ZMAT(K,J)*VLAG(K)
          DO 590 K=1,NPT
  590     W(K)=W(K)+SUM*ZMAT(K,J)
          DO 610 K=1,NPT
          SUM=ZERO
          DO 600 J=1,N
  600     SUM=SUM+XPT(K,J)*XOPT(J)
          W(K+NPT)=W(K)
  610     W(K)=SUM*W(K)
          GQSQ=ZERO
          GISQ=ZERO
          DO 630 I=1,N
          SUM=ZERO
          DO 620 K=1,NPT
  620     SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
          IF (XOPT(I) .EQ. SL(I)) THEN
              GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
              GISQ=GISQ+DMIN1(ZERO,SUM)**2
          ELSE IF (XOPT(I) .EQ. SU(I)) THEN
              GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
              GISQ=GISQ+DMAX1(ZERO,SUM)**2
          ELSE
              GQSQ=GQSQ+GOPT(I)**2
              GISQ=GISQ+SUM*SUM
          END IF
  630     VLAG(NPT+I)=SUM
C
C     Test whether to replace the new quadratic model by the least Frobenius
C     norm interpolant, making the replacement if the test is satisfied.
C
          ITEST=ITEST+1
          IF (GQSQ .LT. TEN*GISQ) ITEST=0
          IF (ITEST .GE. 3) THEN
              DO 640 I=1,MAX0(NPT,NH)
              IF (I .LE. N) GOPT(I)=VLAG(NPT+I)
              IF (I .LE. NPT) PQ(I)=W(NPT+I)
              IF (I .LE. NH) HQ(I)=ZERO
              ITEST=0
  640         CONTINUE
          END IF
      END IF
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case NTRITS=0 occurs
C     when the new interpolation point was reached by an alternative step.
C
      IF (NTRITS .EQ. 0) GOTO 60
      IF (F .LE. FOPT+TENTH*VQUAD) GOTO 60
C
C     Alternatively, find out if the interpolation points are close enough
C       to the best point so far.
C
      DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
      DO 670 K=1,NPT
      SUM=ZERO
      DO 660 J=1,N
  660 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  670 CONTINUE
C
C     If KNEW is positive, then ALTMOV finds alternative new positions for
C     the KNEW-th interpolation point within distance ADELT of XOPT. It is
C     reached via label 90. Otherwise, there is a branch to label 60 for
C     another trust region iteration, unless the calculations with the
C     current RHO are complete.
C
      IF (KNEW .GT. 0) THEN
          DIST=DSQRT(DISTSQ)
          IF (NTRITS .EQ. -1) THEN
              DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
              IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          END IF
          NTRITS=0
          ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          GOTO 90
      END IF
      IF (NTRITS .EQ. -1) GOTO 680
      IF (RATIO .GT. ZERO) GOTO 60
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 60
C
C     The calculations with the current value of RHO are complete. Pick the
C       next values of RHO and DELTA.
C
  680 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 690
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
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C       it is too short to have been tried before.
C
      IF (NTRITS .EQ. -1) GOTO 360
  720 IF (FVAL(KOPT) .LE. FSAVE) THEN
          DO 730 I=1,N
          X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
          IF (XOPT(I) .EQ. SL(I)) X(I)=XL(I)
          IF (XOPT(I) .EQ. SU(I)) X(I)=XU(I)
  730     CONTINUE
          F=FVAL(KOPT)
      END IF
      IF (IPRINT .GE. 1) THEN
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X,
     1      'Number of function values =',I6)
          PRINT 710, F,(X(I),I=1,N)
      END IF
      RETURN
      END
