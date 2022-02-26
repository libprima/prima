      SUBROUTINE LINCOB (N,NPT,M,AMAT,B,X,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XSAV,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
     2  STEP,SP,XNEW,IACT,RESCON,QFAC,RFAC,PQW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
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
     2  STEP,PQW,W)
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
      DO 30 I=1,N
   30 XOPTSQ=XOPTSQ+XOPT(I)**2
      IF (XOPTSQ .GE. 1.0D4*DELTA*DELTA) THEN
          QOPTSQ=0.25D0*XOPTSQ
          DO 50 K=1,NPT
          SUM=ZERO
          DO 40 I=1,N
   40     SUM=SUM+XPT(K,I)*XOPT(I)
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          SP(K)=ZERO
          DO 50 I=1,N
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          STEP(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+QOPTSQ*XOPT(I)
          IP=NPT+I
          DO 50 J=1,I
   50     BMAT(IP,J)=BMAT(IP,J)+STEP(I)*W(J)+W(I)*STEP(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 90 K=1,NPTM
          SUMZ=ZERO
          DO 60 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
   60     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 80 J=1,N
          SUM=QOPTSQ*SUMZ*XOPT(J)
          DO 70 I=1,NPT
   70     SUM=SUM+W(I)*XPT(I,J)
          STEP(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 80 I=1,NPT
   80     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 90 I=1,N
          IP=I+NPT
          TEMP=STEP(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 90 J=1,I
   90     BMAT(IP,J)=BMAT(IP,J)+TEMP*STEP(J)
C
C     Update the right hand sides of the constraints.
C
          IF (M .GT. 0) THEN
              DO 110 J=1,M
              TEMP=ZERO
              DO 100 I=1,N
  100         TEMP=TEMP+AMAT(I,J)*XOPT(I)
  110         B(J)=B(J)-TEMP
          END IF
C
C     The following instructions complete the shift of XBASE, including the
C       changes to the parameters of the quadratic model.
C
          IH=0
          DO 130 J=1,N
          W(J)=ZERO
          DO 120 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  120     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 130 I=1,J
          IH=IH+1
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  130     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 140 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
          XOPT(J)=ZERO
  140     XPT(KOPT,J)=ZERO
      END IF
C
C     In the case KNEW=0, generate the next trust region step by calling
C       TRSTEP, where SNORM is the current trust region radius initially.
C       The final value of SNORM is the length of the calculated step,
C       except that SNORM is zero on return if the projected gradient is
C       unsuitable for starting the conjugate gradient iterations.
C
      DELSAV=DELTA
      KSAVE=KNEW
      IF (KNEW .EQ. 0) THEN
          SNORM=DELTA
          DO 150 I=1,N
  150     XNEW(I)=GOPT(I)
          CALL TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
     1      QFAC,RFAC,SNORM,STEP,XNEW,W,W(M+1),PQW,PQW(NP),W(M+NP))
C
C     A trust region step is applied whenever its length, namely SNORM, is at
C       least HALF*DELTA. It is also applied if its length is at least 0.1999
C       times DELTA and if a line search of TRSTEP has caused a change to the
C       active set. Otherwise there is a branch below to label 530 or 560.
C
          TEMP=HALF*DELTA
          IF (XNEW(1) .GE. HALF) TEMP=0.1999D0*DELTA
          IF (SNORM .LE. TEMP) THEN
              DELTA=HALF*DELTA
              IF (DELTA .LE. 1.4D0*RHO) DELTA=RHO
              NVALA=NVALA+1
              NVALB=NVALB+1
              TEMP=SNORM/RHO
              IF (DELSAV .GT. RHO) TEMP=ONE
              IF (TEMP .GE. HALF) NVALA=ZERO
              IF (TEMP .GE. TENTH) NVALB=ZERO
              IF (DELSAV .GT. RHO) GOTO 530
              IF (NVALA .LT. 5 .AND. NVALB .LT. 3) GOTO 530
              IF (SNORM .GT. ZERO) KSAVE=-1
              GOTO 560
          END IF
          NVALA=ZERO
          NVALB=ZERO
C
C     Alternatively, KNEW is positive. Then the model step is calculated
C       within a trust region of radius DEL, after setting the gradient at
C       XBASE and the second derivative parameters of the KNEW-th Lagrange
C       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
C
      ELSE
          DEL=DMAX1(TENTH*DELTA,RHO)
          DO 160 I=1,N
  160     W(I)=BMAT(KNEW,I)
          DO 170 K=1,NPT
  170     PQW(K)=ZERO
          DO 180 J=1,NPTM
          TEMP=ZMAT(KNEW,J)
          IF (J .LT. IDZ) TEMP=-TEMP
          DO 180 K=1,NPT
  180     PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
          CALL QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,RESCON,
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
      DO 190 J=1,N
      VQUAD=VQUAD+STEP(J)*GOPT(J)
      DO 190 I=1,J
      IH=IH+1
      TEMP=STEP(I)*STEP(J)
      IF (I .EQ. J) TEMP=HALF*TEMP
  190 VQUAD=VQUAD+TEMP*HQ(IH)
      DO 210 K=1,NPT
      TEMP=ZERO
      DO 200 J=1,N
      TEMP=TEMP+XPT(K,J)*STEP(J)
  200 SP(NPT+K)=TEMP
  210 VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
      IF (KSAVE .EQ. 0 .AND. VQUAD .GE. ZERO) GOTO 530
C
C     Calculate the next value of the objective function. The difference
C       between the actual new value of F and the value predicted by the
C       model is recorded in DIFF.
C
  220 NF=NF+1
      IF (NF .GT. MAXFUN) THEN
          NF=NF-1
          IF (IPRINT .GT. 0) PRINT 230
  230     FORMAT (/4X,'Return from LINCOA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 600
      END IF
      XDIFF=ZERO
      DO 240 I=1,N
      XNEW(I)=XOPT(I)+STEP(I)
      X(I)=XBASE(I)+XNEW(I)
  240 XDIFF=XDIFF+(X(I)-XSAV(I))**2
      XDIFF=DSQRT(XDIFF)
      IF (KSAVE .EQ. -1) XDIFF=RHO
      IF (XDIFF .LE. TENTH*RHO .OR. XDIFF .GE. DELTA+DELTA) THEN
          IFEAS=0
          IF (IPRINT .GT. 0) PRINT 250
  250     FORMAT (/4X,'Return from LINCOA because rounding errors',
     1      ' prevent reasonable changes to X.')
          GOTO 600
      END IF
      IF (KSAVE .LE. 0) IFEAS=1
      F=DFLOAT(IFEAS)
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 260, NF,F,(X(I),I=1,N)
  260     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (KSAVE .EQ. -1) GOTO 600
      DIFF=F-FOPT-VQUAD
C
C     If X is feasible, then set DFFALT to the difference between the new
C       value of F and the value predicted by the alternative model.
C
      IF (IFEAS .EQ. 1 .AND. ITEST .LT. 3) THEN
          DO 270 K=1,NPT
          PQW(K)=ZERO
  270     W(K)=FVAL(K)-FVAL(KOPT)
          DO 290 J=1,NPTM
          SUM=ZERO
          DO 280 I=1,NPT
  280     SUM=SUM+W(I)*ZMAT(I,J)
          IF (J .LT. IDZ) SUM=-SUM
          DO 290 K=1,NPT
  290     PQW(K)=PQW(K)+SUM*ZMAT(K,J)
          VQALT=ZERO
          DO 310 K=1,NPT
          SUM=ZERO
          DO 300 J=1,N
  300     SUM=SUM+BMAT(K,J)*STEP(J)
          VQALT=VQALT+SUM*W(K)
  310     VQALT=VQALT+PQW(K)*SP(NPT+K)*(HALF*SP(NPT+K)+SP(K))
          DFFALT=F-FOPT-VQALT
      END IF
      IF (ITEST .EQ. 3) THEN
          DFFALT=DIFF
          ITEST=0
      END IF
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (KSAVE .EQ. 0) THEN
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO .LE. TENTH) THEN
              DELTA=HALF*DELTA
          ELSE IF (RATIO .LE. 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,SNORM)
          ELSE 
              TEMP=DSQRT(2.0D0)*DELTA
              DELTA=DMAX1(HALF*DELTA,SNORM+SNORM)
              DELTA=DMIN1(DELTA,TEMP)
          END IF
          IF (DELTA .LE. 1.4D0*RHO) DELTA=RHO
      END IF
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C       can be moved. If STEP is a trust region step, then KNEW is zero at
C       present, but a positive value is picked by subroutine UPDATE.
C
      CALL UPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,KOPT,
     1  KNEW,PQW,W)
      IF (KNEW .EQ. 0) THEN
          IF (IPRINT .GT. 0) PRINT 320
  320     FORMAT (/4X,'Return from LINCOA because the denominator'
     1      ' of the updating formula is zero.')
          GOTO 600
      END IF
C
C     If ITEST is increased to 3, then the next quadratic model is the
C       one whose second derivative matrix is least subject to the new
C       interpolation conditions. Otherwise the new model is constructed
C       by the symmetric Broyden method in the usual way.
C
      IF (IFEAS .EQ. 1) THEN
          ITEST=ITEST+1
          IF (DABS(DFFALT) .GE. TENTH*DABS(DIFF)) ITEST=0
      END IF
C
C     Update the second derivatives of the model by the symmetric Broyden
C       method, using PQW for the second derivative parameters of the new
C       KNEW-th Lagrange function. The contribution from the old parameter
C       PQ(KNEW) is included in the second derivative matrix HQ. W is used
C       later for the gradient of the new KNEW-th Lagrange function.       
C
      IF (ITEST .LT. 3) THEN
          DO 330 K=1,NPT
  330     PQW(K)=ZERO
          DO 350 J=1,NPTM
          TEMP=ZMAT(KNEW,J)
          IF (TEMP .NE. ZERO) THEN
              IF (J .LT. IDZ) TEMP=-TEMP
              DO 340 K=1,NPT
  340         PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
          END IF
  350     CONTINUE
          IH=0
          DO 360 I=1,N
          W(I)=BMAT(KNEW,I)
          TEMP=PQ(KNEW)*XPT(KNEW,I)
          DO 360 J=1,I
          IH=IH+1
  360     HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
          PQ(KNEW)=ZERO
          DO 370 K=1,NPT
  370     PQ(K)=PQ(K)+DIFF*PQW(K)
      END IF
C
C     Include the new interpolation point with the corresponding updates of
C       SP. Also make the changes of the symmetric Broyden method to GOPT at
C       the old XOPT if ITEST is less than 3.
C
      FVAL(KNEW)=F
      SP(KNEW)=SP(KOPT)+SP(NPT+KOPT)
      SSQ=ZERO
      DO 380 I=1,N
      XPT(KNEW,I)=XNEW(I)
  380 SSQ=SSQ+STEP(I)**2
      SP(NPT+KNEW)=SP(NPT+KOPT)+SSQ
      IF (ITEST .LT. 3) THEN
          DO 390 K=1,NPT
          TEMP=PQW(K)*SP(K)
          DO 390 I=1,N
  390     W(I)=W(I)+TEMP*XPT(K,I)
          DO 400 I=1,N
  400     GOPT(I)=GOPT(I)+DIFF*W(I)
      END IF
C
C     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
C       least calculated value so far with a feasible vector of variables.
C
      IF (F .LT. FOPT .AND. IFEAS .EQ. 1) THEN
          FOPT=F
          DO 410 J=1,N
          XSAV(J)=X(J)
  410     XOPT(J)=XNEW(J)
          KOPT=KNEW
          SNORM=DSQRT(SSQ)
          DO 430 J=1,M
          IF (RESCON(J) .GE. DELTA+SNORM) THEN
              RESCON(J)=SNORM-RESCON(J)
          ELSE
              RESCON(J)=RESCON(J)+SNORM
              IF (RESCON(J)+DELTA .GT. ZERO) THEN
                  TEMP=B(J)
                  DO 420 I=1,N
  420             TEMP=TEMP-XOPT(I)*AMAT(I,J)
                  TEMP=DMAX1(TEMP,ZERO)
                  IF (TEMP .GE. DELTA) TEMP=-TEMP
                  RESCON(J)=TEMP
              END IF
          END IF
  430     CONTINUE
          DO 440 K=1,NPT
  440     SP(K)=SP(K)+SP(NPT+K)
C
C     Also revise GOPT when symmetric Broyden updating is applied.
C
          IF (ITEST .LT. 3) THEN
              IH=0
              DO 450 J=1,N
              DO 450 I=1,J
              IH=IH+1
              IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*STEP(I)
  450         GOPT(I)=GOPT(I)+HQ(IH)*STEP(J)
              DO 460 K=1,NPT
              TEMP=PQ(K)*SP(NPT+K)
              DO 460 I=1,N
  460         GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
          END IF
      END IF
C
C     Replace the current model by the least Frobenius norm interpolant if
C       this interpolant gives substantial reductions in the predictions
C       of values of F at feasible points.
C
      IF (ITEST .EQ. 3) THEN
          DO 470 K=1,NPT
          PQ(K)=ZERO
  470     W(K)=FVAL(K)-FVAL(KOPT)
          DO 490 J=1,NPTM
          SUM=ZERO
          DO 480 I=1,NPT
  480     SUM=SUM+W(I)*ZMAT(I,J)
          IF (J .LT. IDZ) SUM=-SUM
          DO 490 K=1,NPT
  490     PQ(K)=PQ(K)+SUM*ZMAT(K,J)
          DO 500 J=1,N
          GOPT(J)=ZERO
          DO 500 I=1,NPT
  500     GOPT(J)=GOPT(J)+W(I)*BMAT(I,J)
          DO 510 K=1,NPT
          TEMP=PQ(K)*SP(K)
          DO 510 I=1,N
  510     GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
          DO 520 IH=1,NH
  520     HQ(IH)=ZERO
      END IF
C
C     If a trust region step has provided a sufficient decrease in F, then
C       branch for another trust region calculation. Every iteration that
C       takes a model step is followed by an attempt to take a trust region
C       step.
C
      KNEW=0
      IF (KSAVE .GT. 0) GOTO 20
      IF (RATIO .GE. TENTH) GOTO 20
C
C     Alternatively, find out if the interpolation points are close enough
C       to the best point so far.
C
  530 DISTSQ=DMAX1(DELTA*DELTA,4.0D0*RHO*RHO)
      DO 550 K=1,NPT
      SUM=ZERO
      DO 540 J=1,N
  540 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  550 CONTINUE
C
C     If KNEW is positive, then branch back for the next iteration, which
C       will generate a "model step". Otherwise, if the current iteration
C       has reduced F, or if DELTA was above its lower bound when the last
C       trust region step was calculated, then try a "trust region" step
C       instead.
C
      IF (KNEW .GT. 0) GOTO 20
      KNEW=0
      IF (FOPT .LT. FSAVE) GOTO 20
      IF (DELSAV .GT. RHO) GOTO 20
C
C     The calculations with the current value of RHO are complete.
C       Pick the next value of RHO.
C
  560 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          IF (RHO .GT. 250.0D0*RHOEND) THEN
              RHO=TENTH*RHO
          ELSE IF (RHO .LE. 16.0D0*RHOEND) THEN
              RHO=RHOEND
          ELSE
              RHO=DSQRT(RHO*RHOEND)
          END IF 
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 570
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
      IF (KSAVE .EQ. -1) GOTO 220
  600 IF (FOPT .LE. F .OR. IFEAS .EQ. 0) THEN
          DO 610 I=1,N
  610     X(I)=XSAV(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
          PRINT 620, NF
  620     FORMAT (/4X,'At the return from LINCOA',5X,
     1      'Number of function values =',I6)
          PRINT 590, F,(X(I),I=1,N)
      END IF
      W(1)=F
      W(2)=DFLOAT(NF)+HALF
      RETURN
      END
