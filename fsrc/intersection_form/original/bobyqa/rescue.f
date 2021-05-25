      SUBROUTINE RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,
     2  KOPT,VLAG,PTSAUX,PTSID,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),
     1  GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),
     2  VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
C
C     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
C       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
C       the corresponding arguments of BOBYQB on the entry to RESCUE.
C     NF is maintained as the number of calls of CALFUN so far, except that
C       NF is set to -1 if the value of MAXFUN prevents further progress.
C     KOPT is maintained so that FVAL(KOPT) is the least calculated function
C       value. Its correct value must be given on entry. It is updated if a
C       new least function value is found, but the corresponding changes to
C       XOPT and GOPT have to be made later by the calling program.
C     DELTA is the current trust region radius.
C     VLAG is a working space vector that will be used for the values of the
C       provisional Lagrange functions at each of the interpolation points.
C       They are part of a product that requires VLAG to be of length NDIM.
C     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
C       PTSAUX(2,J) specify the two positions of provisional interpolation
C       points when a nonzero step is taken along e_J (the J-th coordinate
C       direction) through XBASE+XOPT, as specified below. Usually these
C       steps have length DELTA, but other lengths are chosen if necessary
C       in order to satisfy the given bounds on the variables.
C     PTSID is also a working space array. It has NPT components that denote
C       provisional new positions of the original interpolation points, in
C       case changes are needed to restore the linear independence of the
C       interpolation conditions. The K-th point is a candidate for change
C       if and only if PTSID(K) is nonzero. In this case let p and q be the
C       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
C       and q are both positive, the step from XBASE+XOPT to the new K-th
C       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
C       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
C       p=0, respectively.
C     The first NDIM+NPT elements of the array W are used for working space. 
C     The final elements of BMAT and ZMAT are set in a well-conditioned way
C       to the values that are appropriate for the new interpolation points.
C     The elements of GOPT, HQ and PQ are also revised to the values that are
C       appropriate to the final quadratic model.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NP=N+1
      SFRAC=HALF/DFLOAT(NP)
      NPTM=NPT-NP
C
C     Shift the interpolation points so that XOPT becomes the origin, and set
C     the elements of ZMAT to zero. The value of SUMPQ is required in the
C     updating of HQ below. The squares of the distances from XOPT to the
C     other interpolation points are set at the end of W. Increments of WINC
C     may be added later to these squares to balance the consideration of
C     the choice of point that is going to become current.
C
      SUMPQ=ZERO
      WINC=ZERO
      DO 20 K=1,NPT
      DISTSQ=ZERO
      DO 10 J=1,N
      XPT(K,J)=XPT(K,J)-XOPT(J)
   10 DISTSQ=DISTSQ+XPT(K,J)**2
      SUMPQ=SUMPQ+PQ(K)
      W(NDIM+K)=DISTSQ
      WINC=DMAX1(WINC,DISTSQ)
      DO 20 J=1,NPTM
   20 ZMAT(K,J)=ZERO
C
C     Update HQ so that HQ and PQ define the second derivatives of the model
C     after XBASE has been shifted to the trust region centre.
C
      IH=0
      DO 40 J=1,N
      W(J)=HALF*SUMPQ*XOPT(J)
      DO 30 K=1,NPT
   30 W(J)=W(J)+PQ(K)*XPT(K,J)
      DO 40 I=1,J
      IH=IH+1
   40 HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
C
C     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
C     also set the elements of PTSAUX.
C
      DO 50 J=1,N
      XBASE(J)=XBASE(J)+XOPT(J)
      SL(J)=SL(J)-XOPT(J)
      SU(J)=SU(J)-XOPT(J)
      XOPT(J)=ZERO
      PTSAUX(1,J)=DMIN1(DELTA,SU(J))
      PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
      IF (PTSAUX(1,J)+PTSAUX(2,J) .LT. ZERO) THEN
          TEMP=PTSAUX(1,J)
          PTSAUX(1,J)=PTSAUX(2,J)
          PTSAUX(2,J)=TEMP
      END IF
      IF (DABS(PTSAUX(2,J)) .LT. HALF*DABS(PTSAUX(1,J))) THEN
          PTSAUX(2,J)=HALF*PTSAUX(1,J)
      END IF
      DO 50 I=1,NDIM
   50 BMAT(I,J)=ZERO
      FBASE=FVAL(KOPT)
C
C     Set the identifiers of the artificial interpolation points that are
C     along a coordinate direction from XOPT, and set the corresponding
C     nonzero elements of BMAT and ZMAT.
C
      PTSID(1)=SFRAC
      DO 60 J=1,N
      JP=J+1
      JPN=JP+N
      PTSID(JP)=DFLOAT(J)+SFRAC
      IF (JPN .LE. NPT) THEN
          PTSID(JPN)=DFLOAT(J)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
          BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
          BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
          BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
          ZMAT(1,J)=DSQRT(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
          ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
          ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
      ELSE
          BMAT(1,J)=-ONE/PTSAUX(1,J)
          BMAT(JP,J)=ONE/PTSAUX(1,J)
          BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
      END IF
   60 CONTINUE
C
C     Set any remaining identifiers with their nonzero elements of ZMAT.
C
      IF (NPT .GE. N+NP) THEN
          DO 70 K=2*NP,NPT
          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
          IP=K-NP-IW*N
          IQ=IP+IW
          IF (IQ .GT. N) IQ=IQ-N
          PTSID(K)=DFLOAT(IP)+DFLOAT(IQ)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
          ZMAT(1,K-NP)=TEMP
          ZMAT(IP+1,K-NP)=-TEMP
          ZMAT(IQ+1,K-NP)=-TEMP
   70     ZMAT(K,K-NP)=TEMP
      END IF
      NREM=NPT
      KOLD=1
      KNEW=KOPT
C
C     Reorder the provisional points in the way that exchanges PTSID(KOLD)
C     with PTSID(KNEW).
C
   80 DO 90 J=1,N
      TEMP=BMAT(KOLD,J)
      BMAT(KOLD,J)=BMAT(KNEW,J)
   90 BMAT(KNEW,J)=TEMP
      DO 100 J=1,NPTM
      TEMP=ZMAT(KOLD,J)
      ZMAT(KOLD,J)=ZMAT(KNEW,J)
  100 ZMAT(KNEW,J)=TEMP
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      IF (KNEW .NE. KOPT) THEN
          TEMP=VLAG(KOLD)
          VLAG(KOLD)=VLAG(KNEW)
          VLAG(KNEW)=TEMP
C
C     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
C     interpolation point can be changed from provisional to original. The
C     branch to label 350 occurs if all the original points are reinstated.
C     The nonnegative values of W(NDIM+K) are required in the search below.
C
          CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
          IF (NREM .EQ. 0) GOTO 350
          DO 110 K=1,NPT
  110     W(NDIM+K)=DABS(W(NDIM+K))
      END IF
C
C     Pick the index KNEW of an original interpolation point that has not
C     yet replaced one of the provisional interpolation points, giving
C     attention to the closeness to XOPT and to previous tries with KNEW.
C
  120 DSQMIN=ZERO
      DO 130 K=1,NPT
      IF (W(NDIM+K) .GT. ZERO) THEN
          IF (DSQMIN .EQ. ZERO .OR. W(NDIM+K) .LT. DSQMIN) THEN
              KNEW=K
              DSQMIN=W(NDIM+K)
          END IF
      END IF
  130 CONTINUE
      IF (DSQMIN .EQ. ZERO) GOTO 260
C
C     Form the W-vector of the chosen original interpolation point.
C
      DO 140 J=1,N
  140 W(NPT+J)=XPT(KNEW,J)
      DO 160 K=1,NPT
      SUM=ZERO
      IF (K .EQ. KOPT) THEN
          CONTINUE
      ELSE IF (PTSID(K) .EQ. ZERO) THEN
          DO 150 J=1,N
  150     SUM=SUM+W(NPT+J)*XPT(K,J)
      ELSE
          IP=PTSID(K)
          IF (IP .GT. 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IF (IQ .GT. 0) THEN
              IW=1
              IF (IP .EQ. 0) IW=2
              SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
          END IF
      END IF
  160 W(K)=HALF*SUM*SUM
C
C     Calculate VLAG and BETA for the required updating of the H matrix if
C     XPT(KNEW,.) is reinstated in the set of interpolation points.
C
      DO 180 K=1,NPT
      SUM=ZERO
      DO 170 J=1,N
  170 SUM=SUM+BMAT(K,J)*W(NPT+J)
  180 VLAG(K)=SUM
      BETA=ZERO
      DO 200 J=1,NPTM
      SUM=ZERO
      DO 190 K=1,NPT
  190 SUM=SUM+ZMAT(K,J)*W(K)
      BETA=BETA-SUM*SUM
      DO 200 K=1,NPT
  200 VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
      BSUM=ZERO
      DISTSQ=ZERO
      DO 230 J=1,N
      SUM=ZERO
      DO 210 K=1,NPT
  210 SUM=SUM+BMAT(K,J)*W(K)
      JP=J+NPT
      BSUM=BSUM+SUM*W(JP)
      DO 220 IP=NPT+1,NDIM
  220 SUM=SUM+BMAT(IP,J)*W(IP)
      BSUM=BSUM+SUM*W(JP)
      VLAG(JP)=SUM
  230 DISTSQ=DISTSQ+XPT(KNEW,J)**2
      BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     KOLD is set to the index of the provisional interpolation point that is
C     going to be deleted to make way for the KNEW-th original interpolation
C     point. The choice of KOLD is governed by the avoidance of a small value
C     of the denominator in the updating calculation of UPDATE.
C
      DENOM=ZERO
      VLMXSQ=ZERO
      DO 250 K=1,NPT
      IF (PTSID(K) .NE. ZERO) THEN
          HDIAG=ZERO
          DO 240 J=1,NPTM
  240     HDIAG=HDIAG+ZMAT(K,J)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          IF (DEN .GT. DENOM) THEN
              KOLD=K
              DENOM=DEN
          END IF
      END IF
  250 VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
      IF (DENOM .LE. 1.0D-2*VLMXSQ) THEN
          W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
          GOTO 120
      END IF
      GOTO 80
C
C     When label 260 is reached, all the final positions of the interpolation
C     points have been chosen although any changes have not been included yet
C     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
C     from the shift of XBASE, the updating of the quadratic model remains to
C     be done. The following cycle through the new interpolation points begins
C     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
C     except that a RETURN occurs if MAXFUN prohibits another value of F.
C
  260 DO 340 KPT=1,NPT
      IF (PTSID(KPT) .EQ. ZERO) GOTO 340
      IF (NF .GE. MAXFUN) THEN
          NF=-1
          GOTO 350
      END IF
      IH=0
      DO 270 J=1,N
      W(J)=XPT(KPT,J)
      XPT(KPT,J)=ZERO
      TEMP=PQ(KPT)*W(J)
      DO 270 I=1,J
      IH=IH+1
  270 HQ(IH)=HQ(IH)+TEMP*W(I)
      PQ(KPT)=ZERO
      IP=PTSID(KPT)
      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
      IF (IP .GT. 0) THEN
          XP=PTSAUX(1,IP)
          XPT(KPT,IP)=XP
      END IF
      IF (IQ .GT. 0) THEN
          XQ=PTSAUX(1,IQ)
          IF (IP .EQ. 0) XQ=PTSAUX(2,IQ)
          XPT(KPT,IQ)=XQ
      END IF
C
C     Set VQUAD to the value of the current model at the new point.
C
      VQUAD=FBASE
      IF (IP .GT. 0) THEN
          IHP=(IP+IP*IP)/2
          VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
      END IF
      IF (IQ .GT. 0) THEN
          IHQ=(IQ+IQ*IQ)/2
          VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
          IF (IP .GT. 0) THEN
              IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
              VQUAD=VQUAD+XP*XQ*HQ(IW)
          END IF
      END IF
      DO 280 K=1,NPT
      TEMP=ZERO
      IF (IP .GT. 0) TEMP=TEMP+XP*XPT(K,IP)
      IF (IQ .GT. 0) TEMP=TEMP+XQ*XPT(K,IQ)
  280 VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
C
C     Calculate F at the new interpolation point, and set DIFF to the factor
C     that is going to multiply the KPT-th Lagrange function when the model
C     is updated to provide interpolation to the new function value.
C
      DO 290 I=1,N
      W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
      IF (XPT(KPT,I) .EQ. SL(I)) W(I)=XL(I)
      IF (XPT(KPT,I) .EQ. SU(I)) W(I)=XU(I)
  290 CONTINUE
      NF=NF+1
      CALL CALFUN (N,W,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 300, NF,F,(W(I),I=1,N)
  300     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(KPT)=F
      IF (F .LT. FVAL(KOPT)) KOPT=KPT
      DIFF=F-VQUAD
C
C     Update the quadratic model. The RETURN from the subroutine occurs when
C     all the new interpolation points are included in the model.
C
      DO 310 I=1,N
  310 GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
      DO 330 K=1,NPT
      SUM=ZERO
      DO 320 J=1,NPTM
  320 SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
      TEMP=DIFF*SUM
      IF (PTSID(K) .EQ. ZERO) THEN
          PQ(K)=PQ(K)+TEMP
      ELSE
          IP=PTSID(K)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IHQ=(IQ*IQ+IQ)/2
          IF (IP .EQ. 0) THEN
              HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
          ELSE
              IHP=(IP*IP+IP)/2
              HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
              IF (IQ .GT. 0) THEN
                  HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                  IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                  HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
              END IF
          END IF
      END IF
  330 CONTINUE
      PTSID(KPT)=ZERO
  340 CONTINUE
  350 RETURN
      END
