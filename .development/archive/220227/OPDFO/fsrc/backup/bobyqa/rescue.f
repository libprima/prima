      SUBROUTINE RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  KOPT,VLAG,PTSAUX,PTSID,W)
     2  KOPT,VLAG,PTSAUX,PTSID,W,F,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ALMOST_INFINITY=HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      DO K=1,NPT
          DISTSQ=ZERO
          DO J=1,N
              XPT(K,J)=XPT(K,J)-XOPT(J)
              DISTSQ=DISTSQ+XPT(K,J)**2
          END DO
          SUMPQ=SUMPQ+PQ(K)
          W(NDIM+K)=DISTSQ
          WINC=DMAX1(WINC,DISTSQ)
          DO J=1,NPTM
              ZMAT(K,J)=ZERO
          END DO
      END DO
C
C     Update HQ so that HQ and PQ define the second derivatives of the model
C     after XBASE has been shifted to the trust region centre.
C
      IH=0
      DO J=1,N
          W(J)=HALF*SUMPQ*XOPT(J)
          DO K=1,NPT
              W(J)=W(J)+PQ(K)*XPT(K,J)
          END DO
          DO I=1,J
              IH=IH+1
              HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
          END DO
      END DO
C
C     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
C     also set the elements of PTSAUX.
C
      DO J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
          SL(J)=SL(J)-XOPT(J)
          SU(J)=SU(J)-XOPT(J)
          XOPT(J)=ZERO
          PTSAUX(1,J)=DMIN1(DELTA,SU(J))
          PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
          IF (PTSAUX(1,J)+PTSAUX(2,J) < ZERO) THEN
              TEMP=PTSAUX(1,J)
              PTSAUX(1,J)=PTSAUX(2,J)
              PTSAUX(2,J)=TEMP
          END IF
          IF (DABS(PTSAUX(2,J)) < HALF*DABS(PTSAUX(1,J))) THEN
              PTSAUX(2,J)=HALF*PTSAUX(1,J)
          END IF
          DO I=1,NDIM
              BMAT(I,J)=ZERO
          END DO
      END DO
      FBASE=FVAL(KOPT)
C
C     Set the identifiers of the artificial interpolation points that are
C     along a coordinate direction from XOPT, and set the corresponding
C     nonzero elements of BMAT and ZMAT.
C
      PTSID(1)=SFRAC
      DO J=1,N
          JP=J+1
          JPN=JP+N
          PTSID(JP)=DFLOAT(J)+SFRAC
          IF (JPN <= NPT) THEN
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
      END DO
C
C     Set any remaining identifiers with their nonzero elements of ZMAT.
C
      IF (NPT >= N+NP) THEN
          DO K=2*NP,NPT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
              IW=INT((DFLOAT(K-NP)-HALF)/DFLOAT(N))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IP=K-NP-IW*N
              IQ=IP+IW
              IF (IQ > N) IQ=IQ-N
              PTSID(K)=DFLOAT(IP)+DFLOAT(IQ)/DFLOAT(NP)+SFRAC
              TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
              ZMAT(1,K-NP)=TEMP
              ZMAT(IP+1,K-NP)=-TEMP
              ZMAT(IQ+1,K-NP)=-TEMP
              ZMAT(K,K-NP)=TEMP
          END DO
      END IF
      NREM=NPT
      KOLD=1
      KNEW=KOPT
C
C     Reorder the provisional points in the way that exchanges PTSID(KOLD)
C     with PTSID(KNEW).
C
   80 DO J=1,N
          TEMP=BMAT(KOLD,J)
          BMAT(KOLD,J)=BMAT(KNEW,J)
          BMAT(KNEW,J)=TEMP
      END DO
      DO J=1,NPTM
          TEMP=ZMAT(KOLD,J)
          ZMAT(KOLD,J)=ZMAT(KNEW,J)
          ZMAT(KNEW,J)=TEMP
      END DO
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      IF (KNEW /= KOPT) THEN
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
          IF (NREM == 0) GOTO 350
          DO K=1,NPT
              W(NDIM+K)=DABS(W(NDIM+K))
          END DO
      END IF
C
C     Pick the index KNEW of an original interpolation point that has not
C     yet replaced one of the provisional interpolation points, giving
C     attention to the closeness to XOPT and to previous tries with KNEW.
C
  120 DSQMIN=ZERO
      DO K=1,NPT
          IF (W(NDIM+K) > ZERO) THEN
              IF (DSQMIN == ZERO .OR. W(NDIM+K) < DSQMIN) THEN
                  KNEW=K
                  DSQMIN=W(NDIM+K)
              END IF
          END IF
      END DO
      IF (DSQMIN == ZERO) GOTO 260
C
C     Form the W-vector of the chosen original interpolation point.
C
      DO J=1,N
          W(NPT+J)=XPT(KNEW,J)
      END DO
      DO K=1,NPT
          SUM=ZERO
          IF (K == KOPT) THEN
              CONTINUE
          ELSE IF (PTSID(K) == ZERO) THEN
              DO J=1,N
                  SUM=SUM+W(NPT+J)*XPT(K,J)
              END DO
          ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IP=PTSID(K)
              IP=INT(PTSID(K))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (IP > 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
              IQ=INT(DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IF (IQ > 0) THEN
                  IW=1
                  IF (IP == 0) IW=2
                  SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
              END IF
          END IF
          W(K)=HALF*SUM*SUM
      END DO
C
C     Calculate VLAG and BETA for the required updating of the H matrix if
C     XPT(KNEW,.) is reinstated in the set of interpolation points.
C
      DO K=1,NPT
          SUM=ZERO
          DO J=1,N
              SUM=SUM+BMAT(K,J)*W(NPT+J)
          END DO
          VLAG(K)=SUM
      END DO
      BETA=ZERO
      DO J=1,NPTM
          SUM=ZERO
          DO K=1,NPT
              SUM=SUM+ZMAT(K,J)*W(K)
          END DO
          BETA=BETA-SUM*SUM
          DO K=1,NPT
              VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
          END DO
      END DO
      BSUM=ZERO
      DISTSQ=ZERO
      DO J=1,N
          SUM=ZERO
          DO K=1,NPT
              SUM=SUM+BMAT(K,J)*W(K)
          END DO
          JP=J+NPT
          BSUM=BSUM+SUM*W(JP)
          DO IP=NPT+1,NDIM
              SUM=SUM+BMAT(IP,J)*W(IP)
          END DO
          BSUM=BSUM+SUM*W(JP)
          VLAG(JP)=SUM
          DISTSQ=DISTSQ+XPT(KNEW,J)**2
      END DO
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
      DO K=1,NPT
          IF (PTSID(K) /= ZERO) THEN
              HDIAG=ZERO
              DO J=1,NPTM
                  HDIAG=HDIAG+ZMAT(K,J)**2
              END DO
              DEN=BETA*HDIAG+VLAG(K)**2
              IF (DEN > DENOM) THEN
                  KOLD=K
                  DENOM=DEN
              END IF
          END IF
          VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
      END DO
      IF (DENOM <= 1.0D-2*VLMXSQ) THEN
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
  260 DO KPT=1,NPT
          IF (PTSID(KPT) == ZERO) CYCLE 
          IF (NF >= MAXFUN) THEN
              NF=-1
              GOTO 350
          END IF
          IH=0
          DO J=1,N
              W(J)=XPT(KPT,J)
              XPT(KPT,J)=ZERO
              TEMP=PQ(KPT)*W(J)
              DO I=1,J
                  IH=IH+1
                  HQ(IH)=HQ(IH)+TEMP*W(I)
              END DO
          END DO
          PQ(KPT)=ZERO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IP=PTSID(KPT)
C      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
          IP=INT(PTSID(KPT))
          IQ=INT(DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (IP > 0) THEN
              XP=PTSAUX(1,IP)
              XPT(KPT,IP)=XP
          END IF
          IF (IQ > 0) THEN
              XQ=PTSAUX(1,IQ)
              IF (IP == 0) XQ=PTSAUX(2,IQ)
              XPT(KPT,IQ)=XQ
          END IF
C
C     Set VQUAD to the value of the current model at the new point.
C
          VQUAD=FBASE
          IF (IP > 0) THEN
              IHP=(IP+IP*IP)/2
              VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
          END IF
          IF (IQ > 0) THEN
              IHQ=(IQ+IQ*IQ)/2
              VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
              IF (IP > 0) THEN
                  IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
                  VQUAD=VQUAD+XP*XQ*HQ(IW)
              END IF
          END IF
          DO K=1,NPT
              TEMP=ZERO
              IF (IP > 0) TEMP=TEMP+XP*XPT(K,IP)
              IF (IQ > 0) TEMP=TEMP+XQ*XPT(K,IQ)
              VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
          END DO
C
C     Calculate F at the new interpolation point, and set DIFF to the factor
C     that is going to multiply the KPT-th Lagrange function when the model
C     is updated to provide interpolation to the new function value.
C
          DO I=1,N
              W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
              IF (XPT(KPT,I) == SL(I)) W(I)=XL(I)
              IF (XPT(KPT,I) == SU(I)) W(I)=XU(I)
          END DO
          NF=NF+1
          CALL CALFUN (N,W,F)
          IF (IPRINT == 3) THEN
              PRINT 300, NF,F,(W(I),I=1,N)
  300         FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1          '    The corresponding X is:'/(2X,5D15.6))
          END IF
          FVAL(KPT)=F
          IF (F < FVAL(KOPT)) KOPT=KPT
          DIFF=F-VQUAD
C
C     Update the quadratic model. The RETURN from the subroutine occurs when
C     all the new interpolation points are included in the model.
C
          DO I=1,N
              GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
          END DO
          DO K=1,NPT
              SUM=ZERO
              DO J=1,NPTM
                  SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
              END DO
              TEMP=DIFF*SUM
              IF (PTSID(K) == ZERO) THEN
                  PQ(K)=PQ(K)+TEMP
              ELSE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C          IP=PTSID(K)
C          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
                  IP=INT(PTSID(K))
                  IQ=INT(DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IHQ=(IQ*IQ+IQ)/2
                  IF (IP == 0) THEN
                      HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
                  ELSE
                      IHP=(IP*IP+IP)/2
                      HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
                      IF (IQ > 0) THEN
                          HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                          IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                          HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
                      END IF
                  END IF
              END IF
          END DO
          PTSID(KPT)=ZERO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     By Tom (on 03-06-2019):
C     If a NaN or an infinite value has been reached during the
C     evaluation of the objective function, the loop exit after setting
C     all the parameters, not to raise an exception. KOPT is set to KPT
C     to check in BOBYQB weather FVAL(KOPT) is NaN or infinite value or
C     not.
          IF (F /= F .OR. F > ALMOST_INFINITY) THEN
              EXIT
          END IF
C     By Tom (on 04-06-2019):
C     If the target function value is reached, the loop exit and KOPT is
C     set to KPT to check in BOBYQB weather FVAL(KOPT) .LE. FTARGET
          IF (F <= FTARGET) THEN
              EXIT
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END DO
  350 RETURN
      END
