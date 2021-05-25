!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of rescue.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 25-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==rescue.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  KOPT,VLAG,PTSAUX,PTSID,W)
            SUBROUTINE RESCUE(N,Npt,Xl,Xu,Iprint,Maxfun,Xbase,Xpt,Fval,X&
     &opt, Gopt,Hq,Pq,Bmat,Zmat,Ndim,Sl,Su,Nf,Delta,Kopt, Vlag,Ptsaux,Pt&
     &sid,W,F,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            USE F77KINDS
            USE S_CALFUN
            USE S_UPDATE
            IMPLICIT NONE
!*--RESCUE14
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER :: N
            INTEGER :: Npt
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Xl
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Xu
            INTEGER , INTENT(IN) :: Iprint
            INTEGER , INTENT(IN) :: Maxfun
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xbase
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Xpt
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Fval
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xopt
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Gopt
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Hq
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Pq
            REAL*8 , INTENT(INOUT) , DIMENSION(Ndim,*) :: Bmat
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Zmat
            INTEGER :: Ndim
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Sl
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Su
            INTEGER , INTENT(INOUT) :: Nf
            REAL*8 , INTENT(IN) :: Delta
            INTEGER , INTENT(INOUT) :: Kopt
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Vlag
            REAL*8 , INTENT(INOUT) , DIMENSION(2,*) :: Ptsaux
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Ptsid
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
            REAL*8 :: F
            REAL*8 , INTENT(IN) :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: almost_infinity , beta , bsum , den , denom , diff&
     & , distsq , dsqmin , fbase , half , hdiag , one , sfrac , sum , su&
     &mpq , temp , vlmxsq , vquad , winc , xp , xq , zero
            INTEGER :: i , ih , ihp , ihq , ip , iq , iw , j , jp , jpn &
     &, k , knew , kold , kpt , np , nptm , nrem
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonzero. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            zero = 0.0D0
            np = N + 1
            sfrac = half/DFLOAT(np)
            nptm = Npt - np
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            almost_infinity = HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
            sumpq = zero
            winc = zero
            DO k = 1 , Npt
               distsq = zero
               DO j = 1 , N
                  Xpt(k,j) = Xpt(k,j) - Xopt(j)
                  distsq = distsq + Xpt(k,j)**2
               ENDDO
               sumpq = sumpq + Pq(k)
               W(Ndim+k) = distsq
               winc = DMAX1(winc,distsq)
               DO j = 1 , nptm
                  Zmat(k,j) = zero
               ENDDO
            ENDDO
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
            ih = 0
            DO j = 1 , N
               W(j) = half*sumpq*Xopt(j)
               DO k = 1 , Npt
                  W(j) = W(j) + Pq(k)*Xpt(k,j)
               ENDDO
               DO i = 1 , j
                  ih = ih + 1
                  Hq(ih) = Hq(ih) + W(i)*Xopt(j) + W(j)*Xopt(i)
               ENDDO
            ENDDO
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
            DO j = 1 , N
               Xbase(j) = Xbase(j) + Xopt(j)
               Sl(j) = Sl(j) - Xopt(j)
               Su(j) = Su(j) - Xopt(j)
               Xopt(j) = zero
               Ptsaux(1,j) = DMIN1(Delta,Su(j))
               Ptsaux(2,j) = DMAX1(-Delta,Sl(j))
               IF ( Ptsaux(1,j)+Ptsaux(2,j)<zero ) THEN
                  temp = Ptsaux(1,j)
                  Ptsaux(1,j) = Ptsaux(2,j)
                  Ptsaux(2,j) = temp
               ENDIF
               IF ( DABS(Ptsaux(2,j))<half*DABS(Ptsaux(1,j)) ) Ptsaux(2,&
     &j) = half*Ptsaux(1,j)
               DO i = 1 , Ndim
                  Bmat(i,j) = zero
               ENDDO
            ENDDO
            fbase = Fval(Kopt)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
            Ptsid(1) = sfrac
            DO j = 1 , N
               jp = j + 1
               jpn = jp + N
               Ptsid(jp) = DFLOAT(j) + sfrac
               IF ( jpn<=Npt ) THEN
                  Ptsid(jpn) = DFLOAT(j)/DFLOAT(np) + sfrac
                  temp = one/(Ptsaux(1,j)-Ptsaux(2,j))
                  Bmat(jp,j) = -temp + one/Ptsaux(1,j)
                  Bmat(jpn,j) = temp + one/Ptsaux(2,j)
                  Bmat(1,j) = -Bmat(jp,j) - Bmat(jpn,j)
                  Zmat(1,j) = DSQRT(2.0D0)/DABS(Ptsaux(1,j)*Ptsaux(2,j))
                  Zmat(jp,j) = Zmat(1,j)*Ptsaux(2,j)*temp
                  Zmat(jpn,j) = -Zmat(1,j)*Ptsaux(1,j)*temp
               ELSE
                  Bmat(1,j) = -one/Ptsaux(1,j)
                  Bmat(jp,j) = one/Ptsaux(1,j)
                  Bmat(j+Npt,j) = -half*Ptsaux(1,j)**2
               ENDIF
            ENDDO
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
            IF ( Npt>=N+np ) THEN
               DO k = 2*np , Npt
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
                  iw = INT((DFLOAT(k-np)-half)/DFLOAT(N))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ip = k - np - iw*N
                  iq = ip + iw
                  IF ( iq>N ) iq = iq - N
                  Ptsid(k) = DFLOAT(ip) + DFLOAT(iq)/DFLOAT(np) + sfrac
                  temp = one/(Ptsaux(1,ip)*Ptsaux(1,iq))
                  Zmat(1,k-np) = temp
                  Zmat(ip+1,k-np) = -temp
                  Zmat(iq+1,k-np) = -temp
                  Zmat(k,k-np) = temp
               ENDDO
            ENDIF
            nrem = Npt
            kold = 1
            knew = Kopt
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
       100 DO j = 1 , N
               temp = Bmat(kold,j)
               Bmat(kold,j) = Bmat(knew,j)
               Bmat(knew,j) = temp
            ENDDO
            DO j = 1 , nptm
               temp = Zmat(kold,j)
               Zmat(kold,j) = Zmat(knew,j)
               Zmat(knew,j) = temp
            ENDDO
            Ptsid(kold) = Ptsid(knew)
            Ptsid(knew) = zero
            W(Ndim+knew) = zero
            nrem = nrem - 1
            IF ( knew/=Kopt ) THEN
               temp = Vlag(kold)
               Vlag(kold) = Vlag(knew)
               Vlag(knew) = temp
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
               CALL UPDATE(N,Npt,Bmat,Zmat,Ndim,Vlag,beta,denom,knew,W)
               IF ( nrem==0 ) GOTO 99999
               DO k = 1 , Npt
                  W(Ndim+k) = DABS(W(Ndim+k))
               ENDDO
            ENDIF
            DO
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
               dsqmin = zero
               DO k = 1 , Npt
                  IF ( W(Ndim+k)>zero ) THEN
                     IF ( dsqmin==zero .OR. W(Ndim+k)<dsqmin ) THEN
                        knew = k
                        dsqmin = W(Ndim+k)
                     ENDIF
                  ENDIF
               ENDDO
               IF ( dsqmin==zero ) THEN
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
                  DO kpt = 1 , Npt
                     IF ( Ptsid(kpt)==zero ) CYCLE
                     IF ( Nf>=Maxfun ) THEN
                        Nf = -1
                        EXIT
                     ENDIF
                     ih = 0
                     DO j = 1 , N
                        W(j) = Xpt(kpt,j)
                        Xpt(kpt,j) = zero
                        temp = Pq(kpt)*W(j)
                        DO i = 1 , j
                           ih = ih + 1
                           Hq(ih) = Hq(ih) + temp*W(i)
                        ENDDO
                     ENDDO
                     Pq(kpt) = zero
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IP=PTSID(KPT)
!      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
                     ip = INT(Ptsid(kpt))
                     iq = INT(DFLOAT(np)*Ptsid(kpt)-DFLOAT(ip*np))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     IF ( ip>0 ) THEN
                        xp = Ptsaux(1,ip)
                        Xpt(kpt,ip) = xp
                     ENDIF
                     IF ( iq>0 ) THEN
                        xq = Ptsaux(1,iq)
                        IF ( ip==0 ) xq = Ptsaux(2,iq)
                        Xpt(kpt,iq) = xq
                     ENDIF
!
!     Set VQUAD to the value of the current model at the new point.
!
                     vquad = fbase
                     IF ( ip>0 ) THEN
                        ihp = (ip+ip*ip)/2
                        vquad = vquad + xp*(Gopt(ip)+half*xp*Hq(ihp))
                     ENDIF
                     IF ( iq>0 ) THEN
                        ihq = (iq+iq*iq)/2
                        vquad = vquad + xq*(Gopt(iq)+half*xq*Hq(ihq))
                        IF ( ip>0 ) THEN
                           iw = MAX0(ihp,ihq) - IABS(ip-iq)
                           vquad = vquad + xp*xq*Hq(iw)
                        ENDIF
                     ENDIF
                     DO k = 1 , Npt
                        temp = zero
                        IF ( ip>0 ) temp = temp + xp*Xpt(k,ip)
                        IF ( iq>0 ) temp = temp + xq*Xpt(k,iq)
                        vquad = vquad + half*Pq(k)*temp*temp
                     ENDDO
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
                     DO i = 1 , N
                        W(i) = DMIN1(DMAX1(Xl(i),Xbase(i)+Xpt(kpt,i)),Xu&
     &(i))
                        IF ( Xpt(kpt,i)==Sl(i) ) W(i) = Xl(i)
                        IF ( Xpt(kpt,i)==Su(i) ) W(i) = Xu(i)
                     ENDDO
                     Nf = Nf + 1
                     CALL CALFUN(N,W,F)
                     IF ( Iprint==3 ) THEN
                        PRINT 99001 , Nf , F , (W(i),i=1,N)
      99001 FORMAT (/4X,'Function number',I6,' F =',1PD18.10, ' The corr&
     &esponding X is:'/(2X,5D15.6))
                     ENDIF
                     Fval(kpt) = F
                     IF ( F<Fval(Kopt) ) Kopt = kpt
                     diff = F - vquad
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
                     DO i = 1 , N
                        Gopt(i) = Gopt(i) + diff*Bmat(kpt,i)
                     ENDDO
                     DO k = 1 , Npt
                        sum = zero
                        DO j = 1 , nptm
                           sum = sum + Zmat(k,j)*Zmat(kpt,j)
                        ENDDO
                        temp = diff*sum
                        IF ( Ptsid(k)==zero ) THEN
                           Pq(k) = Pq(k) + temp
                        ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IP=PTSID(K)
!          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
                           ip = INT(Ptsid(k))
                           iq = INT(DFLOAT(np)*Ptsid(k)-DFLOAT(ip*np))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           ihq = (iq*iq+iq)/2
                           IF ( ip==0 ) THEN
                              Hq(ihq) = Hq(ihq) + temp*Ptsaux(2,iq)**2
                           ELSE
                              ihp = (ip*ip+ip)/2
                              Hq(ihp) = Hq(ihp) + temp*Ptsaux(1,ip)**2
                              IF ( iq>0 ) THEN
                                 Hq(ihq) = Hq(ihq) + temp*Ptsaux(1,iq)**&
     &2
                                 iw = MAX0(ihp,ihq) - IABS(iq-ip)
                                 Hq(iw) = Hq(iw) + temp*Ptsaux(1,ip) *Pt&
     &saux(1,iq)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                     Ptsid(kpt) = zero
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 03-06-2019):
!     If a NaN or an infinite value has been reached during the
!     evaluation of the objective function, the loop exit after setting
!     all the parameters, not to raise an exception. KOPT is set to KPT
!     to check in BOBYQB weather FVAL(KOPT) is NaN or infinite value or
!     not.
                     IF ( F/=F .OR. F>almost_infinity ) EXIT
!     By Tom (on 04-06-2019):
!     If the target function value is reached, the loop exit and KOPT is
!     set to KPT to check in BOBYQB weather FVAL(KOPT) .LE. FTARGET
                     IF ( F<=Ftarget ) EXIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ENDDO
                  EXIT
               ELSE
!
!     Form the W-vector of the chosen original interpolation point.
!
                  DO j = 1 , N
                     W(Npt+j) = Xpt(knew,j)
                  ENDDO
                  DO k = 1 , Npt
                     sum = zero
                     IF ( k==Kopt ) THEN
                     ELSEIF ( Ptsid(k)==zero ) THEN
                        DO j = 1 , N
                           sum = sum + W(Npt+j)*Xpt(k,j)
                        ENDDO
                     ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IP=PTSID(K)
                        ip = INT(Ptsid(k))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        IF ( ip>0 ) sum = W(Npt+ip)*Ptsaux(1,ip)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
                        iq = INT(DFLOAT(np)*Ptsid(k)-DFLOAT(ip*np))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        IF ( iq>0 ) THEN
                           iw = 1
                           IF ( ip==0 ) iw = 2
                           sum = sum + W(Npt+iq)*Ptsaux(iw,iq)
                        ENDIF
                     ENDIF
                     W(k) = half*sum*sum
                  ENDDO
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
                  DO k = 1 , Npt
                     sum = zero
                     DO j = 1 , N
                        sum = sum + Bmat(k,j)*W(Npt+j)
                     ENDDO
                     Vlag(k) = sum
                  ENDDO
                  beta = zero
                  DO j = 1 , nptm
                     sum = zero
                     DO k = 1 , Npt
                        sum = sum + Zmat(k,j)*W(k)
                     ENDDO
                     beta = beta - sum*sum
                     DO k = 1 , Npt
                        Vlag(k) = Vlag(k) + sum*Zmat(k,j)
                     ENDDO
                  ENDDO
                  bsum = zero
                  distsq = zero
                  DO j = 1 , N
                     sum = zero
                     DO k = 1 , Npt
                        sum = sum + Bmat(k,j)*W(k)
                     ENDDO
                     jp = j + Npt
                     bsum = bsum + sum*W(jp)
                     DO ip = Npt + 1 , Ndim
                        sum = sum + Bmat(ip,j)*W(ip)
                     ENDDO
                     bsum = bsum + sum*W(jp)
                     Vlag(jp) = sum
                     distsq = distsq + Xpt(knew,j)**2
                  ENDDO
                  beta = half*distsq*distsq + beta - bsum
                  Vlag(Kopt) = Vlag(Kopt) + one
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
                  denom = zero
                  vlmxsq = zero
                  DO k = 1 , Npt
                     IF ( Ptsid(k)/=zero ) THEN
                        hdiag = zero
                        DO j = 1 , nptm
                           hdiag = hdiag + Zmat(k,j)**2
                        ENDDO
                        den = beta*hdiag + Vlag(k)**2
                        IF ( den>denom ) THEN
                           kold = k
                           denom = den
                        ENDIF
                     ENDIF
                     vlmxsq = DMAX1(vlmxsq,Vlag(k)**2)
                  ENDDO
                  IF ( denom<=1.0D-2*vlmxsq ) THEN
                     W(Ndim+knew) = -W(Ndim+knew) - winc
                     CYCLE
                  ENDIF
                  GOTO 100
               ENDIF
            ENDDO
      99999 END SUBROUTINE RESCUE