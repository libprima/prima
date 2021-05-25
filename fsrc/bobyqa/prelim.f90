!*==prelim.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
      SUBROUTINE PRELIM(N,Npt,X,Xl,Xu,Rhobeg,Iprint,Maxfun,Xbase,Xpt,   &
     &                  Fval,Gopt,Hq,Pq,Bmat,Zmat,Ndim,Sl,Su,Nf,Kopt,F, &
     &                  Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
      USE F77KINDS                        
      USE S_CALFUN
      IMPLICIT NONE
!*--PRELIM13
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      INTEGER :: N
      INTEGER , INTENT(IN) :: Npt
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: X
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xl
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Xu
      REAL*8 , INTENT(IN) :: Rhobeg
      INTEGER , INTENT(IN) :: Iprint
      INTEGER , INTENT(IN) :: Maxfun
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xbase
      REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Xpt
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Fval
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Gopt
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: Hq
      REAL*8 , INTENT(OUT) , DIMENSION(*) :: Pq
      REAL*8 , INTENT(INOUT) , DIMENSION(Ndim,*) :: Bmat
      REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Zmat
      INTEGER , INTENT(IN) :: Ndim
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Sl
      REAL*8 , INTENT(IN) , DIMENSION(*) :: Su
      INTEGER , INTENT(INOUT) :: Nf
      INTEGER , INTENT(INOUT) :: Kopt
      REAL*8 :: F
      REAL*8 , INTENT(IN) :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      REAL*8 :: almost_infinity , diff , fbeg , half , one , recip ,      &
     &        rhosq , stepa , stepb , temp , two , zero
      INTEGER :: i , ih , ipt , itemp , j , jpt , k , nfm , nfx , np
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
      half = 0.5D0
      one = 1.0D0
      two = 2.0D0
      zero = 0.0D0
      rhosq = Rhobeg*Rhobeg
      recip = one/rhosq
      np = N + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      almost_infinity = HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
      DO j = 1 , N
         Xbase(j) = X(j)
         DO k = 1 , Npt
            Xpt(k,j) = zero
         ENDDO
         DO i = 1 , Ndim
            Bmat(i,j) = zero
         ENDDO
      ENDDO
      DO ih = 1 , (N*np)/2
         Hq(ih) = zero
      ENDDO
      DO k = 1 , Npt
         Pq(k) = zero
         DO j = 1 , Npt - np
            Zmat(k,j) = zero
         ENDDO
      ENDDO
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
      Nf = 0
      DO
         nfm = Nf
         nfx = Nf - N
         Nf = Nf + 1
         IF ( nfm>2*N ) THEN
            itemp = (nfm-np)/N
            jpt = nfm - itemp*N - N
            ipt = jpt + itemp
            IF ( ipt>N ) THEN
               itemp = jpt
               jpt = ipt - N
               ipt = itemp
            ENDIF
            Xpt(Nf,ipt) = Xpt(ipt+1,ipt)
            Xpt(Nf,jpt) = Xpt(jpt+1,jpt)
         ELSEIF ( nfm>=1 .AND. nfm<=N ) THEN
            stepa = Rhobeg
            IF ( Su(nfm)==zero ) stepa = -stepa
            Xpt(Nf,nfm) = stepa
         ELSEIF ( nfm>N ) THEN
            stepa = Xpt(Nf-N,nfx)
            stepb = -Rhobeg
            IF ( Sl(nfx)==zero ) stepb = DMIN1(two*Rhobeg,Su(nfx))
            IF ( Su(nfx)==zero ) stepb = DMAX1(-two*Rhobeg,Sl(nfx))
            Xpt(Nf,nfx) = stepb
         ENDIF
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
         DO j = 1 , N
            X(j) = DMIN1(DMAX1(Xl(j),Xbase(j)+Xpt(Nf,j)),Xu(j))
            IF ( Xpt(Nf,j)==Sl(j) ) X(j) = Xl(j)
            IF ( Xpt(Nf,j)==Su(j) ) X(j) = Xu(j)
         ENDDO
         CALL CALFUN(N,X,F)
         IF ( Iprint==3 ) THEN
            PRINT 99001 , Nf , F , (X(i),i=1,N)
99001       FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,        &
     &              '    The corresponding X is:'/(2X,5D15.6))
         ENDIF
         Fval(Nf) = F
         IF ( Nf==1 ) THEN
            fbeg = F
            Kopt = 1
         ELSEIF ( F<Fval(Kopt) ) THEN
            Kopt = Nf
         ENDIF
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
         IF ( Nf>2*N+1 ) THEN
            ih = (ipt*(ipt-1))/2 + jpt
            Zmat(1,nfx) = recip
            Zmat(Nf,nfx) = recip
            Zmat(ipt+1,nfx) = -recip
            Zmat(jpt+1,nfx) = -recip
            temp = Xpt(Nf,ipt)*Xpt(Nf,jpt)
            Hq(ih) = (fbeg-Fval(ipt+1)-Fval(jpt+1)+F)/temp
         ELSEIF ( Nf>=2 .AND. Nf<=N+1 ) THEN
            Gopt(nfm) = (F-fbeg)/stepa
            IF ( Npt<Nf+N ) THEN
               Bmat(1,nfm) = -one/stepa
               Bmat(Nf,nfm) = one/stepa
               Bmat(Npt+nfm,nfm) = -half*rhosq
            ENDIF
         ELSEIF ( Nf>=N+2 ) THEN
            ih = (nfx*(nfx+1))/2
            temp = (F-fbeg)/stepb
            diff = stepb - stepa
            Hq(ih) = two*(temp-Gopt(nfx))/diff
            Gopt(nfx) = (Gopt(nfx)*stepb-temp*stepa)/diff
            IF ( stepa*stepb<zero ) THEN
               IF ( F<Fval(Nf-N) ) THEN
                  Fval(Nf) = Fval(Nf-N)
                  Fval(Nf-N) = F
                  IF ( Kopt==Nf ) Kopt = Nf - N
                  Xpt(Nf-N,nfx) = stepb
                  Xpt(Nf,nfx) = stepa
               ENDIF
            ENDIF
            Bmat(1,nfx) = -(stepa+stepb)/(stepa*stepb)
            Bmat(Nf,nfx) = -half/Xpt(Nf-N,nfx)
            Bmat(Nf-N,nfx) = -Bmat(1,nfx) - Bmat(Nf,nfx)
            Zmat(1,nfx) = DSQRT(two)/(stepa*stepb)
            Zmat(Nf,nfx) = DSQRT(half)/rhosq
            Zmat(Nf-N,nfx) = -Zmat(1,nfx) - Zmat(Nf,nfx)
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
         ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
!     If the evaluation returns an NaN or an infinity value, this
!     subroutine is stopped.
         IF ( F/=F .OR. F>almost_infinity ) EXIT
!     By Tom (on 04-06-2019):
!     If the target value is reached, stop the algorithm.
         IF ( F<=Ftarget ) EXIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
         IF ( Nf>=Npt .OR. Nf>=Maxfun ) EXIT
      ENDDO
      END SUBROUTINE PRELIM
