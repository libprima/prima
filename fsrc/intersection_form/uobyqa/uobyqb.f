!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of uobyqb.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 31-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==uobyqb.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  XOPT,XNEW,XPT,PQ,PL,H,G,D,VLAG,W)
            SUBROUTINE UOBYQB(N,X,Rhobeg,Rhoend,Iprint,Maxfun,Npt,Xbase,&
     &Xopt, Xnew,Xpt,Pq,Pl,H,G,D,Vlag,W,F,Info,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--UOBYQB14
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER :: N
            REAL*8 , INTENT(INOUT) , DIMENSION(N) :: X
            REAL*8 , INTENT(IN) :: Rhobeg
            REAL*8 , INTENT(IN) :: Rhoend
            INTEGER , INTENT(IN) :: Iprint
            INTEGER , INTENT(IN) :: Maxfun
            INTEGER , INTENT(IN) :: Npt
            REAL*8 , INTENT(INOUT) , DIMENSION(N) :: Xbase
            REAL*8 , INTENT(INOUT) , DIMENSION(N) :: Xopt
            REAL*8 , INTENT(INOUT) , DIMENSION(N) :: Xnew
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,N) :: Xpt
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt-1) :: Pq
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,Npt-1) :: Pl
            REAL*8 , INTENT(INOUT) , DIMENSION(N,N*N) :: H
            REAL*8 , INTENT(INOUT) , DIMENSION(N) :: G
            REAL*8 , INTENT(INOUT) , DIMENSION(N) :: D
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt) :: Vlag
            REAL*8 , INTENT(INOUT) , DIMENSION(MAX(6*N,(N**2+3*N+2)/2)) &
     &:: W
            REAL*8 , INTENT(INOUT) :: F
            INTEGER , INTENT(OUT) :: Info
            REAL*8 , INTENT(IN) :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: almost_infinity , ddknew , delta , detrat , diff ,&
     & distest , dnorm , errtol , estim , evalue , fbase , fopt , fsave &
     &, half , one , ratio , rho , rhosq , sixthm , sum , sumg , sumh , &
     &temp , tempa , tol , two , tworsq , vmax , vquad , wmult , zero
            INTEGER :: i , ih , ip , iq , iw , j , jswitch , k , knew , &
     &kopt , ksave , ktemp , nf , nftest , nnp , nptm
!*++
!*++ End of declarations rewritten by SPAG
!*++
!      DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),PQ(*),
!     1  PL(NPT,*),H(N,*),G(*),D(*),VLAG(*),W(*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
!       the corresponding arguments in SUBROUTINE UOBYQA.
!     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
!     XBASE will contain a shift of origin that reduces the contributions from
!       rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     PQ will contain the parameters of the quadratic model.
!     PL will contain the parameters of the Lagrange functions.
!     H will provide the second derivatives that TRSTEP and LAGMAX require.
!     G will provide the first derivatives that TRSTEP and LAGMAX require.
!     D is reserved for trial steps from XOPT, except that it will contain
!       diagonal second derivatives during the initialization procedure.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!     The array W will be used for working space. Its length must be at least
!     max [ 6*N, ( N**2 + 3*N + 2 ) / 2 ].
!
!     Set some constants.
!
            one = 1.0D0
            two = 2.0D0
            zero = 0.0D0
            half = 0.5D0
            tol = 0.01D0
            nnp = N + N + 1
            nptm = Npt - 1
            nftest = MAX0(Maxfun,1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            almost_infinity = HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Initialization. NF is the number of function calculations so far.
!
            rho = Rhobeg
            rhosq = rho*rho
            nf = 0
            DO i = 1 , N
               Xbase(i) = X(i)
               DO k = 1 , Npt
                  Xpt(k,i) = zero
               ENDDO
            ENDDO
            DO k = 1 , Npt
               DO j = 1 , nptm
                  Pl(k,j) = zero
               ENDDO
            ENDDO
!
!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.
!
       100 DO i = 1 , N
               X(i) = Xbase(i) + Xpt(nf+1,i)
            ENDDO
            GOTO 500
       200 tworsq = (two*rho)**2
            rhosq = rho*rho
!
!     Form the gradient of the quadratic model at the trust region centre.
!
       300 knew = 0
            ih = N
            DO j = 1 , N
               Xopt(j) = Xpt(kopt,j)
               G(j) = Pq(j)
               DO i = 1 , j
                  ih = ih + 1
                  G(i) = G(i) + Pq(ih)*Xopt(j)
                  IF ( i<j ) G(j) = G(j) + Pq(ih)*Xopt(i)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTEM or LAGMAX is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended
! behavior of the code, including uninitialized indices.
!   80 H(I,J)=PQ(IH)
                  H(i,j) = Pq(ih)
                  IF ( H(i,j)/=H(i,j) ) THEN
                     Info = -3
                     GOTO 700
                  ENDIF
               ENDDO
            ENDDO
            DO i = 1 , N
               IF ( G(i)/=G(i) ) THEN
                  Info = -3
                  GOTO 700
               ENDIF
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
!
            CALL TRSTEP(N,G,H,delta,tol,D,W(1),W(N+1),W(2*N+1),W(3*N+1),&
     & W(4*N+1),W(5*N+1),evalue)
            temp = zero
            DO i = 1 , N
               temp = temp + D(i)**2
            ENDDO
            dnorm = DMIN1(delta,DSQRT(temp))
            errtol = -one
            IF ( dnorm<half*rho ) THEN
               knew = -1
               errtol = half*evalue*rho*rho
               IF ( nf<=Npt+9 ) errtol = zero
               GOTO 600
            ENDIF
!
!     Calculate the next value of the objective function.
!
       400 DO i = 1 , N
               Xnew(i) = Xopt(i) + D(i)
               X(i) = Xbase(i) + Xnew(i)
            ENDDO
       500 IF ( nf>=nftest ) THEN
               IF ( Iprint>0 ) PRINT 99001
      99001 FORMAT (/4X,'Return from UOBYQA because CALFUN has been', ' &
     &called MAXFUN times')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               Info = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               GOTO 700
            ENDIF
            nf = nf + 1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            DO i = 1 , N
               IF ( X(i)/=X(i) ) THEN
                  F = X(i) ! Set F to NaN
                  IF ( nf==1 ) THEN
                     fopt = F
                     DO j = 1 , N
                        Xopt(j) = zero
                     ENDDO
                  ENDIF
                  Info = -1
                  GOTO 700
               ENDIF
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            CALL CALFUN(N,X,F)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Zaikun (commented on 02-06-2019; implemented in 2016):
!     Exit if F has an NaN or almost infinite value.
!     If this happends at the very first function evaluation (i.e.,
!     NF=1), then it is necessary to set FOPT and XOPT before going to
!     530, because these two variables have not been set yet.
            IF ( F/=F .OR. F>almost_infinity ) THEN
               IF ( nf==1 ) THEN
                  fopt = F
                  DO i = 1 , N
                     Xopt(i) = zero
                  ENDDO
               ENDIF
               Info = -2
               GOTO 700
            ENDIF
!     By Zaikun (commented on 02-06-2019; implemented in 2016):
!     Exit if F .LE. FTARGET.
            IF ( F<=Ftarget ) THEN
               Info = 1
               GOTO 800
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            IF ( Iprint==3 ) THEN
               PRINT 99002 , nf , F , (X(i),i=1,N)
      99002 FORMAT (/4X,'Function number',I6,' F =',1PD18.10, ' The corr&
     &esponding X is:'/(2X,5D15.6))
            ENDIF
            IF ( nf<=Npt ) THEN
               IF ( nf==1 ) THEN
                  fopt = F
                  kopt = nf
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  DO i = 1 , N
                     Xopt(i) = Xpt(1,i)
                  ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  fbase = F
                  j = 0
                  jswitch = -1
                  ih = N
               ELSEIF ( F<fopt ) THEN
                  fopt = F
                  kopt = nf
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  DO i = 1 , N
                     Xopt(i) = Xpt(nf,i)
                  ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ENDIF
!
!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.
!
               IF ( nf<=nnp ) THEN
                  jswitch = -jswitch
                  IF ( jswitch>0 ) THEN
                     IF ( j>=1 ) THEN
                        ih = ih + j
                        IF ( W(j)<zero ) THEN
                           D(j) = (fsave+F-two*fbase)/rhosq
                           Pq(j) = (fsave-F)/(two*rho)
                           Pl(1,ih) = -two/rhosq
                           Pl(nf-1,j) = half/rho
                           Pl(nf-1,ih) = one/rhosq
                        ELSE
                           Pq(j) = (4.0D0*fsave-3.0D0*fbase-F)/(two*rho)
                           D(j) = (fbase+F-two*fsave)/rhosq
                           Pl(1,j) = -1.5D0/rho
                           Pl(1,ih) = one/rhosq
                           Pl(nf-1,j) = two/rho
                           Pl(nf-1,ih) = -two/rhosq
                        ENDIF
                        Pq(ih) = D(j)
                        Pl(nf,j) = -half/rho
                        Pl(nf,ih) = one/rhosq
                     ENDIF
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
!
                     IF ( j<N ) THEN
                        j = j + 1
                        Xpt(nf+1,j) = rho
                     ENDIF
                  ELSE
                     fsave = F
                     IF ( F<fbase ) THEN
                        W(j) = rho
                        Xpt(nf+1,j) = two*rho
                     ELSE
                        W(j) = -rho
                        Xpt(nf+1,j) = -rho
                     ENDIF
                  ENDIF
                  IF ( nf<nnp ) GOTO 100
!
!     Form the off-diagonal second derivatives of the initial quadratic model.
!
                  ih = N
                  ip = 1
                  iq = 2
               ENDIF
               ih = ih + 1
               IF ( nf>nnp ) THEN
                  temp = one/(W(ip)*W(iq))
                  tempa = F - fbase - W(ip)*Pq(ip) - W(iq)*Pq(iq)
                  Pq(ih) = (tempa-half*rhosq*(D(ip)+D(iq)))*temp
                  Pl(1,ih) = temp
                  iw = ip + ip
                  IF ( W(ip)<zero ) iw = iw + 1
                  Pl(iw,ih) = -temp
                  iw = iq + iq
                  IF ( W(iq)<zero ) iw = iw + 1
                  Pl(iw,ih) = -temp
                  Pl(nf,ih) = temp
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
!
                  ip = ip + 1
               ENDIF
               IF ( ip==iq ) THEN
                  ih = ih + 1
                  ip = 1
                  iq = iq + 1
               ENDIF
               IF ( nf<Npt ) THEN
                  Xpt(nf+1,ip) = W(ip)
                  Xpt(nf+1,iq) = W(iq)
                  GOTO 100
               ENDIF
!
!     Set parameters to begin the iterations for the current RHO.
!
               sixthm = zero
               delta = rho
               GOTO 200
            ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (KNEW .EQ. -1) GOTO 420
               IF ( knew==-1 ) THEN
                  Info = 0
                  GOTO 700
               ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.
!
               vquad = zero
               ih = N
               DO j = 1 , N
                  W(j) = D(j)
                  vquad = vquad + W(j)*Pq(j)
                  DO i = 1 , j
                     ih = ih + 1
                     W(ih) = D(i)*Xnew(j) + D(j)*Xopt(i)
                     IF ( i==j ) W(ih) = half*W(ih)
                     vquad = vquad + W(ih)*Pq(ih)
                  ENDDO
               ENDDO
               DO k = 1 , Npt
                  temp = zero
                  DO j = 1 , nptm
                     temp = temp + W(j)*Pl(k,j)
                  ENDDO
                  Vlag(k) = temp
               ENDDO
               Vlag(kopt) = Vlag(kopt) + one
!
!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.
!
               diff = F - fopt - vquad
               sum = zero
               DO k = 1 , Npt
                  temp = zero
                  DO i = 1 , N
                     temp = temp + (Xpt(k,i)-Xnew(i))**2
                  ENDDO
                  temp = DSQRT(temp)
                  sum = sum + DABS(temp*temp*temp*Vlag(k))
               ENDDO
               sixthm = DMAX1(sixthm,DABS(diff)/sum)
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. Then branch if D is not a trust region step.
!
               fsave = fopt
               IF ( F<fopt ) THEN
                  fopt = F
                  DO i = 1 , N
                     Xopt(i) = Xnew(i)
                  ENDDO
               ENDIF
               ksave = knew
               IF ( knew<=0 ) THEN
!
!     Pick the next value of DELTA after a trust region step.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (VQUAD .GE. ZERO) THEN
                  IF ( vquad>=zero ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     IF ( Iprint>0 ) PRINT 99003
      99003 FORMAT (/4X,'Return from UOBYQA because a trust', ' region s&
     &tep has failed to reduce Q')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                     Info = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     GOTO 700
                  ENDIF
                  ratio = (F-fsave)/vquad
                  IF ( ratio<=0.1D0 ) THEN
                     delta = half*dnorm
                  ELSEIF ( ratio<=0.7D0 ) THEN
                     delta = DMAX1(half*delta,dnorm)
                  ELSE
                     delta = DMAX1(delta,1.25D0*dnorm,dnorm+rho)
                  ENDIF
                  IF ( delta<=1.5D0*rho ) delta = rho
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
                  ktemp = 0
                  detrat = zero
                  IF ( F>=fsave ) THEN
                     ktemp = kopt
                     detrat = one
                  ENDIF
                  DO k = 1 , Npt
                     sum = zero
                     DO i = 1 , N
                        sum = sum + (Xpt(k,i)-Xopt(i))**2
                     ENDDO
                     temp = DABS(Vlag(k))
                     IF ( sum>rhosq ) temp = temp*(sum/rhosq)**1.5D0
                     IF ( temp>detrat .AND. k/=ktemp ) THEN
                        detrat = temp
                        ddknew = sum
                        knew = k
                     ENDIF
                  ENDDO
                  IF ( knew==0 ) GOTO 600
               ENDIF
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!
               DO i = 1 , N
                  Xpt(knew,i) = Xnew(i)
               ENDDO
               temp = one/Vlag(knew)
               DO j = 1 , nptm
                  Pl(knew,j) = temp*Pl(knew,j)
                  Pq(j) = Pq(j) + diff*Pl(knew,j)
               ENDDO
               DO k = 1 , Npt
                  IF ( k/=knew ) THEN
                     temp = Vlag(k)
                     DO j = 1 , nptm
                        Pl(k,j) = Pl(k,j) - temp*Pl(knew,j)
                     ENDDO
                  ENDIF
               ENDDO
!
!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
               IF ( F<fsave ) THEN
                  kopt = knew
                  GOTO 300
               ENDIF
               IF ( ksave>0 ) GOTO 300
               IF ( dnorm>two*rho ) GOTO 300
               IF ( ddknew>tworsq ) GOTO 300
            ENDIF
!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
       600 DO k = 1 , Npt
               W(k) = zero
               DO i = 1 , N
                  W(k) = W(k) + (Xpt(k,i)-Xopt(i))**2
               ENDDO
            ENDDO
            DO
               knew = -1
               distest = tworsq
               DO k = 1 , Npt
                  IF ( W(k)>distest ) THEN
                     knew = k
                     distest = W(k)
                  ENDIF
               ENDDO
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.
!
               IF ( knew>0 ) THEN
                  ih = N
                  sumh = zero
                  DO j = 1 , N
                     G(j) = Pl(knew,j)
                     DO i = 1 , j
                        ih = ih + 1
                        temp = Pl(knew,ih)
                        G(j) = G(j) + temp*Xopt(i)
                        IF ( i<j ) THEN
                           G(i) = G(i) + temp*Xopt(j)
                           sumh = sumh + temp*temp
                        ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: See the comments below line number 70
!  330     H(I,J)=TEMP
                        H(i,j) = temp
                        IF ( H(i,j)/=H(i,j) ) THEN
                           Info = -3
                           GOTO 700
                        ENDIF
                     ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     sumh = sumh + half*temp*temp
                  ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: See the comments below line number 70
                  DO i = 1 , N
                     IF ( G(i)/=G(i) ) THEN
                        Info = -3
                        GOTO 700
                     ENDIF
                  ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
                  IF ( errtol>zero ) THEN
                     W(knew) = zero
                     sumg = zero
                     DO i = 1 , N
                        sumg = sumg + G(i)**2
                     ENDDO
                     estim = rho*(DSQRT(sumg)+rho*DSQRT(half*sumh))
                     wmult = sixthm*distest**1.5D0
                     IF ( wmult*estim<=errtol ) CYCLE
                  ENDIF
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
                  CALL LAGMAX(N,G,H,rho,D,Xnew,vmax)
                  IF ( errtol>zero ) THEN
                     IF ( wmult*vmax<=errtol ) CYCLE
                  ENDIF
                  GOTO 400
               ENDIF
               IF ( dnorm>rho ) GOTO 300
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
               IF ( rho>Rhoend ) THEN
                  ih = N
                  DO j = 1 , N
                     Xbase(j) = Xbase(j) + Xopt(j)
                     DO k = 1 , Npt
                        Xpt(k,j) = Xpt(k,j) - Xopt(j)
                     ENDDO
                     DO i = 1 , j
                        ih = ih + 1
                        Pq(i) = Pq(i) + Pq(ih)*Xopt(j)
                        IF ( i<j ) THEN
                           Pq(j) = Pq(j) + Pq(ih)*Xopt(i)
                           DO k = 1 , Npt
                              Pl(k,j) = Pl(k,j) + Pl(k,ih)*Xopt(i)
                           ENDDO
                        ENDIF
                        DO k = 1 , Npt
                           Pl(k,i) = Pl(k,i) + Pl(k,ih)*Xopt(j)
                        ENDDO
                     ENDDO
                  ENDDO
!
!     Pick the next values of RHO and DELTA.
!
                  delta = half*rho
                  ratio = rho/Rhoend
                  IF ( ratio<=16.0D0 ) THEN
                     rho = Rhoend
                  ELSEIF ( ratio<=250.0D0 ) THEN
                     rho = DSQRT(ratio)*Rhoend
                  ELSE
                     rho = 0.1D0*rho
                  ENDIF
                  delta = DMAX1(delta,rho)
                  IF ( Iprint>=2 ) THEN
                     IF ( Iprint>=3 ) PRINT 99004
      99004 FORMAT (5X)
                     PRINT 99005 , rho , nf
      99005 FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of', ' function v&
     &alues =',I6)
                     PRINT 99007 , fopt , (Xbase(i),i=1,N)
                  ENDIF
                  GOTO 200
               ENDIF
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               Info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF ( errtol<zero ) EXIT
               GOTO 400
            ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  420 IF (FOPT .LE. F) THEN
       700 IF ( fopt<=F .OR. F/=F ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  X(i) = Xbase(i) + Xopt(i)
               ENDDO
               F = fopt
            ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (IPRINT .GE. 1) THEN
       800 IF ( Iprint>=1 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               PRINT 99006 , nf
      99006 FORMAT (/4X,'At the return from UOBYQA',5X, 'Number of funct&
     &ion values =',I6)
               PRINT 99007 , F , (X(i),i=1,N)
            ENDIF
      99007 FORMAT (4X,'Least value of F =',1PD23.15,9X, 'The correspond&
     &ing X is:'/(2X,5D15.6))
            END SUBROUTINE UOBYQB