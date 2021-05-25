!*==lincob.f90  processed by SPAG 7.50RE at 23:18 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  STEP, SP, XNEW, IACT, RESCON, QFAC, RFAC, PQW, W)
      SUBROUTINE LINCOB(N,Npt,M,Amat,B,X,Rhobeg,Rhoend,Iprint,Maxfun,   &
     &                  Xbase,Xpt,Fval,Xsav,Xopt,Gopt,Hq,Pq,Bmat,Zmat,  &
     &                  Ndim,Step,Sp,Xnew,Iact,Rescon,Qfac,Rfac,Pqw,W,F,&
     &                  Info,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H, O-Z)
      IMPLICIT NONE
!*--LINCOB12
!*** Start of declarations inserted by SPAG
      REAL*8 almost_infinity , Amat , B , Bmat , del , delsav , delta ,   &
     &     dffalt , diff , distsq , F , fopt , fsave , Ftarget , Fval , &
     &     Gopt , half , Hq , one , Pq
      REAL*8 Pqw , Qfac , qoptsq , ratio , Rescon , Rfac , rho , Rhobeg , &
     &     Rhoend , snorm , Sp , ssq , Step , sum , sumz , temp ,       &
     &     tenth , vqalt , vquad , W
      REAL*8 X , Xbase , xdiff , Xnew , Xopt , xoptsq , Xpt , Xsav ,      &
     &     zero , Zmat
      INTEGER i , Iact , idz , ifeas , ih , imprv , Info , ip , Iprint ,&
     &        itest , j , k , knew , kopt , ksave , M , Maxfun , N ,    &
     &        nact , Ndim
      INTEGER nf , nh , np , Npt , nptm , nvala , nvalb
!*** End of declarations inserted by SPAG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION Amat(N,*) , B(*) , X(*) , Xbase(*) , Xpt(Npt,*) ,       &
     &          Fval(*) , Xsav(*) , Xopt(*) , Gopt(*) , Hq(*) , Pq(*) , &
     &          Bmat(Ndim,*) , Zmat(Npt,*) , Step(*) , Sp(*) , Xnew(*) ,&
     &          Iact(*) , Rescon(*) , Qfac(N,*) , Rfac(*) , Pqw(*) ,    &
     &          W(*)
!
!     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
!       identical to the corresponding arguments in SUBROUTINE LINCOA.
!     AMAT is a matrix whose columns are the constraint gradients, scaled
!       so that they have unit length.
!     B contains on entry the right hand sides of the constraints, scaled
!       as above, but later B is modified for variables relative to XBASE.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT contains the interpolation point coordinates relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XSAV holds the best feasible vector of variables so far, without any
!       shift of origin.
!     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
!       the feasible vector of variables that provides the least calculated
!       F so far, this vector being the current trust region centre.
!     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of the big inverse matrix H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix
!       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
!       where the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     STEP is employed for trial steps from XOPT. It is also used for working
!       space when XBASE is shifted and in PRELIM.
!     SP is reserved for the scalar products XOPT^T XPT(K, .), K = 1, 2, ...,NPT,
!       followed by STEP^T XPT(K, .), K = 1, 2, ...,NPT.
!     XNEW is the displacement from XBASE of the vector of variables for
!       the current calculation of F, except that SUBROUTINE TRSTEP uses it
!       for working space.
!     IACT is an integer array for the indices of the active constraints.
!     RESCON holds useful information about the constraint residuals. Every
!       nonnegative RESCON(J) is the residual of the J-th constraint at the
!       current trust region centre. Otherwise, if RESCON(J) is negative, the
!       J-th constraint holds as a strict inequality at the trust region
!       centre, its residual being at least |RESCON(J)|; further, the value
!       of |RESCON(J)| is at least the current trust region radius DELTA.
!     QFAC is the orthogonal part of the QR factorization of the matrix of
!       active constraint gradients, these gradients being ordered in
!       accordance with IACT. When NACT is less than N, columns are added
!       to QFAC to complete an N by N orthogonal matrix, which is important
!       for keeping calculated steps sufficiently close to the boundaries
!       of the active constraints.
!     RFAC is the upper triangular part of this QR factorization, beginning
!       with the first diagonal element, followed by the two elements in the
!       upper triangular part of the second column and so on.
!     PQW is used for working space, mainly for storing second derivative
!       coefficients of quadratic functions. Its length is NPT+N.
!     The array W is also used for working space. The required number of
!       elements, namely MAX[M+3*N, 2*M+N, 2*NPT], is set in LINCOA.
!
!     Set some constants.
!
      half = 0.5D0
      one = 1.0D0
      tenth = 0.1D0
      zero = 0.0D0
      np = N + 1
      nh = (N*np)/2
      nptm = Npt - np
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      almost_infinity = HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 15-08-2019
! See the comments below line number 210
      imprv = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
!       ZMAT and SP for the first iteration. An important feature is that,
!       if the interpolation point XPT(K, .) is not feasible, where K is any
!       integer from [1, NPT], then a change is made to XPT(K, .) if necessary
!       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
!       is set so that XPT(KOPT, .) is the initial trust region centre.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  STEP, PQW, W)
      CALL PRELIM(N,Npt,M,Amat,B,X,Rhobeg,Iprint,Xbase,Xpt,Fval,Xsav,   &
     &            Xopt,Gopt,kopt,Hq,Pq,Bmat,Zmat,idz,Ndim,Sp,Rescon,    &
     &            Step,Pqw,W,F,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
      IF ( F/=F .OR. F>almost_infinity ) THEN
         fopt = Fval(kopt)
         Info = -2
         GOTO 600
      ENDIF
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
!     Note that we should NOT compare F and FTARGET, because X may not
!     be feasible at the exit of PRELIM.
      IF ( Fval(kopt)<=Ftarget ) THEN
         F = Fval(kopt)
         X(1:N) = Xsav(1:N)
         Info = 1
         GOTO 700
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Begin the iterative procedure.
!
      nf = Npt
      fopt = Fval(kopt)
      rho = Rhobeg
      delta = rho
      ifeas = 0
      nact = 0
      itest = 3
 100  knew = 0
      nvala = 0
      nvalb = 0
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!       to BMAT that do not depend on ZMAT.
!
 200  fsave = fopt
      xoptsq = zero
      DO i = 1 , N
         xoptsq = xoptsq + Xopt(i)**2
      ENDDO
      IF ( xoptsq>=1.0D4*delta*delta ) THEN
         qoptsq = 0.25D0*xoptsq
         DO k = 1 , Npt
            sum = zero
            DO i = 1 , N
               sum = sum + Xpt(k,i)*Xopt(i)
            ENDDO
            sum = sum - half*xoptsq
            W(Npt+k) = sum
            Sp(k) = zero
            DO i = 1 , N
               Xpt(k,i) = Xpt(k,i) - half*Xopt(i)
               Step(i) = Bmat(k,i)
               W(i) = sum*Xpt(k,i) + qoptsq*Xopt(i)
               ip = Npt + i
               DO j = 1 , i
                  Bmat(ip,j) = Bmat(ip,j) + Step(i)*W(j) + W(i)*Step(j)
               ENDDO
            ENDDO
         ENDDO
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
         DO k = 1 , nptm
            sumz = zero
            DO i = 1 , Npt
               sumz = sumz + Zmat(i,k)
               W(i) = W(Npt+i)*Zmat(i,k)
            ENDDO
            DO j = 1 , N
               sum = qoptsq*sumz*Xopt(j)
               DO i = 1 , Npt
                  sum = sum + W(i)*Xpt(i,j)
               ENDDO
               Step(j) = sum
               IF ( k<idz ) sum = -sum
               DO i = 1 , Npt
                  Bmat(i,j) = Bmat(i,j) + sum*Zmat(i,k)
               ENDDO
            ENDDO
            DO i = 1 , N
               ip = i + Npt
               temp = Step(i)
               IF ( k<idz ) temp = -temp
               DO j = 1 , i
                  Bmat(ip,j) = Bmat(ip,j) + temp*Step(j)
               ENDDO
            ENDDO
         ENDDO
!
!     Update the right hand sides of the constraints.
!
         IF ( M>0 ) THEN
            DO j = 1 , M
               temp = zero
               DO i = 1 , N
                  temp = temp + Amat(i,j)*Xopt(i)
               ENDDO
               B(j) = B(j) - temp
            ENDDO
         ENDIF
!
!     The following instructions complete the shift of XBASE, including the
!       changes to the parameters of the quadratic model.
!
         ih = 0
         DO j = 1 , N
            W(j) = zero
            DO k = 1 , Npt
               W(j) = W(j) + Pq(k)*Xpt(k,j)
               Xpt(k,j) = Xpt(k,j) - half*Xopt(j)
            ENDDO
            DO i = 1 , j
               ih = ih + 1
               Hq(ih) = Hq(ih) + W(i)*Xopt(j) + Xopt(i)*W(j)
               Bmat(Npt+i,j) = Bmat(Npt+j,i)
            ENDDO
         ENDDO
         DO j = 1 , N
            Xbase(j) = Xbase(j) + Xopt(j)
            Xopt(j) = zero
            Xpt(kopt,j) = zero
         ENDDO
      ENDIF
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 21-03-2020
! Exit if BMAT or ZMAT contians NaN
      DO j = 1 , N
         DO i = 1 , Ndim
            IF ( Bmat(i,j)/=Bmat(i,j) ) THEN
               Info = -3
               GOTO 600
            ENDIF
         ENDDO
      ENDDO
      DO j = 1 , nptm
         DO i = 1 , Npt
            IF ( Zmat(i,j)/=Zmat(i,j) ) THEN
               Info = -3
               GOTO 600
            ENDIF
         ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!
!     In the case KNEW = 0, generate the next trust region step by calling
!       TRSTEP, where SNORM is the current trust region radius initially.
!       The final value of SNORM is the length of the calculated step,
!       except that SNORM is zero on return if the projected gradient is
!       unsuitable for starting the conjugate gradient iterations.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTEP or QMSTEP is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended
! behavior of the code, including uninitialized indices, which can lead
! to segmentation faults.
      DO j = 1 , N
         IF ( Gopt(j)/=Gopt(j) ) THEN
            Info = -3
            GOTO 600
         ENDIF
      ENDDO
      DO i = 1 , nh
         IF ( Hq(i)/=Hq(i) ) THEN
            Info = -3
            GOTO 600
         ENDIF
      ENDDO
      DO i = 1 , Npt
         IF ( Pq(i)/=Pq(i) ) THEN
            Info = -3
            GOTO 600
         ENDIF
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      delsav = delta
      ksave = knew
      IF ( knew==0 ) THEN
         snorm = delta
         DO i = 1 , N
            Xnew(i) = Gopt(i)
         ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 19-03-2020: B is never used in TRSTEP
!          CALL TRSTEP (N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CALL TRSTEP(N,Npt,M,Amat,Xpt,Hq,Pq,nact,Iact,Rescon,Qfac,Rfac, &
     &               snorm,Step,Xnew,W,W(M+1),Pqw,Pqw(np),W(M+np))
!
!     A trust region step is applied whenever its length, namely SNORM, is at
!       least HALF*DELTA. It is also applied if its length is at least 0.1999
!       times DELTA and if a line search of TRSTEP has caused a change to the
!       active set. Otherwise there is a branch below to label 530 or 560.
!
         temp = half*delta
         IF ( Xnew(1)>=half ) temp = 0.1999D0*delta
         IF ( snorm<=temp ) THEN
            delta = half*delta
            IF ( delta<=1.4D0*rho ) delta = rho
            nvala = nvala + 1
            nvalb = nvalb + 1
            temp = snorm/rho
            IF ( delsav>rho ) temp = one
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 24-07-2019
!              IF (TEMP .GE. HALF) NVALA = ZERO
!              IF (TEMP .GE. TENTH) NVALB = ZERO
            IF ( temp>=half ) nvala = 0
            IF ( temp>=tenth ) nvalb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( delsav>rho ) GOTO 400
            IF ( nvala<5 .AND. nvalb<3 ) GOTO 400
            IF ( snorm>zero ) ksave = -1
            GOTO 500
         ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 24-07-2019
!          NVALA = ZERO
!          NVALB = ZERO
         nvala = 0
         nvalb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Alternatively, KNEW is positive. Then the model step is calculated
!       within a trust region of radius DEL, after setting the gradient at
!       XBASE and the second derivative parameters of the KNEW-th Lagrange
!       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
!
      ELSE
         del = DMAX1(tenth*delta,rho)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: See the comments below line number 140
!          DO 160 I = 1, N
!  160     W(I)=BMAT(KNEW, I)
         DO i = 1 , N
            W(i) = Bmat(knew,i)
            IF ( W(i)/=W(i) ) THEN
               Info = -3
               GOTO 600
            ENDIF
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO k = 1 , Npt
            Pqw(k) = zero
         ENDDO
         DO j = 1 , nptm
            temp = Zmat(knew,j)
            IF ( j<idz ) temp = -temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: See the comments below line number 140
! Note that the data in PQW is used in QMSTEP below
!          DO 180 K = 1, NPT
!  180     PQW(K)=PQW(K)+TEMP*ZMAT(K, J)
            DO k = 1 , Npt
               Pqw(k) = Pqw(k) + temp*Zmat(k,j)
               IF ( Pqw(k)/=Pqw(k) ) THEN
                  Info = -3
                  GOTO 600
               ENDIF
            ENDDO
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used in QMSTEP
!          CALL QMSTEP (N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CALL QMSTEP(N,Npt,M,Amat,Xpt,Xopt,nact,Iact,Rescon,Qfac,kopt,  &
     &               knew,del,Step,W,Pqw,W(np),W(np+M),ifeas)
      ENDIF
!
!     Set VQUAD to the change to the quadratic model when the move STEP is
!       made from XOPT. If STEP is a trust region step, then VQUAD should be
!       negative. If it is nonnegative due to rounding errors in this case,
!       there is a branch to label 530 to try to improve the model.
!
      vquad = zero
      ih = 0
      DO j = 1 , N
         vquad = vquad + Step(j)*Gopt(j)
         DO i = 1 , j
            ih = ih + 1
            temp = Step(i)*Step(j)
            IF ( i==j ) temp = half*temp
            vquad = vquad + temp*Hq(ih)
         ENDDO
      ENDDO
      DO k = 1 , Npt
         temp = zero
         DO j = 1 , N
            temp = temp + Xpt(k,j)*Step(j)
            Sp(Npt+k) = temp
         ENDDO
         vquad = vquad + half*Pq(k)*temp*temp
      ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 15-08-2019
! Although very rarely, with the original code, an infinite loop can occur
! in the following scenario.
! Suppose that, at an certain iteration,
! KNEW = 0, SNORM > 0.5*DELTA > RHO, VQUAD >= 0, and
! sum_{K = 1}^NPT ||XPT(K, :)-XOPT(:)||^2 < DELTA^2
! (i.e., DELTA is large and SNORM is not small, yet VQUAD >= 0 due to
! rounding errors and XPT are not far from XOPT).
! Then the program will goto 530 and then goto 20, where XBASE may be
! shifted to the current best point, in the hope of reducing rounding
! errors and 'improve' the model. Afterwards, another trust region step
! is produced by the 'improved' model. Note that DELTA remains unchanged
! in this process. If the new trust region step turns out to satisfy
! SNORM > 0.5*DELTA and VQUAD >= 0 again (i.e., the 'improved' model
! still suffers from rounding errors), then the program will goto 530
! and then goto 20, where shifting will not happen because either XBASE
! was already shifted to the current best point in last step, or XBASE
! is close to the current best point. Consequently, the model will
! remain unchanged, and produce the same trust region step again. This
! leads to an infinite loop.
! The infinite loop did happen when the MATLAB interface was applied to
! min atan(x+100) s.t. x <= -99 (x0 = -99, npt = 3, rhobeg = 1, rhoend = 1e-6).
! The problem does not exist in NEWUOA or BOBYQA, where the program will
! exit immediately when VQUAD >= 0.
! To prevent such a loop, here we use IMPRV to record whether the path
! 530--> 20 has already happened for last trust region step. IMPRV = 1
! implies that last trust region step satisfies VQUAD >= 0 and followed
! 530--> 20. With IMPRV = 1, if VQUAD is again nonnegative for the new trust
! region step, we should not goto 530 but goto 560, where IMPRV will be
! set to 0 and DELTA will be reduced. Otherwise, an infinite loop would happen.
!      IF (KSAVE .EQ. 0 .AND. VQUAD .GE. ZERO) GOTO 530
      IF ( ksave==0 .AND. .NOT.(vquad<zero) ) THEN
         IF ( imprv==1 ) GOTO 500
         imprv = 1
         GOTO 400
      ELSE
         imprv = 0
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the next value of the objective function. The difference
!       between the actual new value of F and the value predicted by the
!       model is recorded in DIFF.
!
 300  nf = nf + 1
      IF ( nf>Maxfun ) THEN
         nf = nf - 1
         IF ( Iprint>0 ) PRINT 99001
99001    FORMAT (/4X,'Return from LINCOA because CALFUN has been',      &
     &           ' called MAXFUN times.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         Info = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         GOTO 600
      ENDIF
      xdiff = zero
      DO i = 1 , N
         Xnew(i) = Xopt(i) + Step(i)
         X(i) = Xbase(i) + Xnew(i)
         xdiff = xdiff + (X(i)-Xsav(i))**2
      ENDDO
      xdiff = DSQRT(xdiff)
      IF ( ksave==-1 ) xdiff = rho
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (XDIFF .LE. TENTH*RHO .OR. XDIFF .GE. DELTA+DELTA) THEN
      IF ( xdiff<=tenth*rho .OR. xdiff>=delta+delta ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ifeas = 0
         IF ( Iprint>0 ) PRINT 99002
99002    FORMAT (/4X,'Return from LINCOA because rounding errors',      &
     &           ' prevent reasonable changes to X.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         Info = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         GOTO 600
      ENDIF
      IF ( ksave<=0 ) ifeas = 1
      F = DFLOAT(ifeas)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO i = 1 , N
         IF ( X(i)/=X(i) ) THEN
            F = X(i)    ! Set F to NaN
            IF ( nf==1 ) THEN
               fopt = F
               Xopt(1:N) = zero
            ENDIF
            Info = -1
            GOTO 600
         ENDIF
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL CALFUN(N,X,F)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
      IF ( F/=F .OR. F>almost_infinity ) THEN
         IF ( nf==1 ) THEN
            fopt = F
            Xopt(1:N) = zero
         ENDIF
         Info = -2
         GOTO 600
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( Iprint==3 ) THEN
         PRINT 99003 , nf , F , (X(i),i=1,N)
99003    FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,           &
     &           '    The corresponding X is:'/(2X,5D15.6))
      ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (KSAVE .EQ. -1) GOTO 600
      IF ( ksave==-1 ) THEN
         Info = 0
         GOTO 600
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff = F - fopt - vquad
!
!     If X is feasible, then set DFFALT to the difference between the new
!       value of F and the value predicted by the alternative model.
!
      IF ( ifeas==1 .AND. itest<3 ) THEN
         DO k = 1 , Npt
            Pqw(k) = zero
            W(k) = Fval(k) - Fval(kopt)
         ENDDO
         DO j = 1 , nptm
            sum = zero
            DO i = 1 , Npt
               sum = sum + W(i)*Zmat(i,j)
            ENDDO
            IF ( j<idz ) sum = -sum
            DO k = 1 , Npt
               Pqw(k) = Pqw(k) + sum*Zmat(k,j)
            ENDDO
         ENDDO
         vqalt = zero
         DO k = 1 , Npt
            sum = zero
            DO j = 1 , N
               sum = sum + Bmat(k,j)*Step(j)
            ENDDO
            vqalt = vqalt + sum*W(k)
            vqalt = vqalt + Pqw(k)*Sp(Npt+k)*(half*Sp(Npt+k)+Sp(k))
         ENDDO
         dffalt = F - fopt - vqalt
      ENDIF
      IF ( itest==3 ) THEN
         dffalt = diff
         itest = 0
      ENDIF
!
!     Pick the next value of DELTA after a trust region step.
!
      IF ( ksave==0 ) THEN
         ratio = (F-fopt)/vquad
         IF ( ratio<=tenth ) THEN
            delta = half*delta
         ELSEIF ( ratio<=0.7D0 ) THEN
            delta = DMAX1(half*delta,snorm)
         ELSE
            temp = DSQRT(2.0D0)*delta
            delta = DMAX1(half*delta,snorm+snorm)
            delta = DMIN1(delta,temp)
         ENDIF
         IF ( delta<=1.4D0*rho ) delta = rho
      ENDIF
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!       can be moved. If STEP is a trust region step, then KNEW is zero at
!       present, but a positive value is picked by subroutine UPDATE.
!
      CALL UPDATE(N,Npt,Xpt,Bmat,Zmat,idz,Ndim,Sp,Step,kopt,knew,Pqw,W)
      IF ( knew==0 ) THEN
         IF ( Iprint>0 ) PRINT 99004
99004    FORMAT (/4X,                                                   &
     &'Return from LINCOA because the denominator'' of the updating form&
     &ula is zero.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         Info = 9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         GOTO 600
      ENDIF
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 19-03-2020
! Exit if BMAT or ZMAT contians NaN
      DO j = 1 , N
         DO i = 1 , Ndim
            IF ( Bmat(i,j)/=Bmat(i,j) ) THEN
               Info = -3
               GOTO 600
            ENDIF
         ENDDO
      ENDDO
      DO j = 1 , nptm
         DO i = 1 , Npt
            IF ( Zmat(i,j)/=Zmat(i,j) ) THEN
               Info = -3
               GOTO 600
            ENDIF
         ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!
!     If ITEST is increased to 3, then the next quadratic model is the
!       one whose second derivative matrix is least subject to the new
!       interpolation conditions. Otherwise the new model is constructed
!       by the symmetric Broyden method in the usual way.
!
      IF ( ifeas==1 ) THEN
         itest = itest + 1
         IF ( DABS(dffalt)>=tenth*DABS(diff) ) itest = 0
      ENDIF
!
!     Update the second derivatives of the model by the symmetric Broyden
!       method, using PQW for the second derivative parameters of the new
!       KNEW-th Lagrange function. The contribution from the old parameter
!       PQ(KNEW) is included in the second derivative matrix HQ. W is used
!       later for the gradient of the new KNEW-th Lagrange function.
!
      IF ( itest<3 ) THEN
         DO k = 1 , Npt
            Pqw(k) = zero
         ENDDO
         DO j = 1 , nptm
            temp = Zmat(knew,j)
            IF ( temp/=zero ) THEN
               IF ( j<idz ) temp = -temp
               DO k = 1 , Npt
                  Pqw(k) = Pqw(k) + temp*Zmat(k,j)
               ENDDO
            ENDIF
         ENDDO
         ih = 0
         DO i = 1 , N
            W(i) = Bmat(knew,i)
            temp = Pq(knew)*Xpt(knew,i)
            DO j = 1 , i
               ih = ih + 1
               Hq(ih) = Hq(ih) + temp*Xpt(knew,j)
            ENDDO
         ENDDO
         Pq(knew) = zero
         DO k = 1 , Npt
            Pq(k) = Pq(k) + diff*Pqw(k)
         ENDDO
      ENDIF
!
!     Include the new interpolation point with the corresponding updates of
!       SP. Also make the changes of the symmetric Broyden method to GOPT at
!       the old XOPT if ITEST is less than 3.
!
      Fval(knew) = F
      Sp(knew) = Sp(kopt) + Sp(Npt+kopt)
      ssq = zero
      DO i = 1 , N
         Xpt(knew,i) = Xnew(i)
         ssq = ssq + Step(i)**2
      ENDDO
      Sp(Npt+knew) = Sp(Npt+kopt) + ssq
      IF ( itest<3 ) THEN
         DO k = 1 , Npt
            temp = Pqw(k)*Sp(k)
            DO i = 1 , N
               W(i) = W(i) + temp*Xpt(k,i)
            ENDDO
         ENDDO
         DO i = 1 , N
            Gopt(i) = Gopt(i) + diff*W(i)
         ENDDO
      ENDIF
!
!     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
!       least calculated value so far with a feasible vector of variables.
!
      IF ( F<fopt .AND. ifeas==1 ) THEN
         fopt = F
         DO j = 1 , N
            Xsav(j) = X(j)
            Xopt(j) = Xnew(j)
         ENDDO
         kopt = knew
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         By Tom (on 04-06-2019):
         IF ( fopt<=Ftarget ) THEN
            Info = 1
            GOTO 700
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         snorm = DSQRT(ssq)
         DO j = 1 , M
            IF ( Rescon(j)>=delta+snorm ) THEN
               Rescon(j) = snorm - Rescon(j)
            ELSE
               Rescon(j) = Rescon(j) + snorm
               IF ( Rescon(j)+delta>zero ) THEN
                  temp = B(j)
                  DO i = 1 , N
                     temp = temp - Xopt(i)*Amat(i,j)
                  ENDDO
                  temp = DMAX1(temp,zero)
                  IF ( temp>=delta ) temp = -temp
                  Rescon(j) = temp
               ENDIF
            ENDIF
         ENDDO
         DO k = 1 , Npt
            Sp(k) = Sp(k) + Sp(Npt+k)
         ENDDO
!
!     Also revise GOPT when symmetric Broyden updating is applied.
!
         IF ( itest<3 ) THEN
            ih = 0
            DO j = 1 , N
               DO i = 1 , j
                  ih = ih + 1
                  IF ( i<j ) Gopt(j) = Gopt(j) + Hq(ih)*Step(i)
                  Gopt(i) = Gopt(i) + Hq(ih)*Step(j)
               ENDDO
            ENDDO
            DO k = 1 , Npt
               temp = Pq(k)*Sp(Npt+k)
               DO i = 1 , N
                  Gopt(i) = Gopt(i) + temp*Xpt(k,i)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!     Replace the current model by the least Frobenius norm interpolant if
!       this interpolant gives substantial reductions in the predictions
!       of values of F at feasible points.
!
      IF ( itest==3 ) THEN
         DO k = 1 , Npt
            Pq(k) = zero
            W(k) = Fval(k) - Fval(kopt)
         ENDDO
         DO j = 1 , nptm
            sum = zero
            DO i = 1 , Npt
               sum = sum + W(i)*Zmat(i,j)
            ENDDO
            IF ( j<idz ) sum = -sum
            DO k = 1 , Npt
               Pq(k) = Pq(k) + sum*Zmat(k,j)
            ENDDO
         ENDDO
         DO j = 1 , N
            Gopt(j) = zero
            DO i = 1 , Npt
               Gopt(j) = Gopt(j) + W(i)*Bmat(i,j)
            ENDDO
         ENDDO
         DO k = 1 , Npt
            temp = Pq(k)*Sp(k)
            DO i = 1 , N
               Gopt(i) = Gopt(i) + temp*Xpt(k,i)
            ENDDO
         ENDDO
         DO ih = 1 , nh
            Hq(ih) = zero
         ENDDO
      ENDIF
!
!     If a trust region step has provided a sufficient decrease in F, then
!       branch for another trust region calculation. Every iteration that
!       takes a model step is followed by an attempt to take a trust region
!       step.
!
      knew = 0
      IF ( ksave>0 ) GOTO 200
      IF ( ratio>=tenth ) GOTO 200
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
 400  distsq = DMAX1(delta*delta,4.0D0*rho*rho)
      DO k = 1 , Npt
         sum = zero
         DO j = 1 , N
            sum = sum + (Xpt(k,j)-Xopt(j))**2
         ENDDO
         IF ( sum>distsq ) THEN
            knew = k
            distsq = sum
         ENDIF
      ENDDO
!
!     If KNEW is positive, then branch back for the next iteration, which
!       will generate a "model step". Otherwise, if the current iteration
!       has reduced F, or if DELTA was above its lower bound when the last
!       trust region step was calculated, then try a "trust region" step
!       instead.
!
      IF ( knew>0 ) GOTO 200
      knew = 0
      IF ( fopt<fsave ) GOTO 200
      IF ( delsav>rho ) GOTO 200
!
!     The calculations with the current value of RHO are complete.
!       Pick the next value of RHO.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 15-08-2019
! See the comments below line number 210
!  560 IF (RHO .GT. RHOEND) THEN
 500  imprv = 0
      IF ( rho>Rhoend ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         delta = half*rho
         IF ( rho>250.0D0*Rhoend ) THEN
            rho = tenth*rho
         ELSEIF ( rho<=16.0D0*Rhoend ) THEN
            rho = Rhoend
         ELSE
            rho = DSQRT(rho*Rhoend)
         ENDIF
         delta = DMAX1(delta,rho)
         IF ( Iprint>=2 ) THEN
            IF ( Iprint>=3 ) PRINT 99005
99005       FORMAT (5X)
            PRINT 99006 , rho , nf
99006       FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',             &
     &              ' function values =',I6)
            PRINT 99008 , fopt , (Xbase(i)+Xopt(i),i=1,N)
         ENDIF
         GOTO 100
      ENDIF
!
!     Return from the calculation, after branching to label 220 for another
!       Newton-Raphson step if it has not been tried before.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( ksave==-1 ) GOTO 300
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  600 IF (FOPT .LE. F .OR. IFEAS .EQ. 0) THEN
 600  IF ( fopt<=F .OR. ifeas==0 .OR. F/=F ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = 1 , N
            X(i) = Xsav(i)
         ENDDO
         F = fopt
      ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (IPRINT .GE. 1) THEN
 700  IF ( Iprint>=1 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         PRINT 99007 , nf
99007    FORMAT (/4X,'At the return from LINCOA',5X,                    &
     &           'Number of function values =',I6)
         PRINT 99008 , F , (X(i),i=1,N)
      ENDIF
      W(1) = F
      W(2) = DFLOAT(nf) + half
99008 FORMAT (4X,'Least value of F =',1PD23.15,9X,                      &
     &        'The corresponding X is:'/(2X,5D15.6))
      END SUBROUTINE LINCOB
