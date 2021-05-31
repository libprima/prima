!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of prelim.f90.
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


!*==prelim.f90  processed by SPAG 7.50RE at 17:53 on 31 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  SP,RESCON,STEP,PQW,W)
            SUBROUTINE PRELIM(N,Npt,M,Amat,B,X,Rhobeg,Iprint,Xbase,Xpt,F&
     &val, Xsav,Xopt,Gopt,Kopt,Hq,Pq,Bmat,Zmat,Idz,Ndim,Sp, Rescon,Step,&
     &Pqw,W,F,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--PRELIM14
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER :: N
            INTEGER :: Npt
            INTEGER , INTENT(IN) :: M
            REAL*8 , INTENT(IN) , DIMENSION(N,*) :: Amat
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: B
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: X
            REAL*8 , INTENT(IN) :: Rhobeg
            INTEGER , INTENT(IN) :: Iprint
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xbase
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Xpt
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Fval
            REAL*8 , INTENT(OUT) , DIMENSION(*) :: Xsav
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Xopt
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Gopt
            INTEGER , INTENT(INOUT) :: Kopt
            REAL*8 , INTENT(OUT) , DIMENSION(*) :: Hq
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Pq
            REAL*8 , INTENT(INOUT) , DIMENSION(Ndim,*) :: Bmat
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Zmat
            INTEGER :: Idz
            INTEGER :: Ndim
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Sp
            REAL*8 , INTENT(OUT) , DIMENSION(*) :: Rescon
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Step
            REAL*8 , DIMENSION(*) :: Pqw
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
            REAL*8 , INTENT(INOUT) :: F
            REAL*8 , INTENT(IN) :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: almost_infinity , bigv , feas , half , one , recip&
     & , reciq , resid , rhosq , temp , test , zero
            INTEGER :: i , ipt , itemp , j , jp , jpt , jsav , k , kbase&
     & , nf , nptm
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
!       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
!       same as the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
!       region centre.
!     IDZ is going to be set to one, so that every element of Diag(DZ) is
!       one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
!       factorization of the leading NPT by NPT submatrix of H.
!     STEP, PQW and W are used for working space, the arrays STEP and PQW
!       being taken from LINCOB. The length of W must be at least N+NPT.
!
!     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
!       for the first iteration, an important feature being that, if any of
!       of the columns of XPT is an infeasible point, then the largest of
!       the constraint violations there is at least 0.2*RHOBEG. It also sets
!       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            zero = 0.0D0
            nptm = Npt - N - 1
            rhosq = Rhobeg*Rhobeg
            recip = one/rhosq
            reciq = DSQRT(half)/rhosq
            test = 0.2D0*Rhobeg
            Idz = 1
            kbase = 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            almost_infinity = HUGE(0.0D0)/2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the initial elements of XPT, BMAT, SP and ZMAT to zero.
!
            DO j = 1 , N
               Xbase(j) = X(j)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               Xopt(j) = zero
               Xsav(j) = Xbase(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO k = 1 , Npt
                  Xpt(k,j) = zero
               ENDDO
               DO i = 1 , Ndim
                  Bmat(i,j) = zero
               ENDDO
            ENDDO
            DO k = 1 , Npt
               Sp(k) = zero
               DO j = 1 , Npt - N - 1
                  Zmat(k,j) = zero
               ENDDO
            ENDDO
!
!     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
!       but they may be altered later to make a constraint violation
!       sufficiently large. The initial nonzero elements of BMAT and of
!       the first min[N,NPT-N-1] columns of ZMAT are set also.
!
            DO j = 1 , N
               Xpt(j+1,j) = Rhobeg
               IF ( j<Npt-N ) THEN
                  jp = N + j + 1
                  Xpt(jp,j) = -Rhobeg
                  Bmat(j+1,j) = half/Rhobeg
                  Bmat(jp,j) = -half/Rhobeg
                  Zmat(1,j) = -reciq - reciq
                  Zmat(j+1,j) = reciq
                  Zmat(jp,j) = reciq
               ELSE
                  Bmat(1,j) = -one/Rhobeg
                  Bmat(j+1,j) = one/Rhobeg
                  Bmat(Npt+j,j) = -half*rhosq
               ENDIF
            ENDDO
!
!     Set the remaining initial nonzero elements of XPT and ZMAT when the
!       number of interpolation points exceeds 2*N+1.
!
            IF ( Npt>2*N+1 ) THEN
               DO k = N + 1 , Npt - N - 1
                  itemp = (k-1)/N
                  ipt = k - itemp*N
                  jpt = ipt + itemp
                  IF ( jpt>N ) jpt = jpt - N
                  Xpt(N+k+1,ipt) = Rhobeg
                  Xpt(N+k+1,jpt) = Rhobeg
                  Zmat(1,k) = recip
                  Zmat(ipt+1,k) = -recip
                  Zmat(jpt+1,k) = -recip
                  Zmat(N+k+1,k) = recip
               ENDDO
            ENDIF
!
!     Update the constraint right hand sides to allow for the shift XBASE.
!
            IF ( M>0 ) THEN
               DO j = 1 , M
                  temp = zero
                  DO i = 1 , N
                     temp = temp + Amat(i,j)*Xbase(i)
                  ENDDO
                  B(j) = B(j) - temp
               ENDDO
            ENDIF
!
!     Go through the initial points, shifting every infeasible point if
!       necessary so that its constraint violation is at least 0.2*RHOBEG.
!
            DO nf = 1 , Npt
               feas = one
               bigv = zero
               j = 0
               DO
                  j = j + 1
                  IF ( j<=M .AND. nf>=2 ) THEN
                     resid = -B(j)
                     DO i = 1 , N
                        resid = resid + Xpt(nf,i)*Amat(i,j)
                     ENDDO
                     IF ( resid<=bigv ) CYCLE
                     bigv = resid
                     jsav = j
                     IF ( resid<=test ) THEN
                        feas = -one
                        CYCLE
                     ENDIF
                     feas = zero
                  ENDIF
                  IF ( feas<zero ) THEN
                     DO i = 1 , N
                        Step(i) = Xpt(nf,i) + (test-bigv)*Amat(i,jsav)
                     ENDDO
                     DO k = 1 , Npt
                        Sp(Npt+k) = zero
                        DO j = 1 , N
                           Sp(Npt+k) = Sp(Npt+k) + Xpt(k,j)*Step(j)
                        ENDDO
                     ENDDO
                     CALL UPDATE(N,Npt,Xpt,Bmat,Zmat,Idz,Ndim,Sp,Step,kb&
     &ase, nf,Pqw,W)
                     DO i = 1 , N
                        Xpt(nf,i) = Step(i)
                     ENDDO
                  ENDIF
!
!     Calculate the objective function at the current interpolation point,
!       and set KOPT to the index of the first trust region centre.
!
                  DO j = 1 , N
                     X(j) = Xbase(j) + Xpt(nf,j)
                  ENDDO
                  F = feas
                  CALL CALFUN(N,X,F)
                  IF ( Iprint==3 ) THEN
                     PRINT 99001 , nf , F , (X(i),i=1,N)
      99001 FORMAT (/4X,'Function number',I6,' F =',1PD18.10, ' The corr&
     &esponding X is:'/(2X,5D15.6))
                  ENDIF
                  IF ( nf==1 ) THEN
                     Kopt = 1
                  ELSEIF ( F<Fval(Kopt) .AND. feas>zero ) THEN
                     Kopt = nf
                  ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  150 FVAL(NF)=F
                  Fval(nf) = F
!     By Tom/Zaikun (on 03-06-2019/07-06-2019):
!     If the objective function reached a NaN or infinite value, or if
!     the value is under the target value, the algorithm go back to
!     LINCOB with updated KOPT and XSAV.
!     Note that we should NOT compare F and FTARGET, because X may not
!     be feasible.
                  IF ( F/=F .OR. F>almost_infinity .OR. Fval(Kopt)<=Ftar&
     &get ) GOTO 100
                  EXIT
               ENDDO
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set PQ for the first quadratic model.
!
       100 DO j = 1 , nptm
               W(j) = zero
               DO k = 1 , Npt
                  W(j) = W(j) + Zmat(k,j)*Fval(k)
               ENDDO
            ENDDO
            DO k = 1 , Npt
               Pq(k) = zero
               DO j = 1 , nptm
                  Pq(k) = Pq(k) + Zmat(k,j)*W(j)
               ENDDO
            ENDDO
!
!     Set XOPT, SP, GOPT and HQ for the first quadratic model.
!
            DO j = 1 , N
               Xopt(j) = Xpt(Kopt,j)
               Xsav(j) = Xbase(j) + Xopt(j)
               Gopt(j) = zero
            ENDDO
            DO k = 1 , Npt
               Sp(k) = zero
               DO j = 1 , N
                  Sp(k) = Sp(k) + Xpt(k,j)*Xopt(j)
               ENDDO
               temp = Pq(k)*Sp(k)
               DO j = 1 , N
                  Gopt(j) = Gopt(j) + Fval(k)*Bmat(k,j) + temp*Xpt(k,j)
               ENDDO
            ENDDO
            DO i = 1 , (N*N+N)/2
               Hq(i) = zero
            ENDDO
!
!     Set the initial elements of RESCON.
!
            DO j = 1 , M
               temp = B(j)
               DO i = 1 , N
                  temp = temp - Xopt(i)*Amat(i,j)
               ENDDO
               temp = DMAX1(temp,zero)
               IF ( temp>=Rhobeg ) temp = -temp
               Rescon(j) = temp
            ENDDO
            END SUBROUTINE PRELIM