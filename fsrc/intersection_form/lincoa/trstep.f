!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trstep.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 27-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==trstep.f90  processed by SPAG 7.50RE at 17:53 on 31 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 19-03-2020: B is never used
!      SUBROUTINE TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            SUBROUTINE TRSTEP(N,Npt,M,Amat,Xpt,Hq,Pq,Nact,Iact,Rescon,Qf&
     &ac, Rfac,Snorm,Step,G,Resnew,Resact,D,Dw,W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--TRSTEP13
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER :: N
            INTEGER , INTENT(IN) :: Npt
            INTEGER :: M
            REAL*8 , DIMENSION(N,*) :: Amat
            REAL*8 , INTENT(IN) , DIMENSION(Npt,*) :: Xpt
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Hq
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Pq
            INTEGER :: Nact
            INTEGER , DIMENSION(*) :: Iact
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Rescon
            REAL*8 , DIMENSION(N,*) :: Qfac
            REAL*8 , DIMENSION(*) :: Rfac
            REAL*8 , INTENT(INOUT) :: Snorm
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Step
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: G
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Resnew
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Resact
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: D
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Dw
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: ad , adw , alpbd , alpha , alphm , alpht , beta , &
     &ctest , dd , dg , dgd , ds , gamma , half , one , reduct , resmax &
     &, rhs , scale , snsq , ss , sum , temp , tiny , wgd , zero
            INTEGER :: i , icount , ih , ir , j , jsav , k , ncall
!*++
!*++ End of declarations rewritten by SPAG
!*++
! B is never used
!      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),HQ(*),PQ(*),IACT(*),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
!       are the same as the terms with these names in LINCOB. If RESCON(J)
!       is negative, then |RESCON(J)| must be no less than the trust region
!       radius, so that the J-th constraint can be ignored.
!     SNORM is set to the trust region radius DELTA initially. On the
!       return, however, it is the length of the calculated STEP, which is
!       set to zero if the constraints do not allow a long enough step.
!     STEP is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current STEP, except that
!       on return the value of G(1) is set to ONE instead of to ZERO if
!       and only if GETACT is called more than once.
!     RESNEW, RESACT, D, DW and W are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       zero value of RESNEW(J) indicates that the J-th constraint is active,
!       and otherwise RESNEW(J) is set to the greater of TINY and the actual
!       residual of the J-th constraint for the current STEP. RESACT holds
!       the residuals of the active constraints, which may be positive.
!       D is the search direction of each line search. DW is either another
!       search direction or the change in gradient along D. The length of W
!       must be at least MAX[M,2*N].
!
!     Set some numbers for the conjugate gradient iterations.
!
            half = 0.5D0
            one = 1.0D0
            tiny = 1.0D-60
            zero = 0.0D0
            ctest = 0.01D0
            snsq = Snorm*Snorm
!
!     Set the initial elements of RESNEW, RESACT and STEP.
!
            IF ( M>0 ) THEN
               DO j = 1 , M
                  Resnew(j) = Rescon(j)
                  IF ( Rescon(j)>=Snorm ) THEN
                     Resnew(j) = -one
                  ELSEIF ( Rescon(j)>=zero ) THEN
                     Resnew(j) = DMAX1(Resnew(j),tiny)
                  ENDIF
               ENDDO
               IF ( Nact>0 ) THEN
                  DO k = 1 , Nact
                     Resact(k) = Rescon(Iact(k))
                     Resnew(Iact(k)) = zero
                  ENDDO
               ENDIF
            ENDIF
            DO i = 1 , N
               Step(i) = zero
            ENDDO
            ss = zero
            reduct = zero
            ncall = 0
            DO
!
!     GETACT picks the active set for the current STEP. It also sets DW to
!       the vector closest to -G that is orthogonal to the normals of the
!       active constraints. DW is scaled to have length 0.2*SNORM, as then
!       a move of DW from STEP is allowed by the linear constraints.
!
               ncall = ncall + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used in GETACT
!      CALL GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,RESNEW,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               CALL GETACT(N,M,Amat,Nact,Iact,Qfac,Rfac,Snorm,Resnew,Res&
     &act,G, Dw,W,W(N+1))
               IF ( W(N+1)==zero ) EXIT
               scale = 0.2D0*Snorm/DSQRT(W(N+1))
               DO i = 1 , N
                  Dw(i) = scale*Dw(i)
               ENDDO
!
!     If the modulus of the residual of an active constraint is substantial,
!       then set D to the shortest move from STEP to the boundaries of the
!       active constraints.
!
               resmax = zero
               IF ( Nact>0 ) THEN
                  DO k = 1 , Nact
                     resmax = DMAX1(resmax,Resact(k))
                  ENDDO
               ENDIF
               gamma = zero
               IF ( resmax>1.0D-4*Snorm ) THEN
                  ir = 0
                  DO k = 1 , Nact
                     temp = Resact(k)
                     IF ( k>=2 ) THEN
                        DO i = 1 , k - 1
                           ir = ir + 1
                           temp = temp - Rfac(ir)*W(i)
                        ENDDO
                     ENDIF
                     ir = ir + 1
                     W(k) = temp/Rfac(ir)
                  ENDDO
                  DO i = 1 , N
                     D(i) = zero
                     DO k = 1 , Nact
                        D(i) = D(i) + W(k)*Qfac(i,k)
                     ENDDO
                  ENDDO
!
!     The vector D that has just been calculated is also the shortest move
!       from STEP+DW to the boundaries of the active constraints. Set GAMMA
!       to the greatest steplength of this move that satisfies the trust
!       region bound.
!
                  rhs = snsq
                  ds = zero
                  dd = zero
                  DO i = 1 , N
                     sum = Step(i) + Dw(i)
                     rhs = rhs - sum*sum
                     ds = ds + D(i)*sum
                     dd = dd + D(i)**2
                  ENDDO
                  IF ( rhs>zero ) THEN
                     temp = DSQRT(ds*ds+dd*rhs)
                     IF ( ds<=zero ) THEN
                        gamma = (temp-ds)/dd
                     ELSE
                        gamma = rhs/(temp+ds)
                     ENDIF
                  ENDIF
!
!     Reduce the steplength GAMMA if necessary so that the move along D
!       also satisfies the linear constraints.
!
                  j = 0
                  DO WHILE ( gamma>zero )
                     j = j + 1
                     IF ( Resnew(j)>zero ) THEN
                        ad = zero
                        adw = zero
                        DO i = 1 , N
                           ad = ad + Amat(i,j)*D(i)
                           adw = adw + Amat(i,j)*Dw(i)
                        ENDDO
                        IF ( ad>zero ) THEN
                           temp = DMAX1((Resnew(j)-adw)/ad,zero)
                           gamma = DMIN1(gamma,temp)
                        ENDIF
                     ENDIF
                     IF ( j>=M ) EXIT
                  ENDDO
                  gamma = DMIN1(gamma,one)
               ENDIF
!
!     Set the next direction for seeking a reduction in the model function
!       subject to the trust region bound and the linear constraints.
!
               IF ( gamma<=zero ) THEN
                  DO i = 1 , N
                     D(i) = Dw(i)
                  ENDDO
                  icount = Nact
               ELSE
                  DO i = 1 , N
                     D(i) = Dw(i) + gamma*D(i)
                  ENDDO
                  icount = Nact - 1
               ENDIF
               alpbd = one
               DO
!
!     Set ALPHA to the steplength from STEP along D to the trust region
!       boundary. Return if the first derivative term of this step is
!       sufficiently small or if no further progress is possible.
!
                  icount = icount + 1
                  rhs = snsq - ss
                  IF ( rhs<=zero ) GOTO 100
                  dg = zero
                  ds = zero
                  dd = zero
                  DO i = 1 , N
                     dg = dg + D(i)*G(i)
                     ds = ds + D(i)*Step(i)
                     dd = dd + D(i)**2
                  ENDDO
                  IF ( dg>=zero ) GOTO 100
                  temp = DSQRT(rhs*dd+ds*ds)
                  IF ( ds<=zero ) THEN
                     alpha = (temp-ds)/dd
                  ELSE
                     alpha = rhs/(temp+ds)
                  ENDIF
                  IF ( -alpha*dg<=ctest*reduct ) GOTO 100
!
!     Set DW to the change in gradient along D.
!
                  ih = 0
                  DO j = 1 , N
                     Dw(j) = zero
                     DO i = 1 , j
                        ih = ih + 1
                        IF ( i<j ) Dw(j) = Dw(j) + Hq(ih)*D(i)
                        Dw(i) = Dw(i) + Hq(ih)*D(j)
                     ENDDO
                  ENDDO
                  DO k = 1 , Npt
                     temp = zero
                     DO j = 1 , N
                        temp = temp + Xpt(k,j)*D(j)
                     ENDDO
                     temp = Pq(k)*temp
                     DO i = 1 , N
                        Dw(i) = Dw(i) + temp*Xpt(k,i)
                     ENDDO
                  ENDDO
!
!     Set DGD to the curvature of the model along D. Then reduce ALPHA if
!       necessary to the value that minimizes the model.
!
                  dgd = zero
                  DO i = 1 , N
                     dgd = dgd + D(i)*Dw(i)
                  ENDDO
                  alpht = alpha
                  IF ( dg+alpha*dgd>zero ) alpha = -dg/dgd
!
!     Make a further reduction in ALPHA if necessary to preserve feasibility,
!       and put some scalar products of D with constraint gradients in W.
!
                  alphm = alpha
                  jsav = 0
                  IF ( M>0 ) THEN
                     DO j = 1 , M
                        ad = zero
                        IF ( Resnew(j)>zero ) THEN
                           DO i = 1 , N
                              ad = ad + Amat(i,j)*D(i)
                           ENDDO
                           IF ( alpha*ad>Resnew(j) ) THEN
                              alpha = Resnew(j)/ad
                              jsav = j
                           ENDIF
                        ENDIF
                        W(j) = ad
                     ENDDO
                  ENDIF
                  alpha = DMAX1(alpha,alpbd)
                  alpha = DMIN1(alpha,alphm)
                  IF ( icount==Nact ) alpha = DMIN1(alpha,one)
!
!     Update STEP, G, RESNEW, RESACT and REDUCT.
!
                  ss = zero
                  DO i = 1 , N
                     Step(i) = Step(i) + alpha*D(i)
                     ss = ss + Step(i)**2
                     G(i) = G(i) + alpha*Dw(i)
                  ENDDO
                  IF ( M>0 ) THEN
                     DO j = 1 , M
                        IF ( Resnew(j)>zero ) Resnew(j) = DMAX1(Resnew(j&
     &)-alpha*W(j),tiny)
                     ENDDO
                  ENDIF
                  IF ( icount==Nact .AND. Nact>0 ) THEN
                     DO k = 1 , Nact
                        Resact(k) = (one-gamma)*Resact(k)
                     ENDDO
                  ENDIF
                  reduct = reduct - alpha*(dg+half*alpha*dgd)
!
!     Test for termination. Branch to label 40 if there is a new active
!       constraint and if the distance from STEP to the trust region
!       boundary is at least 0.2*SNORM.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: the code can encounter infinite cycling due to NaN
! values. Exit when NCALL is large or NaN detected.
                  IF ( ncall>MIN(10000,100*(M+1)*N) .OR. alpha/=alpha .O&
     &R. alpht/=alpht .OR. alphm/=alphm .OR. dgd/=dgd .OR. dg/=dg .OR. s&
     &s/=ss .OR. snsq/=snsq .OR. reduct/=reduct ) GOTO 100
                  ! Note: M can be 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IF ( alpha==alpht ) GOTO 100
                  temp = -alphm*(dg+half*alphm*dgd)
                  IF ( temp<=ctest*reduct ) GOTO 100
                  IF ( jsav>0 ) THEN
                     IF ( ss>0.64D0*snsq ) GOTO 100
                     EXIT
                  ENDIF
                  IF ( icount==N ) GOTO 100
!
!     Calculate the next search direction, which is conjugate to the
!       previous one except in the case ICOUNT=NACT.
!
                  IF ( Nact>0 ) THEN
                     DO j = Nact + 1 , N
                        W(j) = zero
                        DO i = 1 , N
                           W(j) = W(j) + G(i)*Qfac(i,j)
                        ENDDO
                     ENDDO
                     DO i = 1 , N
                        temp = zero
                        DO j = Nact + 1 , N
                           temp = temp + Qfac(i,j)*W(j)
                        ENDDO
                        W(N+i) = temp
                     ENDDO
                  ELSE
                     DO i = 1 , N
                        W(N+i) = G(i)
                     ENDDO
                  ENDIF
                  IF ( icount==Nact ) THEN
                     beta = zero
                  ELSE
                     wgd = zero
                     DO i = 1 , N
                        wgd = wgd + W(N+i)*Dw(i)
                     ENDDO
                     beta = wgd/dgd
                  ENDIF
                  DO i = 1 , N
                     D(i) = -W(N+i) + beta*D(i)
                  ENDDO
                  alpbd = zero
               ENDDO
            ENDDO
!
!     Return from the subroutine.
!
100    Snorm = zero
            IF ( reduct>zero ) Snorm = DSQRT(ss)
            G(1) = zero
            IF ( ncall>1 ) G(1) = one
            END SUBROUTINE TRSTEP