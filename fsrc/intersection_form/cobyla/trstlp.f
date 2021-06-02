!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trstlp.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 02-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==trstlp.f90  processed by SPAG 7.50RE at 00:16 on 26 May 2021
            SUBROUTINE TRSTLP(N,M,A,B,Rho,Dx,Ifull,Iact,Z,Zdota,Vmultc,S&
     &dirn, Dxnew,Vmultd)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            IMPLICIT NONE
!*--TRSTLP7
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER , INTENT(IN) :: N
            INTEGER , INTENT(IN) :: M
            REAL*8 , INTENT(IN) , DIMENSION(N,*) :: A
            REAL*8 , INTENT(IN) , DIMENSION(*) :: B
            REAL*8 , INTENT(IN) :: Rho
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Dx
            INTEGER , INTENT(OUT) :: Ifull
            INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iact
            REAL*8 , INTENT(INOUT) , DIMENSION(N,*) :: Z
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Zdota
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Vmultc
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Sdirn
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Dxnew
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Vmultd
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: acca , accb , alpha , beta , dd , optnew , optold &
     &, ratio , resmax , resold , sd , sp , spabs , ss , step , stpful ,&
     & sum , sumabs , sumd , temp , tempa , tot , vsave , zdotv , zdotw &
     &, zdvabs , zdwabs
            REAL*8 , DIMENSION(N) :: dsav
            INTEGER :: i , icon , icount , isave , iterc , j , k , kk , &
     &kl , kp , kw , mcon , nact , nactx
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This subroutine calculates an N-component vector DX by applying the
!     following two stages. In the first stage, DX is set to the shortest
!     vector that minimizes the greatest violation of the constraints
!       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
!     subject to the Euclidean length of DX being at most RHO. If its length is
!     strictly less than RHO, then we use the resultant freedom in DX to
!     minimize the objective function
!              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
!     subject to no increase in any greatest constraint violation. This
!     notation allows the gradient of the objective function to be regarded as
!     the gradient of a constraint. Therefore the two stages are distinguished
!     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
!     degeneracy may prevent DX from attaining the target length RHO. Then the
!     value IFULL=0 would be set, but usually IFULL=1 on return.
!
!     In general NACT is the number of constraints in the active set and
!     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
!     contains a permutation of the remaining constraint indices. Further, Z is
!     an orthogonal matrix whose first NACT columns can be regarded as the
!     result of Gram-Schmidt applied to the active constraint gradients. For
!     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
!     column of Z with the gradient of the J-th active constraint. DX is the
!     current vector of variables and here the residuals of the active
!     constraints should be zero. Further, the active constraints have
!     nonnegative Lagrange multipliers that are held at the beginning of
!     VMULTC. The remainder of this vector holds the residuals of the inactive
!     constraints at DX, the ordering of the components of VMULTC being in
!     agreement with the permutation of the indices of the constraints that is
!     in IACT. All these residuals are nonnegative, which is achieved by the
!     shift RESMAX that makes the least residual zero.
!
!     Initialize Z and some other variables. The value of RESMAX will be
!     appropriate to DX=0, while ICON will be the index of a most violated
!     constraint if RESMAX is positive. Usually during the first stage the
!     vector SDIRN gives a search direction that reduces all the active
!     constraint violations by one simultaneously.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
            iterc = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Ifull = 1
            mcon = M
            nact = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESMAX=0.0
            resmax = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO i = 1 , N
               DO j = 1 , N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   10 Z(I,J)=0.0
!      Z(I,I)=1.0
!   20 DX(I)=0.0
                  Z(i,j) = 0.0D0
               ENDDO
               Z(i,i) = 1.0D0
               Dx(i) = 0.0D0
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( M>=1 ) THEN
               DO k = 1 , M
                  IF ( B(k)>resmax ) THEN
                     resmax = B(k)
                     icon = k
                  ENDIF
               ENDDO
               DO k = 1 , M
                  Iact(k) = k
                  Vmultc(k) = resmax - B(k)
               ENDDO
            ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RESMAX .EQ. 0.0) GOTO 480
            IF ( resmax==0.0D0 ) GOTO 600
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO i = 1 , N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   50 SDIRN(I)=0.0
               Sdirn(i) = 0.0D0
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
            DO i = 1 , N
               dsav(i) = Dx(i)
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     End the current stage of the calculation if 3 consecutive iterations
!     have either failed to reduce the best calculated value of the objective
!     function or to increase the number of active constraints since the best
!     value was calculated. This strategy prevents cycling, but there is a
!     remote possibility that it will cause premature termination.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   60 OPTOLD=0.0
       100 optold = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            icount = 0
       200 IF ( mcon==M ) THEN
               optnew = resmax
            ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          OPTNEW=0.0
               optnew = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  optnew = optnew - Dx(i)*A(i,mcon)
               ENDDO
            ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019
! The original code can still encounter infinite cycling,
! which did happen when testing the following CUTEst problems:
! DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in A (up to 10^219).
! To avoid wasting energy, we do the following.
            sumd = 0.0D0
            DO i = 1 , N
               sumd = sumd + DABS(Dx(i))
            ENDDO
            IF ( sumd>=1.0D100 .OR. sumd/=sumd ) THEN
               DO i = 1 , N
                  Dx(i) = dsav(i)
               ENDDO
               GOTO 99999
            ELSE
               DO i = 1 , N
                  dsav(i) = Dx(i)
               ENDDO
            ENDIF
            iterc = iterc + 1
            IF ( iterc>MIN(10000,100*N) ) GOTO 99999
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( icount==0 .OR. optnew<optold ) THEN
               optold = optnew
               nactx = nact
               icount = 3
            ELSEIF ( nact>nactx ) THEN
               nactx = nact
               icount = 3
            ELSE
               icount = icount - 1
               IF ( icount==0 ) GOTO 700
            ENDIF
!
!     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
!     the active set. Apply Givens rotations so that the last N-NACT-1 columns
!     of Z are orthogonal to the gradient of the new constraint, a scalar
!     product being set to zero if its nonzero value could be due to computer
!     rounding errors. The array DXNEW is used for working space.
!
            IF ( icon<=nact ) THEN
!
!     Delete the constraint that has the index IACT(ICON) from the active set.
!
               IF ( icon<nact ) THEN
                  isave = Iact(icon)
                  vsave = Vmultc(icon)
                  k = icon
                  DO
                     kp = k + 1
                     kk = Iact(kp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
                     sp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     DO i = 1 , N
                        sp = sp + Z(i,k)*A(i,kk)
                     ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
                     temp = DSQRT(sp*sp+Zdota(kp)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     alpha = Zdota(kp)/temp
                     beta = sp/temp
                     Zdota(kp) = alpha*Zdota(k)
                     Zdota(k) = temp
                     DO i = 1 , N
                        temp = alpha*Z(i,kp) + beta*Z(i,k)
                        Z(i,kp) = alpha*Z(i,k) - beta*Z(i,kp)
                        Z(i,k) = temp
                     ENDDO
                     Iact(k) = kk
                     Vmultc(k) = Vmultc(kp)
                     k = kp
                     IF ( k>=nact ) THEN
                        Iact(k) = isave
                        Vmultc(k) = vsave
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
               nact = nact - 1
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
               IF ( mcon>M ) GOTO 400
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
               temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  temp = temp + Sdirn(i)*Z(i,nact+1)
               ENDDO
               DO i = 1 , N
                  Sdirn(i) = Sdirn(i) - temp*Z(i,nact+1)
               ENDDO
               GOTO 500
            ELSE
               kk = Iact(icon)
               DO i = 1 , N
                  Dxnew(i) = A(i,kk)
               ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TOT=0.0
               tot = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               k = N
               DO
                  IF ( k>nact ) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
!          SPABS=0.0
                     sp = 0.0D0
                     spabs = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     DO i = 1 , N
                        temp = Z(i,k)*Dxnew(i)
                        sp = sp + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  110     SPABS=SPABS+ABS(TEMP)
!          ACCA=SPABS+0.1*ABS(SP)
!          ACCB=SPABS+0.2*ABS(SP)
!          IF (SPABS .GE. ACCA .OR. ACCA .GE. ACCB) SP=0.0
!          IF (TOT .EQ. 0.0) THEN
                        spabs = spabs + DABS(temp)
                     ENDDO
                     acca = spabs + 0.1D0*DABS(sp)
                     accb = spabs + 0.2D0*DABS(sp)
                     IF ( spabs>=acca .OR. acca>=accb ) sp = 0.0D0
                     IF ( tot==0.0D0 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        tot = sp
                     ELSE
                        kp = k + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=SQRT(SP*SP+TOT*TOT)
                        temp = DSQRT(sp*sp+tot*tot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        alpha = sp/temp
                        beta = tot/temp
                        tot = temp
                        DO i = 1 , N
                           temp = alpha*Z(i,k) + beta*Z(i,kp)
                           Z(i,kp) = alpha*Z(i,kp) - beta*Z(i,k)
                           Z(i,k) = temp
                        ENDDO
                     ENDIF
                     k = k - 1
                     CYCLE
                  ENDIF
!
!     Add the new constraint if this can be done without a deletion from the
!     active set.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TOT .NE. 0.0) THEN
                  IF ( tot/=0.0D0 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     nact = nact + 1
                     Zdota(nact) = tot
                     Vmultc(icon) = Vmultc(nact)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTC(NACT)=0.0
                     Vmultc(nact) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     GOTO 300
                  ENDIF
!
!     The next instruction is reached if a deletion has to be made from the
!     active set in order to make room for the new active constraint, because
!     the new constraint gradient is a linear combination of the gradients of
!     the old active constraints. Set the elements of VMULTD to the multipliers
!     of the linear combination. Further, set IOUT to the index of the
!     constraint to be deleted, but branch if no suitable index can be found.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=-1.0
                  ratio = -1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  k = nact
                  EXIT
               ENDDO
               DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  130 ZDOTV=0.0
!      ZDVABS=0.0
                  zdotv = 0.0D0
                  zdvabs = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  DO i = 1 , N
                     temp = Z(i,k)*Dxnew(i)
                     zdotv = zdotv + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  140 ZDVABS=ZDVABS+ABS(TEMP)
!      ACCA=ZDVABS+0.1*ABS(ZDOTV)
!      ACCB=ZDVABS+0.2*ABS(ZDOTV)
                     zdvabs = zdvabs + DABS(temp)
                  ENDDO
                  acca = zdvabs + 0.1D0*DABS(zdotv)
                  accb = zdvabs + 0.2D0*DABS(zdotv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  IF ( zdvabs<acca .AND. acca<accb ) THEN
                     temp = zdotv/Zdota(k)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .GT. 0.0 .AND. IACT(K) .LE. M) THEN
                     IF ( temp>0.0D0 .AND. Iact(k)<=M ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        tempa = Vmultc(k)/temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              IF (RATIO .LT. 0.0 .OR. TEMPA .LT. RATIO) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: IOUT is never used
!                  IOUT=K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        IF ( ratio<0.0D0 .OR. tempa<ratio ) ratio = temp&
     &a
                     ENDIF
                     IF ( k>=2 ) THEN
                        kw = Iact(k)
                        DO i = 1 , N
                           Dxnew(i) = Dxnew(i) - temp*A(i,kw)
                        ENDDO
                     ENDIF
                     Vmultd(k) = temp
                  ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTD(K)=0.0
                     Vmultd(k) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ENDIF
                  k = k - 1
                  IF ( k<=0 ) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RATIO .LT. 0.0) GOTO 490
                     IF ( ratio<0.0D0 ) GOTO 700
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Revise the Lagrange multipliers and reorder the active constraints so
!     that the one to be replaced is at the end of the list. Also calculate the
!     new value of ZDOTA(NACT) and branch if it is not acceptable.
!
                     DO k = 1 , nact
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  160 VMULTC(K)=AMAX1(0.0,VMULTC(K)-RATIO*VMULTD(K))
                        Vmultc(k) = DMAX1(0.0D0,Vmultc(k)-ratio*Vmultd(k&
     &))
                     ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     IF ( icon<nact ) THEN
                        isave = Iact(icon)
                        vsave = Vmultc(icon)
                        k = icon
                        DO
                           kp = k + 1
                           kw = Iact(kp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
                           sp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           DO i = 1 , N
                              sp = sp + Z(i,k)*A(i,kw)
                           ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
                           temp = DSQRT(sp*sp+Zdota(kp)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           alpha = Zdota(kp)/temp
                           beta = sp/temp
                           Zdota(kp) = alpha*Zdota(k)
                           Zdota(k) = temp
                           DO i = 1 , N
                              temp = alpha*Z(i,kp) + beta*Z(i,k)
                              Z(i,kp) = alpha*Z(i,k) - beta*Z(i,kp)
                              Z(i,k) = temp
                           ENDDO
                           Iact(k) = kw
                           Vmultc(k) = Vmultc(kp)
                           k = kp
                           IF ( k>=nact ) THEN
                              Iact(k) = isave
                              Vmultc(k) = vsave
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
                     temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     DO i = 1 , N
                        temp = temp + Z(i,nact)*A(i,kk)
                     ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TEMP .EQ. 0.0) GOTO 490
                     IF ( temp==0.0D0 ) GOTO 700
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     Zdota(nact) = temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(ICON)=0.0
                     Vmultc(icon) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     Vmultc(nact) = ratio
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
!
!     Update IACT and ensure that the objective function continues to be
!     treated as the last active constraint when MCON>M.
!
       300 Iact(icon) = Iact(nact)
            Iact(nact) = kk
            IF ( mcon>M .AND. kk/=mcon ) THEN
               k = nact - 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
               sp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  sp = sp + Z(i,k)*A(i,kk)
               ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
               temp = DSQRT(sp*sp+Zdota(nact)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               alpha = Zdota(nact)/temp
               beta = sp/temp
               Zdota(nact) = alpha*Zdota(k)
               Zdota(k) = temp
               DO i = 1 , N
                  temp = alpha*Z(i,nact) + beta*Z(i,k)
                  Z(i,nact) = alpha*Z(i,k) - beta*Z(i,nact)
                  Z(i,k) = temp
               ENDDO
               Iact(nact) = Iact(k)
               Iact(k) = kk
               temp = Vmultc(k)
               Vmultc(k) = Vmultc(nact)
               Vmultc(nact) = temp
            ENDIF
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
            IF ( mcon<=M ) THEN
               kk = Iact(nact)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
               temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  temp = temp + Sdirn(i)*A(i,kk)
               ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=TEMP-1.0
               temp = temp - 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               temp = temp/Zdota(nact)
               DO i = 1 , N
                  Sdirn(i) = Sdirn(i) - temp*Z(i,nact)
               ENDDO
               GOTO 500
            ENDIF
!
!     Pick the next search direction of stage two.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  320 TEMP=1.0/ZDOTA(NACT)
       400 temp = 1.0D0/Zdota(nact)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO i = 1 , N
               Sdirn(i) = temp*Z(i,nact)
            ENDDO
!
!     Calculate the step to the boundary of the trust region or take the step
!     that reduces RESMAX to zero. The two statements below that include the
!     factor 1.0E-6 prevent some harmless underflows that occurred in a test
!     calculation. Further, we skip the step if it could be zero within a
!     reasonable tolerance for computer rounding errors.
!
       500 dd = Rho*Rho
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SD=0.0
!      SS=0.0
            sd = 0.0D0
            ss = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO i = 1 , N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
               IF ( DABS(Dx(i))>=1.0D-6*Rho ) dd = dd - Dx(i)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               sd = sd + Dx(i)*Sdirn(i)
               ss = ss + Sdirn(i)**2
            ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (DD .LE. 0.0) GOTO 490
!      TEMP=SQRT(SS*DD)
!      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
            IF ( dd<=0.0D0 ) GOTO 700
            temp = DSQRT(ss*dd)
            IF ( DABS(sd)>=1.0D-6*temp ) temp = DSQRT(ss*dd+sd*sd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            stpful = dd/(temp+sd)
            step = stpful
            IF ( mcon==M ) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          ACCA=STEP+0.1*RESMAX
!          ACCB=STEP+0.2*RESMAX
               acca = step + 0.1D0*resmax
               accb = step + 0.2D0*resmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF ( step>=acca .OR. acca>=accb ) GOTO 600
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          STEP=AMIN1(STEP,RESMAX)
               step = DMIN1(step,resmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ENDIF
!
!     Set DXNEW to the new variables if STEP is the steplength, and reduce
!     RESMAX to the corresponding maximum residual if stage one is being done.
!     Because DXNEW will be changed during the calculation of some Lagrange
!     multipliers, it will be restored to the following value later.
!
            DO i = 1 , N
               Dxnew(i) = Dx(i) + step*Sdirn(i)
            ENDDO
            IF ( mcon==M ) THEN
               resold = resmax
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=0.0
               resmax = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO k = 1 , nact
                  kk = Iact(k)
                  temp = B(kk)
                  DO i = 1 , N
                     temp = temp - A(i,kk)*Dxnew(i)
                  ENDDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=AMAX1(RESMAX,TEMP)
                  resmax = DMAX1(resmax,temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ENDDO
            ENDIF
!
!     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
!     device is included to force VMULTD(K)=0.0 if deviations from this value
!     can be attributed to computer rounding errors. First calculate the new
!     Lagrange multipliers.
!
            k = nact
            DO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  390 ZDOTW=0.0
!      ZDWABS=0.0
               zdotw = 0.0D0
               zdwabs = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  temp = Z(i,k)*Dxnew(i)
                  zdotw = zdotw + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  400 ZDWABS=ZDWABS+ABS(TEMP)
!      ACCA=ZDWABS+0.1*ABS(ZDOTW)
!      ACCB=ZDWABS+0.2*ABS(ZDOTW)
!      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0
                  zdwabs = zdwabs + DABS(temp)
               ENDDO
               acca = zdwabs + 0.1D0*DABS(zdotw)
               accb = zdwabs + 0.2D0*DABS(zdotw)
               IF ( zdwabs>=acca .OR. acca>=accb ) zdotw = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               Vmultd(k) = zdotw/Zdota(k)
               IF ( k>=2 ) THEN
                  kk = Iact(k)
                  DO i = 1 , N
                     Dxnew(i) = Dxnew(i) - Vmultd(k)*A(i,kk)
                  ENDDO
                  k = k - 1
                  CYCLE
               ENDIF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (MCON .GT. M) VMULTD(NACT)=AMAX1(0.0,VMULTD(NACT))
               IF ( mcon>M ) Vmultd(nact) = DMAX1(0.0D0,Vmultd(nact))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Complete VMULTC by finding the new constraint residuals.
!
               DO i = 1 , N
                  Dxnew(i) = Dx(i) + step*Sdirn(i)
               ENDDO
               IF ( mcon>nact ) THEN
                  kl = nact + 1
                  DO k = kl , mcon
                     kk = Iact(k)
                     sum = resmax - B(kk)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SUMABS=RESMAX+ABS(B(KK))
                     sumabs = resmax + DABS(B(kk))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     DO i = 1 , N
                        temp = A(i,kk)*Dxnew(i)
                        sum = sum + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  430     SUMABS=SUMABS+ABS(TEMP)
!          ACCA=SUMABS+0.1*ABS(SUM)
!          ACCB=SUMABS+0.2*ABS(SUM)
!          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
                        sumabs = sumabs + DABS(temp)
                     ENDDO
                     acca = sumabs + 0.1D0*DABS(sum)
                     accb = sumabs + 0.2D0*DABS(sum)
                     IF ( sumabs>=acca .OR. acca>=accb ) sum = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     Vmultd(k) = sum
                  ENDDO
               ENDIF
!
!     Calculate the fraction of the step from DX to DXNEW that will be taken.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=1.0
               ratio = 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               icon = 0
               DO k = 1 , mcon
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (VMULTD(K) .LT. 0.0) THEN
                  IF ( Vmultd(k)<0.0D0 ) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     temp = Vmultc(k)/(Vmultc(k)-Vmultd(k))
                     IF ( temp<ratio ) THEN
                        ratio = temp
                        icon = k
                     ENDIF
                  ENDIF
               ENDDO
!
!     Update DX, VMULTC and RESMAX.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0-RATIO
               temp = 1.0D0 - ratio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               DO i = 1 , N
                  Dx(i) = temp*Dx(i) + ratio*Dxnew(i)
               ENDDO
               DO k = 1 , mcon
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
                  Vmultc(k) = DMAX1(0.0D0,temp*Vmultc(k)+ratio*Vmultd(k)&
     &)
               ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               IF ( mcon==M ) resmax = resold + ratio*(resmax-resold)
!
!     If the full step is not acceptable then begin another iteration.
!     Otherwise switch to stage two or end the calculation.
!
               IF ( icon>0 ) GOTO 200
               IF ( step/=stpful ) EXIT
               GOTO 99999
            ENDDO
       600 mcon = M + 1
            icon = mcon
            Iact(mcon) = mcon
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(MCON)=0.0
            Vmultc(mcon) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            GOTO 100
!
!     We employ any freedom that may be available to reduce the objective
!     function before returning a DX whose length is less than RHO.
!
       700 IF ( mcon==M ) GOTO 600
            Ifull = 0
      99999 END SUBROUTINE TRSTLP