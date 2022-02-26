!*==trstep.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
      SUBROUTINE TRSTEP(N,G,H,Delta,Tol,D,Gg,Td,Tn,W,Piv,Z,Evalue)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
      IMPLICIT NONE
!*--TRSTEP7
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      INTEGER , INTENT(IN) :: N
      REAL*8 , INTENT(IN) , DIMENSION(*) :: G
      REAL*8 , INTENT(INOUT) , DIMENSION(N,*) :: H
      REAL*8 , INTENT(IN) :: Delta
      REAL*8 , INTENT(IN) :: Tol
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: D
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Gg
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Td
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Tn
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Piv
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Z
      REAL*8 , INTENT(OUT) :: Evalue
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      REAL*8 :: delsq , dhd , dnorm , dsq , dtg , dtz , gam , gnorm ,     &
     &        gsq , hnorm , one , par , parl , parlest , paru ,         &
     &        paruest , phi , phil , phiu , pivksv , pivot , posdef ,   &
     &        scale , shfmax , shfmin , shift , slope , sum , sumd ,    &
     &        tdmin , temp , tempa , tempb , two , wsq , wwsq , wz ,    &
     &        zero , zsq
      REAL*8 , DIMENSION(N) :: dsav
      INTEGER :: i , iterc , j , jp , k , kp , kpp , ksav , ksave , nm
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      DIMENSION G(*),H(N,*),D(*),GG(*),TD(*),TN(*),W(*),PIV(*),Z(*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q. Only the upper triangular and diagonal
!       parts need be set. The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.
!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
!     EVALUE will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step. Then EVALUE will be positive, but otherwise it
!     will be set to zero.
!
!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| .LEQ. DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties
!
!             ||D|| .LEQ. DELTA  and  ACTRED .LEQ. (1-TOL)*MAXRED.
!
!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.
!
!     Initialization.
!
      one = 1.0D0
      two = 2.0D0
      zero = 0.0D0
      delsq = Delta*Delta
      Evalue = zero
      nm = N - 1
      DO i = 1 , N
         D(i) = zero
         Td(i) = H(i,i)
         DO j = 1 , i
            H(i,j) = H(j,i)
         ENDDO
      ENDDO
!
!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H. Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.
!
      DO k = 1 , nm
         kp = k + 1
         sum = zero
         IF ( kp<N ) THEN
            kpp = kp + 1
            DO i = kpp , N
               sum = sum + H(i,k)**2
            ENDDO
         ENDIF
         IF ( sum==zero ) THEN
            Tn(k) = H(kp,k)
            H(kp,k) = zero
         ELSE
            temp = H(kp,k)
            Tn(k) = DSIGN(DSQRT(sum+temp*temp),temp)
            H(kp,k) = -sum/(temp+Tn(k))
            temp = DSQRT(two/(sum+H(kp,k)**2))
            DO i = kp , N
               W(i) = temp*H(i,k)
               H(i,k) = W(i)
               Z(i) = Td(i)*W(i)
            ENDDO
            wz = zero
            DO j = kp , nm
               jp = j + 1
               DO i = jp , N
                  Z(i) = Z(i) + H(i,j)*W(j)
                  Z(j) = Z(j) + H(i,j)*W(i)
               ENDDO
               wz = wz + W(j)*Z(j)
            ENDDO
            wz = wz + W(N)*Z(N)
            DO j = kp , N
               Td(j) = Td(j) + W(j)*(wz*W(j)-two*Z(j))
               IF ( j<N ) THEN
                  jp = j + 1
                  DO i = jp , N
                     H(i,j) = H(i,j) - W(i)*Z(j) - W(j)*(Z(i)-wz*W(i))
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
!     Form GG by applying the similarity transformation to G.
!
      gsq = zero
      DO i = 1 , N
         Gg(i) = G(i)
         gsq = gsq + G(i)**2
      ENDDO
      gnorm = DSQRT(gsq)
      DO k = 1 , nm
         kp = k + 1
         sum = zero
         DO i = kp , N
            sum = sum + Gg(i)*H(i,k)
         ENDDO
         DO i = kp , N
            Gg(i) = Gg(i) - sum*H(i,k)
         ENDDO
      ENDDO
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.
!
      hnorm = DABS(Td(1)) + DABS(Tn(1))
      tdmin = Td(1)
      Tn(N) = zero
      DO i = 2 , N
         temp = DABS(Tn(i-1)) + DABS(Td(i)) + DABS(Tn(i))
         hnorm = DMAX1(hnorm,temp)
         tdmin = DMIN1(tdmin,Td(i))
      ENDDO
      IF ( hnorm==zero ) THEN
         IF ( gnorm==zero ) GOTO 99999
         scale = Delta/gnorm
         DO i = 1 , N
            D(i) = -scale*Gg(i)
         ENDDO
         GOTO 500
      ENDIF
!
!     Set the initial values of PAR and its bounds.
!
      parl = DMAX1(zero,-tdmin,gnorm/Delta-hnorm)
      parlest = parl
      par = parl
      paru = zero
      paruest = zero
      posdef = zero
      iterc = 0
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the lines below line number 140
      DO i = 1 , N
         dsav(i) = D(i)
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
 100  iterc = iterc + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019
! The original code can encounter infinite cycling, which did happen
! when testing the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in H (up to 10^219).
! To avoid wasting energy, we do the following
      sumd = zero
      DO i = 1 , N
         sumd = sumd + DABS(D(i))
      ENDDO
      IF ( sumd>=1.0D100 .OR. sumd/=sumd ) THEN
         DO i = 1 , N
            D(i) = dsav(i)
         ENDDO
         GOTO 500
      ELSE
         DO i = 1 , N
            dsav(i) = D(i)
         ENDDO
      ENDIF
      IF ( iterc>MIN(10000,100*N) ) GOTO 500
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ksav = 0
      Piv(1) = Td(1) + par
      k = 1
      DO
         IF ( Piv(k)>zero ) THEN
            Piv(k+1) = Td(k+1) + par - Tn(k)**2/Piv(k)
         ELSE
            IF ( Piv(k)<zero .OR. Tn(k)/=zero ) EXIT
            ksav = k
            Piv(k+1) = Td(k+1) + par
         ENDIF
         k = k + 1
         IF ( k>=N ) THEN
            IF ( Piv(k)<zero ) EXIT
            IF ( Piv(k)==zero ) ksav = k
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
!
            IF ( ksav==0 .AND. gsq>zero ) THEN
!
!     Calculate D for the current PAR in the positive definite case.
!
               W(1) = -Gg(1)/Piv(1)
               DO i = 2 , N
                  W(i) = (-Gg(i)-Tn(i-1)*W(i-1))/Piv(i)
               ENDDO
               D(N) = W(N)
               DO i = nm , 1 , -1
                  D(i) = W(i) - Tn(i)*D(i+1)/Piv(i)
               ENDDO
!
!     Branch if a Newton-Raphson step is acceptable.
!
               dsq = zero
               wsq = zero
               DO i = 1 , N
                  dsq = dsq + D(i)**2
                  wsq = wsq + Piv(i)*W(i)**2
               ENDDO
               IF ( par==zero .AND. dsq<=delsq ) THEN
!
!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
!
                  shfmin = zero
                  pivot = Td(1)
                  shfmax = pivot
                  DO k = 2 , N
                     pivot = Td(k) - Tn(k-1)**2/pivot
                     shfmax = DMIN1(shfmax,pivot)
                  ENDDO
!
!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
                  ksave = 0
                  GOTO 400
               ELSE
!
!     Make the usual test for acceptability of a full trust region step.
!
                  dnorm = DSQRT(dsq)
                  phi = one/dnorm - one/Delta
                  temp = Tol*(one+par*dsq/wsq) - dsq*phi*phi
                  IF ( temp>=zero ) THEN
                     scale = Delta/dnorm
                     DO i = 1 , N
                        D(i) = scale*D(i)
                     ENDDO
                     GOTO 500
                  ENDIF
                  IF ( iterc>=2 .AND. par<=parl ) GOTO 500
                  IF ( paru>zero .AND. par>=paru ) GOTO 500
!
!     Complete the iteration when PHI is negative.
!
                  IF ( phi<zero ) THEN
                     parlest = par
                     IF ( posdef==one ) THEN
                        IF ( phi<=phil ) GOTO 500
                        slope = (phi-phil)/(par-parl)
                        parlest = par - phi/slope
                     ENDIF
                     slope = one/gnorm
                     IF ( paru>zero ) slope = (phiu-phi)/(paru-par)
                     temp = par - phi/slope
                     IF ( paruest>zero ) temp = DMIN1(temp,paruest)
                     paruest = temp
                     posdef = one
                     parl = par
                     phil = phi
                     GOTO 300
                  ENDIF
!
!     If required, calculate Z for the alternative test for convergence.
!
                  IF ( posdef==zero ) THEN
                     W(1) = one/Piv(1)
                     DO i = 2 , N
                        temp = -Tn(i-1)*W(i-1)
                        W(i) = (DSIGN(one,temp)+temp)/Piv(i)
                     ENDDO
                     Z(N) = W(N)
                     DO i = nm , 1 , -1
                        Z(i) = W(i) - Tn(i)*Z(i+1)/Piv(i)
                     ENDDO
                     wwsq = zero
                     zsq = zero
                     dtz = zero
                     DO i = 1 , N
                        wwsq = wwsq + Piv(i)*W(i)**2
                        zsq = zsq + Z(i)**2
                        dtz = dtz + D(i)*Z(i)
                     ENDDO
!
!     Apply the alternative test for convergence.
!
                     tempa = DABS(delsq-dsq)
                     tempb = DSQRT(dtz*dtz+tempa*zsq)
                     gam = tempa/(DSIGN(tempb,dtz)+dtz)
                     temp = Tol*(wsq+par*delsq) - gam*gam*wwsq
                     IF ( temp>=zero ) THEN
                        DO i = 1 , N
                           D(i) = D(i) + gam*Z(i)
                        ENDDO
                        GOTO 500
                     ENDIF
                     parlest = DMAX1(parlest,par-wwsq/zsq)
                  ENDIF
!
!     Complete the iteration when PHI is positive.
!
                  slope = one/gnorm
                  IF ( paru>zero ) THEN
                     IF ( phi>=phiu ) GOTO 500
                     slope = (phiu-phi)/(paru-par)
                  ENDIF
                  parlest = DMAX1(parlest,par-phi/slope)
                  paruest = par
                  IF ( posdef==one ) THEN
                     slope = (phi-phil)/(par-parl)
                     paruest = par - phi/slope
                  ENDIF
                  paru = par
                  phiu = phi
                  GOTO 300
               ENDIF
            ELSE
               IF ( gsq==zero ) THEN
                  IF ( par==zero ) GOTO 500
                  paru = par
                  paruest = par
                  IF ( ksav==0 ) GOTO 200
               ENDIF
               k = ksav
               EXIT
            ENDIF
         ENDIF
      ENDDO
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
      D(k) = one
      IF ( DABS(Tn(k))<=DABS(Piv(k)) ) THEN
         dsq = one
         dhd = Piv(k)
      ELSE
         temp = Td(k+1) + par
         IF ( temp<=DABS(Piv(k)) ) THEN
            D(k+1) = DSIGN(one,-Tn(k))
            dhd = Piv(k) + temp - two*DABS(Tn(k))
         ELSE
            D(k+1) = -Tn(k)/temp
            dhd = Piv(k) + Tn(k)*D(k+1)
         ENDIF
         dsq = one + D(k+1)**2
      ENDIF
      DO
         IF ( k>1 ) THEN
            k = k - 1
            IF ( Tn(k)/=zero ) THEN
               D(k) = -Tn(k)*D(k+1)/Piv(k)
               dsq = dsq + D(k)**2
               CYCLE
            ENDIF
            DO i = 1 , k
               D(i) = zero
            ENDDO
         ENDIF
         parl = par
         parlest = par - dhd/dsq
         EXIT
      ENDDO
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
 200  temp = paruest
      IF ( gsq==zero ) temp = temp*(one-Tol)
      IF ( paruest>zero .AND. parlest>=temp ) THEN
         dtg = zero
         DO i = 1 , N
            dtg = dtg + D(i)*Gg(i)
         ENDDO
         scale = -DSIGN(Delta/DSQRT(dsq),dtg)
         DO i = 1 , N
            D(i) = scale*D(i)
         ENDDO
         GOTO 500
      ENDIF
!
!     Pick the value of PAR for the next iteration.
!
 300  IF ( paru==zero ) THEN
         par = two*parlest + gnorm/Delta
      ELSE
         par = 0.5D0*(parl+paru)
         par = DMAX1(par,parlest)
      ENDIF
      IF ( paruest>zero ) par = DMIN1(par,paruest)
      GOTO 100
 400  shift = 0.5D0*(shfmin+shfmax)
      k = 1
      temp = Td(1) - shift
      DO
         IF ( temp>zero ) THEN
            Piv(k) = temp
            IF ( k<N ) THEN
               temp = Td(k+1) - shift - Tn(k)**2/temp
               k = k + 1
               CYCLE
            ENDIF
            shfmin = shift
         ELSEIF ( k<ksave ) THEN
            Evalue = shfmin
            EXIT
         ELSEIF ( k==ksave ) THEN
            IF ( pivksv==zero ) THEN
               Evalue = shfmin
               EXIT
            ELSEIF ( Piv(k)-temp<temp-pivksv ) THEN
               pivksv = temp
               shfmax = shift
            ELSE
               pivksv = zero
               shfmax = (shift*Piv(k)-shfmin*temp)/(Piv(k)-temp)
            ENDIF
         ELSE
            ksave = k
            pivksv = temp
            shfmax = shift
         ENDIF
         IF ( shfmin<=0.99D0*shfmax ) GOTO 400
         Evalue = shfmin
         EXIT
      ENDDO
!
!     Apply the inverse Householder transformations to D.
!
 500  nm = N - 1
      DO k = nm , 1 , -1
         kp = k + 1
         sum = zero
         DO i = kp , N
            sum = sum + D(i)*H(i,k)
         ENDDO
         DO i = kp , N
            D(i) = D(i) - sum*H(i,k)
         ENDDO
      ENDDO
!
!     Return from the subroutine.
!
99999 END SUBROUTINE TRSTEP
