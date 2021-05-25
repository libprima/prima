!*==getact.f90  processed by SPAG 7.50RE at 00:12 on 26 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used
!      SUBROUTINE GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GETACT(N,M,Amat,Nact,Iact,Qfac,Rfac,Snorm,Resnew,      &
     &                  Resact,G,Dw,Vlam,W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
      IMPLICIT NONE
!*--GETACT12
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: M
      REAL*8 , INTENT(IN) , DIMENSION(N,*) :: Amat
      INTEGER , INTENT(INOUT) :: Nact
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iact
      REAL*8 , INTENT(INOUT) , DIMENSION(N,*) :: Qfac
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Rfac
      REAL*8 , INTENT(IN) :: Snorm
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Resnew
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Resact
      REAL*8 , INTENT(IN) , DIMENSION(*) :: G
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Dw
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Vlam
      REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      REAL*8 :: cosv , ctol , cval , dd , ddsav , dnorm , one , rdiag ,   &
     &        sinv , sprod , sum , sval , tdel , temp , test , tiny ,   &
     &        violmx , vmult , zero
      INTEGER :: i , ic , idiag , iflag , j , jc , jcp , jdiag , jw ,   &
     &           k , l , nactp
!*++
!*++ End of declarations rewritten by SPAG
!*++
!      DIMENSION AMAT(N,*),B(*),IACT(*),QFAC(N,*),RFAC(*),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms
!       with these names in SUBROUTINE LINCOB. The current values must be
!       set on entry. NACT, IACT, QFAC and RFAC are kept up to date when
!       GETACT changes the current active set.
!     SNORM, RESNEW, RESACT, G and DW are the same as the terms with these
!       names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are
!       also kept up to date.
!     VLAM and W are used for working space, the vector VLAM being reserved
!       for the Lagrange multipliers of the calculation. Their lengths must
!       be at least N.
!     The main purpose of GETACT is to pick the current active set. It is
!       defined by the property that the projection of -G into the space
!       orthogonal to the active constraint normals is as large as possible,
!       subject to this projected steepest descent direction moving no closer
!       to the boundary of every constraint whose current residual is at most
!       0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are
!       all appropriate to this choice of active set.
!     Occasionally this projected direction is zero, and then the final value
!       of W(1) is set to zero. Otherwise, the direction itself is returned
!       in DW, and W(1) is set to the square of the length of the direction.
!
!     Set some constants and a temporary VLAM.
!
      one = 1.0D0
      tiny = 1.0D-60
      zero = 0.0D0
      tdel = 0.2D0*Snorm
      ddsav = zero
      DO i = 1 , N
         ddsav = ddsav + G(i)**2
         Vlam(i) = zero
      ENDDO
      ddsav = ddsav + ddsav
!
!     Set the initial QFAC to the identity matrix in the case NACT=0.
!
      IF ( Nact==0 ) THEN
         DO i = 1 , N
            DO j = 1 , N
               Qfac(i,j) = zero
            ENDDO
            Qfac(i,i) = one
         ENDDO
         GOTO 400
      ENDIF
!
!     Remove any constraints from the initial active set whose residuals
!       exceed TDEL.
!
      iflag = 1
      ic = Nact
 100  IF ( Resact(ic)>tdel ) GOTO 900
 200  ic = ic - 1
      IF ( ic>0 ) GOTO 100
!
!     Remove any constraints from the initial active set whose Lagrange
!       multipliers are nonnegative, and set the surviving multipliers.
!
      iflag = 2
 300  IF ( Nact/=0 ) THEN
         ic = Nact
         DO
            temp = zero
            DO i = 1 , N
               temp = temp + Qfac(i,ic)*G(i)
            ENDDO
            idiag = (ic*ic+ic)/2
            IF ( ic<Nact ) THEN
               jw = idiag + ic
               DO j = ic + 1 , Nact
                  temp = temp - Rfac(jw)*Vlam(j)
                  jw = jw + j
               ENDDO
            ENDIF
            IF ( temp>=zero ) GOTO 900
            Vlam(ic) = temp/Rfac(idiag)
            ic = ic - 1
            IF ( ic<=0 ) EXIT
         ENDDO
      ENDIF
!
!     Set the new search direction D. Terminate if the 2-norm of D is zero
!       or does not decrease, or if NACT=N holds. The situation NACT=N
!       occurs for sufficiently large SNORM if the origin is in the convex
!       hull of the constraint gradients.
!
 400  IF ( Nact==N ) THEN
         dd = zero
         GOTO 800
      ELSE
         DO j = Nact + 1 , N
            W(j) = zero
            DO i = 1 , N
               W(j) = W(j) + Qfac(i,j)*G(i)
            ENDDO
         ENDDO
         dd = zero
         DO i = 1 , N
            Dw(i) = zero
            DO j = Nact + 1 , N
               Dw(i) = Dw(i) - W(j)*Qfac(i,j)
            ENDDO
            dd = dd + Dw(i)**2
         ENDDO
         IF ( dd>=ddsav ) THEN
            dd = zero
            GOTO 800
         ELSE
            IF ( dd==zero ) GOTO 800
            ddsav = dd
            dnorm = DSQRT(dd)
!
!     Pick the next integer L or terminate, a positive value of L being
!       the index of the most violated constraint. The purpose of CTOL
!       below is to estimate whether a positive value of VIOLMX may be
!       due to computer rounding errors.
!
            l = 0
            IF ( M>0 ) THEN
               test = dnorm/Snorm
               violmx = zero
               DO j = 1 , M
                  IF ( Resnew(j)>zero .AND. Resnew(j)<=tdel ) THEN
                     sum = zero
                     DO i = 1 , N
                        sum = sum + Amat(i,j)*Dw(i)
                     ENDDO
                     IF ( sum>test*Resnew(j) ) THEN
                        IF ( sum>violmx ) THEN
                           l = j
                           violmx = sum
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               ctol = zero
               temp = 0.01D0*dnorm
               IF ( violmx>zero .AND. violmx<temp ) THEN
                  IF ( Nact>0 ) THEN
                     DO k = 1 , Nact
                        j = Iact(k)
                        sum = zero
                        DO i = 1 , N
                           sum = sum + Dw(i)*Amat(i,j)
                        ENDDO
                        ctol = DMAX1(ctol,DABS(sum))
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
            W(1) = one
            IF ( l==0 ) GOTO 800
            IF ( violmx<=10.0D0*ctol ) GOTO 800
!
!     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
!       the first (NACT+1) columns of QFAC are the ones required for the
!       addition of the L-th constraint, and add the appropriate column
!       to RFAC.
!
            nactp = Nact + 1
            idiag = (nactp*nactp-nactp)/2
            rdiag = zero
            DO j = N , 1 , -1
               sprod = zero
               DO i = 1 , N
                  sprod = sprod + Qfac(i,j)*Amat(i,l)
               ENDDO
               IF ( j<=Nact ) THEN
                  Rfac(idiag+j) = sprod
               ELSEIF ( DABS(rdiag)<=1.0D-20*DABS(sprod) ) THEN
                  rdiag = sprod
               ELSE
                  temp = DSQRT(sprod*sprod+rdiag*rdiag)
                  cosv = sprod/temp
                  sinv = rdiag/temp
                  rdiag = temp
                  DO i = 1 , N
                     temp = cosv*Qfac(i,j) + sinv*Qfac(i,j+1)
                     Qfac(i,j+1) = -sinv*Qfac(i,j) + cosv*Qfac(i,j+1)
                     Qfac(i,j) = temp
                  ENDDO
               ENDIF
            ENDDO
            IF ( rdiag<zero ) THEN
               DO i = 1 , N
                  Qfac(i,nactp) = -Qfac(i,nactp)
               ENDDO
            ENDIF
            Rfac(idiag+nactp) = DABS(rdiag)
            Nact = nactp
            Iact(Nact) = l
            Resact(Nact) = Resnew(l)
            Vlam(Nact) = zero
            Resnew(l) = zero
         ENDIF
      ENDIF
!
!     Set the components of the vector VMU in W.
!
 500  W(Nact) = one/Rfac((Nact*Nact+Nact)/2)**2
      IF ( Nact>1 ) THEN
         DO i = Nact - 1 , 1 , -1
            idiag = (i*i+i)/2
            jw = idiag + i
            sum = zero
            DO j = i + 1 , Nact
               sum = sum - Rfac(jw)*W(j)
               jw = jw + j
            ENDDO
            W(i) = sum/Rfac(idiag)
         ENDDO
      ENDIF
!
!     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
!
      vmult = violmx
      ic = 0
      j = 1
      DO
         IF ( j<Nact ) THEN
            IF ( Vlam(j)>=vmult*W(j) ) THEN
               ic = j
               vmult = Vlam(j)/W(j)
            ENDIF
            j = j + 1
            CYCLE
         ENDIF
         DO j = 1 , Nact
            Vlam(j) = Vlam(j) - vmult*W(j)
         ENDDO
         IF ( ic>0 ) Vlam(ic) = zero
         violmx = DMAX1(violmx-vmult,zero)
         IF ( ic==0 ) violmx = zero
!
!     Reduce the active set if necessary, so that all components of the
!       new VLAM are negative, with resetting of the residuals of the
!       constraints that become inactive.
!
         iflag = 3
         ic = Nact
         EXIT
      ENDDO
 600  IF ( Vlam(ic)>=zero ) THEN
         Resnew(Iact(ic)) = DMAX1(Resact(ic),tiny)
         GOTO 900
      ENDIF
 700  ic = ic - 1
      IF ( ic>0 ) GOTO 600
!
!     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
!       as then the active constraints imply D=0. Otherwise, go to label
!       100, to calculate the new D and to test for termination.
!
      IF ( violmx>zero ) GOTO 500
      IF ( Nact<N ) GOTO 400
      dd = zero
 800  W(1) = dd
      RETURN
!
!     These instructions rearrange the active constraints so that the new
!       value of IACT(NACT) is the old value of IACT(IC). A sequence of
!       Givens rotations is applied to the current QFAC and RFAC. Then NACT
!       is reduced by one.
!
 900  Resnew(Iact(ic)) = DMAX1(Resact(ic),tiny)
      jc = ic
      DO
         IF ( jc<Nact ) THEN
            jcp = jc + 1
            idiag = jc*jcp/2
            jw = idiag + jcp
            temp = DSQRT(Rfac(jw-1)**2+Rfac(jw)**2)
            cval = Rfac(jw)/temp
            sval = Rfac(jw-1)/temp
            Rfac(jw-1) = sval*Rfac(idiag)
            Rfac(jw) = cval*Rfac(idiag)
            Rfac(idiag) = temp
            IF ( jcp<Nact ) THEN
               DO j = jcp + 1 , Nact
                  temp = sval*Rfac(jw+jc) + cval*Rfac(jw+jcp)
                  Rfac(jw+jcp) = cval*Rfac(jw+jc) - sval*Rfac(jw+jcp)
                  Rfac(jw+jc) = temp
                  jw = jw + j
               ENDDO
            ENDIF
            jdiag = idiag - jc
            DO i = 1 , N
               IF ( i<jc ) THEN
                  temp = Rfac(idiag+i)
                  Rfac(idiag+i) = Rfac(jdiag+i)
                  Rfac(jdiag+i) = temp
               ENDIF
               temp = sval*Qfac(i,jc) + cval*Qfac(i,jcp)
               Qfac(i,jcp) = cval*Qfac(i,jc) - sval*Qfac(i,jcp)
               Qfac(i,jc) = temp
            ENDDO
            Iact(jc) = Iact(jcp)
            Resact(jc) = Resact(jcp)
            Vlam(jc) = Vlam(jcp)
            jc = jcp
            CYCLE
         ENDIF
         Nact = Nact - 1
         IF ( iflag==1 ) GOTO 200
         IF ( iflag==2 ) GOTO 300
         IF ( iflag==3 ) GOTO 700
         EXIT
      ENDDO
      END SUBROUTINE GETACT
