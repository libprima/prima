!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of lagmax.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 10-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==lagmax.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
            SUBROUTINE LAGMAX(N,G,H,Rho,D,V,Vmax)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--LAGMAX7
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER , INTENT(IN) :: N
            REAL*8 , INTENT(IN) , DIMENSION(*) :: G
            REAL*8 , INTENT(INOUT) , DIMENSION(N,*) :: H
            REAL*8 , INTENT(IN) :: Rho
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: D
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: V
            REAL*8 , INTENT(OUT) :: Vmax
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: dd , dhd , dlin , dsq , gd , gg , ghg , gnorm , ha&
     &lf , halfrt , hmax , one , ratio , scale , sum , sumv , temp , tem&
     &pa , tempb , tempc , tempd , tempv , vhg , vhv , vhw , vlin , vmu &
     &, vnorm , vsq , vv , wcos , whw , wsin , wsq , zero
            INTEGER :: i , j , k
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     RHO is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.
!
!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. RHO
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.
!
!     Preliminary calculations.
!
            half = 0.5D0
            halfrt = DSQRT(half)
            one = 1.0D0
            zero = 0.0D0
!
!     Pick V such that ||HV|| / ||V|| is large.
!
            hmax = zero
            DO i = 1 , N
               sum = zero
               DO j = 1 , N
                  H(j,i) = H(i,j)
                  sum = sum + H(i,j)**2
               ENDDO
               IF ( sum>hmax ) THEN
                  hmax = sum
                  k = i
               ENDIF
            ENDDO
            DO j = 1 , N
               V(j) = H(k,j)
            ENDDO
!
!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.
!
            vsq = zero
            vhv = zero
            dsq = zero
            DO i = 1 , N
               vsq = vsq + V(i)**2
               D(i) = zero
               DO j = 1 , N
                  D(i) = D(i) + H(i,j)*V(j)
               ENDDO
               vhv = vhv + V(i)*D(i)
               dsq = dsq + D(i)**2
            ENDDO
            IF ( vhv*vhv<=0.9999D0*dsq*vsq ) THEN
               temp = vhv/vsq
               wsq = zero
               DO i = 1 , N
                  D(i) = D(i) - temp*V(i)
                  wsq = wsq + D(i)**2
               ENDDO
               whw = zero
               ratio = DSQRT(wsq/vsq)
               DO i = 1 , N
                  temp = zero
                  DO j = 1 , N
                     temp = temp + H(i,j)*D(j)
                  ENDDO
                  whw = whw + temp*D(i)
                  V(i) = ratio*V(i)
               ENDDO
               vhv = ratio*ratio*vhv
               vhw = ratio*wsq
               temp = half*(whw-vhv)
               temp = temp + DSIGN(DSQRT(temp**2+vhw**2),whw+vhv)
               DO i = 1 , N
                  D(i) = vhw*V(i) + temp*D(i)
               ENDDO
            ENDIF
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
            gg = zero
            gd = zero
            dd = zero
            dhd = zero
            DO i = 1 , N
               gg = gg + G(i)**2
               gd = gd + G(i)*D(i)
               dd = dd + D(i)**2
               sum = zero
               DO j = 1 , N
                  sum = sum + H(i,j)*D(j)
               ENDDO
               dhd = dhd + sum*D(i)
            ENDDO
            temp = gd/gg
            vv = zero
            scale = DSIGN(Rho/DSQRT(dd),gd*dhd)
            DO i = 1 , N
               V(i) = D(i) - temp*G(i)
               vv = vv + V(i)**2
               D(i) = scale*D(i)
            ENDDO
            gnorm = DSQRT(gg)
            IF ( gnorm*dd<=0.5D-2*Rho*DABS(dhd) .OR. vv/dd<=1.0D-4 ) THE&
     &N
               Vmax = DABS(scale*(gd+half*scale*dhd))
               GOTO 99999
            ENDIF
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.
!
            ghg = zero
            vhg = zero
            vhv = zero
            DO i = 1 , N
               sum = zero
               sumv = zero
               DO j = 1 , N
                  sum = sum + H(i,j)*G(j)
                  sumv = sumv + H(i,j)*V(j)
               ENDDO
               ghg = ghg + sum*G(i)
               vhg = vhg + sumv*G(i)
               vhv = vhv + sumv*V(i)
            ENDDO
            vnorm = DSQRT(vv)
            ghg = ghg/gg
            vhg = vhg/(vnorm*gnorm)
            vhv = vhv/vv
            IF ( DABS(vhg)<=0.01D0*DMAX1(DABS(ghg),DABS(vhv)) ) THEN
               vmu = ghg - vhv
               wcos = one
               wsin = zero
            ELSE
               temp = half*(ghg-vhv)
               vmu = temp + DSIGN(DSQRT(temp**2+vhg**2),temp)
               temp = DSQRT(vmu**2+vhg**2)
               wcos = vmu/temp
               wsin = vhg/temp
            ENDIF
            tempa = wcos/gnorm
            tempb = wsin/vnorm
            tempc = wcos/vnorm
            tempd = wsin/gnorm
            DO i = 1 , N
               D(i) = tempa*G(i) + tempb*V(i)
               V(i) = tempc*V(i) - tempd*G(i)
            ENDDO
!
!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.
!
            dlin = wcos*gnorm/Rho
            vlin = -wsin*gnorm/Rho
            tempa = DABS(dlin) + half*DABS(vmu+vhv)
            tempb = DABS(vlin) + half*DABS(ghg-vmu)
            tempc = halfrt*(DABS(dlin)+DABS(vlin)) + 0.25D0*DABS(ghg+vhv&
     &)
            IF ( tempa>=tempb .AND. tempa>=tempc ) THEN
               tempd = DSIGN(Rho,dlin*(vmu+vhv))
               tempv = zero
            ELSEIF ( tempb>=tempc ) THEN
               tempd = zero
               tempv = DSIGN(Rho,vlin*(ghg-vmu))
            ELSE
               tempd = DSIGN(halfrt*Rho,dlin*(ghg+vhv))
               tempv = DSIGN(halfrt*Rho,vlin*(ghg+vhv))
            ENDIF
            DO i = 1 , N
               D(i) = tempd*D(i) + tempv*V(i)
            ENDDO
            Vmax = Rho*Rho*DMAX1(tempa,tempb,tempc)
      99999 END SUBROUTINE LAGMAX