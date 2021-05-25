!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of update.f90.
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


!*==update.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
            SUBROUTINE UPDATE(N,Npt,Bmat,Zmat,Ndim,Vlag,Beta,Denom,Knew,&
     &W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            USE F77KINDS
            IMPLICIT NONE
!*--UPDATE7
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER , INTENT(IN) :: N
            INTEGER , INTENT(IN) :: Npt
            REAL*8 , INTENT(INOUT) , DIMENSION(Ndim,*) :: Bmat
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Zmat
            INTEGER , INTENT(IN) :: Ndim
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Vlag
            REAL*8 , INTENT(IN) :: Beta
            REAL*8 , INTENT(IN) :: Denom
            INTEGER , INTENT(IN) :: Knew
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: alpha , one , tau , temp , tempa , tempb , zero , &
     &ztest
            INTEGER :: i , j , jp , k , nptm
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.
!
!     Set some constants.
!
            one = 1.0D0
            zero = 0.0D0
            nptm = Npt - N - 1
            ztest = zero
            DO k = 1 , Npt
               DO j = 1 , nptm
                  ztest = DMAX1(ztest,DABS(Zmat(k,j)))
               ENDDO
            ENDDO
            ztest = 1.0D-20*ztest
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: JL is never used
!      JL=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO j = 2 , nptm
               IF ( DABS(Zmat(Knew,j))>ztest ) THEN
                  temp = DSQRT(Zmat(Knew,1)**2+Zmat(Knew,j)**2)
                  tempa = Zmat(Knew,1)/temp
                  tempb = Zmat(Knew,j)/temp
                  DO i = 1 , Npt
                     temp = tempa*Zmat(i,1) + tempb*Zmat(i,j)
                     Zmat(i,j) = tempa*Zmat(i,j) - tempb*Zmat(i,1)
                     Zmat(i,1) = temp
                  ENDDO
               ENDIF
               Zmat(Knew,j) = zero
            ENDDO
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
            DO i = 1 , Npt
               W(i) = Zmat(Knew,1)*Zmat(i,1)
            ENDDO
            alpha = W(Knew)
            tau = Vlag(Knew)
            Vlag(Knew) = Vlag(Knew) - one
!
!     Complete the updating of ZMAT.
!
            temp = DSQRT(Denom)
            tempb = Zmat(Knew,1)/temp
            tempa = tau/temp
            DO i = 1 , Npt
               Zmat(i,1) = tempa*Zmat(i,1) - tempb*Vlag(i)
            ENDDO
!
!     Finally, update the matrix BMAT.
!
            DO j = 1 , N
               jp = Npt + j
               W(jp) = Bmat(Knew,j)
               tempa = (alpha*Vlag(jp)-tau*W(jp))/Denom
               tempb = (-Beta*W(jp)-tau*Vlag(jp))/Denom
               DO i = 1 , jp
                  Bmat(i,j) = Bmat(i,j) + tempa*Vlag(i) + tempb*W(i)
                  IF ( i>Npt ) Bmat(jp,i-Npt) = Bmat(i,j)
               ENDDO
            ENDDO
            END SUBROUTINE UPDATE