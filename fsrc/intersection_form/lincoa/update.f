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
! on 26-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==update.f90  processed by SPAG 7.50RE at 00:12 on 26 May 2021
            SUBROUTINE UPDATE(N,Npt,Xpt,Bmat,Zmat,Idz,Ndim,Sp,Step,Kopt,&
     &Knew, Vlag,W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--UPDATE8
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER , INTENT(IN) :: N
            INTEGER , INTENT(IN) :: Npt
            REAL*8 , INTENT(IN) , DIMENSION(Npt,*) :: Xpt
            REAL*8 , INTENT(INOUT) , DIMENSION(Ndim,*) :: Bmat
            REAL*8 , INTENT(INOUT) , DIMENSION(Npt,*) :: Zmat
            INTEGER , INTENT(INOUT) :: Idz
            INTEGER , INTENT(IN) :: Ndim
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Sp
            REAL*8 , INTENT(IN) , DIMENSION(*) :: Step
            INTEGER , INTENT(IN) :: Kopt
            INTEGER , INTENT(INOUT) :: Knew
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: Vlag
            REAL*8 , INTENT(INOUT) , DIMENSION(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            REAL*8 :: alpha , beta , bsum , denabs , denmax , denom , di&
     &stsq , dx , half , hdiag , one , scala , scalb , sqrtdn , ssq , su&
     &m , tau , tausq , temp , tempa , tempb , zero
            INTEGER :: i , iflag , j , ja , jb , jl , jp , k , nptm
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
!       identical to the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is such that XPT(KOPT,.) is the current trust region centre.
!     KNEW on exit is usually positive, and then it is the index of an
!       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
!       It is set on entry either to its final value or to 0. In the latter
!       case, the final value of KNEW is chosen to maximize the denominator
!       of the matrix updating formula times a weighting factor.
!     VLAG and W are used for working space, the first NPT+N elements of
!       both of these vectors being required.
!
!     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
!       the ones that are suitable after the shift of the KNEW-th point to
!       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero
!       occurs if the calculation fails due to a zero denominator in the
!       updating formula, which should never happen.
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            zero = 0.0D0
            nptm = Npt - N - 1
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!
            DO k = 1 , Npt
               W(k) = Sp(Npt+k)*(half*Sp(Npt+k)+Sp(k))
               sum = zero
               DO j = 1 , N
                  sum = sum + Bmat(k,j)*Step(j)
               ENDDO
               Vlag(k) = sum
            ENDDO
            beta = zero
            DO k = 1 , nptm
               sum = zero
               DO i = 1 , Npt
                  sum = sum + Zmat(i,k)*W(i)
               ENDDO
               IF ( k<Idz ) THEN
                  beta = beta + sum*sum
                  sum = -sum
               ELSE
                  beta = beta - sum*sum
               ENDIF
               DO i = 1 , Npt
                  Vlag(i) = Vlag(i) + sum*Zmat(i,k)
               ENDDO
            ENDDO
            bsum = zero
            dx = zero
            ssq = zero
            DO j = 1 , N
               sum = zero
               DO i = 1 , Npt
                  sum = sum + W(i)*Bmat(i,j)
               ENDDO
               bsum = bsum + sum*Step(j)
               jp = Npt + j
               DO k = 1 , N
                  sum = sum + Bmat(jp,k)*Step(k)
               ENDDO
               Vlag(jp) = sum
               bsum = bsum + sum*Step(j)
               dx = dx + Step(j)*Xpt(Kopt,j)
               ssq = ssq + Step(j)**2
            ENDDO
            beta = dx*dx + ssq*(Sp(Kopt)+dx+dx+half*ssq) + beta - bsum
            Vlag(Kopt) = Vlag(Kopt) + one
!
!     If KNEW is zero initially, then pick the index of the interpolation
!       point to be deleted, by maximizing the absolute value of the
!       denominator of the updating formula times a weighting factor.
!
!
            IF ( Knew==0 ) THEN
               denmax = zero
               DO k = 1 , Npt
                  hdiag = zero
                  DO j = 1 , nptm
                     temp = one
                     IF ( j<Idz ) temp = -one
                     hdiag = hdiag + temp*Zmat(k,j)**2
                  ENDDO
                  denabs = DABS(beta*hdiag+Vlag(k)**2)
                  distsq = zero
                  DO j = 1 , N
                     distsq = distsq + (Xpt(k,j)-Xpt(Kopt,j))**2
                  ENDDO
                  temp = denabs*distsq*distsq
                  IF ( temp>denmax ) THEN
                     denmax = temp
                     Knew = k
                  ENDIF
               ENDDO
            ENDIF
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
            jl = 1
            IF ( nptm>=2 ) THEN
               DO j = 2 , nptm
                  IF ( j==Idz ) THEN
                     jl = Idz
                  ELSEIF ( Zmat(Knew,j)/=zero ) THEN
                     temp = DSQRT(Zmat(Knew,jl)**2+Zmat(Knew,j)**2)
                     tempa = Zmat(Knew,jl)/temp
                     tempb = Zmat(Knew,j)/temp
                     DO i = 1 , Npt
                        temp = tempa*Zmat(i,jl) + tempb*Zmat(i,j)
                        Zmat(i,j) = tempa*Zmat(i,j) - tempb*Zmat(i,jl)
                        Zmat(i,jl) = temp
                     ENDDO
                     Zmat(Knew,j) = zero
                  ENDIF
               ENDDO
            ENDIF
!
!     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
!       into W, and calculate the parameters of the updating formula.
!
            tempa = Zmat(Knew,1)
            IF ( Idz>=2 ) tempa = -tempa
            IF ( jl>1 ) tempb = Zmat(Knew,jl)
            DO i = 1 , Npt
               W(i) = tempa*Zmat(i,1)
               IF ( jl>1 ) W(i) = W(i) + tempb*Zmat(i,jl)
            ENDDO
            alpha = W(Knew)
            tau = Vlag(Knew)
            tausq = tau*tau
            denom = alpha*beta + tausq
            Vlag(Knew) = Vlag(Knew) - one
            IF ( denom==zero ) THEN
               Knew = 0
               GOTO 99999
            ENDIF
            sqrtdn = DSQRT(DABS(denom))
!
!     Complete the updating of ZMAT when there is only one nonzero element
!       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
!       the value of IDZ is going to be reduced.
!
            iflag = 0
            IF ( jl==1 ) THEN
               tempa = tau/sqrtdn
               tempb = Zmat(Knew,1)/sqrtdn
               DO i = 1 , Npt
                  Zmat(i,1) = tempa*Zmat(i,1) - tempb*Vlag(i)
               ENDDO
               IF ( denom<zero ) THEN
                  IF ( Idz==1 ) THEN
                     Idz = 2
                  ELSE
                     iflag = 1
                  ENDIF
               ENDIF
            ELSE
!
!     Complete the updating of ZMAT in the alternative case.
!
               ja = 1
               IF ( beta>=zero ) ja = jl
               jb = jl + 1 - ja
               temp = Zmat(Knew,jb)/denom
               tempa = temp*beta
               tempb = temp*tau
               temp = Zmat(Knew,ja)
               scala = one/DSQRT(DABS(beta)*temp*temp+tausq)
               scalb = scala*sqrtdn
               DO i = 1 , Npt
                  Zmat(i,ja) = scala*(tau*Zmat(i,ja)-temp*Vlag(i))
                  Zmat(i,jb) = scalb*(Zmat(i,jb)-tempa*W(i)-tempb*Vlag(i&
     &))
               ENDDO
               IF ( denom<=zero ) THEN
                  IF ( beta<zero ) THEN
                     Idz = Idz + 1
                  ELSE
                     iflag = 1
                  ENDIF
               ENDIF
            ENDIF
!
!     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
!       ZMAT^T factorization gains another positive element. Then exchange
!       the first and IDZ-th columns of ZMAT.
!
            IF ( iflag==1 ) THEN
! Zaikun 2020-06-28: I came here when reading the corresponding part of
! the NEWUOA code, which seems to have a bug. In NEWUOA, IDZ is redued
! only if IDZ >= 2, which is reasonable. Here there seems no such
! restriction. Why? Did Powell use a different definition for IDZ? In
! NEWUOA, IDZ is an intger used to represent the leading NPT sub-matrix
! of H, which is called OMEGA in the paper and represented in the code
! as
!
! OMEGA = sum_{K = 1}^{NPT-N-1} S_K*ZMAT(:, K)*ZMAT(:, K)',
! where S(1:IDZ-1) = -1 and S(IDZ:NPT-N-1) = 1.
!
! Indeed, theoretically, OMEGA is positive semidefinite, and S should be
! all positive. The negative entries of S result the rounding errors
! that cause OMEGA to lose the postive semidefiniteness. Therefore, in
! most cases, IDZ is small (e.g., IDZ=1, meaning that OMEGA has not lost
! the positive semidefiniteness), but it cannot be nonpositive in the
! NEWUOA code. Is it different in the LINCOA code??? To be studied.
! Unfortunately, Powell did not write a LINCOA paper!!!
!
! The BOBYQA code does not have this part --- it does not use IDZ at
! all. Why?

               Idz = Idz - 1
               DO i = 1 , Npt
                  temp = Zmat(i,1)
                  Zmat(i,1) = Zmat(i,Idz)
                  Zmat(i,Idz) = temp
               ENDDO
            ENDIF
!
!     Finally, update the matrix BMAT.
!
            DO j = 1 , N
               jp = Npt + j
               W(jp) = Bmat(Knew,j)
               tempa = (alpha*Vlag(jp)-tau*W(jp))/denom
               tempb = (-beta*W(jp)-tau*Vlag(jp))/denom
               DO i = 1 , jp
                  Bmat(i,j) = Bmat(i,j) + tempa*Vlag(i) + tempb*W(i)
                  IF ( i>Npt ) Bmat(jp,i-Npt) = Bmat(i,j)
               ENDDO
            ENDDO
      99999 END SUBROUTINE UPDATE