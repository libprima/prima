!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of calfun.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 01-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==calfun.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
            SUBROUTINE CALFUN(N,X,F)
            IMPLICIT NONE
!*--CALFUN5
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER , INTENT(IN) :: N
            REAL*8(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
            REAL*8(R8KIND) , INTENT(INOUT) :: F
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            INTEGER :: i , iw , j , np
            REAL*8(R8KIND) :: sum
            REAL*8(R8KIND) , DIMENSION(10,10) :: y
!*++
!*++ End of declarations rewritten by SPAG
!*++
            DO j = 1 , N
               y(1,j) = 1.0D0
               y(2,j) = 2.0D0*X(j) - 1.0D0
            ENDDO
            DO i = 2 , N
               DO j = 1 , N
                  y(i+1,j) = 2.0D0*y(2,j)*y(i,j) - y(i-1,j)
               ENDDO
            ENDDO
            F = 0.0D0
            np = N + 1
            iw = 1
            DO i = 1 , np
               sum = 0.0D0
               DO j = 1 , N
                  sum = sum + y(i,j)
               ENDDO
               sum = sum/DFLOAT(N)
               IF ( iw>0 ) sum = sum + 1.0/DFLOAT(i*i-2*i)
               iw = -iw
               F = F + sum*sum
            ENDDO
            END SUBROUTINE CALFUN