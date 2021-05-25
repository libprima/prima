!*==aa0001.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.
!
      IMPLICIT NONE
!*--AA00018
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      INTEGER :: i , iprint , maxfun , n
      REAL*8(R8KIND) :: rhobeg , rhoend
      REAL*8(R8KIND) , DIMENSION(10000) :: w
      REAL*8(R8KIND) , DIMENSION(10) :: x
!*++
!*++ End of declarations rewritten by SPAG
!*++
      iprint = 2
      maxfun = 5000
      rhoend = 1.0D-8
      DO n = 2 , 8 , 2
         DO i = 1 , n
            x(i) = DFLOAT(i)/DFLOAT(n+1)
         ENDDO
         rhobeg = 0.2D0*x(1)
         PRINT 99001 , n
99001    FORMAT (//5X,'******************'/5X,'Results with N =',I2,/5X,&
     &           '******************')
         CALL UOBYQA(n,x,rhobeg,rhoend,iprint,maxfun,w)
      ENDDO
      END PROGRAM AA0001
