!*==calfun.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
      subroutine CALFUN(N, X, F)
      implicit none
!*--CALFUN5
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      integer, intent(IN) :: N
      real*8, intent(IN), dimension(*) :: X
      real*8, intent(INOUT) :: F
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      integer :: i, j
      real*8 :: temp
!*++
!*++ End of declarations rewritten by SPAG
!*++
      F = 0.0D0
      do i = 4, N, 2
          do j = 2, i - 2, 2
              temp = (X(i - 1) - X(j - 1))**2 + (X(i) - X(j))**2
              temp = DMAX1(temp, 1.0D-6)
              F = F + 1.0D0 / DSQRT(temp)
          end do
      end do
      end subroutine CALFUN
