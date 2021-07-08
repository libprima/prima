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
! on 08-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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