!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of calcfc.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 25-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==calcfc.f90  processed by SPAG 7.50RE at 00:16 on 26 May 2021
            subroutine CALCFC(N, M, X, F, Con)
            implicit none
!*--CALCFC5
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer :: N, NProb
            integer :: M
            real*8, intent(IN), dimension(*) :: X
            real*8, intent(OUT) :: F
            real*8, intent(OUT), dimension(*) :: Con
!*++
!*++ End of declarations rewritten by SPAG
!*++
            if (NPRob == 1) then
!
!     Test problem 1 (Simple quadratic)
!
                F = 10.0 * (X(1) + 1.0)**2 + X(2)**2
            elseif (NPRob == 2) then
!
!    Test problem 2 (2D unit circle calculation)
!
                F = X(1) * X(2)
                Con(1) = 1.0 - X(1)**2 - X(2)**2
            elseif (NPRob == 3) then
!
!     Test problem 3 (3D ellipsoid calculation)
!
                F = X(1) * X(2) * X(3)
                Con(1) = 1.0 - X(1)**2 - 2.0 * X(2)**2 - 3.0 * X(3)**2
            elseif (NPRob == 4) then
!
!     Test problem 4 (Weak Rosenbrock)
!
                F = (X(1)**2 - X(2))**2 + (1.0 + X(1))**2
            elseif (NPRob == 5) then
!
!     Test problem 5 (Intermediate Rosenbrock)
!
                F = 10.0 * (X(1)**2 - X(2))**2 + (1.0 + X(1))**2
            elseif (NPRob == 6) then
!
!     Test problem 6 (Equation (9.1.15) in Fletcher's book)
!
                F = -X(1) - X(2)
                Con(1) = X(2) - X(1)**2
                Con(2) = 1.0 - X(1)**2 - X(2)**2
            elseif (NPRob == 7) then
!
!     Test problem 7 (Equation (14.4.2) in Fletcher's book)
!
                F = X(3)
                Con(1) = 5.0 * X(1) - X(2) + X(3)
                Con(2) = X(3) - X(1)**2 - X(2)**2 - 4.0 * X(2)
                Con(3) = X(3) - 5.0 * X(1) - X(2)
            elseif (NPRob == 8) then
!
!     Test problem 8 (Rosen-Suzuki)
!
                F = X(1)**2 + X(2)**2 + 2.0 * X(3)**2 + X(4)**2 - 5.0 * &
     &X(1) - 5.0 * X(2) - 21.0 * X(3) + 7.0 * X(4)
                Con(1) = 8.0 - X(1)**2 - X(2)**2 - X(3)**2 - X(4)**2 - X&
     &(1) + X(2) - X(3) + X(4)
                Con(2) = 10.0 - X(1)**2 - 2.0 * X(2)**2 - X(3)**2 - 2.0 &
     &* X(4) **2 + X(1) + X(4)
                Con(3) = 5.0 - 2.0 * X(1)**2 - X(2)**2 - X(3)**2 - 2.0 *&
     & X(1) + X(2) + X(4)
            elseif (NPRob == 9) then
!
!     Test problem 9 (Hock and Schittkowski 100)
!
                F = (X(1) - 10.0)**2 + 5.0 * (X(2) - 12.0)**2 + X(3) **4&
     & + 3.0 * (X(4) - 11.0)**2 + 10.0 * X(5)**6 + 7.0 * X(6) **2 + X(7)&
     &**4 - 4.0 * X(6) * X(7) - 10.0 * X(6) - 8.0 * X(7)
                Con(1) = 127.0 - 2.0 * X(1)**2 - 3.0 * X(2)**4 - X(3) - &
     &4.0 * X(4) **2 - 5.0 * X(5)
                Con(2) = 282.0 - 7.0 * X(1) - 3.0 * X(2) - 10.0 * X(3)**&
     &2 - X(4) + X(5)
                Con(3) = 196.0 - 23.0 * X(1) - X(2)**2 - 6.0 * X(6)**2 +&
     & 8.0 * X(7)
                Con(4) = -4.0 * X(1)**2 - X(2)**2 + 3.0 * X(1) * X(2) - &
     &2.0 * X(3) **2 - 5.0 * X(6) + 11.0 * X(7)
            elseif (NPRob == 10) then
!
!     Test problem 10 (Hexagon area)
!
                F = -0.5 * (X(1) * X(4) - X(2) * X(3) + X(3) * X(9) - X(&
     &5) * X(9) + X(5) * X(8) - X(6) * X(7))
                Con(1) = 1.0 - X(3)**2 - X(4)**2
                Con(2) = 1.0 - X(9)**2
                Con(3) = 1.0 - X(5)**2 - X(6)**2
                Con(4) = 1.0 - X(1)**2 - (X(2) - X(9))**2
                Con(5) = 1.0 - (X(1) - X(5))**2 - (X(2) - X(6))**2
                Con(6) = 1.0 - (X(1) - X(7))**2 - (X(2) - X(8))**2
                Con(7) = 1.0 - (X(3) - X(5))**2 - (X(4) - X(6))**2
                Con(8) = 1.0 - (X(3) - X(7))**2 - (X(4) - X(8))**2
                Con(9) = 1.0 - X(7)**2 - (X(8) - X(9))**2
                Con(10) = X(1) * X(4) - X(2) * X(3)
                Con(11) = X(3) * X(9)
                Con(12) = -X(5) * X(9)
                Con(13) = X(5) * X(8) - X(6) * X(7)
                Con(14) = X(9)
            end if
            end subroutine CALCFC