!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of main.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 16-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      program AA0001
!*==aa0001.f90  processed by SPAG 7.50RE at 00:16 on 26 May 2021
      implicit none
!*--AA00015
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      integer :: i, icase, iprint, m, maxfun, n, NPROB, info
      integer, dimension(51) :: iact
      real*8 :: rhobeg, rhoend, temp, tempa, tempb, tempc, tempd, f, fta&
     &rget, resmax
      real*8, dimension(3000) :: w, con
      real*8, dimension(10) :: x, xopt
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!     Main program of test problems in Report DAMTP 1992/NA5.
!------------------------------------------------------------------------------
      do NPRob = 1, 10
          if (NPRob == 1) then
!
!     Minimization of a simple quadratic function of two variables.
!
              print 99001
99001 format(/7X, 'Output from test problem 1 (Simple quadratic)')
              n = 2
              m = 0
              xopt(1) = -1.0
              xopt(2) = 0.0
          elseif (NPRob == 2) then
!
!     Easy two dimensional minimization in unit circle.
!
              print 99002
99002 format(/7X, 'Output from test problem 2 (2D unit circle ', 'calcul&
     &ation)')
              n = 2
              m = 1
              xopt(1) = sqrt(0.5)
              xopt(2) = -xopt(1)
          elseif (NPRob == 3) then
!
!     Easy three dimensional minimization in ellipsoid.
!
              print 99003
99003 format(/7X, 'Output from test problem 3 (3D ellipsoid ', 'calculat&
     &ion)')
              n = 3
              m = 1
              xopt(1) = 1.0 / sqrt(3.0)
              xopt(2) = 1.0 / sqrt(6.0)
              xopt(3) = -1.0 / 3.0
          elseif (NPRob == 4) then
!
!     Weak version of Rosenbrock's problem.
!
              print 99004
99004 format(/7X, 'Output from test problem 4 (Weak Rosenbrock)')
              n = 2
              m = 0
              xopt(1) = -1.0
              xopt(2) = 1.0
          elseif (NPRob == 5) then
!
!     Intermediate version of Rosenbrock's problem.
!
              print 99005
99005 format(/7X, 'Output from test problem 5 (Intermediate ', 'Rosenbro&
     &ck)')
              n = 2
              m = 0
              xopt(1) = -1.0
              xopt(2) = 1.0
          elseif (NPRob == 6) then
!
!     This problem is taken from Fletcher's book Practical Methods of
!     Optimization and has the equation number (9.1.15).
!
              print 99006
99006 format(/7X, 'Output from test problem 6 (Equation ', '(9.1.15) in &
     &Fletcher)')
              n = 2
              m = 2
              xopt(1) = sqrt(0.5)
              xopt(2) = xopt(1)
          elseif (NPRob == 7) then
!
!     This problem is taken from Fletcher's book Practical Methods of
!     Optimization and has the equation number (14.4.2).
!
              print 99007
99007 format(/7X, 'Output from test problem 7 (Equation ', '(14.4.2) in &
     &Fletcher)')
              n = 3
              m = 3
              xopt(1) = 0.0
              xopt(2) = -3.0
              xopt(3) = -3.0
          elseif (NPRob == 8) then
!
!     This problem is taken from page 66 of Hock and Schittkowski's book Test
!     Examples for Nonlinear Programming Codes. It is their test problem Number
!     43, and has the name Rosen-Suzuki.
!
              print 99008
99008 format(/7X, 'Output from test problem 8 (Rosen-Suzuki)')
              n = 4
              m = 3
              xopt(1) = 0.0
              xopt(2) = 1.0
              xopt(3) = 2.0
              xopt(4) = -1.0
          elseif (NPRob == 9) then
!
!     This problem is taken from page 111 of Hock and Schittkowski's
!     book Test Examples for Nonlinear Programming Codes. It is their
!     test problem Number 100.
!
              print 99009
99009 format(/7X, 'Output from test problem 9 (Hock and ', 'Schittkowski&
     & 100)')
              n = 7
              m = 4
              xopt(1) = 2.330499
              xopt(2) = 1.951372
              xopt(3) = -0.4775414
              xopt(4) = 4.365726
              xopt(5) = -0.624487
              xopt(6) = 1.038131
              xopt(7) = 1.594227
          elseif (NPRob == 10) then
!
!     This problem is taken from page 415 of Luenberger's book Applied
!     Nonlinear Programming. It is to maximize the area of a hexagon of
!     unit diameter.
!
              print 99010
99010 format(/7X, 'Output from test problem 10 (Hexagon area)')
              n = 9
              m = 14
          end if
          do icase = 1, 2
              do i = 1, n
                  x(i) = 1.0
              end do
              rhobeg = 0.5
              rhoend = 0.001
              if (icase == 2) rhoend = 0.0001
              iprint = 1
              maxfun = 2000
              call COBYLA(n, m, x, rhobeg, rhoend, iprint, maxfun, w, ia&
     &ct, f, info, ftarget, resmax, con)
              if (NPRob == 10) then
                  tempa = x(1) + x(3) + x(5) + x(7)
                  tempb = x(2) + x(4) + x(6) + x(8)
                  tempc = 0.5 / sqrt(tempa * tempa + tempb * tempb)
                  tempd = tempc * sqrt(3.0)
                  xopt(1) = tempd * tempa + tempc * tempb
                  xopt(2) = tempd * tempb - tempc * tempa
                  xopt(3) = tempd * tempa - tempc * tempb
                  xopt(4) = tempd * tempb + tempc * tempa
                  do i = 1, 4
                      xopt(i + 4) = xopt(i)
                  end do
              end if
              temp = 0.0
              do i = 1, n
                  temp = temp + (x(i) - xopt(i))**2
              end do
              print 99011, sqrt(temp)
99011 format(/5X, 'Least squares error in variables =', 1PE16.6)
          end do
          print 99012
99012 format(2X, '----------------------------------------------', '----&
     &----------------')
      end do
      end program AA0001