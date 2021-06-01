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
! on 01-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==aa0001.f90  processed by SPAG 7.50RE at 17:58 on 25 May 2021
!
!     Test problem for BOBYQA, the objective function being the sum of
!     the reciprocals of all pairwise distances between the points P_I,
!     I=1,2,...,M in two dimensions, where M=N/2 and where the components
!     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
!     defines the M points P_I. The initial X gives equally spaced points
!     on a circle. Four different choices of the pairs (N,NPT) are tried,
!     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
!     minimum that is not global occurs in both the N=10 cases. The details
!     of the results are highly sensitive to computer rounding errors. The
!     choice IPRINT=2 provides the current X and optimal F so far whenever
!     RHO is reduced. The bound constraints of the problem require every
!     component of X to be in the interval [-1,1].
      program AA000119
!
      implicit none
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      real*8 :: bdl, bdu, rhobeg, rhoend, temp, twopi, f, ftarget
      integer :: i, iprint, j, jcase, m, maxfun, n, npt, info
      real*8, dimension(500000) :: w
      real*8, dimension(100) :: x, xl, xu
!*++
!*++ End of declarations rewritten by SPAG
!*++
      twopi = 8.0D0 * DATAN(1.0D0)
      bdl = -1.0D0
      bdu = 1.0D0
      iprint = 2
      maxfun = 500000
      rhobeg = 1.0D-1
      rhoend = 1.0D-6
      m = 5
      do
          n = 2 * m
          do i = 1, n
              xl(i) = bdl
              xu(i) = bdu
          end do
          do jcase = 1, 2
              npt = n + 6
              if (jcase == 2) npt = 2 * n + 1
              print 99001, m, n, npt
      99001 format(//5X, '2D output with M =', I4, ', N =', I4, ' and NP&
     &T =', I4)
              do j = 1, m
                  temp = DFLOAT(j) * twopi / DFLOAT(m)
                  x(2 * j - 1) = DCOS(temp)
                  x(2 * j) = DSIN(temp)
              end do
              call BOBYQA(n, npt, x, xl, xu, rhobeg, rhoend, iprint, max&
     &fun, w, f, info, ftarget)
          end do
          m = m + m
          if (m > 10) exit
      end do
      end program AA000119