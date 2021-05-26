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
! on 26-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      program AA0001
!*==aa0001.f90  processed by SPAG 7.50RE at 23:21 on 25 May 2021
!     Calculate the tetrahedron of least volume that encloses the points
!       (XP(J),YP(J),ZP(J)), J=1,2,...,NP. Our method requires the origin
!       to be strictly inside the convex hull of these points. There are
!       twelve variables that define the four faces of each tetrahedron
!       that is considered. Each face has the form ALPHA*X+BETA*Y+GAMMA*Z=1,
!       the variables X(3K-2), X(3K-1) and X(3K) being the values of ALPHA,
!       BETA and GAMMA for the K-th face, K=1,2,3,4. Let the set T contain
!       all points in three dimensions that can be reached from the origin
!       without crossing a face. Because the volume of T may be infinite,
!       the objective function is the smaller of FMAX and the volume of T,
!       where FMAX is set to an upper bound on the final volume initially.
!       There are 4*NP linear constraints on the variables, namely that each
!       of the given points (XP(J),YP(J),ZP(J)) shall be in T. Let XS = min
!       XP(J), YS = min YP(J), ZS = min ZP(J) and SS = max XP(J)+YP(J)+ZP(J),
!       where J runs from 1 to NP. The initial values of the variables are
!       X(1)=1/XS, X(5)=1/YS, X(9)=1/ZS, X(2)=X(3)=X(4)=X(6)=X(7) =X(8)=0
!       and X(10)=X(11)=X(12)=1/SS, which satisfy the linear constraints,
!       and which provide the bound FMAX=(SS-XS-YS-ZS)**3/6. Other details
!       of the test calculation are given below, including the choice of
!       the data points (XP(J),YP(J),ZP(J)), J=1,2,...,NP. The smaller final
!       value of the objective function in the case NPT=35 shows that the
!       problem has local minima.
!
      implicit none
!*--AA000126
!*** Start of declarations inserted by SPAG
      real * 8 a, b, FMAx, one, pi, rhobeg, rhoend, ss, sumx, sumy, sumz&
     &, theta, two, w, x, xp, xs, yp, ys, zero
      real * 8 zp, zs, f, ftarget
      integer i, ia, iprint, iw, j, jcase, k, m, maxfun, n, np, npt, inf&
     &o
!*** End of declarations inserted by SPAG
      common FMAx
      dimension xp(50), yp(50), zp(50), a(12, 200), b(200), x(12), w(500&
     &000)
!
!     Set some constants.
!
      one = 1.0D0
      two = 2.0D0
      zero = 0.0D0
      pi = 4.0D0 * DATAN(one)
      ia = 12
      n = 12
!
!     Set the data points.
!
      np = 50
      sumx = zero
      sumy = zero
      sumz = zero
      do j = 1, np
          theta = DFLOAT(j - 1) * pi / DFLOAT(np - 1)
          xp(j) = DCOS(theta) * DCOS(two * theta)
          sumx = sumx + xp(j)
          yp(j) = DSIN(theta) * DCOS(two * theta)
          sumy = sumy + yp(j)
          zp(j) = DSIN(two * theta)
          sumz = sumz + zp(j)
      end do
      sumx = sumx / DFLOAT(np)
      sumy = sumy / DFLOAT(np)
      sumz = sumz / DFLOAT(np)
      do j = 1, np
          xp(j) = xp(j) - sumx
          yp(j) = yp(j) - sumy
          zp(j) = zp(j) - sumz
      end do
!
!     Set the linear constraints.
!
      m = 4 * np
      do k = 1, m
          b(k) = one
          do i = 1, n
              a(i, k) = zero
          end do
      end do
      do j = 1, np
          do i = 1, 4
              k = 4 * j + i - 4
              iw = 3 * i
              a(iw - 2, k) = xp(j)
              a(iw - 1, k) = yp(j)
              a(iw, k) = zp(j)
          end do
      end do
!
!     Set the initial vector of variables. The JCASE=1,6 loop gives six
!       different choices of NPT when LINCOA is called.
!
      xs = zero
      ys = zero
      zs = zero
      ss = zero
      do j = 1, np
          xs = DMIN1(xs, xp(j))
          ys = DMIN1(ys, yp(j))
          zs = DMIN1(zs, zp(j))
          ss = DMAX1(ss, xp(j) + yp(j) + zp(j))
      end do
      FMAx = (ss - xs - ys - zs)**3 / 6.0D0
      do jcase = 1, 6
          do i = 2, 8
              x(i) = zero
          end do
          x(1) = one / xs
          x(5) = one / ys
          x(9) = one / zs
          x(10) = one / ss
          x(11) = one / ss
          x(12) = one / ss
!
!     Call of LINCOA, which provides the printing given at the end of this
!       note.
!
          npt = 5 * jcase + 10
          rhobeg = 1.0D0
          rhoend = 1.0D-6
          iprint = 1
          maxfun = 10000
          print 99001, npt, rhoend
      99001 format(//4X, 'Output from LINCOA with NPT =', I4, ' and RHOE&
     &ND =', 1PD12.4)
          call LINCOA(n, npt, m, a, ia, b, x, rhobeg, rhoend, iprint, ma&
     &xfun, w, f, info, ftarget)
      end do
      end program AA0001