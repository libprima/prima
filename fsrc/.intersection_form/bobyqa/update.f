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
! Generated using the interform.m script by Zaikun ZHANG (www.zhangzk.net)
! on 11-Feb-2022.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine update(n, npt, bmat, zmat, ndim, vlag, beta, denom, kne&
     &w, w)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
      implicit real(kind(0.0D0)) (a - h, o - z)
      implicit integer(i - n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension bmat(ndim, *), zmat(npt, *), vlag(*), w(*)
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
      nptm = npt - n - 1
      ztest = zero
      do k = 1, npt
          do j = 1, nptm
              ztest = dmax1(ztest, dabs(zmat(k, j)))
          end do
      end do
      ztest = 1.0D-20 * ztest
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
! Zaikun 2019-08-15: JL is never used
!      JL=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 2, nptm
          if (dabs(zmat(knew, j)) > ztest) then
              temp = dsqrt(zmat(knew, 1)**2 + zmat(knew, j)**2)
              tempa = zmat(knew, 1) / temp
              tempb = zmat(knew, j) / temp
              do i = 1, npt
                  temp = tempa * zmat(i, 1) + tempb * zmat(i, j)
                  zmat(i, j) = tempa * zmat(i, j) - tempb * zmat(i, 1)
                  zmat(i, 1) = temp
              end do
          end if
          zmat(knew, j) = zero
      end do
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
      do i = 1, npt
          w(i) = zmat(knew, 1) * zmat(i, 1)
      end do
      alpha = w(knew)
      tau = vlag(knew)
      vlag(knew) = vlag(knew) - one
!
!     Complete the updating of ZMAT.
!
      temp = dsqrt(denom)
      tempb = zmat(knew, 1) / temp
      tempa = tau / temp
      do i = 1, npt
          zmat(i, 1) = tempa * zmat(i, 1) - tempb * vlag(i)
      end do
!
!     Finally, update the matrix BMAT.
!
      do j = 1, n
          jp = npt + j
          w(jp) = bmat(knew, j)
          tempa = (alpha * vlag(jp) - tau * w(jp)) / denom
          tempb = (-beta * w(jp) - tau * vlag(jp)) / denom
          do i = 1, jp
              bmat(i, j) = bmat(i, j) + tempa * vlag(i) + tempb * w(i)
              if (i > npt) bmat(jp, i - npt) = bmat(i, j)
          end do
      end do
      return
      end