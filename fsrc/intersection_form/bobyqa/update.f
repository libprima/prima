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
! on 06-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==update.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
            subroutine UPDATE(N, Npt, Bmat, Zmat, Ndim, Vlag, Beta, Deno&
     &m, Knew, W)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--UPDATE7
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer, intent(IN) :: N
            integer, intent(IN) :: Npt
            real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
            real*8, intent(INOUT), dimension(Npt, *) :: Zmat
            integer, intent(IN) :: Ndim
            real*8, intent(INOUT), dimension(*) :: Vlag
            real*8, intent(IN) :: Beta
            real*8, intent(IN) :: Denom
            integer, intent(IN) :: Knew
            real*8, intent(INOUT), dimension(*) :: W
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: alpha, one, tau, temp, tempa, tempb, zero, ztest
            integer :: i, j, jp, k, nptm
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            nptm = Npt - N - 1
            ztest = zero
            do k = 1, Npt
                do j = 1, nptm
                    ztest = DMAX1(ztest, DABS(Zmat(k, j)))
                end do
            end do
            ztest = 1.0D-20 * ztest
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: JL is never used
!      JL=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do j = 2, nptm
                if (DABS(Zmat(Knew, j)) > ztest) then
                    temp = DSQRT(Zmat(Knew, 1)**2 + Zmat(Knew, j)**2)
                    tempa = Zmat(Knew, 1) / temp
                    tempb = Zmat(Knew, j) / temp
                    do i = 1, Npt
                        temp = tempa * Zmat(i, 1) + tempb * Zmat(i, j)
                        Zmat(i, j) = tempa * Zmat(i, j) - tempb * Zmat(i&
     &, 1)
                        Zmat(i, 1) = temp
                    end do
                end if
                Zmat(Knew, j) = zero
            end do
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
            do i = 1, Npt
                W(i) = Zmat(Knew, 1) * Zmat(i, 1)
            end do
            alpha = W(Knew)
            tau = Vlag(Knew)
            Vlag(Knew) = Vlag(Knew) - one
!
!     Complete the updating of ZMAT.
!
            temp = DSQRT(Denom)
            tempb = Zmat(Knew, 1) / temp
            tempa = tau / temp
            do i = 1, Npt
                Zmat(i, 1) = tempa * Zmat(i, 1) - tempb * Vlag(i)
            end do
!
!     Finally, update the matrix BMAT.
!
            do j = 1, N
                jp = Npt + j
                W(jp) = Bmat(Knew, j)
                tempa = (alpha * Vlag(jp) - tau * W(jp)) / Denom
                tempb = (-Beta * W(jp) - tau * Vlag(jp)) / Denom
                do i = 1, jp
                    Bmat(i, j) = Bmat(i, j) + tempa * Vlag(i) + tempb * &
     &W(i)
                    if (i > Npt) Bmat(jp, i - Npt) = Bmat(i, j)
                end do
            end do
            end subroutine UPDATE