      module vlagbeta_mod

      implicit none
      private
      public :: vlagbeta


      contains

!      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
!     & vlag, beta, wcheck)
      subroutine vlagbeta(idz, kopt, bmat, zmat, xpt, xopt, d, vlag,    &
     & beta, wcheck, dsq, xoptsq)

      use consts_mod, only : RP, IK, ONE, HALF, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      integer(IK), intent(in) :: idz, kopt
      real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
      real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
      real(RP), intent(in) :: xopt(:)  ! XOPT(N)
      real(RP), intent(in) :: d(:)  ! D(N)
      real(RP), intent(out) :: vlag(:)  ! VLAG(NPT+N)
      real(RP), intent(out) :: beta
      real(RP), intent(out) :: wcheck(:)  ! WCHECK(NPT)

      integer(IK) :: k, j, n, npt
      real(RP) :: bw(size(bmat, 1)), bwvd
      real(RP) :: wz(size(zmat, 2)), wzsave(size(wz))
      real(RP) :: dx, dsq, xoptsq
      character(len = SRNLEN), parameter :: srname = 'VLAGBETA'
      

      ! Get and verify the sizes
      n = size(xpt, 1)
      npt = size(xpt, 2)

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          if (size(bmat, 1) /= n .or. size(bmat, 2) /= npt + n) then
              call errstop(srname, 'SIZE(BMAT) is invalid')
          end if
          if (size(zmat, 1) /= npt .or. size(zmat, 2) /= npt-n-1) then
              call errstop(srname, 'SIZE(ZMAT) is invalid')
          end if
          if (size(xopt) /= n) then
              call errstop(srname, 'SIZE(XOPT) /= N')
          end if
          if (size(d) /= n) then
              call errstop(srname, 'SIZE(D) /= N')
          end if
          if (size(vlag) /= n + npt) then
              call errstop(srname, 'SIZE(VLAG) /= N + NPT')
          end if
          if (size(wcheck) /= npt) then
              call errstop(srname, 'SIZE(WCHECK) /= NPT')
          end if
      end if
      

!----------------------------------------------------------------------!
      ! This is the one of the two places where WCHECK is calculated,
      ! the other being BIGDEN. 
      ! WCHECK contains the first NPT entries of (w-v) for the vectors 
      ! w and v defined in eq(4.10) and eq(4.24) of the NEWUOA paper,
      ! and also \hat{w} in eq(6.5) of 
      !
      ! M. J. D. Powell, Least Frobenius norm updating of quadratic
      ! models that satisfy interpolation conditions. Math. Program.,
      ! 100:183--215, 2004
      !
      ! WCHECK is used ONLY in CALQUAD, which evaluates the qudratic
      ! model. Indeed, CALQUAD can be implemented without WCHECK.
      wcheck = matmul(d, xpt)
      wcheck = wcheck*(HALF*wcheck + matmul(xopt, xpt))
!----------------------------------------------------------------------!

      vlag(1 : npt) = matmul(d, bmat(:, 1 : npt))

      wz = matmul(wcheck, zmat)
      wzsave = wz
      wz(1 : idz - 1) = -wz(1 : idz - 1)
      beta = -dot_product(wzsave, wz)
!----------------------------------------------------------------------!
      ! The following DO LOOP implements the update below. The results
      ! will not be identical due to the non-associativity of
      ! floating point arithmetic addition.
!-----!vlag(1 : npt) = vlag(1 : npt) + matmul(zmat, wz) !--------------!
      do k = 1, npt - n - 1
          vlag(1 : npt) = vlag(1 : npt) + wz(k)*zmat(:, k)
      end do
!----------------------------------------------------------------------!

      bw = matmul(bmat(:, 1 : npt), wcheck)
!----------------------------------------------------------------------!
      ! The following DO LOOP implements the update below. The results
      ! will not be identical due to the non-associativity of
      ! floating point arithmetic addition.
!-----!vlag(npt + 1: npt + n) = bw + matmul(d, bmat(:, npt + 1 : npt + n))
      vlag(npt + 1 : npt + n) = bw
      do k = 1, n
          vlag(npt + 1 : npt + n) = vlag(npt + 1 : npt + n) +           &
     &     bmat(k, npt + 1 : npt + n)*d(k)
      end do
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
      ! The following DO LOOP implements the dot product below. The
      ! results will not be identical due to the non-associativity of
      ! floating point arithmetic addition.
!-----!bwvd = dot_product(bw + vlag(npt+1 : npt+n), d) !---------------!
      bwvd = ZERO
      do j = 1, n
          bwvd = bwvd + bw(j)*d(j) + vlag(npt + j)*d(j)
      end do
!----------------------------------------------------------------------!

      dx = dot_product(d, xopt)

      !dsq = dot_product(d, d)
      !xoptsq = dot_product(xopt, xopt)

      beta = dx*dx + dsq*(xoptsq + dx + dx + HALF*dsq) + beta - bwvd
      vlag(kopt) = vlag(kopt) + ONE

      return

      end subroutine vlagbeta

      end module vlagbeta_mod
