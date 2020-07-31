      module vlagbeta_mod

      implicit none
      private
      public :: vlagbeta


      contains

      subroutine vlagbeta(idz, kopt, bmat, d, xopt, xpt, zmat,beta,vlag)

      use consts_mod, only : RP, IK, ONE, HALF, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      ! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
      real(RP), intent(in) :: d(:)  ! D(N)
      real(RP), intent(in) :: xopt(:)  ! XOPT(N)
      real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

      ! Outputs
      real(RP), intent(out) :: beta
      real(RP), intent(out) :: vlag(:)  ! VLAG(NPT+N)

      ! Intermediate variables
      integer(IK) :: n, npt
      real(RP) :: bw(size(bmat, 1)), bwvd
      real(RP) :: wcheck(size(zmat, 1)) 
      real(RP) :: wz(size(zmat, 2)), wzsave(size(wz))
      real(RP) :: dx, dsq, xoptsq
      character(len = SRNLEN), parameter :: srname = 'VLAGBETA'
      

      ! Get and verify the sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(bmat, n, npt + n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(xopt, n)
          call verisize(d, n)
          call verisize(vlag, n + npt)
      end if
      

!----------------------------------------------------------------------!
      ! This is the one of the two places where WCHECK is calculated,
      ! the other one being BIGDEN (now removed). 
      ! WCHECK contains the first NPT entries of (w-v) for the vectors 
      ! w and v defined in eq(4.10) and eq(4.24) of the NEWUOA paper,
      ! and also \hat{w} in eq(6.5) of 
      !
      ! M. J. D. Powell, Least Frobenius norm updating of quadratic
      ! models that satisfy interpolation conditions. Math. Program.,
      ! 100:183--215, 2004
      !
      wcheck = matmul(d, xpt)
      wcheck = wcheck*(HALF*wcheck + matmul(xopt, xpt))
!----------------------------------------------------------------------!

      vlag(1 : npt) = matmul(d, bmat(:, 1 : npt))

      wz = matmul(wcheck, zmat)
      wzsave = wz
      wz(1 : idz - 1) = -wz(1 : idz - 1)
      beta = -dot_product(wzsave, wz)
!----------------------------------------------------------------------!
!-----!vlag(1 : npt) = vlag(1 : npt) + matmul(zmat, wz) !--------------!
       vlag(1 : npt) = Ax_plus_y(zmat, wz, vlag(1 : npt))
!----------------------------------------------------------------------!

      bw = matmul(bmat(:, 1 : npt), wcheck)
!----------------------------------------------------------------------!
!-----!vlag(npt + 1: npt + n) = bw + matmul(d, bmat(:, npt + 1 : npt + n))
      vlag(npt+1: npt+n) = xA_plus_y(bmat(:, npt+1 : npt+n), d, bw)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!-----!bwvd = dot_product(bw + vlag(npt+1 : npt+n), d) !---------------!
      bwvd = xpy_dot_z(bw, vlag(npt+1 : npt+n), d)
!----------------------------------------------------------------------!

      dx = dot_product(d, xopt)

      dsq = dot_product(d, d)
      xoptsq = dot_product(xopt, xopt)

      beta = dx*dx + dsq*(xoptsq + dx + dx + HALF*dsq) + beta - bwvd
      vlag(kopt) = vlag(kopt) + ONE

      return

      end subroutine vlagbeta

      end module vlagbeta_mod
