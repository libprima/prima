!      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
!     & vlag, beta, wcheck)
      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
     & vlag, beta, wcheck, dsq, xoptsq)

      use consts, only : RP, IK, ONE, HALF, ZERO
      use lina
      implicit none

      integer(IK), intent(in) :: n, npt, idz, kopt
      real(RP), intent(in) :: bmat(n, npt+n), zmat(npt, npt - n - 1)
      real(RP), intent(in) :: xpt(n, npt), xopt(n), d(n)
      real(RP), intent(out) :: vlag(npt+n), beta, wcheck(npt)

      integer(IK) :: k, j
      real(RP) :: bw(n), bwvd, wz(npt - n - 1), wzsave(npt - n - 1) 
      real(RP) :: dx, dsq, xoptsq


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
      ! model. Indeed, we may calculate WCHECK internally in CALQUAD.
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
