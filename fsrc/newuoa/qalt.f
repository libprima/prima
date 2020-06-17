      subroutine qalt(gq, hq, pq, fval, smat, zmat, n, npt, kopt, idz)
!     QALT calculates the alternative model, namely the model that
!     minimizes the F-norm of the Hessian subject to the interpolation
!     conditions. 
!     Note that SMAT = BMAT(1:NPT, 1:N)

          use pdfomod, only : rp, zero
          implicit none

          integer, intent(in) :: n, npt, kopt, idz
          real(kind = rp), intent(in) :: fval(npt), smat(npt, n),       &
     &    zmat(npt, npt - n - 1)
          real(kind = rp), intent(out) :: gq(n), hq((n*(n+1))/2),pq(npt)

          real(kind = rp) :: vlag(npt), vz(npt - n - 1)
          integer :: i, k 

!          vlag = fval - fval(kopt)
!          gq = matmul(vlag, smat)
!          hq = zero
!          vz = matmul(vlag, zmat)
!          vz(1 : idz - 1) = - vz(1 : idz - 1)
!          pq = matmul(zmat, vz)

          vlag = fval - fval(kopt)
          gq = zero
          do k = 1, npt
              gq = gq + vlag(k)*smat(k, :)
          end do

          hq = zero

          vz = zero
          do k = 1, npt
              vz = vz + vlag(k)*zmat(k, :)
          end do
          vz(1 : idz - 1) = -vz(1 : idz - 1)
          pq = zero
          do i = 1, npt - n - 1
               pq = pq + zmat(:, i)*vz(i)
          end do

      end subroutine qalt
