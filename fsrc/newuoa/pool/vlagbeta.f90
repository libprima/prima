!      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
!     & vlag, beta, wcheck)
      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
     & vlag, beta, wcheck, dsq, xoptsq)

      use consts, only : rp, one, half, zero
      implicit none
      
      integer, intent(in) :: n, npt, idz, kopt
      real(kind = rp), intent(in) :: bmat(n, npt+n), zmat(npt, npt-n-1),&
     & xpt(n, npt), xopt(n), d(n)
      real(kind = rp), intent(out) :: vlag(npt+n), beta, wcheck(npt)

      integer :: i, j, jp, k
      real(kind = rp) :: summationa, summationb, summation, bsummation, &
     & dx, dsq, xoptsq


          real(kind = rp) :: bmar(npt+n, n), xpr(npt, n)
          bmar = transpose(bmat)
          xpr = transpose(xpt)

      do k = 1, npt
          summationa = zero
          summationb = zero
          summation = zero
          do j = 1, n
              summationa = summationa + xpr(k, j)*d(j)
              summationb = summationb + xpr(k, j)*xopt(j)
              summation = summation + bmar(k, j)*d(j)
          end do
          wcheck(k) = summationa*(half*summationa + summationb)
          vlag(k) = summation
      end do

      beta = zero
      do k = 1, npt - n - 1
          summation = zero
          do i = 1, npt
              summation = summation + zmat(i, k)*wcheck(i)
          end do
          if (k < idz) then
              beta = beta + summation*summation
              summation = -summation
          else
              beta = beta - summation*summation
          end if
          do i = 1, npt
              vlag(i) = vlag(i) + summation*zmat(i, k)
          end do
      end do

      bsummation = zero
      dx = zero
      do j = 1, n
          summation = zero
          do i = 1, npt
              summation = summation + wcheck(i)*bmar(i, j)
          end do
          bsummation = bsummation + summation*d(j)
          jp = npt + j
          do k = 1, n
              summation = summation + bmar(jp, k)*d(k)
          end do
          vlag(jp) = summation
          bsummation = bsummation + summation*d(j)
          dx = dx + d(j)*xopt(j)
      end do

              
      beta = dx*dx + dsq*(xoptsq + dx + dx + half*dsq) + beta-bsummation
      vlag(kopt) = vlag(kopt) + one

      end subroutine vlagbeta
