!      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
!     & vlag, beta, wcheck)
      subroutine vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,  &
     & vlag, beta, wcheck, dsq, xoptsq)

      use pdfomod, only : rp, one, half, zero
      implicit none
      
      integer, intent(in) :: n, npt, idz, kopt
      real(kind = rp), intent(in) :: bmat(npt+n, n), zmat(npt, npt-n-1),&
     & xpt(npt, n), xopt(n), d(n)
      real(kind = rp), intent(out) :: vlag(npt+n), beta, wcheck(npt)

      integer :: i, j, jp, k
      real(kind = rp) :: summationa, summationb, summation, bsummation, &
     & dx, dsq, xoptsq

      do k = 1, npt
          summationa = zero
          summationb = zero
          summation = zero
          do j = 1, n
              summationa = summationa + xpt(k, j)*d(j)
              summationb = summationb + xpt(k, j)*xopt(j)
              summation = summation + bmat(k, j)*d(j)
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
              summation = summation + wcheck(i)*bmat(i, j)
          end do
          bsummation = bsummation + summation*d(j)
          jp = npt + j
          do k = 1, n
              summation = summation + bmat(jp, k)*d(k)
          end do
          vlag(jp) = summation
          bsummation = bsummation + summation*d(j)
          dx = dx + d(j)*xopt(j)
      end do
       
      !real(kind = rp) :: wz(npt - n - 1), wb(n), bd(n)
      !wcheck = matmul(bmat, d)*(half*matmul(xpt, d)+matmul(xpt, xopt))
      !vlag(1 : npt) = matmul(bmat(1 : npt, :), d)
      !wz = matmul(wcheck, zmat)
      !beta = sum(wz(1 : idz - 1)**2) - sum(wz(idz : npt -n -1)**2)
      !wz(1 : idz - 1) = -wz(1 : idz - 1)
      !vlag(1 : npt) = vlag(1 : npt) + matmul(zmat, wz(1:npt - n - 1))
      !wb = matmul(wcheck, bmat(1 : npt, :)
      !bsummation = dot_product(wb, d)
      !bd = matmul(bmat(npt + 1 : npt + n, :), d)
      !summation = summation + matmul(bmat(npt + 1 : npt + n, :), d)
      !vlag(npt+1 : npt+n) = wb + bd 
      !bsummation = dot_product(2*wb + bd, d) !????
      !dx = dot_product(d, xopt)
      !dsq = dot_product(d, d)
      !xoptsq = dot_product(xopt, xopt)
        
      beta = dx*dx + dsq*(xoptsq + dx + dx + half*dsq) + beta-bsummation
      vlag(kopt) = vlag(kopt) + one

      end subroutine vlagbeta
