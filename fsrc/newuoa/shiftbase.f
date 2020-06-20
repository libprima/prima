      subroutine shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq,hq,xpt)

      use pdfomod, only : rp, zero, half
      implicit none

      integer, intent(in) :: idz, n, npt

      real(kind = rp), intent(in) :: xopt(n), pq(npt)
      real(kind = rp), intent(inout) :: bmat(npt + n, n),               &
     & zmat(npt, npt - n - 1), gq(n), hq((n*(n + 1))/2), xpt(npt, n)

      integer :: i, ih, j, k
      real(kind = rp) :: summation, summationz, temp, tempq, vlag(n),   &
     & w(2*npt), xoptsq

      ! First make the changes to BMAT that do not depend on ZMAT.
      !xoptsq = dot_product(xopt, xopt)
      xoptsq = zero
      do i = 1, n
          xoptsq = xoptsq + xopt(i)**2
      end do
      tempq = 0.25_rp * xoptsq
      do k = 1, npt
          !summation = dot_product(xpt(k, :), xopt)
          summation = zero
          do i = 1, n
              summation = summation + xpt(k, i)*xopt(i)
          end do
          temp = pq(k)*summation
          summation = summation - half*xoptsq
          w(npt + k) = summation
          do i = 1, n
              gq(i) = gq(i) + temp*xpt(k, i)
              xpt(k, i) = xpt(k, i) - half*xopt(i)
              vlag(i) = bmat(k, i)
              w(i) = summation*xpt(k, i) + tempq*xopt(i)
              do j = 1, i
                  bmat(npt + i, j) = bmat(npt + i, j) + vlag(i)*w(j) +  &
     &             w(i)*vlag(j)
              end do
          end do
      end do

      ! Then the revisions of BMAT that depend on ZMAT are calculated.
      do k = 1, npt - n - 1
          summationz = zero
          do i = 1, npt
              summationz = summationz + zmat(i, k)
              w(i) = w(npt + i)*zmat(i, k)
          end do
          do j = 1, n
              summation = tempq*summationz*xopt(j)
              do i = 1, npt
                  summation = summation + w(i)*xpt(i, j)
              end do
              vlag(j) = summation
              if (k < idz) summation = -summation
              do i = 1, npt
                  bmat(i, j) = bmat(i, j) + summation*zmat(i, k)
              end do
          end do
          do i = 1, n
              temp = vlag(i)
              if (k < idz) temp = -temp
              do j = 1, i
                  bmat(npt + i, j) = bmat(npt + i, j) + temp*vlag(j)
              end do
          end do
      end do

      do j = 1, n
          do i = 1, j
              bmat(npt + i, j) = bmat(npt + j, i)
          end do
      end do

      ! The following instructions complete the shift of XBASE,
      ! including the changes to the parameters of the quadratic model.
      ih = 0
      do j = 1, n
          w(j) = zero
          do k = 1, npt
              w(j) = w(j) + pq(k)*xpt(k, j)
              xpt(k, j) = xpt(k, j) - half*xopt(j)
          end do
          do i = 1, j
              ih = ih + 1
              if (i < j) gq(j) = gq(j) + hq(ih)*xopt(i)
              gq(i) = gq(i) + hq(ih)*xopt(j)
              hq(ih) = hq(ih) + w(i)*xopt(j) + xopt(i)*w(j)
          end do
      end do

      return 

      end subroutine shiftbase
