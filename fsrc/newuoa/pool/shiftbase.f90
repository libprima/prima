      subroutine shiftbase(n, npt, idz, xopt, pq, bmat, zmat,gq,hq,xpt)

      implicit none
      integer, parameter :: dp = kind(0.0d0)
      real(kind = dp), parameter :: almost_infinity = huge(0.0d0)/2.0d0,&
     & half = 0.5d0, zero = 0.0d0

      integer, intent(in) :: idz, n, npt
      integer :: info

      real(kind = dp), intent(in) :: xopt(n), pq(npt)
      real(kind = dp), intent(inout) :: bmat(n, npt + n),               &
     & zmat(npt, npt - n - 1), gq(n), hq((n*(n + 1))/2), xpt(n, npt)

      integer :: i, ih, j, k
      real(kind = dp) :: summation, summationz, temp, tempq, vlag(n),   &
     & w(2*npt), xoptsq

          real(kind = dp) :: xpr(npt, n), bmar(npt + n, n)
          bmar = transpose(bmat)
          xpr = transpose(xpt)

      ! First make the changes to BMAT that do not depend on ZMAT.
      !xoptsq = dot_product(xopt, xopt)
      xoptsq = zero
      do i = 1, n
          xoptsq = xoptsq + xopt(i)**2
      end do
      tempq = 0.25d0 * xoptsq
      do k = 1, npt
          !summation = dot_product(xpr(k, :), xopt)
          summation = zero
          do i = 1, n
              summation = summation + xpr(k, i)*xopt(i)
          end do
          temp = pq(k)*summation
          summation = summation - half*xoptsq
          w(npt + k) = summation
          do i = 1, n
              gq(i) = gq(i) + temp*xpr(k, i)
              xpr(k, i) = xpr(k, i) - half*xopt(i)
              vlag(i) = bmar(k, i)
              w(i) = summation*xpr(k, i) + tempq*xopt(i)
              do j = 1, i
                  bmar(npt + i, j) = bmar(npt + i, j) + vlag(i)*w(j) +  &
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
                  summation = summation + w(i)*xpr(i, j)
              end do
              vlag(j) = summation
              if (k < idz) summation = -summation
              do i = 1, npt
                  bmar(i, j) = bmar(i, j) + summation*zmat(i, k)
              end do
          end do
          do i = 1, n
              temp = vlag(i)
              if (k < idz) temp = -temp
              do j = 1, i
                  bmar(npt + i, j) = bmar(npt + i, j) + temp*vlag(j)
              end do
          end do
      end do

      do j = 1, n
          do i = 1, j
              bmar(npt + i, j) = bmar(npt + j, i)
          end do
      end do

      ! The following instructions complete the shift of XBASE,
      ! including the changes to the parameters of the quadratic model.
      ih = 0
      do j = 1, n
          w(j) = zero
          do k = 1, npt
              w(j) = w(j) + pq(k)*xpr(k, j)
              xpr(k, j) = xpr(k, j) - half*xopt(j)
          end do
          do i = 1, j
              ih = ih + 1
              if (i < j) gq(j) = gq(j) + hq(ih)*xopt(i)
              gq(i) = gq(i) + hq(ih)*xopt(j)
              hq(ih) = hq(ih) + w(i)*xopt(j) + xopt(i)*w(j)
          end do
      end do

      xpt = transpose(xpr)
      bmat = transpose(bmar)

      info = -100

      end subroutine shiftbase
