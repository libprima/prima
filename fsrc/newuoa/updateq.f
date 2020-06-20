      subroutine updateq(n, npt, idz, knew, fqdiff, xptknew, bmatknew,  &
     & zmat, gq, hq, pq)

      use pdfomod, only : rp, zero
      implicit none

      integer, intent(in) :: n, npt, idz, knew

      real(kind = rp), intent(in) :: fqdiff, xptknew(n), bmatknew(n),   &
     & zmat(npt, npt - n - 1)
      real(kind = rp), intent(inout) :: gq(n), hq(n), pq(n) 

      integer :: i, ih, j, k
      real(kind = rp) :: temp

      ! Update the explicit part of second derivatives.
      ih = 0
      do i = 1, n
          temp = pq(knew)*xptknew(i)
          do j = 1, i
              ih = ih + 1
              hq(ih) = hq(ih) + temp*xptknew(j)
          end do
      end do
      
      ! Update the implicit part of second derivatives.
      pq(knew) = zero
      do j = 1, npt - n - 1
          temp = fqdiff*zmat(knew, j)
          if (j < idz) temp = -temp
          do k = 1, npt
              pq(k) = pq(k) + temp*zmat(k, j)
          end do
      end do

      ! Update the gradient.
      gq = gq + fqdiff*bmatknew

      return 

      end subroutine updateq
