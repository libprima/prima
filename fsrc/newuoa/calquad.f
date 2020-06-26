      !subroutine calquad(vquad, d, x, xpt, gq, hq, pq, n, npt)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CALQUAD can be implemented without WCHECK. For the moment, to 
      ! produce exactly the same result as Powell's code, we use WCHECK.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calquad(vquad, d, x, xpt, gq, hq, pq, n, npt, wcheck)
      ! CALQUAD calculates VQUAD = Q(X + D) - Q(X), where Q is the
      ! quadratic function defined by (GQ, HQ, PQ).

      use consts, only : rp, half, zero
      implicit none

      integer, intent(in) :: n, npt

      real(kind = rp), intent(in) :: d(n), x(n), xpt(npt, n), gq(n),    &
     & hq((n*(n+1))/2), pq(npt)
      real(kind = rp), intent(out) :: vquad 

      integer :: i, ih, j, k 
      real(kind = rp) :: s(n), temp!,sd
      real(kind = rp) :: wcheck(npt)

      s = x + d  ! Does NOT apply to the version below

      vquad = zero

      ih = 0
      do j = 1, n
          vquad = vquad + d(j)*gq(j)
          do i = 1, j
              ih = ih + 1
              temp = d(i)*s(j) + d(j)*x(i)
              if (i == j) then 
                  temp = half*temp
              end if
              vquad = vquad + temp*hq(ih)
          end do
      end do

      !wcheck = matmul(xpt, d)
      !wcheck = wcheck*(half*wcheck + matmul(xpt, x))

      do k = 1, npt
          vquad = vquad + pq(k)*wcheck(k)
      end do

      !!!!!!!!!!!!!!!!!!IMPLEMENTATON WITHOUT WCHECK!!!!!!!!!!!!!!!!!!!!
!      s = d + x + x  ! Different from the above version.
!
!      ! 1st order term plus implicit 2nd-order term
!      vquad = dot_product(d,gq)+half*sum(pq*matmul(xpt,s)*matmul(xpt,d)) 
!
!      ! explicit 2nd-order term
!      ih = 0 
!      do i = 1, n
!          do j = 1, i 
!              if (i == j) then
!                  sd = s(i)*d(i)
!              else
!                  sd = s(i)*d(j) + s(j)*d(i)
!              end if
!              ih = ih + 1
!              vquad = vquad + half * hq(ih) * sd 
!          end do
!      end do
      !!!!!!!!!!!!!!IMPLEMENTATON WITHOUT WCHECK ENDS!!!!!!!!!!!!!!!!!!!
      return

      end subroutine calquad
