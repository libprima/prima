      module calquad_mod

      implicit none
      private
      public calquad


      contains

      !subroutine calquad(vquad, d, x, xpt, gq, hq, pq)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CALQUAD is the only place where WCHECK is used.
      ! CALQUAD can be implemented without WCHECK. For the moment, to 
      ! produce exactly the same result as Powell's code, we use WCHECK.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calquad(vquad, d, x, xpt, gq, hq, pq, wcheck)
      ! CALQUAD calculates VQUAD = Q(X + D) - Q(X), where Q is the
      ! quadratic function defined by (GQ, HQ, PQ).

      use consts_mod, only : RP, IK, HALF, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none


      real(RP), intent(in) :: d(:)  ! D(N)
      real(RP), intent(in) :: x(:)  ! X(N)
      real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
      real(RP), intent(in) :: gq(:)  ! GQ(N)
      real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
      real(RP), intent(in) :: pq(:)  ! PQ(NPT)
      real(RP), intent(out) :: vquad 

      integer(IK) :: i, ih, j, k, n, npt
      real(RP) :: s(size(x)), temp!,sd
      real(RP) :: wcheck(size(pq))
      character(len = SRNLEN), parameter :: srname = 'CALQUAD'


      ! Get and verify the sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(d, n)
          call verisize(x, n)
          call verisize(gq, n)
          call verisize(hq, n, n)
          call verisize(pq, npt)
      end if


      s = x + d  

      ! First order term and explicit second order term
      vquad = ZERO
      ih = 0
      do j = 1, n
          vquad = vquad + d(j)*gq(j)
          do i = 1, j
              ih = int(ih + 1, kind(ih))
              temp = d(i)*s(j) + d(j)*x(i)
              if (i == j) then 
                  temp = HALF*temp
              end if
              vquad = vquad + temp*hq(i, j)
          end do
      end do

!----------------------------------------------------------------------!
      ! WCHECK is the following vector in theory. We really do not need
      ! to pass it as an input. The only place that can make WCHECK
      ! different from this value is in BIGDEN, where WCHECK is
      ! calculated in another way (that is mathematically equivalent). 
!-----!wcheck = matmul(d, xpt) !---------------------------------------!
!-----!wcheck = wcheck*(half*wcheck + matmul(x, xpt)) !----------------!
!----------------------------------------------------------------------!

      ! Implicit second order term
      do k = 1, npt
          vquad = vquad + pq(k)*wcheck(k)
      end do

      !!!!!!!!!!!!!!!!!!IMPLEMENTATON WITHOUT WCHECK!!!!!!!!!!!!!!!!!!!!
!      s = d + x + x  ! Different from the above version.
!
!      ! 1st order term 
!      vquad = dot_product(d, gq)
!      ! implicit 2nd-order term
!      vquad = vquad + half*sum(pq*(matmul(s, xpt)*matmul(d, xpt)))
!      ! explicit 2nd-order term
!      vquad = vquad + half*dot_product(s, matmul(hq, d))
!
      !!!!!!!!!!!!!!IMPLEMENTATON WITHOUT WCHECK ENDS!!!!!!!!!!!!!!!!!!!
      return

      end subroutine calquad

      end module calquad_mod
