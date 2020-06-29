      subroutine updateq(n, npt, idz, knew, fqdiff, xptknew, bmatknew,  &
     & zmat, gq, hq, pq)

      use consts, only : rp, zero
      use lina
      implicit none

      integer, intent(in) :: n, npt, idz, knew

      real(kind = rp), intent(in) :: fqdiff, xptknew(n), bmatknew(n),   &
     & zmat(npt, npt - n - 1)
      real(kind = rp), intent(inout) :: gq(n), hq(n, n), pq(npt) 

      integer :: i, ih, j
      real(kind = rp) :: fqdz(npt - n - 1)


      !----------------------------------------------------------------!
      ! This update does NOT preserve symmetry. Symmetrization needed! 
      hq = hq + outprod(xptknew, pq(knew)*xptknew)
      do i = 1, n
          hq(i, 1:i-1) = hq(1:i-1, i)
      end do
      !---------- A probably better implementation: -------------------!
      ! This is better because it preserves symmetry even with floating
      ! point arithmetic.
!-----!hq = hq + pq(knew)*outprod(xptknew, xptknew) !------------------!
      !----------------------------------------------------------------!

      ! Update the implicit part of second derivatives.
      fqdz = fqdiff*zmat(knew, :)
      fqdz(1 : idz - 1) = -fqdz(1 : idz - 1)
      pq(knew) = zero
      !----------------------------------------------------------------!
      !!! The following DO LOOP implements the update given below, yet
      !!! the result will not be exactly the same due to the
      !!! non-associativity of floating point arithmetic addition.
      !!! In future versions, we will use MATMUL instead of DO LOOP.
!----!pq = pq + matmul(zmat, fqdz) !-----------------------------------!
      !----------------------------------------------------------------!
      do j = 1, npt - n - 1
          pq = pq + fqdz(j)*zmat(:, j)
      end do

      ! Update the gradient.
      gq = gq + fqdiff*bmatknew

      return 

      end subroutine updateq
