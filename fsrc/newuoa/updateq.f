      subroutine updateq(n, npt, idz, knew, fqdiff, xptknew, bmatknew,  &
     & zmat, gq, hq, pq)

      use consts, only : rp, zero
      use lina
      implicit none

      integer, intent(in) :: n, npt, idz, knew

      real(kind = rp), intent(in) :: fqdiff, xptknew(n), bmatknew(n),   &
     & zmat(npt, npt - n - 1)
      real(kind = rp), intent(inout) :: gq(n), hq((n*(n+1))/2), pq(npt) 

      integer :: i, ih, j
      real(kind = rp) :: fqdz(npt - n - 1)

      ! Update the explicit part of second derivatives. It adds 
      ! PQKNEW*outprod(XPTKNEW, XPTKNEW) to the explicit HESSIAN.
      ih = 0
      do i = 1, n
          hq(ih+1:ih+i)=hq(ih+1:ih+i)+(pq(knew)*xptknew(i))*xptknew(1:i)
          ih = ih + i
      end do
      
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
