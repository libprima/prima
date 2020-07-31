      module shiftbase_mod

      implicit none
      private
      public :: shiftbase


      contains 

      subroutine shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt,xpt)
      ! SHIFTBASE shifts the base point to XBASE + XOPT and updates GQ,
      ! HQ, and BMAT accordingly. PQ and ZMAT remain the same after the
      ! shifting. See Section 7 of the NEWUOA paper.

      use consts_mod, only : RP, IK, ZERO, ONE, HALF, QUART 
      use consts_mod, only : DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      ! Inputs
      integer(IK), intent(in) :: idz
      real(RP), intent(in) :: pq(:)  ! PQ(NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

      ! In-outputs
      real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
      real(RP), intent(inout) :: gq(:)  ! GQ(N)
      real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)
      real(RP), intent(inout) :: xbase(:)  ! XBASE(N)
      real(RP), intent(inout) :: xopt(:)  ! XOPT(N)
      real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)

      ! Intermediate variables
      integer(IK) :: k, n, npt
      real(RP) :: sumz(size(zmat, 2)), vlag(size(xopt))
      real(RP) :: qxoptq, xoptsq, xpq(size(xopt)), bmatk(size(bmat, 1)) 
      REAL(RP) :: w1(size(pq)), w2(size(xopt))
      character(len = SRNLEN), parameter :: srname = 'SHIFTBASE' 


      ! Get and verify the sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))
      
      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(pq, npt)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(bmat, n, npt + n)
          call verisize(gq, n)
          call verisize(hq, n, n)
          call verisize(xopt, n)
          call verisize(xbase, n)
      end if

            
      xoptsq = dot_product(xopt, xopt)
      qxoptq = QUART * xoptsq

      !----------------------------------------------------------------!   
      ! The update for gq can indeed be done by the following 3 lines:
      !real(RP) :: hxopt(n)
      !call hessmul(n, npt, xpt, hq, pq, xopt, hxopt)
      !gq = gq + hxopt 
      gq = Ax_plus_y(xpt, pq*matmul(xopt, xpt), gq)
      gq = Ax_plus_y(hq, xopt, gq)
      !----------------------------------------------------------------!
      

      w1 = matmul(xopt, xpt) - HALF*xoptsq
      ! W1 equals MATMUL(XPT, XOPT) after XPT is updated as follows.
      xpt = xpt - HALF*spread(xopt, dim = 2, ncopies = npt)
      !do k = 1, npt
      !    xpt(:, k) = xpt(:, k) - HALF*xopt
      !end do

      ! Update HQ. It has to be done after the above revision to XPT!!!
      xpq = matmul(xpt, pq)
      
      !----------------------------------------------------------------!
      ! Implement R2UPDATE properly so that it ensures HQ is symmetric.
      call r2update(hq, ONE, xopt, xpq)
      !----------------------------------------------------------------!


      ! Make the changes to BMAT that do not depend on ZMAT.
      do k = 1, npt
          bmatk = bmat(:, k)
          w2 = w1(k)*xpt(:, k) + qxoptq*xopt
          ! Implement R2UPDATE properly so that it ensures 
          ! bmat(:, npt+1:npt+n) is symmetric.
          call r2update(bmat(:, npt+1 : npt+n), ONE, bmatk, w2)
      end do
      
      ! Then the revisions of BMAT that depend on ZMAT are calculated.
      sumz = sum(zmat, dim = 1)
      do k = 1, int(idz - 1, kind(k))
!----------------------------------------------------------------------!
!---------!vlag = qxoptq*sumz(k)*xopt + matmul(xpt, w1*zmat(:, k)) !---!
          vlag = Ax_plus_y(xpt, w1*zmat(:, k), qxoptq*sumz(k)*xopt)
!----------------------------------------------------------------------!
          call r1update(bmat(:, 1:npt), -ONE, vlag, zmat(:, k))
          ! Implement R1UPDATE properly so that it ensures 
          ! bmat(:, npt+1:npt+n) is symmetric.
          call r1update(bmat(:, npt+1 : npt+n), -ONE, vlag)
      end do
      do k = idz, int(npt - n - 1, kind(k))
!----------------------------------------------------------------------!
!---------!vlag = qxoptq*sumz(k)*xopt + matmul(xpt, w1*zmat(:, k)) !---!
          vlag = Ax_plus_y(xpt, w1*zmat(:, k), qxoptq*sumz(k)*xopt)
!----------------------------------------------------------------------!
          call r1update(bmat(:, 1:npt), ONE, vlag, zmat(:, k))
          ! Implement R1UPDATE properly so that it ensures 
          ! bmat(:, npt+1:npt+n) is symmetric.
          call r1update(bmat(:, npt+1 : npt+n), ONE, vlag)
      end do

      !----------------------------------------------------------------!
      ! The following instructions complete the shift of XBASE.
      ! Recall the we have already subtracted HALF*XOPT from XPT. 
      ! Therefore, overall, the new XPT is XPT - XOPT.
      xpt = xpt - HALF*spread(xopt, dim = 2, ncopies = npt)
      !do k = 1, npt
      !    xpt(:, k) = xpt(:, k) - HALF*xopt
      !end do

      xbase = xbase + xopt
      xopt = ZERO

      return 
      end subroutine shiftbase

      end module shiftbase_mod
