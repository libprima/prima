      module shiftbase_mod

      implicit none
      private
      public :: shiftbase


      contains 

      subroutine shiftbase(idz, xopt, pq, bmat, zmat, gq, hq, xpt)
      ! SHIFTBASE shifts the base point to XBASE + XOPT and updates GQ,
      ! HQ, and BMAT accordingly. PQ and ZMAT remain the same after the
      ! shifting. See Section 7 of the NEWUOA paper.

      use consts_mod, only : RP, IK, ONE, HALF, QUART 
      use consts_mod, only : DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      integer(IK), intent(in) :: idz

      real(RP), intent(in) :: xopt(:)  ! XOPT(N)
      real(RP), intent(in) :: pq(:)  ! PQ(NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
      real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
      real(RP), intent(inout) :: gq(:)  ! GQ(N)
      real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)
      real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)

      integer(IK) :: i, j, k, n, npt
      real(RP) :: sumz(size(zmat, 2)), vlag(size(xopt))
      real(RP) :: qxoptq, xoptsq, xpq(size(xopt)), bmatk(size(bmat, 1)) 
      REAL(RP) :: w1(size(pq)), w2(size(xopt)), w3(size(pq)) 
      character(len = SRNLEN), parameter :: srname = 'SHIFTBASE' 


      ! Get and verify the sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))
      
      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(xopt, n)
          call verisize(pq, npt)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(bmat, n, npt + n)
          call verisize(gq, n)
          call verisize(hq, n, n)
      end if

            
      xoptsq = dot_product(xopt, xopt)
      qxoptq = QUART * xoptsq

      !----------------------------------------------------------------!   
      ! The update for gq can indeed be done by the following 3 lines:
      !real(RP) :: hxopt(n)
      !call hessmul(n, npt, xpt, hq, pq, xopt, hxopt)
      !gq = gq + hxopt 
      !----------------------------------------------------------------!
      do k = 1, npt
      ! Update of GQ due to the implicit part of the HESSIAN:
      ! GQ = GQ + (\sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)') * XOPT
          gq = gq + pq(k)*dot_product(xpt(:, k), xopt)*xpt(:, k)
      end do
      do j = 1, n
          gq = gq + hq(:, j)*xopt(j)
      end do
      
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

!==============================POWELL'S SCHEME=========================!
!      ! Make the changes to BMAT that do not depend on ZMAT.
!      do k = 1, npt
!          bmatk = bmat(:, k)
!          w2 = w1(k)*xpt(:, k) + qxoptq*xopt
!          ! The following DO LOOP performs the update below, but only
!          ! calculates the upper triangular part since it is symmetric.
!!----------------------------------------------------------------------!
!!          bmat(:, npt + 1 : npt + n) = bmat(:, npt + 1 : npt + n) +     &
!!     &     outprod(w2, bmatk) + outprod(bmatk, w2)
!!----------------------------------------------------------------------!
!          do i = 1, n
!              bmat(1 : i, npt + i) = bmat(1 : i, npt + i) +             &
!     &         bmatk(i)*w2(1 : i) + w2(i)*bmatk(1 : i)
!          end do
!      end do
!
!      ! Then the revisions of BMAT that depend on ZMAT are calculated.
!      sumz = sum(zmat, dim = 1)
!      do k = 1, npt - n - 1
!          ! The following W3 and DO LOOP indeed defines the VLAG below.
!          ! The results are not identical due to the non-associtivity 
!          ! of floating point arithmetic addition.
!!----------------------------------------------------------------------!
!!          vlag = qxoptq*sumz(k)*xopt + matmul(xpt, w1*zmat(:, k))
!!----------------------------------------------------------------------!
!          w3 = w1*zmat(:, k)
!          do j = 1, n
!              vlag(j) = qxoptq*sumz(k)*xopt(j)
!              do i = 1, npt
!                  vlag(j) = vlag(j) + w3(i)*xpt(j, i)
!              end do
!          end do
!
!          if (k < idz) then
!               bmat(:, 1:npt) = bmat(:, 1:npt) - outprod(vlag, zmat(:,k))
!          else
!               bmat(:, 1:npt) = bmat(:, 1:npt) + outprod(vlag, zmat(:,k))
!          end if
!
!          if (k < idz) then 
!          ! The following DO LOOP performs the update below, but only
!          ! calculates the upper triangular part since it is symmetric.
!!----------------------------------------------------------------------!          
!!              bmat(:, npt + 1 : npt + n) = bmat(:, npt + 1 : npt + n) - &
!!     &         outprod(vlag, vlag)
!!----------------------------------------------------------------------!          
!              do i = 1, n
!                  bmat(1:i, npt+i) = bmat(1:i, npt+i) - vlag(i)*vlag(1:i)
!              end do
!          else
!          ! The following DO LOOP performs the update below, but only
!          ! calculates the upper triangular part since it is symmetric.
!!----------------------------------------------------------------------!
!!              bmat(:, npt + 1 : npt + n) = bmat(:, npt + 1 : npt + n) + &
!!     &         outprod(vlag, vlag)
!!----------------------------------------------------------------------!
!              do i = 1, n
!                  bmat(1:i, npt+i) = bmat(1:i, npt+i) + vlag(i)*vlag(1:i)
!              end do
!          end if
!      end do
!=========================POWELL'S SCHEME ENDS=========================!


!!!!!!!!!!!!!!!!!!!!!A COMPACTER YET NOT SO ECONOMIC SCHEME!!!!!!!!!!!!!
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
          ! The following W3 and DO LOOP indeed defines the VLAG below.
          ! The results are not identical due to the non-associtivity 
          ! of floating point arithmetic addition.
!----------------------------------------------------------------------!
          !vlag = qxoptq*sumz(k)*xopt + matmul(xpt, w1*zmat(:, k))
!----------------------------------------------------------------------!
          w3 = w1*zmat(:, k)
          vlag = qxoptq*sumz(k)*xopt
          do i = 1, npt
              vlag = vlag + w3(i)*xpt(:, i)
          end do
          call r1update(bmat(:, 1:npt), -ONE, vlag, zmat(:, k))
          ! Implement R1UPDATE properly so that it ensures 
          ! bmat(:, npt+1:npt+n) is symmetric.
          call r1update(bmat(:, npt+1 : npt+n), -ONE, vlag)
      end do
      do k = idz, int(npt - n - 1, kind(k))
          ! The following W3 and DO LOOP indeed defines the VLAG below.
          ! The results are not identical due to the non-associtivity 
          ! of floating point arithmetic addition.
!----------------------------------------------------------------------!
          !vlag = qxoptq*sumz(k)*xopt + matmul(xpt, w1*zmat(:, k))
!----------------------------------------------------------------------!
          w3 = w1*zmat(:, k)
          vlag = qxoptq*sumz(k)*xopt
          do i = 1, npt
              vlag = vlag + w3(i)*xpt(:, i)
          end do
          call r1update(bmat(:, 1:npt), ONE, vlag, zmat(:, k))
          ! Implement R1UPDATE properly so that it ensures 
          ! bmat(:, npt+1:npt+n) is symmetric.
          call r1update(bmat(:, npt+1 : npt+n), ONE, vlag)
      end do
!!!!!!!!!!!!!!!!!!!!!COMPACT SCHEME ENDS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !----------------------------------------------------------------!
      ! The following instructions complete the shift of XBASE.
      ! Recall the we have already subtracted HALF*XOPT from XPT. 
      ! Therefore, overall, the new XPT is XPT - XOPT.
      xpt = xpt - HALF*spread(xopt, dim = 2, ncopies = npt)
      !do k = 1, npt
      !    xpt(:, k) = xpt(:, k) - HALF*xopt
      !end do

      return 

      end subroutine shiftbase

      end module shiftbase_mod
