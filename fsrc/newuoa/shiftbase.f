      subroutine shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq,hq,xpt)
      ! SHIFTBASE shifts the base point to XBASE + XOPT and updates GQ,
      ! HQ, and BMAT accordingly. PQ and ZMAT remain the same after the
      ! shifting. See Section 7 of the NEWUOA paper.

      use consts_mod, only : RP, IK, HALF, QUART
      use lina_mod
      implicit none

      integer(IK), intent(in) :: idz, n, npt

      real(RP), intent(in) :: xopt(n), pq(npt), zmat(npt,npt-n-1)
      real(RP), intent(inout) :: bmat(n, npt + n), gq(n), hq(n, n),     &
     & xpt(n, npt)

      integer(IK) :: i, j, k
      real(RP) :: sumz(npt-n-1), vlag(n), qxoptq, xoptsq, xpq(n)
      real(RP) :: bmatk(n), w1(npt), w2(n), w3(npt) 

            
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
      ! This update does NOT preserve symmetry. Symmetrization needed! 
      hq = hq + outprod(xpq, xopt) + outprod(xopt, xpq)
      do i = 1, n
          hq(i, 1:i-1) = hq(1:i-1, i)
      end do
      !---------- A probably better implementation: -------------------!
      ! This is better because it preserves symmetry even with floating
      ! point arithmetic.
!-----!hq = hq + ( outprod(xpq, xopt) + outprod(xopt, xpq) ) !---------!
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
          ! This is the only place where non-symmetry of
          ! BMAT(:, NPT+1:NPT+N) can come from. It is because of the
          ! non-associtivity of floating point arithmetic addition. To
          ! make BMAT(:, NPT+1:NPT+N) symmetric, use the following two
          ! lines to replace the code below. The only difference is the
          ! parenthsis around the two outter products. It is probably 
          ! a BETTERE implementation so we should take it in future
          ! versions. 
!----------------------------------------------------------------------!
!          bmat(:, npt + 1 : npt + n) = bmat(:, npt + 1 : npt + n) +     &
!     &     ( outprod(w2, bmatk) + outprod(bmatk, w2) )
!----------------------------------------------------------------------!
          bmat(:, npt + 1 : npt + n) = bmat(:, npt + 1 : npt + n) +     &
     &     outprod(w2, bmatk) + outprod(bmatk, w2)
      end do

      ! Then the revisions of BMAT that depend on ZMAT are calculated.
      sumz = sum(zmat, dim = 1)
      do k = 1, idz - 1
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
          bmat(:, 1:npt) = bmat(:, 1:npt) - outprod(vlag, zmat(:,k))
          bmat(:,npt+1:npt+n) = bmat(:,npt+1:npt+n) - outprod(vlag,vlag)
      end do
      do k = idz, npt - n - 1
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
          bmat(:, 1:npt) = bmat(:, 1:npt) + outprod(vlag, zmat(:,k))
          bmat(:,npt+1:npt+n) = bmat(:,npt+1:npt+n) + outprod(vlag,vlag)
      end do
!!!!!!!!!!!!!!!!!!!!!COMPACT SCHEME ENDS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Set the upper triangular part of BMAT(:,NPT+1:NPT+N) by symmetry 
      ! Note that UPDATE sets the upper triangular part by copying
      ! the lower triangular part, but here it does the opposite. There
      ! seems not any particular reason to keep them different. It was
      ! probably an ad-hoc decision that Powell made when coding. 
      ! As mentioned above, this part can be spared if we put a pair of
      ! parenthsis around the two outter products.
      do j = 1, n
          bmat(j, npt + 1 : npt + j - 1) = bmat(1 : j - 1, npt + j)
      end do
      
      ! The following instructions complete the shift of XBASE.
      ! Recall the we have already subtracted HALF*XOPT from XPT. 
      ! Therefore, overall, the new XPT is XPT - XOPT.
      xpt = xpt - HALF*spread(xopt, dim = 2, ncopies = npt)
      !do k = 1, npt
      !    xpt(:, k) = xpt(:, k) - HALF*xopt
      !end do

      return 

      end subroutine shiftbase
