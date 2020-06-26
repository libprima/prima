      subroutine shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq,hq,xpt)
      ! SHIFTBASE shifts the base point to XBASE + XOPT and updates GQ,
      ! HQ, and BMAT accordingly. PQ and ZMAT remain the same after the
      ! shifting. See Section 7 of the NEWUOA paper.

      use consts, only : rp, half, quart
      use lina
      implicit none

      integer, intent(in) :: idz, n, npt

      real(kind = rp), intent(in) :: xopt(n), pq(npt), zmat(npt,npt-n-1)
      real(kind = rp), intent(inout) :: bmat(npt + n, n), gq(n),        &
     & hq((n*(n + 1))/2), xpt(npt, n)

      integer :: i, ih, j, k
      real(kind = rp) :: sumz(npt-n-1), vlag(n), qxoptq, xoptsq, pqx(n),&
     & bmatk(n), w1(npt), w2(n), w3(npt) 


      xoptsq = dot_product(xopt, xopt)
      qxoptq = quart * xoptsq

      !----------------------------------------------------------------!   
      ! The update for gq can indeed be done by the following 3 lines:
      !real(kind = rp) :: hxopt(n)
      !call hessmul(n, npt, xpt, hq, pq, xopt, hxopt)
      !gq = gq + hxopt 
      !----------------------------------------------------------------!
      do k = 1, npt
      ! Update of GQ due to the implicit part of the HESSIAN:
      ! GQ = GQ + (\sum_{K=1}^NPT PQ(K)*XPT(K, :)'*XPT(K, :)) * XOPT
          gq = gq + pq(k)*dot_product(xpt(k, :), xopt)*xpt(k, :)
      end do
      ih = 0
      do j = 1, n
      ! Update of GQ due to the explicit part of the HESSIAN    
         do i = 1, j
             ih = ih + 1
             if (i < j) then 
                 gq(j) = gq(j) + hq(ih)*xopt(i)
             end if
             gq(i) = gq(i) + hq(ih)*xopt(j)
         end do
      end do
      
      w1 = matmul(xpt, xopt) - half*xoptsq
      ! W1 equals MATMUL(XPT, XOPT) after XPT is updated as follows.
      do k = 1, npt
          xpt(k, :) = xpt(k, :) - half*xopt
      end do

      ! Update HQ. It has to be done after the above revision to XPT!!!
      pqx = matmul(pq, xpt)
      ! The followind DO LOOP indeed adds (pqx*xopt' + xopt*pqx') to the 
      ! explicit part of the HESSIAN
      ih = 0
      do j = 1, n
          hq(ih+1:ih+j)=hq(ih+1:ih+j)+pqx(1:j)*xopt(j)+xopt(1:j)*pqx(j)
          ih = ih + j
      end do

!==============================POWELL'S SCHEME=========================!
!      ! Make the changes to BMAT that do not depend on ZMAT.
!      do k = 1, npt
!          bmatk = bmat(k, :)
!          w2 = w1(k)*xpt(k, :) + qxoptq*xopt
!          ! The following DO LOOP performs the update below, but only
!          ! calculates the lower triangular part since it is symmetric.
!!----------------------------------------------------------------------!
!!          bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) +     &
!!     &     outprod(bmatk, w2) + outprod(w2, bmatk)
!!----------------------------------------------------------------------!
!          do i = 1, n
!              bmat(npt + i, 1 : i) = bmat(npt + i, 1 : i) +             &
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
!!          vlag = qxoptq*sumz(k)*xopt(j) + matmul(w1*zmat(:, k), xpt)
!!----------------------------------------------------------------------!
!          w3 = w1*zmat(:, k)
!          do j = 1, n
!              vlag(j) = qxoptq*sumz(k)*xopt(j)
!              do i = 1, npt
!                  vlag(j) = vlag(j) + w3(i)*xpt(i, j)
!              end do
!          end do
!
!          if (k < idz) then
!               bmat(1:npt, :) = bmat(1:npt, :) - outprod(zmat(:,k),vlag)
!          else
!               bmat(1:npt, :) = bmat(1:npt, :) + outprod(zmat(:,k),vlag)
!          end if
!
!          if (k < idz) then 
!          ! The following DO LOOP performs the update below, but only
!          ! calculates the lower triangular part since it is symmetric.
!!----------------------------------------------------------------------!          
!!              bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) - &
!!     &         outprod(vlag, vlag)
!!----------------------------------------------------------------------!          
!              do i = 1, n
!                  bmat(npt+i, 1:i) = bmat(npt+i,1:i) - vlag(i)*vlag(1:i)
!              end do
!          else
!          ! The following DO LOOP performs the update below, but only
!          ! calculates the lower triangular part since it is symmetric.
!!----------------------------------------------------------------------!
!!              bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) + &
!!     &         outprod(vlag, vlag)
!!----------------------------------------------------------------------!
!              do i = 1, n
!                  bmat(npt+i, 1:i) = bmat(npt+i,1:i) + vlag(i)*vlag(1:i)
!              end do
!          end if
!      end do
!=========================POWELL'S SCHEME ENDS=========================!


!!!!!!!!!!!!!!!!!!!!!A COMPACTER YET NOT SO ECONOMIC SCHEME!!!!!!!!!!!!!
      ! Make the changes to BMAT that do not depend on ZMAT.
      do k = 1, npt
          bmatk = bmat(k, :)
          w2 = w1(k)*xpt(k, :) + qxoptq*xopt
          ! This is the only place where non-symmetry of
          ! BMAT(NPT+1:NPT+N, :) can come from. It is because of the
          ! non-associtivity of floating point arithmetic addition. To
          ! make BMAT(NPT+1:NPT+N,:) symmetric, use the following two
          ! lines to replace the code below. The only difference is the
          ! parenthsis around the two outter products. It is probably 
          ! a BETTERE implementation so we should take it in future
          ! versions. 
!----------------------------------------------------------------------!
!          bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) +     &
!     &     ( outprod(bmatk, w2) + outprod(w2, bmatk) )
!----------------------------------------------------------------------!
          bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) +     &
     &     outprod(bmatk, w2) + outprod(w2, bmatk)
      end do

      ! Then the revisions of BMAT that depend on ZMAT are calculated.
      sumz = sum(zmat, dim = 1)
      do k = 1, npt - n - 1
          ! The following W3 and DO LOOP indeed defines the VLAG below.
          ! The results are not identical due to the non-associtivity 
          ! of floating point arithmetic addition.
!----------------------------------------------------------------------!
          !vlag = qxoptq*sumz(k)*xopt(j) + matmul(w1*zmat(:, k), xpt)
!----------------------------------------------------------------------!
          w3 = w1*zmat(:, k)
          do j = 1, n
              vlag(j) = qxoptq*sumz(k)*xopt(j)
              do i = 1, npt
                  vlag(j) = vlag(j) + w3(i)*xpt(i, j)
              end do
          end do

          if (k < idz) then
               bmat(1:npt, :) = bmat(1:npt, :) - outprod(zmat(:,k),vlag)
              bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) - &
     &         outprod(vlag, vlag)
          else
               bmat(1:npt, :) = bmat(1:npt, :) + outprod(zmat(:,k),vlag)
              bmat(npt + 1 : npt + n, :) = bmat(npt + 1 : npt + n, :) + &
     &         outprod(vlag, vlag)
          end if
      end do
!!!!!!!!!!!!!!!!!!!!!COMPACT SCHEME ENDS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Set the upper triangular part of BMAT(NPT+1:NPT+N,:) by symmetry 
      ! Note that UPDATE sets the lower triangular part by copying
      ! the upper triangular part, but here it does the opposite. There
      ! seems not any particular reason to keep them different. It was
      ! probably an ad-hoc decision that Powell made when coding. 
      do j = 1, n
          bmat(npt + 1 : npt + j - 1, j) = bmat(npt + j, 1 : j - 1)
      end do
      
      ! The following instructions complete the shift of XBASE.
      ! Recall the we have already subtracted HALF*XOPT from XPT. 
      ! Therefore, overall, the new XPT is XPT - XOPT.
      do k = 1, npt
          xpt(k, :) = xpt(k, :) - half*xopt
      end do

      return 

      end subroutine shiftbase
