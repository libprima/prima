      subroutine biglag(n, npt, xopt, xpt, bmat, zmat, idz, knew, delta,&
     & d, alpha)
      ! BIGLAG calculates a D by approximately solving
      !
      ! max |LFUNC(XOPT + D)|, subject to ||D|| <= DELTA, 
      !
      ! where LFUNC is the KNEW-th Lagrange function.

      use consts, only : rp, one, two, half, pi, zero
      use lina
      implicit none

      integer, intent(in) :: n, npt, knew, idz

      real(kind=rp), intent(in) :: xopt(n), xpt(npt, n), bmat(npt+n, n),&
     & zmat(npt, npt-n-1), delta 
      real(kind=rp), intent(out) :: alpha, d(n)

      integer :: i, isave, iterc, iu, k
      real(kind=rp) :: hcol(npt), gc(n), gd(n), s(n), w(n), zk(npt-n-1),&
     & angle, cf(5), cth, dd, denom, dhd, gg, scaling, sp, ss, step,    &
     & sth, tau, taubeg, tauold, taumax, temp, tempa, tempb

       
      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! XOPT is the best interpolation point so far.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N columns of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! submatrix of H.
      ! KNEW is the index of the interpolation point to be removed.
      ! DELTA is the current trust region bound.
      ! D will be set to the step from XOPT to the new point.
      ! ALPHA will be set to the KNEW-th diagonal element of matrix H.
      ! HCOL, GC, GD, S and W will be used for working space.

      
      ! Set HCOL to the leading NPT elements of the KNEW-th column of H.
      iterc = 0
      zk = zmat(knew, :)
      zk(1 : idz - 1) = -zk(1 : idz - 1)
      hcol = matmul(zmat, zk)
      alpha = hcol(knew)

!      do j = 1, npt - n - 1 
!         temp = zmat(knew, j)
!         if (j < idz) temp = - temp
!         do k = 1, npt
!             hcol(k) = hcol(k) + temp*zmat(k, j)
!         end do
!      end do

      ! Set the unscaled initial direction D. Form the gradient of LFUNC
      ! atXOPT, and multiply D by the second derivative matrix of LFUNC.
      d = xpt(knew, :) - xopt
      dd = dot_product(d, d)
      gc = bmat(knew, :)
      gd = zero
!      do i = 1, n
!         d(i) = xpt(knew, i) - xopt(i)
!         gc(i) = bmat(knew, i)
!         gd(i) = zero
!         dd = dd + d(i)**2
!      end do
           
      do k = 1, npt
!         temp = zero
!         sum = zero
!         do j = 1, n
!      temp = temp + xpt(k, j)*xopt(j)
!      summation = summation + xpt(k, j)*d(j)
!         end do
!          temp = hcol(k)*temp
!          summation = hcol(k)*summation
!          gc = gc + temp*xpt(k, :)
!          gd = gd + summation*xpt(k, :)
!          do i = 1, n
!              gc(i) = gc(i) + temp*xpt(k, i)
!              gd(i) = gd(i) + sum*xpt(k, i)
!          end do
!          temp = hcol(k)*(dot_product(xpt(k,:), xopt))
!          summation = hcol(k)*(dot_product(xpt(k, :), d))
          gc = gc + (hcol(k)*dot_product(xpt(k, :), xopt))*xpt(k, :)
          gd = gd + (hcol(k)*dot_product(xpt(k, :), d))*xpt(k, :)
      end do

      ! Scale D and GD, with a sign change if required. Set S to another
      ! vector in the initial two dimensional subspace.
!      gg = zero
!      sp = zero
!      dhd = zero
!      do i = 1, n
!         gg = gg + gc(i)**2
!         sp = sp + d(i)*gc(i)
!         dhd = dhd + d(i)*gd(i)
!      end do
      gg = dot_product(gc, gc)
      sp = dot_product(d, gc)
      dhd = dot_product(d, gd)
      scaling = delta/sqrt(dd)
      if (sp*dhd < zero) then 
          scaling = - scaling
      end if
      temp = zero
      if (sp*sp > 0.99_rp*dd*gg) then 
          temp = one
      end if
      tau = scaling*(abs(sp) + half*scaling*abs(dhd))
      if (gg*(delta*delta) < 0.01_rp*tau*tau) then 
          temp = one
      end if
!      do i = 1, n
!         d(i) = scaling*d(i)
!         gd(i) = scaling*gd(i)
!         s(i) = gc(i) + temp*gd(i)
!      end do
      d = scaling*d
      gd = scaling*gd
      s = gc + temp*gd
      
      ! Begin the iteration by overwriting S with a vector that has the
      ! required length and direction, except that termination occurs if
      ! the given D and S are nearly parallel.
   80 iterc = iterc + 1
!      dd = zero
!      sp = zero
!      ss = zero
!      do i = 1, n
!         dd = dd + d(i)**2
!         sp = sp + d(i)*s(i)
!         ss = ss + s(i)**2
!      end do
      dd = dot_product(d, d)
      sp = dot_product(d, s)
      ss = dot_product(s, s)
      temp = dd*ss - sp*sp
      if (temp <= 1.0e-8_rp*dd*ss) goto 160
      denom = sqrt(temp)
!      do i = 1, n
!         s(i) = (dd*s(i) - sp*d(i))/denom
!         w(i) = zero
!      end do
      s = (dd*s - sp*d)/denom
      w = zero
      
      ! Calculate the coefficients of the objective function on the
      ! circle, beginning with the multiplication of S by the second
      ! derivative matrix.
      do k = 1, npt
!         sum = zero
!         do j = 1, n
!             sum = sum + xpt(k, j)*s(j)
!         end do
!         sum = hcol(k)*sum
!         do i = 1, n
!             w(i) = w(i) + sum*xpt(k, i)
!         end do
!         summation = hcol(k)*dot_product(xpt(k, :), s)
         w = w + (hcol(k)*dot_product(xpt(k, :), s))*xpt(k, :)
      end do

!      cf(1) = zero
!      cf(2) = zero
!      cf(3) = zero
!      cf(4) = zero
!      cf(5) = zero
!      do i = 1, n
!         cf(1) = cf(1) + s(i)*w(i)
!         cf(2) = cf(2) + d(i)*gc(i)
!         cf(3) = cf(3) + s(i)*gc(i)
!         cf(4) = cf(4) + d(i)*gd(i)
!         cf(5) = cf(5) + s(i)*gd(i)
!      end do
      cf(1) = dot_product(s, w)
      cf(2) = dot_product(d, gc)
      cf(3) = dot_product(s, gc)
      cf(4) = dot_product(d, gd)
      cf(5) = dot_product(s, gd)
      cf(1) = half*cf(1)
      cf(4) = half*cf(4) - cf(1)
      
      ! Seek the value of the angle that maximizes the modulus of TAU.
      taubeg = cf(1) + cf(2) + cf(4)
      taumax = taubeg
      tauold = taubeg
      isave = 0
      iu = 49
      temp = (two*pi)/real(iu + 1, rp)
      do i = 1, iu
          angle = real(i, rp)*temp
          cth = cos(angle)
          sth = sin(angle)
          tau = cf(1) + (cf(2) + cf(4)*cth)*cth + (cf(3)+cf(5)*cth)*sth
          if (abs(tau) > abs(taumax)) then
              taumax = tau
              isave = i
              tempa = tauold
          else if (i == isave + 1) then
              tempb = tau
          end if
          tauold = tau
      end do
      if (isave == 0) then 
          tempa = tau
      end if
      if (isave == iu) then 
          tempb = taubeg
      end if
      tempa = tempa - taumax
      tempb = tempb - taumax
      if (abs(tempa - tempb) > zero) then
          step = half*(tempa - tempb)/(tempa + tempb)
      else
          step = zero
      end if
      angle = temp*(real(isave, rp) + step)
      
      ! Calculate the new D and GD. Then test for convergence.
      cth = cos(angle)
      sth = sin(angle)
      tau = cf(1) + (cf(2) + cf(4)*cth)*cth + (cf(3) + cf(5)*cth)*sth
!      do i = 1, n
!          d(i) = cth*d(i) + sth*s(i)
!          gd(i) = cth*gd(i) + sth*w(i)
!          s(i) = gc(i) + gd(i)
!      end do
      d = cth*d + sth*s
      gd = cth*gd + sth*w
      s = gc + gd
      if (abs(tau) <= 1.1_rp*abs(taubeg)) goto 160
      if (iterc < n) goto 80
  160 return

      end subroutine biglag
