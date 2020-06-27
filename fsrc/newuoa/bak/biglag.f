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

      real(kind=rp), intent(in) :: xopt(n), xpt(n, npt), bmat(npt+n, n),&
     & zmat(npt, npt-n-1), delta 
      real(kind=rp), intent(out) :: alpha, d(n)

      integer :: i, isave, iterc, iu, k
      real(kind=rp) :: hcol(npt), gc(n), gd(n), s(n), w(n),             &
     & zknew(npt - n - 1), angle, cf(5), cth, dd, denom, dhd, gg,       &
     & scaling, sp, ss, step, sth, tau, taubeg, tauold, taumax, temp,   &
     & tempa, tempb

          real(kind=rp) :: BMAR(n, npt+n)
          BMAR = transpose(bmat)

       
      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! XOPT is the best interpolation point so far.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N columns of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! sub-matrix of H.
      ! KNEW is the index of the interpolation point to be removed.
      ! DELTA is the current trust region bound.
      ! D will be set to the step from XOPT to the new point.
      ! ALPHA will be set to the KNEW-th diagonal element of matrix H.
      ! HCOL, GC, GD, S and W will be used for working space.

      
      ! Set HCOL to the leading NPT elements of the KNEW-th column of H.
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      hcol = matmul(zmat, zknew)
      alpha = hcol(knew)

      ! Set the unscaled initial direction D. Form the gradient of LFUNC
      ! at XOPT, and multiply D by the Hessian of LFUNC.
      d = xpt(:, knew) - xopt
      dd = dot_product(d, d)

      gd = matmul(xpt, hcol*matmul(d, xpt))

      !----------------------------------------------------------------!
      ! The following DO LOOP calculates the GC below
!-----!gc = BMAR(:, knew) + matmul(xpt, hcol*matmul(xopt, xpt)) !------!
      gc = BMAR(:, knew)
      do k = 1, npt
          gc = gc + (hcol(k)*dot_product(xpt(:, k), xopt))*xpt(:, k)
      end do
      !----------------------------------------------------------------!

      ! Scale D and GD, with a sign change if required. Set S to another
      ! vector in the initial two dimensional subspace.
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
      d = scaling*d
      gd = scaling*gd
      s = gc + temp*gd
      
      ! Begin the iteration by overwriting S with a vector that has the
      ! required length and direction, except that termination occurs if
      ! the given D and S are nearly parallel.
      do iterc = 1, n
          dd = dot_product(d, d)
          sp = dot_product(d, s)
          ss = dot_product(s, s)
          if (dd*ss - sp*sp <= 1.0e-8_rp*dd*ss) then 
              exit
          end if
          denom = sqrt(dd*ss - sp*sp)
          s = (dd*s - sp*d)/denom

          w = matmul(xpt, hcol*matmul(s, xpt))
!          w = zero
!          do k = 1, npt
!             w = w + (hcol(k)*dot_product(xpt(:, k), s))*xpt(:, k)
!          end do
          
          ! Calculate the coefficients of the objective function on the
          ! circle, beginning with the multiplication of S by the second
          ! derivative matrix.
          cf(1) = dot_product(s, w)
          cf(2) = dot_product(d, gc)
          cf(3) = dot_product(s, gc)
          cf(4) = dot_product(d, gd)
          cf(5) = dot_product(s, gd)
          cf(1) = half*cf(1)
          cf(4) = half*cf(4) - cf(1)
          
          ! Seek the value of the angle that maximizes |TAU|.
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
              tau = cf(1) + (cf(2)+cf(4)*cth)*cth+(cf(3)+cf(5)*cth)*sth
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
          if (abs(tempa - tempb) > zero) then
              tempa = tempa - taumax
              tempb = tempb - taumax
              step = half*(tempa - tempb)/(tempa + tempb)
          else
              step = zero
          end if
          angle = temp*(real(isave, rp) + step)
          
          ! Calculate the new D and GD. Then test for convergence.
          cth = cos(angle)
          sth = sin(angle)
          tau = cf(1) + (cf(2) + cf(4)*cth)*cth + (cf(3)+cf(5)*cth)*sth
          d = cth*d + sth*s
          gd = cth*gd + sth*w
          s = gc + gd
          if (abs(tau) <= 1.1_rp*abs(taubeg)) then
              exit
          end if
      end do

      return

      end subroutine biglag
