!!! TODO: Replace the DO LOOP in the update of S with MATMUL.
      subroutine bigden (n, npt, xopt, xpt, bmat, zmat, idz, kopt, knew,&
     & d, wcheck, vlag, beta)
      ! BIGDEN calculates a D by approximately solving
      !
      ! max |LFUNC(XOPT + D)|, subject to ||D|| <= DELTA, 
      !
      ! where LFUNC is the KNEW-th Lagrange function.

      use consts, only : rp, one, two, half, quart, pi, zero
      use lina
      implicit none

      integer, intent(in) :: n, npt, knew, kopt, idz

      real(kind=rp), intent(in) :: xopt(n), xpt(npt, n), bmat(npt+n, n),&
     & zmat(npt, npt - n - 1)
      real(kind=rp), intent(out) :: beta, vlag(npt + n), wcheck(npt + n) 
      real(kind=rp), intent(inout) :: d(n)

      integer :: i, isave, iterc, iu, j, jc, k, nw
      real(kind=rp) :: s(n), wvec(npt + n, 5), prod(npt + n, 5), den(9),&
     & denex(9), par(9), zknew(npt - n - 1), stemp(n), dstemp(npt),     &
     & sstemp(npt), wz(npt - n - 1), angle, dd, denmax, denold, densav, &
     & ds, dtest, ss, ssden, summation, sumold, step, tau, temp, tempa, &
     & tempb, tempc, tempv(npt), xoptd, xopts, xoptsq, alpha

      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! XOPT is the best interpolation point so far.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N columns of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! submatrix of H.
      ! NDIM is the first dimension of BMAT and has the value NPT + N.
      ! KOPT is the index of the optimal interpolation point.
      ! KNEW is the index of the interpolation point to be removed.
      ! D will be set to the step from XOPT to the new point, and on 
      ! entry it should be the D that was calculated by the last call
      ! of BIGLAG. The length of the initial D provides a trust region
      ! bound on the final D.
      ! WCHECK will be set to wcheck for the final choice of D.
      ! VLAG will be set to Theta*WCHECK+e_b for the final choice of D.
      ! BETA will be set to the value that will occur in the updating
      ! formula when the KNEW-th interpolation point is moved to its new
      ! position.
      
      ! D is calculated in a way that should provide a denominator with 
      ! a large modulus in the updating formula when the KNEW-th
      ! interpolation point is shifted to the new position XOPT + D.
      

      ! Store the first NPT elements of the KNEW-th column of H in 
      ! WCHECK(N + 1) to WCHECK(N + NPT).
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      wcheck(n + 1 : n + npt) = matmul(zmat, zknew)
      alpha = wcheck(n + knew)
      
      ! The initial search direction D is taken from the last call of
      ! BIGLAG, and the initial S is set below, usually to the direction
      ! from XOPT to X_KNEW, but a different direction to an 
      ! interpolation point may be chosen, in order to prevent S from
      ! being nearly parallel to D.
!      dd = zero
!      ds = zero
!      ss = zero
!      xoptsq = zero
!      do i = 1, n
!          dd = dd + d(i)**2
!          s(i) = xpt(knew, i) - xopt(i)
!          ds = ds + d(i)*s(i)
!          ss = ss + s(i)**2
!          xoptsq = xoptsq + xopt(i)**2
!      end do
      dd = dot_product(d, d)
      s = xpt(knew, :) - xopt
      ds = dot_product(d, s)
      ss = dot_product(s, s)
      xoptsq = dot_product(xopt, xopt)

      !----------------------------------------------------------------!
      ! Zaikun 2019-08-29: With the original code, if DS, DD, or SS is 
      ! NaN, KSAV will not get a value. This may cause Segmentation 
      ! Fault because XPT(KSAV, :) will later be accessed. 
      !IF (DS*DS .GT. 0.99D0*DD*SS) THEN
      if (.not. (ds*ds <= 0.99_rp*dd*ss)) then
          !ksav = knew
          dtest = ds*ds/ss
          do k = 1, npt
!              if (k /= kopt) then
!                  dstemp = zero
!                  sstemp = zero
!                  do i = 1, n
!                      diff = xpt(k, i) - xopt(i)
!                      dstemp = dstemp + d(i)*diff
!                      sstemp = sstemp + diff*diff
!                  end do
                  stemp = xpt(k, :) - xopt
                  dstemp(k) = dot_product(d, stemp)
                  sstemp(k) = dot_product(stemp, stemp)
!                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  ! Zaikun 2019-08-29: See the comments above.
!                  !IF (DSTEMP*DSTEMP/SSTEMP .LT. DTEST) THEN
!                  if (.not. (dstemp*dstemp/sstemp >= dtest)) then
!                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      ksav = k
!                      dtest = dstemp*dstemp/sstemp
!                      ds = dstemp
!                      ss = sstemp
!                  end if
!              end if
          end do
          
          dstemp(kopt) = two*ds + one  
          sstemp(kopt) = ss     
          k = minloc(dstemp*dstemp/sstemp, dim = 1)
          if ((.not. (dstemp(k)*dstemp(k)/sstemp(k) >= dtest)) .and.    &
     &     k /= kopt) then
          ! Althoguh unlikely, if NaN occurs, it may happen that k=kopt
              s = xpt(k, :) - xopt
              ds = dstemp(k)
              ss = sstemp(k)
          end if
      end if

      ssden = dd*ss - ds*ds
      densav = zero

      ! Begin the iteration by overwriting S with a vector that has the
      ! required length and direction.
      do iterc = 1, n
!      xoptd = zero
!      xopts = zero
!      do i = 1, n
!          s(i) = temp*(dd*s(i) - ds*d(i))
!          xoptd = xoptd + xopt(i)*d(i)
!          xopts = xopts + xopt(i)*s(i)
!      end do
      s = (one/sqrt(ssden))*(dd*s - ds*d)
      xoptd = dot_product(xopt, d)
      xopts = dot_product(xopt, s)

      ! Set the coefficients of the first two terms of BETA.
      tempa = half*xoptd*xoptd
      tempb = half*xopts*xopts
      den(1) = dd*(xoptsq + half*dd) + tempa + tempb
      den(2) = two*xoptd*dd
      den(3) = two*xopts*dd
      den(4) = tempa - tempb
      den(5) = xoptd*xopts
      den(6 : 9) = zero
      
      ! Put the coefficients of WCHECK in WVEC.
      do k = 1, npt
!          tempa = zero
!          tempb = zero
!          tempc = zero
!          do i = 1, n
!              tempa = tempa + xpt(k, i)*d(i)
!              tempb = tempb + xpt(k, i)*s(i)
!              tempc = tempc + xpt(k, i)*xopt(i)
!          end do
          tempa = dot_product(xpt(k, :), d)
          tempb = dot_product(xpt(k, :), s)
          tempc = dot_product(xpt(k, :), xopt)
          wvec(k, 1) = quart*(tempa*tempa + tempb*tempb)
          wvec(k, 2) = tempa*tempc
          wvec(k, 3) = tempb*tempc
          wvec(k, 4) = quart*(tempa*tempa - tempb*tempb)
          wvec(k, 5) = half*tempa*tempb
      end do
!      do i = 1, n
!          ip = i + npt
!          wvec(ip, 1) = zero
!          wvec(ip, 2) = d(i)
!          wvec(ip, 3) = s(i)
!          wvec(ip, 4) = zero
!          wvec(ip, 5) = zero
!      end do
      wvec(npt + 1 : npt + n, 1 : 5) = zero
      wvec(npt + 1 : npt + n, 2) = d
      wvec(npt + 1 : npt + n, 3) = s

      
      ! Put the coefficents of THETA*WCHECK in PROD.
      do jc = 1, 5
!          do k = 1, npt
!              prod(k, jc) = zero
!          end do
!          prod(1 : npt, jc) = zero
!          do j = 1, npt - n - 1
!              sum = zero
!              do k = 1, npt
!                  sum = sum + zmat(k, j)*wvec(k, jc)
!              end do
!              wz(j) = dot_product(zmat(:, j), wvec(:, jc))
!              if (j < idz) sum = - sum
!              do k = 1, npt
!                  prod(k, jc) = prod(k, jc) + sum*zmat(k, j)
!              end do
!          end do

          wz = matmul(wvec(1 : npt, jc), zmat)
          wz(1 : idz - 1) = -wz(1 : idz - 1)
          prod(1 : npt, jc) = matmul(zmat, wz)
          
          
          nw = npt
          if (jc == 2 .or. jc == 3) then
!              nw = npt + n
!          end if
!          if (nw == npt + n) then
!              do k = 1, npt
!                  sum = zero
!                  do j = 1, n
!                      sum = sum + bmat(k, j)*wvec(npt+j, jc)
!                  end do
!                  prod(k, jc) = prod(k, jc) + sum
!              end do
              !bw = matmul(bmat(1:npt, wvec(npt+1 : npt + n, jc))
              prod(1 : npt, jc) = prod(1 : npt, jc) +                   &
     &         matmul(bmat(1 : npt, :), wvec(npt + 1 : npt + n, jc))
              nw = npt + n
          end if
!          do j = 1, n
!              sum = zero
!              do i = 1, nw
!                  sum = sum + bmat(i, j)*wvec(i, jc)
!              end do
!              prod(npt + j, jc) = dot_product(bmat(1:nw, j), wvec(1:nw, jc)) 
!          end do
          prod(npt + 1 : npt + n, jc) = matmul(wvec(1 : nw, jc),        &
     &     bmat(1 : nw, 1 : n))
      end do

      ! Include in DEN the part of BETA that depends on THETA.
      do k = 1, npt + n
!          sum = zero
!          do i = 1, 5
!              par(i) = half*prod(k, i)*wvec(k, i)
!              sum = sum + par(i)
!          end do
          par(1 : 5) = half*prod(k, 1 : 5)*wvec(k, 1 : 5)
          den(1) = den(1) - par(1) - sum(par(1 : 5)) 
          tempa = prod(k, 1)*wvec(k, 2) + prod(k, 2)*wvec(k, 1)
          tempb = prod(k, 2)*wvec(k, 4) + prod(k, 4)*wvec(k, 2)
          tempc = prod(k, 3)*wvec(k, 5) + prod(k, 5)*wvec(k, 3)
          den(2) = den(2) - tempa-half*(tempb + tempc)
          den(6) = den(6) - half*(tempb-tempc)
          tempa = prod(k, 1)*wvec(k, 3) + prod(k, 3)*wvec(k, 1)
          tempb = prod(k, 2)*wvec(k, 5) + prod(k, 5)*wvec(k, 2)
          tempc = prod(k, 3)*wvec(k, 4) + prod(k, 4)*wvec(k, 3)
          den(3) = den(3) - tempa-half*(tempb-tempc)
          den(7) = den(7) - half*(tempb + tempc)
          tempa = prod(k, 1)*wvec(k, 4) + prod(k, 4)*wvec(k, 1)
          den(4) = den(4) - tempa-par(2) + par(3)
          tempa = prod(k, 1)*wvec(k, 5) + prod(k, 5)*wvec(k, 1)
          tempb = prod(k, 2)*wvec(k, 3) + prod(k, 3)*wvec(k, 2)
          den(5) = den(5) - tempa-half*tempb
          den(8) = den(8) - par(4) + par(5)
          tempa = prod(k, 4)*wvec(k, 5) + prod(k, 5)*wvec(k, 4)
          den(9) = den(9) - half*tempa
      end do
      
      ! Extend DEN so that it holds all the coefficients of DENOM.
!      sum = zero
!      do i = 1, 5
!          par(i) = half*prod(knew, i)**2
!          sum = sum + par(i)
!      end do
      par(1 : 5) = half*prod(knew, 1 : 5)**2
      denex(1) = alpha*den(1) + par(1) + sum(par(1 : 5))
      tempa = two*prod(knew, 1)*prod(knew, 2)
      tempb = prod(knew, 2)*prod(knew, 4)
      tempc = prod(knew, 3)*prod(knew, 5)
      denex(2) = alpha*den(2) + tempa + tempb + tempc
      denex(6) = alpha*den(6) + tempb - tempc
      tempa = two*prod(knew, 1)*prod(knew, 3)
      tempb = prod(knew, 2)*prod(knew, 5)
      tempc = prod(knew, 3)*prod(knew, 4)
      denex(3) = alpha*den(3) + tempa + tempb - tempc
      denex(7) = alpha*den(7) + tempb + tempc
      tempa = two*prod(knew, 1)*prod(knew, 4)
      denex(4) = alpha*den(4) + tempa + par(2) - par(3)
      tempa = two*prod(knew, 1)*prod(knew, 5)
      denex(5) = alpha*den(5) + tempa + prod(knew, 2)*prod(knew, 3)
      denex(8) = alpha*den(8) + par(4) - par(5)
      denex(9) = alpha*den(9) + prod(knew, 4)*prod(knew, 5)
      
      ! Seek the value of the angle that maximizes the modulus of DENOM.
      summation = denex(1) + denex(2) + denex(4) + denex(6) + denex(8)
      denold = summation
      denmax = summation
      isave = 0
      iu = 49
      temp = (two*pi)/real(iu + 1, rp)
      par(1) = one
      do i = 1, iu
          angle = real(i, rp)*temp
          par(2) = cos(angle)
          par(3) = sin(angle)
          do j = 4, 8, 2
              par(j) = par(2)*par(j - 2) - par(3)*par(j - 1)
              par(j + 1) = par(2)*par(j - 1) + par(3)*par(j - 2)
          end do
          sumold = summation
!          summation = zero
!          do j = 1, 9
!              summation = summation + denex(j)*par(j)
!          end do
          summation = dot_product(denex(1 : 9), par(1 : 9))
          if (abs(summation) > abs(denmax)) then
              denmax = summation
              isave = i
              tempa = sumold
          else if (i == isave + 1) then
              tempb = summation
          end if
      end do
      if (isave == 0) then
          tempa = summation
      end if
      if (isave == iu) then 
          tempb = denold
      end if
      if (abs(tempa - tempb) > 0) then
          tempa = tempa - denmax
          tempb = tempb - denmax
          step = half*(tempa - tempb)/(tempa + tempb)
      else
          step = zero
      end if
      angle = temp*(real(isave, rp) + step)
      
      ! Calculate the new parameters of the denominator, the new VLAG 
      ! vector and the new D. Then test for convergence.
      par(2) = cos(angle)
      par(3) = sin(angle)
      do j = 4, 8, 2
          par(j) = par(2)*par(j - 2) - par(3)*par(j - 1)
          par(j + 1) = par(2)*par(j - 1) + par(3)*par(j - 2)
      end do

!      beta = zero
!      denmax = zero
!      do j = 1, 9
!          beta = beta + den(j)*par(j)
!          denmax = denmax + denex(j)*par(j)
!      end do

      beta = dot_product(den(1 : 9), par(1 : 9))
      denmax = dot_product(denex(1 : 9), par(1 : 9))

!      do k = 1, npt + n
!          vlag(k) = zero
!          do j = 1, 5
!              vlag(k) = vlag(k) + prod(k, j)*par(j)
!          end do
!      end do

      vlag = matmul(prod(:, 1 : 5), par(1 : 5))

      tau = vlag(knew)

!      dd = zero
!      tempa = zero
!      tempb = zero
!      do i = 1, n
!          d(i) = par(2)*d(i) + par(3)*s(i)
!          wcheck(i) = xopt(i) + d(i)
!          dd = dd + d(i)**2
!          tempa = tempa + d(i)*wcheck(i)
!          tempb = tempb + wcheck(i)*wcheck(i)
!      end do

      d = par(2)*d + par(3)*s


      dd = dot_product(d, d)
      wcheck(1 : n) = xopt + d
      tempa = dot_product(d, wcheck(1 : n))
      tempb = dot_product(wcheck(1 : n), wcheck(1 : n))

      if (iterc > 1) then
          densav = max(densav, denold)
      end if
      if (abs(denmax) <= 1.1_rp*abs(densav)) then 
          exit
      end if
      densav = denmax

      ! Set S to half the gradient of the denominator with respect to D.
      ! Then branch for the next iteration.
!      do i = 1, n
!          temp = tempa*xopt(i) + tempb*d(i) - vlag(npt+i)
!          s(i) = tau*bmat(knew, i) + alpha*temp
!      end do

!!----------------------------------------------------------------------!
!      tempv(1 : n) = tempa*xopt + tempb*d - vlag(npt + 1 : npt + n)
!      s = tau*bmat(knew, :) + alpha*tempv(1 : n)
    
      s = tau*bmat(knew,:)+alpha*(tempa*xopt+tempb*d-vlag(npt+1:npt+n))
!!----------------------------------------------------------------------!

!      
!      do k = 1, npt
!          summation = zero
!          do j = 1, n
!              summation = summation + xpt(k, j)*wcheck(j)
!          end do
!          temp = (tau*wcheck(n + k) - alpha*vlag(k))*summation
!          do i = 1, n
!              s(i) = s(i) + temp*xpt(k, i)
!          end do
!      end do

!----------------------------------------------------------------------!
      tempv = matmul(xpt, wcheck(1 : n))
      tempv = (tau*wcheck(n + 1 : n + npt) - alpha*vlag(1 : npt))*tempv
      !!! In later versiions, the following DO LOOP should be replaced
      !!! by MATMUL as follows:
!-----!s = s + matmul(tempv, xpt) !------------------------------------!
      !!! For the moment, we use a DO LOOP instead of MATMUL to ensure
      !!! that the result is identical to that of Powell's code. If we
      !!! use MATMUL, the result is not always the same, because floating
      !!! point addition (and multiplication) is not associative,
      !!! although it is commutative. The following LOOP calculates 
      !!! s + tempv(1)*xpt(1, :) + tempv(2)*xpt(2, :) + ...
      !!! while s + matmul(tempv, xpt) calculates
      !!! s + ( tempv(1)*xpt(1, :) + tempv(2)*xpt(2, :) + ... )
      do k = 1, npt
          s = s + tempv(k)*xpt(k, :)
      end do 
!----------------------------------------------------------------------!

!      ss = zero
!      ds = zero
!      do i = 1, n
!          ss = ss + s(i)**2
!          ds = ds + d(i)*s(i)
!      end do

      ss = dot_product(s, s)
      ds = dot_product(d, s)
      ssden = dd*ss - ds*ds
      if (ssden < 1.0e-8_rp*dd*ss) then 
          exit
      end if
      end do

      
      ! Set the vector WCHECK before the RETURN from the subroutine.
!      340 do k = 1, npt + n
!          wcheck(k) = zero
!          do j = 1, 5
!              wcheck(k) = wcheck(k) + wvec(k, j)*par(j)
!          end do
!      end do
      wcheck = matmul(wvec(:, 1 : 5), par(1 : 5))
      vlag(kopt) = vlag(kopt) + one

      return

      end subroutine bigden
