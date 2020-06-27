!!! TODO: Replace the DO LOOP in the update of S with MATMUL.
      subroutine bigden (n, npt, xopt, xpt, bmat, zmat, idz, kopt, knew,&
     & d, wcheck, vlag, beta)
      ! BIGDEN calculates a D by approximately solving
      !
      ! max |LFUNC(XOPT + D)|, subject to ||D|| <= DELTA, 
      !
      ! where LFUNC is the KNEW-th Lagrange function.
      ! In addition, it sets VLAG, BETA, and WCHECK for the selected D.

      use consts, only : rp, one, two, half, quart, pi, zero
      use lina
      implicit none

      integer, intent(in) :: n, npt, knew, kopt, idz

      real(kind=rp), intent(in) :: xopt(n), xpt(n, npt), bmat(n, npt+n),&
     & zmat(npt, npt - n - 1)
      real(kind=rp), intent(out) :: beta, vlag(npt + n), wcheck(npt + n) 
      real(kind=rp), intent(inout) :: d(n)

      integer :: i, isave, iterc, iu, j, jc, k, nw
      real(kind=rp) :: s(n), wvec(npt + n, 5), prod(npt + n, 5), den(9),&
     & denex(9), par(9), zknew(npt - n - 1), stemp(n), dstemp(npt),     &
     & sstemp(npt), wz(npt - n - 1), angle, dd, denmax, denold, densav, &
     & ds, dtest, ss, ssden, summation, sumold, step, tau, temp, tempa, &
     & tempb, tempc, tempv(npt), xoptd, xopts, xoptsq, alpha,           &
     & w1(npt), w2(n)
            

      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! XOPT is the best interpolation point so far.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N ROWs of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! sub-matrix of H.
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
      w1 = matmul(zmat, zknew)
      alpha = w1(knew)
      
      ! The initial search direction D is taken from the last call of
      ! BIGLAG, and the initial S is set below, usually to the direction
      ! from XOPT to X_KNEW, but a different direction to an 
      ! interpolation point may be chosen, in order to prevent S from
      ! being nearly parallel to D.
      dd = dot_product(d, d)
      s = xpt(:, knew) - xopt
      ds = dot_product(d, s)
      ss = dot_product(s, s)
      xoptsq = dot_product(xopt, xopt)

      if (.not. (ds*ds <= 0.99_rp*dd*ss)) then
          dtest = ds*ds/ss
          ! The following DO LOOP implements the code below.
          !dstemp = matmul(d, xpt) - dot_product(xopt, d)
          !sstemp = sum((xpt - spread(xopt, dim = 2, ncopies = npt))**2, dim = 1) 
          do k = 1, npt
              stemp = xpt(:, k) - xopt
              dstemp(k) = dot_product(d, stemp)
              sstemp(k) = dot_product(stemp, stemp)
          end do
          
          dstemp(kopt) = two*ds + one  
          sstemp(kopt) = ss     
          k = minloc(dstemp*dstemp/sstemp, dim = 1)
          if ((.not. (dstemp(k)*dstemp(k)/sstemp(k) >= dtest)) .and.    &
     &     k /= kopt) then
          ! Althoguh unlikely, if NaN occurs, it may happen that k=kopt.
              s = xpt(:, k) - xopt
              ds = dstemp(k)
              ss = sstemp(k)
          end if
      end if

      ssden = dd*ss - ds*ds
      densav = zero

      ! Begin the iteration by overwriting S with a vector that has the
      ! required length and direction.
      do iterc = 1, n
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
              tempa = dot_product(xpt(:, k), d)
              tempb = dot_product(xpt(:, k), s)
              tempc = dot_product(xpt(:, k), xopt)
              wvec(k, 1) = quart*(tempa*tempa + tempb*tempb)
              wvec(k, 2) = tempa*tempc
              wvec(k, 3) = tempb*tempc
              wvec(k, 4) = quart*(tempa*tempa - tempb*tempb)
              wvec(k, 5) = half*tempa*tempb
          end do
          wvec(npt + 1 : npt + n, 1 : 5) = zero
          wvec(npt + 1 : npt + n, 2) = d
          wvec(npt + 1 : npt + n, 3) = s
    
          ! Put the coefficents of THETA*WCHECK in PROD.
          do jc = 1, 5
              wz = matmul(wvec(1 : npt, jc), zmat)
              wz(1 : idz - 1) = -wz(1 : idz - 1)
              prod(1 : npt, jc) = matmul(zmat, wz)
              
              nw = npt
              if (jc == 2 .or. jc == 3) then
                  prod(1 : npt, jc) = prod(1 : npt, jc) +               &
     &             matmul(wvec(npt + 1 : npt + n, jc), bmat(:, 1 : npt))
                  nw = npt + n
              end if
              prod(npt + 1 : npt + n, jc) = matmul(bmat(:, 1 : nw),     &
     &         wvec(1 : nw, jc))
          end do
    
          ! Include in DEN the part of BETA that depends on THETA.
          do k = 1, npt + n
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
          
          ! Seek the value of the angle that maximizes the |DENOM|.
          summation = denex(1) + denex(2) + denex(4) + denex(6)+denex(8)
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
          
          ! Calculate the new parameters of the denominator, the new
          ! VLAG vector and the new D. Then test for convergence.
          par(2) = cos(angle)
          par(3) = sin(angle)
          do j = 4, 8, 2
              par(j) = par(2)*par(j - 2) - par(3)*par(j - 1)
              par(j + 1) = par(2)*par(j - 1) + par(3)*par(j - 2)
          end do
    
          beta = dot_product(den(1 : 9), par(1 : 9))
          denmax = dot_product(denex(1 : 9), par(1 : 9))
    
          vlag = matmul(prod(:, 1 : 5), par(1 : 5))
    
          tau = vlag(knew)
    
          d = par(2)*d + par(3)*s
    
    
          dd = dot_product(d, d)
          w2 = xopt + d
          tempa = dot_product(d, w2)
          tempb = dot_product(w2, w2)
          
    
          if (iterc > 1) then
              densav = max(densav, denold)
          end if
          if (abs(denmax) <= 1.1_rp*abs(densav)) then 
              exit
          end if
          densav = denmax
    
          ! Set S to half the gradient of the denominator with respect 
          ! to D. Then branch for the next iteration.
          s = tau*bmat(:, knew) + alpha*(tempa*xopt + tempb*d -         &
     &     vlag(npt+1:npt+n))
          tempv = matmul(w2, xpt)
          tempv = (tau*w1 - alpha*vlag(1 : npt))*tempv
!----------------------------------------------------------------------!
          !!! In later versions, the following DO LOOP should be
          !!! replaced by MATMUL as follows:
!---------!s = s + matmul(xpt, tempv) !--------------------------------!
          !!! For the moment, we use a DO LOOP instead of MATMUL to
          !!! ensure that the result is identical to that of Powell's
          !!! code. If we use MATMUL, the result is not always the same,
          !!! since floating point addition (and multiplication) is not
          !!! associative, although it is commutative. The following 
          !!! LOOP calculates 
          !!! s + tempv(1)*xpt(:, 1) + tempv(2)*xpt(:, 2) + ...
          !!! while s + matmul(xpt, tempv) calculates
          !!! s + ( tempv(1)*xpt(:, 1) + tempv(2)*xpt(:, 2) + ... )
          do k = 1, npt
              s = s + tempv(k)*xpt(:, k)
          end do 
!----------------------------------------------------------------------!
    
          ss = dot_product(s, s)
          ds = dot_product(d, s)
          ssden = dd*ss - ds*ds
          if (ssden < 1.0e-8_rp*dd*ss) then 
              exit
          end if
      end do
      
      ! Set the vector WCHECK before the RETURN from the subroutine.
!----------------------------------------------------------------------!
      ! This is the one of the two places where WCHECK is calculated,
      ! the other being VLAGBETA. 
      ! WCHECK contains the first NPT entries of (w-v) for the vectors 
      ! w and v defined in eq(4.10) and eq(4.24) of the NEWUOA paper,
      ! and also \hat{w} in eq(6.5) of 
      !
      ! M. J. D. Powell, Least Frobenius norm updating of quadratic
      ! models that satisfy interpolation conditions. Math. Program.,
      ! 100:183--215, 2004
      !
      ! WCHECK is used ONLY in CALQUAD, which evaluates the qudratic
      ! model. Indeed, we may calculate WCHECK internally in CALQUAD.
      !
      ! WCHECK is the following vector in theory. 
!-----!wcheck = matmul(d, xpt) !---------------------------------------!
!-----!wcheck = wcheck*(half*wcheck + matmul(xopt, xpt)) !----------------!
      ! The result here is likely different from the theoretical value.
!----------------------------------------------------------------------!
      wcheck = matmul(wvec(:, 1 : 5), par(1 : 5))
      vlag(kopt) = vlag(kopt) + one

      return

      end subroutine bigden
