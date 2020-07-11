      module geometry

      implicit none
      private
      public :: biglag, bigden

      contains

      subroutine biglag(x, xpt, bmat, zmat, idz, knew, delta, d, alpha)
      ! BIGLAG calculates a D by approximately solving
      !
      ! max |LFUNC(X + D)|, subject to ||D|| <= DELTA, 
      !
      ! where LFUNC is the KNEW-th Lagrange function.
      ! In addition, it sets ALPHA for the selected D.

      use consts_mod, only : RP, IK, ONE, TWO, HALF, PI, ZERO,DEBUG_MODE
      use lina_mod
      use warnerror_mod, only : errmssg
      implicit none

      integer(IK), intent(in) ::    idz
      integer(IK), intent(in) ::    knew
      real(RP), intent(in) ::   x(:)        ! X(N)
      real(RP), intent(in) ::   xpt(:, :)   ! XPT(N, NPT)
      real(RP), intent(in) ::   bmat(:, :)  ! BMAT(N, NPT + N)
      real(RP), intent(in) ::   zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
      real(RP), intent(in) ::   delta 
      real(RP), intent(out) ::  alpha
      real(RP), intent(out) ::  d(:)        ! D(N)

      integer(IK) :: i, isave, iterc, iu, k, n, npt
      real(RP) :: hcol(size(xpt, 2)), gc(size(x)), gd(size(x)),         &
     & s(size(x)), w(size(x)), zknew(size(zmat, 2)), angle, cf(5), cth, &
     & dd, denom, dhd, gg, scaling, sp, ss, step, sth, tau, taubeg,     &
     & tauold, taumax, temp, tempa, tempb
      character(len = 100) :: srname

       
      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N ROWs of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! sub-matrix of H.
      ! KNEW is the index of the interpolation point to be removed.
      ! DELTA is the current trust region bound.
      ! D will be set to the step from X to the new point.
      ! ALPHA will be set to the KNEW-th diagonal element of matrix H.
      ! HCOL, GC, GD, S and W will be used for working space.


      srname = 'BIGLAG'  ! Name of the current subroutine.
       
      ! Get and verify the sizes.
      n = size(xpt, 1)
      npt = size(xpt, 2)

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errmssg(srname, 'SIZE(XPT) is invalid')
              stop
          end if
          if (size(x) /= n) then
              call errmssg(srname, 'SIZE(X) /= SIZE(XPT, 1)')
              stop
          end if
          if (size(bmat, 1) /= n .or. size(bmat, 2) /= npt + n) then
              call errmssg(srname, 'SIZE(BMAT) is invalid')
              stop
          end if
          if (size(zmat, 1) /= npt .or. size(zmat, 2) /= npt-n-1) then
              call errmssg(srname, 'SIZE(ZMAT) is invalid')
              stop
          end if
      end if

      ! Set HCOL to the leading NPT elements of the KNEW-th column of H.
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      hcol = matmul(zmat, zknew)
      alpha = hcol(knew)

      ! Set the unscaled initial direction D. Form the gradient of LFUNC
      ! at X, and multiply D by the Hessian of LFUNC.
      d = xpt(:, knew) - x
      dd = dot_product(d, d)

      gd = matmul(xpt, hcol*matmul(d, xpt))

      !----------------------------------------------------------------!
      ! The following DO LOOP calculates the GC below
!-----!gc = bmat(:, knew) + matmul(xpt, hcol*matmul(x, xpt)) !---------!
      gc = bmat(:, knew)
      do k = 1, npt
          gc = gc + (hcol(k)*dot_product(xpt(:, k), x))*xpt(:, k)
      end do
      !----------------------------------------------------------------!

      ! Scale D and GD, with a sign change if required. Set S to another
      ! vector in the initial two dimensional subspace.
      gg = dot_product(gc, gc)
      sp = dot_product(d, gc)
      dhd = dot_product(d, gd)
      scaling = delta/sqrt(dd)
      if (sp*dhd < ZERO) then 
          scaling = - scaling
      end if
      temp = ZERO
      if (sp*sp > 0.99_RP*dd*gg) then 
          temp = ONE 
      end if
      tau = scaling*(abs(sp) + HALF*scaling*abs(dhd))
      if (gg*(delta*delta) < 1.0e-2_RP*tau*tau) then 
          temp = ONE
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
          if (dd*ss - sp*sp <= 1.0e-8_RP*dd*ss) then 
              exit
          end if
          denom = sqrt(dd*ss - sp*sp)
          s = (dd*s - sp*d)/denom

          w = matmul(xpt, hcol*matmul(s, xpt))
!          w = ZERO
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
          cf(1) = HALF*cf(1)
          cf(4) = HALF*cf(4) - cf(1)
          
          ! Seek the value of the angle that maximizes |TAU|.
          taubeg = cf(1) + cf(2) + cf(4)
          taumax = taubeg
          tauold = taubeg
          isave = 0
          iu = 49
          temp = (TWO*PI)/real(iu + 1, RP)
          do i = 1, iu
              angle = real(i, RP)*temp
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
          if (abs(tempa - tempb) > ZERO) then
              tempa = tempa - taumax
              tempb = tempb - taumax
              step = HALF*(tempa - tempb)/(tempa + tempb)
          else
              step = ZERO
          end if
          angle = temp*(real(isave, RP) + step)
          
          ! Calculate the new D and GD. Then test for convergence.
          cth = cos(angle)
          sth = sin(angle)
          tau = cf(1) + (cf(2) + cf(4)*cth)*cth + (cf(3)+cf(5)*cth)*sth
          d = cth*d + sth*s
          gd = cth*gd + sth*w
          s = gc + gd
          if (abs(tau) <= 1.1_RP*abs(taubeg)) then
              exit
          end if
      end do

      return

      end subroutine biglag



      subroutine bigden (x, xpt, bmat, zmat, idz, kopt, knew, d, wcheck,&
     & vlag, beta)
!!! TODO: Replace the DO LOOP in the update of S with MATMUL.
      ! BIGDEN calculates a D by approximately solving
      !
      ! max |SIGMA(X + D)|, subject to ||D|| <= DELTA, 
      !
      ! where SIGMA is the denominator \sigma in the updating formula
      ! (4.11)--(4.12) for H, which is the inverse of the coefficient
      ! matrix for the interplolation system (see (3.12)). Indeed, each
      ! column of H corresponds to a Lagrange basis function of the
      ! interpolation problem. 
      ! In addition, it sets VLAG, BETA, and WCHECK for the selected D.

      use consts_mod, only : RP, IK, ONE, TWO, HALF, QUART, PI, ZERO,       &
     & DEBUG_MODE
      use warnerror_mod, only : errmssg
      use lina_mod
      implicit none

      integer(IK), intent(in) ::        knew
      integer(IK), intent(in) ::        kopt
      integer(IK), intent(in) ::        idz
      real(RP), intent(in) ::       x(:)        ! X(N)
      real(RP), intent(in) ::       xpt(:, :)   ! XPT(N, NPT)
      real(RP), intent(in) ::       bmat(:, :)  ! BMAT(N, NPT+N)
      real(RP), intent(in) ::       zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
      real(RP), intent(out) ::      beta
      real(RP), intent(out) ::      vlag(:)     ! VLAG(NPT + N)
      real(RP), intent(out) ::      wcheck(:)   ! WCHECK(NPT + N) 
      real(RP), intent(inout) ::    d(:)        ! D(N)

      integer(IK) :: i, isave, iterc, iu, j, jc, k, nw, n, npt
      real(RP) :: s(size(x)), wvec(size(xpt, 2) + size(x), 5),          &
     & prod(size(xpt, 2) + size(x), 5), den(9), denex(9), par(9),       
     & zknew(size(zmat, 2)), stemp(size(x)), dstemp(size(xpt, 2)),      &
     & sstemp(size(xpt, 2)), wz(size(zmat, 2)), w1(size(xpt, 2)),       &
     & w2(size(x)), angle, dd, denmax, denold, densav, ds, dtest, ss,   &
     & ssden, summation, sumold, step, tau, temp, tempa, tempb, tempc,  &
     & tempv(size(xpt, 2)), xd, xs, xsq, alpha
      character(len = 100) :: srname

      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! X is the best interpolation point so far.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N ROWs of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! sub-matrix of H.
      ! NDIM is the second dimension of BMAT and has the value NPT + N.
      ! KOPT is the index of the optimal interpolation point.
      ! KNEW is the index of the interpolation point to be removed.
      ! D will be set to the step from X to the new point, and on 
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
      ! interpolation point is shifted to the new position X + D.

      srname = 'BIGDEN'  ! Name of the current subroutine.
       
      ! Get and verify the sizes.
      n = size(xpt, 1)
      npt = size(xpt, 2)

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errmssg(srname, 'SIZE(XPT) is invalid')
              stop
          end if
          if (size(x) /= n) then
              call errmssg(srname, 'SIZE(X) /= SIZE(XPT, 1)')
              stop
          end if
          if (size(bmat, 1) /= n .or. size(bmat, 2) /= npt + n) then
              call errmssg(srname, 'SIZE(BMAT) is invalid')
              stop
          end if
          if (size(zmat, 1) /= npt .or. size(zmat, 2) /= npt-n-1) then
              call errmssg(srname, 'SIZE(ZMAT) is invalid')
              stop
          end if
          if (size(vlag) /= npt + n) then
              call errmssg(srname, 'SIZE(VLAG) is invalid')
              stop
          end if
          if (size(wcheck) /= npt + n) then
              call errmssg(srname, 'SIZE(WCHECK) is invalid')
              stop
          end if
          if (size(d) /= n) then
              call errmssg(srname, 'SIZE(D) /= SIZE(X)')
              stop
          end if
      end if    

      ! Store the first NPT elements of the KNEW-th column of H in 
      ! WCHECK(N + 1) to WCHECK(N + NPT).
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      w1 = matmul(zmat, zknew)
      alpha = w1(knew)
      
      ! The initial search direction D is taken from the last call of
      ! BIGLAG, and the initial S is set below, usually to the direction
      ! from X to X_KNEW, but a different direction to an 
      ! interpolation point may be chosen, in order to prevent S from
      ! being nearly parallel to D.
      dd = dot_product(d, d)
      s = xpt(:, knew) - x
      ds = dot_product(d, s)
      ss = dot_product(s, s)
      xsq = dot_product(x, x)

      if (.not. (ds*ds <= 0.99_RP*dd*ss)) then
          dtest = ds*ds/ss
          ! The following DO LOOP implements the code below.
          !dstemp = matmul(d, xpt) - dot_product(x, d)
          !sstemp = sum((xpt-spread(x, dim=2, ncopies=npt))**2,dim=1) 
          do k = 1, npt
              stemp = xpt(:, k) - x
              dstemp(k) = dot_product(d, stemp)
              sstemp(k) = dot_product(stemp, stemp)
          end do
          
          dstemp(kopt) = TWO*ds + ONE  
          sstemp(kopt) = ss     
          k = minloc(dstemp*dstemp/sstemp, dim = 1)
          if ((.not. (dstemp(k)*dstemp(k)/sstemp(k) >= dtest)) .and.    &
     &     k /= kopt) then
          ! Althoguh unlikely, if NaN occurs, it may happen that k=kopt.
              s = xpt(:, k) - x
              ds = dstemp(k)
              ss = sstemp(k)
          end if
      end if

      ssden = dd*ss - ds*ds
      densav = ZERO

      ! Begin the iteration by overwriting S with a vector that has the
      ! required length and direction.
      do iterc = 1, n
          s = (ONE/sqrt(ssden))*(dd*s - ds*d)
          xd = dot_product(x, d)
          xs = dot_product(x, s)
    
          ! Set the coefficients of the first two terms of BETA.
          tempa = HALF*xd*xd
          tempb = HALF*xs*xs
          den(1) = dd*(xsq + HALF*dd) + tempa + tempb
          den(2) = TWO*xd*dd
          den(3) = TWO*xs*dd
          den(4) = tempa - tempb
          den(5) = xd*xs
          den(6 : 9) = ZERO
          
          ! Put the coefficients of WCHECK in WVEC.
          do k = 1, npt
              tempa = dot_product(xpt(:, k), d)
              tempb = dot_product(xpt(:, k), s)
              tempc = dot_product(xpt(:, k), x)
              wvec(k, 1) = quart*(tempa*tempa + tempb*tempb)
              wvec(k, 2) = tempa*tempc
              wvec(k, 3) = tempb*tempc
              wvec(k, 4) = quart*(tempa*tempa - tempb*tempb)
              wvec(k, 5) = HALF*tempa*tempb
          end do
          wvec(npt + 1 : npt + n, 1 : 5) = ZERO
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
              par(1 : 5) = HALF*prod(k, 1 : 5)*wvec(k, 1 : 5)
              den(1) = den(1) - par(1) - sum(par(1 : 5)) 
              tempa = prod(k, 1)*wvec(k, 2) + prod(k, 2)*wvec(k, 1)
              tempb = prod(k, 2)*wvec(k, 4) + prod(k, 4)*wvec(k, 2)
              tempc = prod(k, 3)*wvec(k, 5) + prod(k, 5)*wvec(k, 3)
              den(2) = den(2) - tempa-HALF*(tempb + tempc)
              den(6) = den(6) - HALF*(tempb-tempc)
              tempa = prod(k, 1)*wvec(k, 3) + prod(k, 3)*wvec(k, 1)
              tempb = prod(k, 2)*wvec(k, 5) + prod(k, 5)*wvec(k, 2)
              tempc = prod(k, 3)*wvec(k, 4) + prod(k, 4)*wvec(k, 3)
              den(3) = den(3) - tempa-HALF*(tempb-tempc)
              den(7) = den(7) - HALF*(tempb + tempc)
              tempa = prod(k, 1)*wvec(k, 4) + prod(k, 4)*wvec(k, 1)
              den(4) = den(4) - tempa-par(2) + par(3)
              tempa = prod(k, 1)*wvec(k, 5) + prod(k, 5)*wvec(k, 1)
              tempb = prod(k, 2)*wvec(k, 3) + prod(k, 3)*wvec(k, 2)
              den(5) = den(5) - tempa-HALF*tempb
              den(8) = den(8) - par(4) + par(5)
              tempa = prod(k, 4)*wvec(k, 5) + prod(k, 5)*wvec(k, 4)
              den(9) = den(9) - HALF*tempa
          end do
          
          par(1 : 5) = HALF*prod(knew, 1 : 5)**2
          denex(1) = alpha*den(1) + par(1) + sum(par(1 : 5))
          tempa = TWO*prod(knew, 1)*prod(knew, 2)
          tempb = prod(knew, 2)*prod(knew, 4)
          tempc = prod(knew, 3)*prod(knew, 5)
          denex(2) = alpha*den(2) + tempa + tempb + tempc
          denex(6) = alpha*den(6) + tempb - tempc
          tempa = TWO*prod(knew, 1)*prod(knew, 3)
          tempb = prod(knew, 2)*prod(knew, 5)
          tempc = prod(knew, 3)*prod(knew, 4)
          denex(3) = alpha*den(3) + tempa + tempb - tempc
          denex(7) = alpha*den(7) + tempb + tempc
          tempa = TWO*prod(knew, 1)*prod(knew, 4)
          denex(4) = alpha*den(4) + tempa + par(2) - par(3)
          tempa = TWO*prod(knew, 1)*prod(knew, 5)
          denex(5) = alpha*den(5) + tempa + prod(knew, 2)*prod(knew, 3)
          denex(8) = alpha*den(8) + par(4) - par(5)
          denex(9) = alpha*den(9) + prod(knew, 4)*prod(knew, 5)
          
          ! Seek the value of the angle that maximizes the |DENOM|.
          summation = denex(1) + denex(2) + denex(4) + denex(6)+denex(8)
          denold = summation
          denmax = summation
          isave = 0
          iu = 49
          temp = (TWO*PI)/real(iu + 1, RP)
          par(1) = ONE
          do i = 1, iu
              angle = real(i, RP)*temp
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
              step = HALF*(tempa - tempb)/(tempa + tempb)
          else
              step = ZERO
          end if
          angle = temp*(real(isave, RP) + step)
          
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
          w2 = x + d
          tempa = dot_product(d, w2)
          tempb = dot_product(w2, w2)
          
    
          if (iterc > 1) then
              densav = max(densav, denold)
          end if
          if (abs(denmax) <= 1.1_RP*abs(densav)) then 
              exit
          end if
          densav = denmax
    
          ! Set S to HALF the gradient of the denominator with respect 
          ! to D. Then branch for the next iteration.
          s = tau*bmat(:,knew)+alpha*(tempa*x+tempb*d-vlag(npt+1:npt+n))
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
          if (ssden < 1.0e-8_RP*dd*ss) then 
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
!-----!wcheck = wcheck*(HALF*wcheck + matmul(x, xpt)) !----------------!
      ! The result here is likely different from the theoretical value.
!----------------------------------------------------------------------!
      wcheck = matmul(wvec(:, 1 : 5), par(1 : 5))
      vlag(kopt) = vlag(kopt) + ONE

      return

      end subroutine bigden

      end module geometry
