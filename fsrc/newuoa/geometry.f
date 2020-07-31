      module geometry_mod

      implicit none
      private
      public :: setremove, ameliorgeo


      contains

      subroutine setremove(idz, kopt, beta, delta, ratio, rho, vlag,    &
     & xopt, xpt, zmat, knew)
      ! SETREMOVE sets KNEW to the index of the interpolation point that
      ! will be deleted AFTER A TRUST REGION STEP.
      ! KNEW will be set in a way ensuring that the geometry of XPT is
      ! "optimal" after XPT(:, KNEW) is replaced by XNEW = XOPT + D,
      ! where D is the trust-region step.  Note that the information of
      ! XNEW is included in VLAG and BETA, which are calculated 
      ! according to D.

      use consts_mod, only : RP, IK, ONE, ZERO, TENTH
      use consts_mod, only : SRNLEN, DEBUG_MODE
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      ! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: beta
      real(RP), intent(in) :: delta 
      real(RP), intent(in) :: ratio 
      real(RP), intent(in) :: rho 
      real(RP), intent(in) :: vlag(:)  ! VLAG(NPT)
      real(RP), intent(in) :: xopt(:)  ! XOPT(N)
      real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

      ! Output
      integer(IK), intent(out) :: knew

      ! Intermediate variables
      integer(IK) :: n, npt
      real(RP) :: hdiag(size(zmat, 1)), xdsq(size(xpt, 2))
      real(RP) :: sigma(size(xpt, 2)), rhosq 
      character(len = SRNLEN), parameter :: srname = 'SETREMOVE'
     

      ! Get and verify the sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(vlag, npt)
          call verisize(xopt, n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
      end if

      rhosq = max(TENTH*delta, rho)**2
      hdiag = -sum(zmat(:, 1 : idz - 1)**2, dim = 2) +                  &
     & sum(zmat(:, idz : npt - n - 1)**2, dim = 2)             
      xdsq = sum((xpt - spread(xopt, dim = 2, ncopies = npt))**2, dim=1)
      sigma = abs(beta*hdiag + vlag(1 : npt)**2)
      sigma = sigma * max(xdsq/rhosq, ONE)**3
      if (ratio <= ZERO) then
      ! When the new F is not better than the current FOPT,
      ! we set SIGMA(KOPT) = -1 to prevent KNEW from being KOPT.
          sigma(kopt) = -ONE
      end if
      if (maxval(sigma) > ONE .or. ratio > ZERO) then
      ! KNEW > 0 unless MAXVAL(SIGMA) <= 1 and RATIO <= ZERO.
      ! If RATIO > ZERO (i.e., the new F is smaller than the current
      ! FOPT), then KNEW > 0, ensuring XNEW to be included into XPT.
          knew = int(maxloc(sigma, dim = 1), kind(knew))
      else
          knew = 0
      end if
      end subroutine setremove

      
      subroutine ameliorgeo(idz, knew, kopt, bmat, delbar, xopt, xpt,   &
     & zmat, d, beta, vlag)
      use consts_mod, only : RP, IK, ONE
      use consts_mod, only : DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      use vlagbeta_mod, only : vlagbeta

      ! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: knew
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
      real(RP), intent(in) :: delbar
      real(RP), intent(in) :: xopt(:)     ! XOPT(N)
      real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

      ! In-output
      real(RP), intent(inout) :: d(:)     ! D(N)

      ! Outputs
      real(RP), intent(out) :: beta
      real(RP), intent(out) :: vlag(:)    ! VLAG(NPT + N)

      ! Intermediate variables
      integer(IK) :: n, npt
      real(RP) :: alpha, zknew(size(zmat, 2))
      character(len = SRNLEN), parameter :: srname = 'AMELIORGEO'


      ! Get and verify the sizes.
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(xopt, n)
          call verisize(bmat, n, npt + n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(vlag, npt + n)
          call verisize(d, n)
      end if    

      call biglag(idz, knew, delbar, bmat, xopt, xpt, zmat, d)

      ! ALPHA is the KNEW-th diagonal entry of H
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      alpha = dot_product(zmat(knew, :), zknew)

      ! Calculate VLAG and BETA for D.
      call vlagbeta(idz, kopt, bmat, d, xopt,xpt,zmat,beta,vlag)

      ! If the cancellation in DENOM is unacceptable, then BIGDEN
      ! calculates an alternative model step D. 
      ! VLAG and BETA for this D are calculated within BIGDEN. 
      if (abs(ONE + alpha*beta/vlag(knew)**2) <= 0.8_RP) then
          call bigden(idz, knew, kopt, bmat, xopt, xpt,zmat,d,beta,vlag)
      end if

      return
      end subroutine ameliorgeo


      subroutine biglag(idz, knew, delbar, bmat, x, xpt, zmat, d)
      ! BIGLAG calculates a D by approximately solving
      !
      ! max |LFUNC(X + D)|, subject to ||D|| <= DELBAR, 
      !
      ! where LFUNC is the KNEW-th Lagrange function.

      use consts_mod, only : RP, IK, ONE, TWO, HALF, PI, ZERO
      use consts_mod, only : DEBUG_MODE, SRNLEN
      use lina_mod
      use warnerror_mod, only : errstop
      implicit none

      ! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: knew
      real(RP), intent(in) :: delbar 
      real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
      real(RP), intent(in) :: x(:)        ! X(N)
      real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

      ! Output
      real(RP), intent(out) :: d(:)       ! D(N)

      ! Intermediate variables
      integer(IK) :: i, isave, iterc, iu, n, npt
      real(RP) :: hcol(size(xpt, 2)), gc(size(x)), gd(size(x))
      real(RP) :: s(size(x)), w(size(x)), zknew(size(zmat, 2))
      real(RP) :: angle, cf(5), cth, dd, denom, dhd, gg, scaling
      real(RP) :: sp, ss, step, sth, tau, taubeg, tauold, taumax
      real(RP) :: unitang, taua, taub, t
      character(len = SRNLEN), parameter :: srname = 'BIGLAG'

       
      ! N is the number of variables.
      ! NPT is the number of interpolation equations.
      ! XPT contains the current interpolation points.
      ! BMAT provides the last N ROWs of H.
      ! ZMAT and IDZ give a factorization of the first NPT by NPT
      ! sub-matrix of H.
      ! KNEW is the index of the interpolation point to be removed.
      ! DELBAR is the trust region bound.
      ! D will be set to the step from X to the new point.
      ! HCOL, GC, GD, S and W will be used for working space.


      ! Get and verify the sizes.
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(x, n)
          call verisize(bmat, n, npt + n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
      end if

      ! Set HCOL to the leading NPT elements of the KNEW-th column of H.
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      hcol = matmul(zmat, zknew)

      ! Set the unscaled initial direction D. Form the gradient of LFUNC
      ! at X, and multiply D by the Hessian of LFUNC.
      d = xpt(:, knew) - x
      dd = dot_product(d, d)

      gd = matmul(xpt, hcol*matmul(d, xpt))

      !----------------------------------------------------------------!
!-----!gc = bmat(:, knew) + matmul(xpt, hcol*matmul(x, xpt)) !---------!
      gc = Ax_plus_y(xpt, hcol*matmul(x, xpt), bmat(:, knew))
      !----------------------------------------------------------------!

      ! Scale D and GD, with a sign change if required. Set S to another
      ! vector in the initial two dimensional subspace.
      gg = dot_product(gc, gc)
      sp = dot_product(d, gc)
      dhd = dot_product(d, gd)
      scaling = delbar/sqrt(dd)
      if (sp*dhd < ZERO) then 
          scaling = - scaling
      end if
      t = ZERO
      if (sp*sp > 0.99_RP*dd*gg) then 
          t = ONE 
      end if
      tau = scaling*(abs(sp) + HALF*scaling*abs(dhd))
      if (gg*(delbar*delbar) < 1.0e-2_RP*tau*tau) then 
          t = ONE
      end if
      d = scaling*d
      gd = scaling*gd
      s = gc + t*gd
      
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
          unitang = (TWO*PI)/real(iu + 1, RP)

          do i = 1, iu
              angle = real(i, RP)*unitang
              cth = cos(angle)
              sth = sin(angle)
              tau = cf(1) + (cf(2)+cf(4)*cth)*cth+(cf(3)+cf(5)*cth)*sth
              if (abs(tau) > abs(taumax)) then
                  taumax = tau
                  isave = i
                  taua = tauold
              else if (i == isave + 1) then
                  taub = tau
              end if
              tauold = tau
          end do

          if (isave == 0) then 
              taua = tau
          end if
          if (isave == iu) then 
              taub = taubeg
          end if
          if (abs(taua - taub) > ZERO) then
              taua = taua - taumax
              taub = taub - taumax
              step = HALF*(taua - taub)/(taua + taub)
          else
              step = ZERO
          end if
          angle = unitang*(real(isave, RP) + step)
          
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


      subroutine bigden(idz, knew, kopt, bmat, x, xpt, zmat,d,beta,vlag)
      ! BIGDEN calculates a D by approximately solving
      !
      ! max |SIGMA(X + D)|, subject to ||D|| <= DELBAR, 
      !
      ! where SIGMA is the denominator \sigma in the updating formula
      ! (4.11)--(4.12) for H, which is the inverse of the coefficient
      ! matrix for the interplolation system (see (3.12)). Indeed, each
      ! column of H corresponds to a Lagrange basis function of the
      ! interpolation problem. 
      ! In addition, it sets VLAG and BETA for the selected D.

      use consts_mod, only : RP, IK, ONE, TWO, HALF, QUART, PI, ZERO
      use consts_mod, only : DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      ! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: knew
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
      real(RP), intent(in) :: x(:)        ! X(N)
      real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

      ! In-output
      real(RP), intent(inout) :: d(:)     ! D(N)

      ! Outputs
      real(RP), intent(out) :: beta
      real(RP), intent(out) :: vlag(:)    ! VLAG(NPT + N)

      ! Intermediate variable
      integer(IK) :: i, isave, iterc, iu, j, jc, k, nw, n, npt
      real(RP) :: s(size(x)) 
      real(RP) :: w(size(xpt, 2) + size(x), 5)
      real(RP) :: prod(size(xpt, 2) + size(x), 5)
      real(RP) :: den(9), denex(9), par(9)       
      real(RP) :: zknew(size(zmat, 2)), dstemp(size(xpt, 2))
      real(RP) :: sstemp(size(xpt, 2)), wz(size(zmat, 2))
      real(RP) :: hcol(size(xpt, 2)), xnew(size(x))
      real(RP) :: alpha, angle, unitang, dd, denmax, denold, densav
      real(RP) :: ds, dtest, ss, ssden, denom, denomold, dena, denb
      real(RP) :: step, tau, tempa, tempb, tempc, v(size(xpt, 2))
      real(RP) :: xd, xs, xsq, dxn, xnsq
      real(RP) :: xptemp(size(xpt, 1), size(xpt, 2))
      character(len = SRNLEN), parameter :: srname = 'BIGDEN'

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
      ! VLAG will be set to Theta*WCHECK+e_b for the final choice of D.
      ! BETA will be set to the value that will occur in the updating
      ! formula when the KNEW-th interpolation point is moved to its new
      ! position.
      
      ! D is calculated in a way that should provide a denominator with 
      ! a large modulus in the updating formula when the KNEW-th
      ! interpolation point is shifted to the new position X + D.

       
      ! Get and verify the sizes.
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(x, n)
          call verisize(bmat, n, npt + n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(vlag, npt + n)
          call verisize(d, n)
      end if    

      ! Store the first NPT elements of the KNEW-th column of H in HCOL.
      zknew = zmat(knew, :)
      zknew(1 : idz - 1) = -zknew(1 : idz - 1)
      hcol = matmul(zmat, zknew)
      alpha = hcol(knew)
      
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
          xptemp = xpt - spread(x, dim = 2, ncopies = npt)
!----------------------------------------------------------------------!
!---------!dstemp = matmul(d, xpt) - dot_product(x, d) !---------------!
          dstemp = matmul(d, xptemp)
!----------------------------------------------------------------------!
          sstemp = sum((xptemp)**2, dim = 1) 
          
          dstemp(kopt) = TWO*ds + ONE  
          sstemp(kopt) = ss     
          k = int(minloc(dstemp*dstemp/sstemp, dim = 1), kind(k))
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
          
          ! Put the coefficients of WCHECK in W.
          do k = 1, npt
              tempa = dot_product(xpt(:, k), d)
              tempb = dot_product(xpt(:, k), s)
              tempc = dot_product(xpt(:, k), x)
              w(k, 1) = QUART*(tempa*tempa + tempb*tempb)
              w(k, 2) = tempa*tempc
              w(k, 3) = tempb*tempc
              w(k, 4) = QUART*(tempa*tempa - tempb*tempb)
              w(k, 5) = HALF*tempa*tempb
          end do
          w(npt + 1 : npt + n, 1 : 5) = ZERO
          w(npt + 1 : npt + n, 2) = d
          w(npt + 1 : npt + n, 3) = s
    
          ! Put the coefficents of THETA*WCHECK in PROD.
          do jc = 1, 5
              wz = matmul(w(1 : npt, jc), zmat)
              wz(1 : idz - 1) = -wz(1 : idz - 1)
              prod(1 : npt, jc) = matmul(zmat, wz)
              
              nw = npt
              if (jc == 2 .or. jc == 3) then
                  prod(1 : npt, jc) = prod(1 : npt, jc) +               &
     &             matmul(w(npt + 1 : npt + n, jc), bmat(:, 1 : npt))
                  nw = npt + n
              end if
              prod(npt + 1 : npt + n, jc) = matmul(bmat(:, 1 : nw),     &
     &         w(1 : nw, jc))
          end do
    
          ! Include in DEN the part of BETA that depends on THETA.
          do k = 1, npt + n
              par(1 : 5) = HALF*prod(k, 1 : 5)*w(k, 1 : 5)
              den(1) = den(1) - par(1) - sum(par(1 : 5)) 
              tempa = prod(k, 1)*w(k, 2) + prod(k, 2)*w(k, 1)
              tempb = prod(k, 2)*w(k, 4) + prod(k, 4)*w(k, 2)
              tempc = prod(k, 3)*w(k, 5) + prod(k, 5)*w(k, 3)
              den(2) = den(2) - tempa - HALF*(tempb + tempc)
              den(6) = den(6) - HALF*(tempb-tempc)
              tempa = prod(k, 1)*w(k, 3) + prod(k, 3)*w(k, 1)
              tempb = prod(k, 2)*w(k, 5) + prod(k, 5)*w(k, 2)
              tempc = prod(k, 3)*w(k, 4) + prod(k, 4)*w(k, 3)
              den(3) = den(3) - tempa - HALF*(tempb - tempc)
              den(7) = den(7) - HALF*(tempb + tempc)
              tempa = prod(k, 1)*w(k, 4) + prod(k, 4)*w(k, 1)
              den(4) = den(4) - tempa - par(2) + par(3)
              tempa = prod(k, 1)*w(k, 5) + prod(k, 5)*w(k, 1)
              tempb = prod(k, 2)*w(k, 3) + prod(k, 3)*w(k, 2)
              den(5) = den(5) - tempa - HALF*tempb
              den(8) = den(8) - par(4) + par(5)
              tempa = prod(k, 4)*w(k, 5) + prod(k, 5)*w(k, 4)
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
          denom = denex(1) + denex(2) + denex(4) + denex(6) + denex(8)
          denold = denom
          denmax = denom
          isave = 0
          iu = 49
          unitang = (TWO*PI)/real(iu + 1, RP)
          par(1) = ONE
          do i = 1, iu
              angle = real(i, RP)*unitang
              par(2) = cos(angle)
              par(3) = sin(angle)
              do j = 4, 8, 2
                  par(j) = par(2)*par(j - 2) - par(3)*par(j - 1)
                  par(j + 1) = par(2)*par(j - 1) + par(3)*par(j - 2)
              end do
              denomold = denom
              denom = dot_product(denex(1 : 9), par(1 : 9))
              if (abs(denom) > abs(denmax)) then
                  denmax = denom
                  isave = i
                  dena = denomold
              else if (i == isave + 1) then
                  denb = denom
              end if
          end do
          if (isave == 0) then
              dena = denom
          end if
          if (isave == iu) then 
              denb = denold
          end if
          if (abs(dena - denb) > 0) then
              dena = dena - denmax
              denb = denb - denmax
              step = HALF*(dena - denb)/(dena + denb)
          else
              step = ZERO
          end if
          angle = unitang*(real(isave, RP) + step)
          
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
          xnew = x + d
          dxn = dot_product(d, xnew)
          xnsq = dot_product(xnew, xnew)
    
          if (iterc > 1) then
              densav = max(densav, denold)
          end if
          if (abs(denmax) <= 1.1_RP*abs(densav)) then 
              exit
          end if
          densav = denmax
    
          ! Set S to HALF the gradient of the denominator with respect 
          ! to D. 
          s = tau*bmat(:,knew) + alpha*(dxn*x+xnsq*d-vlag(npt+1:npt+n))
          v = matmul(xnew, xpt)
          v = (tau*hcol - alpha*vlag(1 : npt))*v
!----------------------------------------------------------------------!
!---------!s = s + matmul(xpt, v) !--------------------------------!
          s = Ax_plus_y(xpt, v, s)
!----------------------------------------------------------------------!
    
          ss = dot_product(s, s)
          ds = dot_product(d, s)
          ssden = dd*ss - ds*ds
          if (ssden < 1.0e-8_RP*dd*ss) then 
              exit
          end if
      end do
      
      vlag(kopt) = vlag(kopt) + ONE

      return

      end subroutine bigden

      end module geometry_mod
