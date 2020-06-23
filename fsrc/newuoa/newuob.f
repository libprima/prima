      subroutine newuob (n, npt, x, rhobeg, rhoend, iprint, maxfun,     &
     & xbase, xopt, xnew, xpt, fval, gq, hq, pq, bmat, zmat,            &
     & d, vlag, w, f, info, ftarget)

      use consts, only : rp, zero, one, half, tenth
      use infnan
      implicit none

      ! inputs
      integer, intent(in) :: n, npt, iprint, maxfun
      integer, intent(out) :: info
      real(kind = rp), intent(in) :: rhobeg, rhoend, ftarget
      real(kind = rp), intent(out) :: f
      real(kind = rp), intent(inout) :: x(n), xbase(n), xopt(n),        &
     & xnew(n), xpt(npt, n), fval(npt), gq(n), hq((n*(n+1))/2), pq(npt)
      real(kind = rp), intent(inout) :: bmat(npt + n, n),               &
     & zmat(npt, npt - n - 1), d(n), vlag(npt + n), w(10*(npt + n))

      ! other variables
      integer :: i, idz, itest, k, knew, kopt, nf, nfsave, subinfo
      integer :: tr, maxtr
      real(kind = rp) :: alpha, beta, crvmin, delta, prederr(3),        &
     & distsq, dnorm, dsq, dstep, xdiff(n), xdsq(npt)
      real(kind = rp) :: fopt, fsave, galt(n), galtsq, gqsq, hdiag(npt),&
     & ratio, rho, rhosq, vquad, xoptsq, wcheck(npt+n), sigma(npt) 

      logical :: model_step, reduce_rho, shortd 

      ! The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
      ! identical to the corresponding arguments in SUBROUTINE NEWUOA.
      ! XBASE will hold a shift of origin that should reduce the
      ! contributions from rounding errors to values of the model and
      ! Lagrange functions.
      ! XOPT will be set to the displacement from XBASE of the vector of
      ! variables that provides the least calculated F so far.
      ! XNEW will be set to the displacement from XBASE of the vector of
      ! variables for the current calculation of F.
      ! XPT will contain the interpolation point coordinates relative to
      ! XBASE.
      ! FVAL will hold the values of F at the interpolation points.
      ! GQ will hold the gradient of the quadratic model at XBASE.
      ! HQ will hold the explicit second order derivatives of the
      ! quadratic model.
      ! PQ will contain the parameters of the implicit second order
      ! derivatives of the quadratic model.
      ! BMAT will hold the last N columns of H. ZMAT will hold the
      ! factorization of the leading NPT by NPT submatrix of H, this
      ! factorization being ZMAT times Diag(DZ) times ZMAT^T, where the
      ! elements of DZ are plus or minus one, as specified by IDZ.
      ! NDIM is the first dimension of BMAT and has the value NPT + N.
      ! D is reserved for trial steps from XOPT.
      ! VLAG will contain the values of the Lagrange functions at a new
      ! point X. They are part of a product that requires VLAG to be of
      ! length NDIM = NPT+N.
      ! The array W will be used for working space. Its length must be
      ! at least 10*NDIM = 10*(NPT + N).

      call initialize(n, npt, rhobeg, x, xbase, xpt, f, fval, xopt,     &
     & fopt, kopt, bmat, zmat, gq, hq, pq, nf, subinfo, ftarget)
      if (subinfo == 1 .or. subinfo == -1 .or. subinfo == -2 .or.       &
     & subinfo == -3) then
          info = subinfo
          x = xbase + xopt
          f = fopt
          return
      end if

      ! Set some more initial values.
      rho = rhobeg
      delta = rho
      idz = 1
      prederr = zero
      itest = 0
      nfsave = nf
      xoptsq = zero
      do i = 1, n
          xoptsq = xoptsq + xopt(i)**2
      end do

      ! Begin the iterative procedure.

      maxtr = maxfun  ! Maximal numer of trust region iterations
      info = 0  ! Exit status

      do tr = 1, maxtr

          ! Is the trust region trial step short?
          shortd = .false.
          ! Will we improve the model after the trust region iteration?
          model_step = .false.
          ! Will we reduce rho after the trust region iteration?
          reduce_rho = .false.
          ! NEWUOA never sets MODEL_STEP = REDUCE_RHO = .TRUE.

          ! Exit if the model contains NaN. Otherwise, the behaviour of
          ! TRSAPP is unpredictable, and Segmentation Fault may occur.
          ! In the future, we may build a new model instead of exiting.
          if (any(is_nan(gq)) .or. any(is_nan(hq)) .or. any(is_nan(pq)))&
     &     then
              info = -3
              exit
          end if

          ! Solve the trust region subproblem.
          ! In Powell's NEWUOA code, VQUAD is not an output of TRSAPP. 
          ! Here we output it but do not use it to align with Powell's
          ! code. VQUAD is later calculated by CALQUAD.
          call trsapp(n, npt, 1.0e-2_rp, xopt, xpt, gq, hq, pq, delta,  &
     &     d, crvmin, vquad, subinfo)
          
          ! Calculate the length of the trial step D.
          !dsq = dot_product(d, d)
          dsq = zero
          do i = 1, n
              dsq = dsq + d(i)**2
          end do
          dnorm = min(delta, sqrt(dsq))

          ! Is the step long enough to invoke a function evaluation?
          if (dnorm < half*rho) then
              shortd = .true.
              if (0.125_rp*crvmin*rho*rho > maxval(abs(prederr)) .and.  &
     &         nf > nfsave + 2) then 
                  ! The 1st possibility (out of 2) that reduce_rho=true
                  reduce_rho = .true.  
              else ! 3 recent values of ||D_k|| and |F−Q| are small.
                  delta = tenth*delta  ! Reduce DELTA by a factor of 10
                  if (delta <= 1.5_rp*rho) then
                      delta = rho  ! Set DELTA to RHO when it is close.
                  end if
                  ! After this, DELTA < DNORM may happen, explaining why 
                  ! we sometimes write MAX(DELTA, DNORM).
              end if
          end if
    
          if (.not. shortd) then
              ! Shift XBASE if XOPT may be too far from XBASE.
              if (dsq <= 1.0e-3_rp*xoptsq) then
                  call shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq, &
     &             hq, xpt)
                  xbase = xbase + xopt
                  xopt = zero
                  xoptsq = zero
              end if
    
              ! Calculate VLAG and BETA for D. The first NPT components
              ! of W_check will be held in WCHECK.
              !call vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,&
!     &         vlag, beta, wcheck)
              call vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,&
     &         vlag, beta, wcheck, dsq, xoptsq)

      !----------------------------------------------------------------!
      
              ! Use the quadratic model to predict the change in F due
              ! to the step D.
              !call calquad(vquad, d, xopt, xpt, gq, hq, pq, n, npt)
              call calquad(vquad, d, xopt, xpt, gq, hq, pq, n, npt,     &
     &         wcheck(1:npt))

              ! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)  ! Set F to NaN. It is necessary.
                  info = -1
                  exit
              end if
              call calfun(n, x, f)
              nf = nf + 1
              ! Record the latest function evaluation with ||D|| > RHO. 
              if (dnorm > rho) then 
                  nfsave = nf 
              end if

              ! PREDERR is the prediction errors of the latest 3 models.
              prederr(2 : size(prederr)) = prederr(1 : size(prederr)-1)
              prederr(1) = f - fopt - vquad
        
              ! Update FOPT and XOPT
              ! FSAVE is needed later when setting MODEL_STEP 
              ! and REDUCE_RHO.
              fsave = fopt  
              if (f < fopt) then
                  fopt = f
                  xopt = xnew
                  xoptsq = zero
                  do i = 1, n
                      xoptsq = xoptsq + xopt(i)**2
                  end do
              end if
              ! Exit if F is NaN or INF. 
              if (is_nan(f) .or. is_posinf(f)) then 
                  info = -2
                  exit
              end if
              ! Exit if F <= FTARGET.
              if (f <= ftarget) then
                  info = 1
                  exit
              end if
              ! Exit if NF >= NFTEST
              if (nf >= maxfun) then
                  info = 3
                  exit
              end if

      !----------------------------------------------------------------!
    
              ! Update DELTA according to RATIO.
              if (is_nan(vquad) .or. vquad >= zero) then
                  info = 2
                  exit
              end if
              ratio = (f - fsave)/vquad
              if (ratio <= tenth) then
                  delta = half*dnorm
              else if (ratio <= 0.7_rp) then
                  delta = max(half*delta, dnorm)
              else
                  delta = max(half*delta, dnorm + dnorm)
              end if
              if (delta <= 1.5_rp*rho) then
                  delta = rho
              end if
              
              ! Set KNEW to the index of the interpolation point that
              ! will be deleted.
              rhosq = max(tenth*delta, rho)**2
              do k = 1, npt
                  hdiag(k) = -sum(zmat(k, 1 : idz - 1)**2) +            &
     &             sum(zmat(k, idz : npt - n - 1)**2)
                  xdiff = xpt(k, :) - xopt
                  xdsq(k) = dot_product(xdiff, xdiff) 
              end do
              sigma = abs(beta*hdiag + vlag(1 : npt)**2)
              sigma = sigma * max(xdsq/rhosq, one)**3
              if (f >= fsave) then
              ! Set SIGMA(KOPT) = -1 to prevent KNEW from being KOPT 
                  sigma(kopt) = -one  
              end if 
              if (maxval(sigma) > one .or. f < fsave) then
              ! KNEW > 0 unless MAXVAL(SIGMA) <= 1 and F >= FSAVE.
              ! If F < FSAVE, then KNEW > 0, ensuring that XNEW 
              ! will be included into XPT.
                  knew = maxloc(sigma, dim = 1)
              else
                  knew = 0
              end if

              if (knew > 0) then
                  ! Update BMAT, ZMAT and IDZ, so that the KNEW-th 
                  ! interpolation point can be removed. 
                  call update(n, npt, bmat, zmat, idz, vlag, beta,      &
     &             knew, w)
                  ! Update the quadratic model
                  call updateq(n, npt, idz, knew, prederr(1),           &
     &             xpt(knew, :), bmat(knew, :), zmat, gq, hq, pq)
    
                  ! Include the new interpolation point. This should be 
                  ! done after updating BMAT, ZMAT, and the model.
                  fval(knew) = f
                  xpt(knew, :) = xnew
    
                  ! Test whether to replace the new quadratic model Q by
                  ! the least Frobenius norm interpolant Q_alt, making
                  ! the replacement if the test is satisfied.
                  ! In the NEWUOA paper, Powell replaces Q with Q_alt
                  ! when RATIO <= 0.01 and ||G_alt|| <= 0.1||GQ||
                  ! hold for 3 consecutive times (equation (8.4)) But
                  ! the original NEWUOA compares ABS(RATIO) with 0.01
                  ! instead of RATIO. Here we use RATIO. It is observed
                  ! in Zhang Zaikun's PhD thesis (Section 3.3.2) that
                  ! it is indeed more efficient to check RATIO rather
                  ! than ABS(RATIO).
                  if (delta <= rho) then  ! Indeed, DELTA == RHO.
                      ! if (abs(ratio) > 1.0e-2_rp) then
                      if (ratio > 1.0e-2_rp) then
                          itest = 0
                      else
                          gqsq = zero
                          do i = 1, n
                              gqsq = gqsq + gq(i)**2
                          end do
                          vlag(1 : npt) = fval - fval(kopt)
                           
                          !galt = matmul(vlag(1 : npt), bmat(1 : npt, 1 : n))
                          !galtsq = dot_product(galt, galt)
                          galt = zero
                          do k = 1, npt
                              galt = galt + vlag(k)*bmat(k, :)
                          end do
                          galtsq = zero
                          do i = 1, n
                              galtsq = galtsq + galt(i)*galt(i) 
                          end do
            
                          if (gqsq < 100.0_rp*galtsq) then
                              itest = 0
                          else
                              itest = itest + 1
                          end if
                      end if
                  end if

                  if (itest >= 3) then
                      call qalt(gq, hq, pq, fval, bmat(1 : npt, :),     &
     &                 zmat, n, npt, kopt, idz)
                      itest = 0
                  end if
        
                  ! Update KOPT to KNEW if RATIO > ZERO 
                  if (f < fsave) then 
                      kopt = knew
                  end if
              end if
          end if
    
!          if (((.not. shortd) .and. (ratio < tenth .or. knew == 0))     &
!     &        .or. (shortd .and. .not. reduce_rho)) then
          if (((.not. shortd) .and. (f>fsave+tenth*vquad .or. knew==0)) &
     &     .or. (shortd .and. .not. reduce_rho)) then
              
              ! Find out if the interpolation points are close enough to
              ! the best point so far.
              distsq = 4.0_rp*delta*delta
              do k = 1, npt
                  xdiff = xpt(k, :) - xopt
                  xdsq(k) = dot_product(xdiff, xdiff)
              end do
              if (maxval(xdsq) > distsq) then
                  knew = maxloc(xdsq, dim = 1)
                  distsq = maxval(xdsq)
              else
                  knew = 0
              end if
            
              ! If KNEW is positive, then a model step will be taken to
              ! improve the geometry of the interpolation set and hence
              ! ameliorate them model.
              if (knew > 0) then
                  ! The only possibility that model_step=true
                  model_step = .true. 
              else if (max(delta, dnorm) <= rho .and. (ratio <= 0 .or.  &
     &         shortd)) then
                  ! The 2nd possibility (out of 2) that reduce_rho=true
                  reduce_rho = .true.
              end if
          end if 

          if (reduce_rho) then
              ! The calculations with the current RHO are complete.
              ! Pick the next values of RHO and DELTA.
              if (rho <= rhoend) then
                  exit
              else
                  delta = half*rho
                  ratio = rho/rhoend
                  if (ratio <= 16.0_rp) then
                      rho = rhoend
                  else if (ratio <= 250.0_rp) then
                      rho = sqrt(ratio)*rhoend
                  else
                      rho = tenth*rho
                  end if
                  delta = max(delta, rho)
                  nfsave = nf  ! Set NFSAVE to NF
              end if
          end if

          if (model_step) then
              ! Set DSTEP, which will be used as the trust region radius
              ! for the model-improving schemes BIGLAG and BIGDEN. We 
              ! also need it to decide whether to shift XBASE or not.
              dstep = max(min(tenth*sqrt(distsq), half*delta), rho)
              dsq = dstep*dstep

              ! Shift XBASE if XOPT may be too far from XBASE.
              if (dsq <= 1.0e-3_rp*xoptsq) then
                  call shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq, &
     &             hq, xpt)
                  xbase = xbase + xopt
                  xopt = zero
                  xoptsq = zero
              end if

              ! Exit if BMAT or ZMAT contains NaN. Otherwise, the 
              ! behaviour of BIGLAG and BIGDEN is unpredictable. 
              ! Segmentation Fault may occur.
              if (any(is_nan(bmat)) .or. any(is_nan(zmat))) then
                  info = -3
                  exit
              end if

              call biglag(n, npt, xopt, xpt, bmat, zmat, idz, knew,     &
     &         dstep, d, alpha)

              ! Calculate VLAG, BETA, and WCHECK for D.
!              call vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,&
!     &         vlag, beta, wcheck)
              call vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d,&
     &         vlag, beta, wcheck, dsq, xoptsq)

              ! If KNEW is positive and if the cancellation in DENOM is
              ! unacceptable, then BIGDEN calculates an alternative 
              ! model step, XNEW being used for working space.
              ! No need to check whether BMAT and ZMAT contain NaN as no
              ! change has been made to them.
              if (abs(one + alpha*beta/vlag(knew)**2) <= 0.8_rp) then
                  call bigden (n, npt, xopt, xpt, bmat, zmat, idz, kopt,&
     &             knew, d, wcheck, vlag, beta)
!                  call bigden (n, npt, xopt, xpt, bmat, zmat, idz,      &
!     &             kopt, knew, d, wcheck, vlag, beta, xnew, w(npt+n+1), &
!     &             w(6*(npt+n)+1))
              end if
 
      !----------------------------------------------------------------!

              ! Use the current quadratic model to predict the change in
              ! F due to the step D. 
              !call calquad(vquad, d, xopt, xpt, gq, hq, pq, n, npt)
              call calquad(vquad, d, xopt, xpt, gq, hq, pq, n, npt,     &
     &         wcheck(1:npt))

              ! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)  ! Set F to NaN. It is necessary.
                  info = -1
                  exit
              end if
              call calfun(n, x, f)
              nf = nf + 1

              ! The following seems different from what is introduced in 
              ! Section 7 (around (7.7)) of the NEUOA paper. Seemingly
              ! we should keep dnorm=||d||.
              if (dnorm > rho) then 
                  nfsave = nf  !? dnorm is from last TR? 
              end if
             
              ! PREDERR is the prediction errors of the latest 3 models.
              prederr(2 : size(prederr)) = prederr(1 : size(prederr)-1)
              prederr(1) = f - fopt - vquad
        
              ! Update FOPT and XOPT if the new F is the least value of 
              ! the objective function so far. 
              fsave = fopt
              if (f < fopt) then
                  fopt = f
                  xopt = xnew
                  xoptsq = zero
                  do i = 1, n
                      xoptsq = xoptsq + xopt(i)**2
                  end do
              end if

              ! Check whether to exit.
              if (is_nan(f) .or. is_posinf(f)) then 
                  info = -2
                  exit
              end if
              if (f <= ftarget) then
                  info = 1
                  exit
              end if
              if (nf >= maxfun) then
                  info = 3
                  exit
              end if

      !----------------------------------------------------------------!
    
              ! Update BMAT, ZMAT and IDZ, so that the KNEW-th
              ! interpolation point can be moved. 
              call update(n, npt, bmat, zmat, idz,vlag,beta,knew,w)
              ! Update the quadratic model.
              call updateq(n, npt, idz, knew, prederr(1), xpt(knew, :), &
     &         bmat(knew, :), zmat, gq, hq, pq)
    
              ! Include the new interpolation point. This should be done
              ! after updating BMAT, ZMAT, and the model.
              fval(knew) = f
              xpt(knew, :) = xnew
              if (f < fsave) then
                  kopt = knew
              end if
          end if
      end do

      ! Return from the calculation, after another Newton-Raphson step,
      ! if it is too short to have been tried before.
      if (shortd .and. nf < maxfun) then
          x = xbase + (xopt + d)
          if (any(is_nan(x))) then
              f = sum(x)  ! Set F to NaN. It is necessary.
              info = -1
          else
              call calfun(n, x, f)
              nf = nf + 1
          end if
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! By Zaikun (commented on 02-06-2019; implemented in 2016):
      ! Note that (FOPT .LE. F) is FALSE if F is NaN; When F is NaN, it
      ! is also necessary to update X and F.
      if (is_nan(f) .or. fopt <= f) then
          x = xbase + xopt
          f = fopt
      end if

      return

      end subroutine newuob
