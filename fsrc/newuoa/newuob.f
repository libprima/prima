      module newuob_mod

      implicit none
      private
      public :: newuob


      contains

      subroutine newuob (iprint, maxfun, npt, ftarget, rhobeg, rhoend,  &
     & x, nf, f, info)

      use consts_mod, only : RP, IK, ZERO, ONE, HALF, TENTH
      use infnan_mod, only : is_nan, is_posinf
      use lina_mod
      use initialize_mod, only : initxf, initq, inith
      use trsubp_mod, only : trsapp
      use geometry_mod, only : biglag, bigden
      use shiftbase_mod, only : shiftbase
      use vlagbeta_mod, only : vlagbeta
      use update_mod, only : updateh, updateq, qalt
      implicit none
            
      ! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: maxfun
      integer(IK), intent(in) :: npt
      real(RP), intent(in) :: ftarget
      real(RP), intent(in) :: rhobeg
      real(RP), intent(in) :: rhoend

      ! In-output
      real(RP), intent(inout) :: x(:)  ! SIZE(X) = N

      ! Outputs 
      integer(IK), intent(out) :: info
      integer(IK), intent(out) :: nf
      real(RP), intent(out) :: f 

      ! Intermediate variables
      integer(IK) :: n, ij(2, npt), idz, itest, knew, kopt, nfsave
      integer(IK) :: tr, maxtr, subinfo
      real(RP) :: alpha, beta, crvmin, delta, prederr(3) 
      real(RP) :: distsq, dnorm, dsq, dstep, xdsq(npt), d(size(x))
      real(RP) :: xbase(size(x)), xopt(size(x)), xnew(size(x))
      real(RP) :: xpt(size(x), npt), fval(npt)
      real(RP) :: gq(size(x)), hq(size(x), size(x)), pq(npt)
      real(RP) :: bmat(size(x), npt+size(x)), zmat(npt, npt-size(x)-1)
      real(RP) :: vlag(npt + size(x)), zknew(npt-size(x)-1)
      real(RP) :: fopt, fsave, galt(size(x)), galtsq, gqsq, hdiag(npt)
      real(RP) :: ratio, rho, rhosq, vquad, xoptsq, sigma(npt), trtol
      logical :: improve_geometry, reduce_rho, shortd

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
      ! XBASE, each COLUMN corresponding to a point.
      ! FVAL will hold the values of F at the interpolation points.
      ! GQ will hold the gradient of the quadratic model at XBASE.
      ! HQ will hold the explicit second order derivatives of the
      ! quadratic model.
      ! PQ will contain the parameters of the implicit second order
      ! derivatives of the quadratic model.
      ! BMAT will hold the last N ROWs of H. ZMAT will hold the
      ! factorization of the leading NPT by NPT sub-matrix of H, this
      ! factorization being ZMAT times Diag(DZ) times ZMAT^T, where the
      ! elements of DZ are plus or minus one, as specified by IDZ.
      ! NDIM is the second dimension of BMAT and has the value NPT + N.
      ! D is reserved for trial steps from XOPT.
      ! VLAG will contain the values of the Lagrange functions at a new
      ! point X. They are part of a product that requires VLAG to be of
      ! length NDIM = NPT+N.

      ! Get size.
      n = int(size(x), kind(n))
      
      maxtr = maxfun  ! Maximal numer of trust region iterations.
      info = 0  ! Exit status. The default value is 0.

      ! Initialize FVAL, XBASE, and XPT
      call initxf(x, rhobeg, ftarget, ij, kopt, nf, fval, xbase, xpt,   &
     & subinfo)
      xopt = xpt(:, kopt)
      fopt = fval(kopt)
      x = xbase + xopt  ! Set X.
      f = fopt  ! Set F.
      if (subinfo == 1 .or. subinfo == -1 .or. subinfo == -2) then
          info = subinfo  ! Set the exit status.
          return
      end if

      ! Initialize GQ, HQ, and PQ
      call initq(ij, fval, xpt, gq, hq, pq, subinfo)

      ! Initialize BMAT and ZMAT
      call inith(ij, xpt, bmat, zmat, subinfo)

      ! Set some more initial values.
      rho = rhobeg
      delta = rho
      idz = 1
      prederr = ZERO
      itest = 0
      nfsave = nf
      xoptsq = dot_product(xopt, xopt)

      ! Begin the iterative procedure.
      do tr = 1, maxtr

          ! Is the trust region trial step short?
          shortd = .false.
          ! Will we improve the model after the trust region iteration?
          improve_geometry = .false.
          ! Will we reduce rho after the trust region iteration?
          reduce_rho = .false.
          ! NEWUOA never sets both MODEL_STEP and REDUCE_RHO to .TRUE. 
          ! at the same time.

          ! Solve the trust region subproblem.
          ! In Powell's NEWUOA code, VQUAD is not an output of TRSAPP.
          ! Here we output it. However, we will not use it; it will sill
          ! be calculated latter by CALQUAD in order to produce the same
          ! results Powell's code. 
          trtol = 1.0e-2_RP  ! Tolerance used in trsapp.
          call trsapp(delta, gq, hq, pq, trtol, xopt, xpt, crvmin,      &
     &     vquad, d, subinfo)

          ! Calculate the length of the trial step D.
          dsq = dot_product(d, d)
          dnorm = min(delta, sqrt(dsq))

          ! Is the step long enough to invoke a function evaluation?
          if (dnorm < HALF*rho) then
              shortd = .true.
              if (0.125_RP*crvmin*rho*rho > maxval(abs(prederr)) .and.  &
     &         nf > nfsave + 2) then
                  ! The 1st possibility (out of 2) that reduce_rho=true
                  reduce_rho = .true.
              else ! 3 recent values of ||D_k|| and |F−Q| are small.
                  delta = TENTH*delta  ! Reduce DELTA by a factor of 10
                  if (delta <= 1.5_RP*rho) then
                      delta = rho  ! Set DELTA to RHO when it is close.
                  end if
                  ! After this, DELTA < DNORM may happen, explaining why
                  ! we sometimes write MAX(DELTA, DNORM).
              end if
          end if

          if (.not. shortd) then
              ! Shift XBASE if XOPT may be too far from XBASE.
              if (dsq <= 1.0e-3_RP*xoptsq) then
                  call shiftbase(idz, pq, xopt, zmat, bmat, gq, hq, xpt)
                  xbase = xbase + xopt
                  xopt = ZERO
                  xoptsq = ZERO
              end if

              ! Calculate VLAG and BETA for D. 
              call vlagbeta(idz, kopt, bmat, d, xopt,xpt,zmat,beta,vlag)

      !----------------------------------------------------------------!

              ! Use the current quadratic model to predict the change in
              ! F due to the step D.
              call calquad(vquad, d, xopt, xpt, gq, hq, pq)

              ! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)  ! Set F to NaN. It is necessary.
                  info = -1
                  exit
              end if
              call calfun(n, x, f)
              nf = int(nf + 1, kind(nf))
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
                  xoptsq = dot_product(xopt, xopt)
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
              if (is_nan(vquad) .or. vquad >= ZERO) then
                  info = 2
                  exit
              end if
              ratio = (f - fsave)/vquad
              if (ratio <= TENTH) then
                  delta = HALF*dnorm
              else if (ratio <= 0.7_RP) then
                  delta = max(HALF*delta, dnorm)
              else
                  delta = max(HALF*delta, dnorm + dnorm)
              end if
              if (delta <= 1.5_RP*rho) then
                  delta = rho
              end if

              ! Set KNEW to the index of the interpolation point that
              ! will be deleted.
              rhosq = max(TENTH*delta, rho)**2
              hdiag = -sum(zmat(:, 1 : idz - 1)**2, dim = 2) +          &
     &         sum(zmat(:, idz : npt - n - 1)**2, dim = 2)             
              xdsq = sum((xpt-spread(xopt,dim=2,ncopies=npt))**2, dim=1)
              sigma = abs(beta*hdiag + vlag(1 : npt)**2)
              sigma = sigma * max(xdsq/rhosq, one)**3
              if (f >= fsave) then
              ! Set SIGMA(KOPT) = -1 to prevent KNEW from being KOPT
                  sigma(kopt) = -one
              end if
              if (maxval(sigma) > ONE .or. f < fsave) then
              ! KNEW > 0 unless MAXVAL(SIGMA) <= 1 and F >= FSAVE.
              ! If F < FSAVE, then KNEW > 0, ensuring that XNEW
              ! will be included into XPT.
                  knew = int(maxloc(sigma, dim = 1), kind(knew))
              else
                  knew = 0
              end if

              if (knew > 0) then
                  ! Update BMAT, ZMAT and IDZ, so that the KNEW-th
                  ! interpolation point can be removed.
                  call updateh(knew, beta, vlag, idz, bmat, zmat)

                  ! Update the quadratic model
                  call updateq(idz, knew, bmat(:, knew), prederr(1),    &
     &             zmat, xpt(:, knew), gq, hq, pq)

                  ! Include the new interpolation point. This should be
                  ! done after updating BMAT, ZMAT, and the model.
                  fval(knew) = f
                  xpt(:, knew) = xnew

                  ! Test whether to replace the new quadratic model Q by
                  ! the least Frobenius norm interpolant Q_alt.
                  ! In the NEWUOA paper, Powell replaces Q with Q_alt
                  ! when RATIO <= 0.01 and ||G_alt|| <= 0.1||GQ||
                  ! hold for 3 consecutive times (equation (8.4)). But
                  ! the original NEWUOA compares ABS(RATIO) with 0.01
                  ! instead of RATIO. Here we use RATIO. It is observed
                  ! in Zhang Zaikun's PhD thesis (Section 3.3.2) that
                  ! it is indeed more efficient to check RATIO rather
                  ! than ABS(RATIO).
                  if (delta <= rho) then  ! Equivalent to  DELTA == RHO.
                      ! if (abs(ratio) > 1.0e-2_RP) then
                      if (ratio > 1.0e-2_RP) then
                          itest = 0
                      else
                          gqsq = dot_product(gq, gq)
                          galt = matmul(bmat(:, 1:npt), fval-fval(kopt))
                          galtsq = dot_product(galt, galt)

                          if (gqsq < 100.0_RP*galtsq) then
                              itest = 0
                          else
                              itest = int(itest + 1, kind(itest))
                          end if
                      end if
                  end if

                  ! Replace Q with Q_alt when ITEST >= 3.
                  if (itest >= 3) then
                      call qalt(idz, kopt, fval, bmat(:, 1:npt), zmat,  &
     &                 gq, hq, pq)
                      itest = 0
                  end if

                  ! Update KOPT to KNEW if RATIO > ZERO
                  if (f < fsave) then
                      kopt = knew
                  end if
              end if
          end if


          if (((.not. shortd) .and. (ratio < TENTH .or. knew == 0))     &
!          if (((.not. shortd) .and. (f>fsave+TENTH*vquad .or. knew==0)) &
     &     .or. (shortd .and. .not. reduce_rho)) then

              ! Find out if the interpolation points are close enough to
              ! the best point so far.
              distsq = 4.0_RP*delta*delta
              xdsq = sum((xpt-spread(xopt,dim=2,ncopies=npt))**2, dim=1)
              if (maxval(xdsq) > distsq) then
                  knew = int(maxloc(xdsq, dim = 1), kind(knew))
                  distsq = maxval(xdsq)
              else
                  knew = 0
              end if

              ! If KNEW is positive, then a model step will be taken to
              ! improve the geometry of the interpolation set and hence
              ! ameliorate them model.
              if (knew > 0) then
                  ! The only possibility that improve_geometry=true
                  improve_geometry = .true.
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
                  delta = HALF*rho
                  ratio = rho/rhoend
                  if (ratio <= 16.0_RP) then
                      rho = rhoend
                  else if (ratio <= 250.0_RP) then
                      rho = sqrt(ratio)*rhoend
                  else
                      rho = TENTH*rho
                  end if
                  delta = max(delta, rho)
                  nfsave = nf  ! Set NFSAVE to NF
              end if
          end if

          if (improve_geometry) then
              ! Set DSTEP, which will be used as the trust region radius
              ! for the model-improving schemes BIGLAG and BIGDEN. We
              ! also need it to decide whether to shift XBASE or not.
              dstep = max(min(TENTH*sqrt(distsq), HALF*delta), rho)
              dsq = dstep*dstep

              ! Shift XBASE if XOPT may be too far from XBASE.
              if (dsq <= 1.0e-3_RP*xoptsq) then
                  call shiftbase(idz, pq, xopt, zmat, bmat, gq, hq, xpt)
                  xbase = xbase + xopt
                  xopt = ZERO
                  xoptsq = ZERO
              end if

              call biglag(idz, knew, dstep, bmat, xopt, xpt, zmat, d)

              ! ALPHA is the KNEW-th diagonal entry of H
              zknew = zmat(knew, :)
              zknew(1 : idz - 1) = -zknew(1 : idz - 1)
              alpha = dot_product(zmat(knew, :), zknew)

              ! Calculate VLAG and BETA for D.
              call vlagbeta(idz, kopt, bmat, d, xopt,xpt,zmat,beta,vlag)

              ! If KNEW is positive and if the cancellation in DENOM is
              ! unacceptable, then BIGDEN calculates an alternative
              ! model step D.
              ! Why don't we call vlagbeta for this D? Because vlag and 
              ! beta are calculated within BIGDEN. 
              ! No need to check whether BMAT and ZMAT contain NaN as no
              ! change has been made to them.
              if (abs(ONE + alpha*beta/vlag(knew)**2) <= 0.8_RP) then
                  call bigden(idz, knew, kopt, bmat, xopt, xpt, zmat, d,&
     &             beta, vlag)
              end if

      !----------------------------------------------------------------!

              ! Use the current quadratic model to predict the change in
              ! F due to the step D.
              call calquad(vquad, d, xopt, xpt, gq, hq, pq)

              ! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)  ! Set F to NaN. It is necessary.
                  info = -1
                  exit
              end if
              call calfun(n, x, f)
              nf = int(nf + 1, kind(nf))

              !--------------------------------------------------------!
              !!!??? Powell's code does not update DNORM. Therefore,
              ! DNORM is the length of last trust-region trial step.
              ! which seems inconsistent with what is introduced in
              ! Section 7 (around (7.7)) of the NEUOA paper. Seemingly
              ! we should keep dnorm=||d|| as we do here. 
              dnorm = min(dstep, sqrt(dot_product(d, d)))
              ! In theory, dnorm=dstep
              !--------------------------------------------------------!

              if (dnorm > rho) then
                  nfsave = nf  !? DNORM is from last TR?
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
                  xoptsq = dot_product(xopt, xopt)
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
              call updateh(knew, beta, vlag, idz, bmat, zmat)

              ! Update the quadratic model.
              call updateq(idz, knew, bmat(:, knew), prederr(1), zmat,  &
     &         xpt(:, knew), gq, hq, pq)

              ! Include the new interpolation point. This should be done
              ! after updating BMAT, ZMAT, and the model.
              fval(knew) = f
              xpt(:, knew) = xnew
              if (f < fsave) then
                  kopt = knew
              end if
          end if
      end do

      ! Return from the calculation, after another Newton-Raphson step,
      ! if it is too short to have been tried before. 
      ! Note that no trust region iteration has been done if maxtr = 0,
      ! and hence we should not check whether shortd = true but return
      ! immediately. 
      if (maxtr > 0 .and. shortd .and. nf < maxfun) then
          x = xbase + (xopt + d)
          if (any(is_nan(x))) then
              f = sum(x)  ! Set F to NaN. It is necessary.
              info = -1
          else
              call calfun(n, x, f)
              nf = int(nf + 1, kind(nf))
          end if
      end if

      ! Note that (FOPT .LE. F) is FALSE if F is NaN; When F is NaN, it
      ! is also necessary to update X and F.
      if (is_nan(f) .or. fopt <= f) then
          x = xbase + xopt
          f = fopt
      end if

      return

      end subroutine newuob

      end module newuob_mod
