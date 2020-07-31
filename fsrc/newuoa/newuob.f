      module newuob_mod

      implicit none
      private
      public :: newuob


      contains

      subroutine newuob(iprint, maxfun, npt, eta1, eta2, ftarget,       &
     & gamma1, gamma2, rhobeg, rhoend, x, nf, f, info)

      use consts_mod, only : RP, IK, ZERO, HALF, TENTH, HUGENUM
      use infos_mod
      use infnan_mod, only : is_nan, is_posinf
      use lina_mod
      use initialize_mod, only : initxf, initq, inith
      use trustregion_mod, only : trsapp, trrad
      use geometry_mod, only : setremove, ameliorgeo
      use shiftbase_mod, only : shiftbase
      use vlagbeta_mod, only : vlagbeta
      use update_mod, only : updateh, updateq, tryqalt
      use output_mod, only : retmssg, rhomssg, fmssg
      implicit none

      ! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: maxfun
      integer(IK), intent(in) :: npt
      real(RP), intent(in) :: eta1
      real(RP), intent(in) :: eta2
      real(RP), intent(in) :: ftarget
      real(RP), intent(in) :: gamma1
      real(RP), intent(in) :: gamma2
      real(RP), intent(in) :: rhobeg
      real(RP), intent(in) :: rhoend

      ! In-output
      real(RP), intent(inout) :: x(:)  ! SIZE(X) = N

      ! Outputs
      integer(IK), intent(out) :: info
      integer(IK), intent(out) :: nf
      real(RP), intent(out) :: f

      ! Intermediate variables
      integer(IK) :: n, ij(2, npt), idz, itest, knew, kopt
      integer(IK) :: tr, maxtr, subinfo
      real(RP) :: beta, crvmin, delta
      real(RP) :: distsq, dnorm, delbar, xdsq(npt), d(size(x))
      real(RP) :: xbase(size(x)), xopt(size(x)), xnew(size(x))
      real(RP) :: xpt(size(x), npt), fval(npt)
      real(RP) :: gq(size(x)), hq(size(x), size(x)), pq(npt)
      real(RP) :: bmat(size(x), npt+size(x)), zmat(npt, npt-size(x)-1)
      real(RP) :: vlag(npt + size(x))
      real(RP) :: fopt, fsave, ratio, rho
      real(RP) :: vquad, trtol, fqdiff
      real(RP) :: moderr(3), dnormsave(size(moderr))
      logical :: terminate, improve_geometry, reduce_rho, shortd
      character(len = 6), parameter :: solver= 'NEWUOA'

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
      terminate = .false. ! Whether to terminate (after initialization).

      ! Initialize FVAL, XBASE, and XPT
      call initxf(iprint, x, rhobeg, ftarget, ij, kopt, nf, fval, xbase,&
     & xpt, subinfo)
      xopt = xpt(:, kopt)
      fopt = fval(kopt)
      x = xbase + xopt  ! Set X.
      f = fopt  ! Set F.

      ! Check whether to return
      if (subinfo == FTARGET_ACHIEVED) then
          terminate = .true.
      else if (subinfo == NAN_X) then
          terminate = .true.
      else if (subinfo == NAN_INF_F) then
          terminate = .true.
      end if

      if (terminate) then
          info = subinfo
          if (iprint >= 1) then
              call retmssg(info, iprint, nf, f, x, solver)
          end if
          return
      end if

      ! Initialize GQ, HQ, and PQ
      call initq(ij, fval, xpt, gq, hq, pq, subinfo)

      ! Initialize BMAT and ZMAT
      call inith(ij, xpt, bmat, zmat, subinfo)

      ! After initializing GQ, HQ, PQ, BMAT, ZMAT, one can also choose
      ! to return if subinfo = NAN_MODEL (NaN occurs in the model). We
      ! do not do it here. If such a modle is harmful, then it will
      ! probably lead to other returns (NaN in X, NaN in F, trust region
      ! subproblem fails, ...); otherwise, the code will continue to run
      ! and possibly get rid of the NaN in the model.

      ! Set some more initial values.
      rho = rhobeg
      delta = rho
      idz = 1
      moderr = HUGENUM
      dnormsave = HUGENUM
      itest = 0

      ! Begin the iterative procedure.
      do tr = 1, maxtr

          ! Is the trust region trial step short?
          shortd = .false.
          ! Will we improve the model after the trust region iteration?
          improve_geometry = .false.
          ! Will we reduce rho after the trust region iteration?
          reduce_rho = .false.
          ! NEWUOA never sets both IMPROVE_GEOMETRY and REDUCE_RHO to
          ! .TRUE. at the same time.

          ! Solve the trust region subproblem.
          ! In Powell's NEWUOA code, VQUAD is not an output of TRSAPP.
          ! Here we output it. However, we will not use it; it will sill
          ! be calculated latter by CALQUAD in order to produce the same
          ! results Powell's code.
          trtol = 1.0e-2_RP  ! Tolerance used in trsapp.
          call trsapp(delta, gq, hq, pq, trtol, xopt, xpt, crvmin,      &
     &     vquad, d, subinfo)

          ! Calculate the length of the trial step D.
          dnorm = min(delta, sqrt(dot_product(d, d)))

          ! Is the step long enough to invoke a function evaluation?
          if (dnorm < HALF*rho) then
              shortd = .true.
              if (0.125_RP*crvmin*rho*rho >= maxval(abs(moderr)) .and.  &
     &         rho >= maxval(dnormsave)) then
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
              ! Save the current FOPT in FSAVE. It is needed later.
              fsave = fopt

              ! Shift XBASE if XOPT may be too far from XBASE.
              !if (dot_product(d, d) <= 1.0e-3_RP*xoptsq) then  ! Powell
              if (dnorm*dnorm <= 1.0e-3_RP*dot_product(xopt, xopt)) then
                  call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase,    &
     &             xopt, xpt)
              end if

              ! Calculate VLAG and BETA for D.
              call vlagbeta(idz, kopt, bmat, d, xopt,xpt,zmat,beta,vlag)

              ! Use the current quadratic model to predict the change in
              ! F due to the step D.
              call calquad(d, gq, hq, pq, xopt, xpt, vquad)

              ! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)  ! Set F to NaN. It is necessary.
                  info = NAN_X
                  exit
              end if
              call calfun(n, x, f)
              nf = int(nf + 1, kind(nf))
              if (iprint >= 3) then
                  call fmssg(iprint, nf, f, x, solver)
              end if

              ! FQDIFF is the error of the current model in predicting
              ! the change in F due to D.
              fqdiff = f - fsave - vquad

              ! Update FOPT and XOPT
              if (f < fopt) then
                  fopt = f
                  xopt = xnew
              end if

              ! Check whether to exit
              if (is_nan(f) .or. is_posinf(f)) then
                  info = NAN_INF_F
                  exit
              end if
              if (f <= ftarget) then
                  info = FTARGET_ACHIEVED
                  exit
              end if
              if (nf >= maxfun) then
                  info = MAXFUN_REACHED
                  exit
              end if

              ! Update DELTA according to RATIO.
              if (is_nan(vquad) .or. vquad >= ZERO) then
                  info = TRSUBP_FAILED
                  exit
              end if
              ratio = (f - fsave)/vquad
              delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2,   &
     &         ratio)
              if (delta <= 1.5_RP*rho) then
                  delta = rho
              end if

              ! Set KNEW to the index of the interpolation point that
              ! will be deleted. KNEW will ensure that the geometry of
              ! XPT is "optimal" after XPT(:, KNEW) is replaced by XNEW.
              ! Note that the information of XNEW is included in VLAG
              ! and BETA, which are calculated according to D=XNEW-XOPT.
              call setremove(idz, kopt, beta, delta, ratio, rho,        &
     & vlag(1:npt), xopt, xpt, zmat, knew)

              if (knew > 0) then
                  ! Update BMAT, ZMAT and IDZ, so that the KNEW-th
                  ! interpolation point can be removed.
                  call updateh(knew, beta, vlag, idz, bmat, zmat)

                  ! Update the quadratic model
                  call updateq(idz, knew, bmat(:, knew), fqdiff, zmat,  &
     &             xpt(:, knew), gq, hq, pq)

                  ! Include the new interpolation point. This should be
                  ! done after updating BMAT, ZMAT, and the model.
                  fval(knew) = f
                  xpt(:, knew) = xnew

                  ! Update KOPT to KNEW if F < FSAVE (i.e., last FOPT)
                  if (f < fsave) then
                      kopt = knew
                  end if

                  if (delta <= rho) then  ! Equivalent to  DELTA == RHO.
                  ! Test whether to replace the new quadratic model Q by
                  ! the least Frobenius norm interpolant Q_alt. Perform
                  ! the replacement if certain ceriteria is satisfied.
                  ! This part is optional, but it is important for
                  ! the performance on a certain class of problems See
                  ! Section 8 of the NEWUOA paper.
                  ! In theory, the FVAL-FSAVE in the following line can
                  ! be replaced by FVAL + C with any constant C. This
                  ! constant will not affect the result in precise
                  ! arithmetic. Powell chose C = -FSAVE.
                      call tryqalt(idz, fval - fsave, ratio,            &
     &                 bmat(:, 1 : npt), zmat, itest, gq, hq, pq)
                  ! Since tryqalt is invoked only when DELTA equals the
                  ! current RHO, why not reset ITEST to 0 when RHO is
                  ! reduced?
                  end if
              end if

              ! DNORMSAVE constains the DNORM corresponding to the
              ! latest 3 function evaluations with the current RHO.
              dnormsave = (/ dnorm, dnormsave(1 : size(dnormsave)-1) /)
              ! MODERR is the prediction errors of the latest 3 models.
              moderr = (/ fqdiff, moderr(1 : size(moderr)-1) /)

          end if


          ! The geometry of XPT probably need improvement if
          ! 1. The trust region step D is not short but RATIO < TENTH or
          !    there is no good member in XPT to be replaced by XOPT+D
          !    (i.e., KNEW=0), or
          ! 2. D is short but the latest three model errors are not
          !    small enough to render REDUCE_RHO = TRUE.
          if (((.not. shortd) .and. (ratio < TENTH .or. knew == 0))     &
     &     .or. (shortd .and. .not. reduce_rho)) then
              ! Find out if the interpolation points are close enough to
              ! the best point so far, i.e., all the points are within
              ! a ball centered at XOPT with a radius of 2*DELTA.
              ! If not, set KNEW to the index of the point that is the
              ! farthest away.
              distsq = 4.0_RP*delta*delta
              xdsq = sum((xpt-spread(xopt,dim=2,ncopies=npt))**2, dim=1)
              if (maxval(xdsq) > distsq) then
                  knew = int(maxloc(xdsq, dim = 1), kind(knew))
                  distsq = maxval(xdsq)
              else
                  knew = 0
              end if

              ! If KNEW is positive (i.e., not all points are close
              ! enough to XOPT), then a model step will be taken to
              ! ameliorgeo the geometry of the interpolation set and
              ! hence improve the model. Otherwise, RHO will be reduced
              ! (if MAX(DELTA, DNORM) <= RHO) or another trust-region
              ! step will be taken.
              if (knew > 0) then
                  ! The only possibility that IMPROVE_GEOMETRY=TRUE
                  improve_geometry = .true.
              else if (max(delta, dnorm) <= rho .and. (ratio <= 0 .or.  &
     &         shortd)) then
                  ! The 2nd possibility (out of 2) that REDUCE_RHO=TRUE
                  reduce_rho = .true.
              end if
          end if

          ! Before the next trust region step, we will either reduce rho
          ! or improve the geometry of XPT according to the values of
          ! REDUCE_RHO and IMPROVE_GEOMETRY.

          if (reduce_rho) then
              ! The calculations with the current RHO are complete.
              ! Pick the next values of RHO and DELTA.
              if (rho <= rhoend) then
                  info = SMALL_TR_RADIUS  
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
                  ! DNORMSAVE constains the DNORM corresponding to the
                  ! latest 3 function evaluations with the current RHO.
                  dnormsave = HUGENUM
                  moderr = HUGENUM
                  if (iprint >= 2) then
                      call rhomssg(iprint, nf, fopt, rho, xbase+xopt,   &
     &                 solver)
                 end if
              end if
          end if

          if (improve_geometry) then

              ! Save the current FOPT in fsave. It is needed later.
              fsave = fopt

              ! Set DELBAR, which will be used as the trust region
              ! radius for the geometry-improving schemes AMELIORGEO.
              ! We also need it to decide whether to shift XBASE or not.
              delbar = max(min(TENTH*sqrt(distsq), HALF*delta), rho)

              ! Shift XBASE if XOPT may be too far from XBASE.
              if (delbar*delbar<=1.0e-3_RP*dot_product(xopt,xopt)) then
                  call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase,    &
     &             xopt, xpt)
              end if

              ! Find a step D so that the geometry of XPT will be
              ! improved when XPT(:, KNEW) is replaced by XOPT + D.
              ! The AMELIORGEO subroutine will call Powell's BIGLAG and
              ! BIGDEN. It will also calculate the VLAG and BETA for
              ! this D.
              call ameliorgeo(idz, knew, kopt, bmat, delbar, xopt, xpt, &
     & zmat, d, beta, vlag)

              ! Use the current quadratic model to predict the change in
              ! F due to the step D.
              call calquad(d, gq, hq, pq, xopt, xpt, vquad)

              ! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)  ! Set F to NaN. It is necessary.
                  info = NAN_X
                  exit
              end if
              call calfun(n, x, f)
              nf = int(nf + 1, kind(nf))
              if (iprint >= 3) then
                  call fmssg(iprint, nf, f, x, solver)
              end if

              ! FQDIFF is the error of the current model in predicting
              ! the change in F due to D.
              fqdiff = f - fsave - vquad

              ! Update FOPT and XOPT
              if (f < fopt) then
                  fopt = f
                  xopt = xnew
              end if

              ! Check whether to exit.
              if (is_nan(f) .or. is_posinf(f)) then
                  info = NAN_INF_F
                  exit
              end if
              if (f <= ftarget) then
                  info = FTARGET_ACHIEVED
                  exit
              end if
              if (nf >= maxfun) then
                  info = MAXFUN_REACHED
                  exit
              end if

              ! Update BMAT, ZMAT and IDZ, so that the KNEW-th
              ! interpolation point can be moved.
              call updateh(knew, beta, vlag, idz, bmat, zmat)

              ! Update the quadratic model.
              call updateq(idz, knew, bmat(:, knew), fqdiff, zmat,      &
     &         xpt(:, knew), gq, hq, pq)

              ! Include the new interpolation point. This should be done
              ! after updating BMAT, ZMAT, and the model.
              fval(knew) = f
              xpt(:, knew) = xnew
              if (f < fsave) then
                  kopt = knew
              end if

              ! DNORMSAVE constains the DNORM corresponding to the
              ! latest 3 function evaluations with the current RHO.
              !--------------------------------------------------------!
              ! Powell's code does not update DNORM. Therefore,
              ! DNORM is the length of last trust-region trial step,
              ! which seems inconsistent with what is described in
              ! Section 7 (around (7.7)) of the NEWUOA paper. Seemingly
              ! we should keep DNORM = ||D|| as we do here. The value of
              ! DNORM will be used when defining NFSAVE.
              dnorm = min(delbar, sqrt(dot_product(d, d)))
              ! In theory, DNORM = DELBAR in this case.
              !--------------------------------------------------------!
              dnormsave = (/ dnorm, dnormsave(1 : size(dnormsave)-1) /)
              ! MODERR is the prediction errors of the latest 3 models.
              moderr = (/ fqdiff, moderr(1 : size(moderr)-1) /)
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
              info = NAN_X
          else
              call calfun(n, x, f)
              nf = int(nf + 1, kind(nf))
              if (iprint >= 3) then
                  call fmssg(iprint, nf, f, x, solver)
              end if
          end if
      end if

      ! Note that (FOPT .LE. F) is FALSE if F is NaN; When F is NaN, it
      ! is also necessary to update X and F.
      if (is_nan(f) .or. fopt <= f) then
          x = xbase + xopt
          f = fopt
      end if

      if (iprint >= 1) then
          call retmssg(info, iprint, nf, f, x, solver)
      end if

      return
      end subroutine newuob

      end module newuob_mod
