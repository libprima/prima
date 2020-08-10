!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of newuob.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 11-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! NEWUOB_MOD is a module that performs the major calculations of NEWUOA.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.


      module newuob_mod

      implicit none
      private
      public :: newuob


      contains


      subroutine newuob(calfun, iprint, maxfun, npt, eta1, eta2, ftarget&
     &, gamma1, gamma2, rhobeg, rhoend, x, nf, f, fhist, xhist, info)
! NEWUOB performs the actual calculations of NEWUOA. The arguments IPRINT,
! MAXFUN, MAXHIST, NPT, ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG,
! RHOEND, X, NF, F, FHIST, XHIST, and INFO are identical to the corresponding
! arguments in subroutine NEWUOA.

! XBASE will hold a shift of origin that should reduce the contributions
! from rounding errors to values of the model and Lagrange functions.
! XOPT will be set to the displacement from XBASE of the vector of variables
! that provides the least calculated F so far.
! XNEW will be set to the displacement from XBASE of the vector of variables
! for the current calculation of F.
! XPT will contain the interpolation point coordinates relative to XBASE,
! each COLUMN corresponding to a point.
! FVAL will hold the values of F at the interpolation points.
! GQ will hold the gradient of the quadratic model at XBASE.
! HQ will hold the explicit second order derivatives of the quadratic model.
! PQ will contain the parameters of the implicit second order derivatives
! of the quadratic model.
! BMAT will hold the last N ROWs of H. ZMAT will hold the factorization
! of the leading NPT by NPT sub-matrix of H, this factorization being
! ZMAT times Diag(DZ) times ZMAT^T, where DZ(1 : IDZ-1) = -1 and
! DZ(IDZ - 1 : NPT) = 1.
! D is reserved for trial steps from XOPT.
! VLAG will contain the values of the Lagrange functions at a new point
! X. They are part of a product that requires VLAG to be of length NPT + N.
!
! See Section 2 of the NEWUOA paper.

! Generic modules
      use consts_mod, only : RP, IK, ZERO, HALF, TENTH, HUGENUM, DEBUGGI&
     &NG, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod

! Solver-specific modules
      use prob_mod, only : funeval
      use initialize_mod, only : initxf, initq, inith
      use trustregion_mod, only : trsapp, trrad
      use geometry_mod, only : setremove, ameliorgeo
      use shiftbase_mod, only : shiftbase
      use vlagbeta_mod, only : vlagbeta
      use update_mod, only : updateh, updateq, tryqalt

      implicit none

! Inputs
      procedure(funeval) :: calfun
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
      real(RP), intent(inout) :: x(:) ! SIZE(X) = N

! Outputs
      integer(IK), intent(out) :: info
      integer(IK), intent(out) :: nf
      real(RP), intent(out) :: f
      real(RP), intent(out) :: fhist(:)
      real(RP), intent(out) :: xhist(:, :)

! Intermediate variables
      integer(IK) :: idz
      integer(IK) :: ij(2, npt)
      integer(IK) :: itest
      integer(IK) :: khist
      integer(IK) :: knew
      integer(IK) :: kopt
      integer(IK) :: maxfhist
      integer(IK) :: maxtr
      integer(IK) :: maxxhist
      integer(IK) :: n
      integer(IK) :: subinfo
      integer(IK) :: tr
      real(RP) :: beta
      real(RP) :: bmat(size(x), npt + size(x))
      real(RP) :: crvmin
      real(RP) :: d(size(x))
      real(RP) :: delbar
      real(RP) :: delta
      real(RP) :: distsq
      real(RP) :: dnorm
      real(RP) :: dnormsave(3)
      real(RP) :: fopt
      real(RP) :: moderr
      real(RP) :: fsave
      real(RP) :: fval(npt)
      real(RP) :: gq(size(x))
      real(RP) :: hq(size(x), size(x))
      real(RP) :: moderrsave(size(dnormsave))
      real(RP) :: pq(npt)
      real(RP) :: ratio
      real(RP) :: rho
      real(RP) :: trtol
      real(RP) :: vlag(npt + size(x))
      real(RP) :: vquad
      real(RP) :: xbase(size(x))
      real(RP) :: xdsq(npt)
      real(RP) :: xnew(size(x))
      real(RP) :: xopt(size(x))
      real(RP) :: xpt(size(x), npt)
      real(RP) :: zmat(npt, npt - size(x) - 1)
      logical :: improve_geometry
      logical :: reduce_rho
      logical :: shortd
      logical :: terminate
      character(len = 6), parameter :: solver = 'NEWUOA'
      character(len = SRNLEN), parameter :: srname = 'NEWUOB'


! Get size.
      n = int(size(x), kind(n))
      maxfhist = int(size(fhist), kind(maxfhist))
      maxxhist = int(size(xhist, 2), kind(maxxhist))

      if (DEBUGGING) then
          if (n == 0) then
              call errstop(srname, 'SIZE(X) is invalid')
          end if
          if (size(xhist, 1) /= n .and. maxxhist > 0) then
              call errstop(srname, 'XHIST is nonempty but SIZE(XHIST, 1)&
     & /= SIZE(X)')
          end if
          if (maxfhist*maxxhist > 0 .and. maxfhist /= maxxhist) then
              call errstop(srname, 'FHIST and XHIST are both nonempty bu&
     &t SIZE(FHIST) /= SIZE(XHIST, 2)')
          end if
      end if

      maxtr = maxfun ! Maximal number of trust region iterations.
      terminate = .false. ! Whether to terminate after initialization.

! Initialize FVAL, XBASE, and XPT.
      call initxf(calfun, iprint, x, rhobeg, ftarget, ij, kopt, nf, fhis&
     &t, fval, xbase, xhist, xpt, subinfo)
      xopt = xpt(:, kopt)
      fopt = fval(kopt)
      x = xbase + xopt ! Set X.
      f = fopt ! Set F.

! Check whether to return.
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
! Rearrange FHIST and XHIST so that they are in the chronological order.
          if (maxfhist >= 1 .and. maxfhist < nf) then
              khist = mod(nf - 1_IK, maxfhist) + 1_IK
              fhist = (/ fhist(khist + 1 : maxfhist), fhist(1 : khist) /&
     &)
          end if
          if (maxxhist >= 1 .and. maxxhist < nf) then
              khist = mod(nf - 1_IK, maxxhist) + 1_IK
              xhist = reshape((/ xhist(:, khist + 1 : maxxhist), xhist(:&
     &, 1 : khist) /), shape(xhist))
          end if
          return
      end if

! Initialize GQ, HQ, and PQ.
      call initq(ij, fval, xpt, gq, hq, pq, subinfo)

! Initialize BMAT and ZMAT, and IDZ.
      call inith(ij, xpt, idz, bmat, zmat, subinfo)

! After initializing GQ, HQ, PQ, BMAT, ZMAT, one can also choose to return
! if subinfo = NAN_MODEL (NaN occurs in the model). We do not do it here.
! If such a modle is harmful, then it will probably lead to other returns
! (NaN in X, NaN in F, trust region subproblem fails, ...); otherwise, the
! code will continue to run and possibly get rid of the NaN in the model.

! Set some more initial values.
      rho = rhobeg
      delta = rho
      moderrsave = HUGENUM
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
! NEWUOA never sets IMPROVE_GEOMETRY and REDUCE_RHO to TRUE simultaneously.

! Solve the trust region subproblem.
! In Powell's NEWUOA code, VQUAD is not an output of TRSAPP. Here we
! output it but will NOT use it; it will sill be calculated latter by
! CALQUAD in order to produce the same results as Powell's code.
          trtol = 1.0e-2_RP ! Tolerance used in trsapp.
          call trsapp(delta, gq, hq, pq, trtol, xopt, xpt, crvmin, vquad&
     &, d, subinfo)

! Calculate the length of the trial step D.
          dnorm = min(delta, sqrt(inprod(d, d)))

! Is the step long enough to invoke a function evaluation?
          if (dnorm < HALF*rho) then
              shortd = .true.
              if (maxval(abs(moderrsave)) <= 0.125_RP*crvmin*rho*rho .an&
     &d. maxval(dnormsave) <= rho) then
! Three recent values of ||D|| and |F−Q| are small.
! The 1st possibility (out of 2) that REDUCE_RHO = TRUE.
                  reduce_rho = .true.
              else
                  delta = TENTH*delta ! Reduce DELTA by a factor of 10.
                  if (delta <= 1.5_RP*rho) then
                      delta = rho ! Set DELTA to RHO when it is close.
                  end if
! After this, DELTA < DNORM may happen, explaining why we
! sometimes write MAX(DELTA, DNORM).
              end if
          end if

          if (.not. shortd) then ! This is the normal case.
! Save the current FOPT in FSAVE. It is needed later.
              fsave = fopt

! Shift XBASE if XOPT may be too far from XBASE.
!if (inprod(d, d) <= 1.0e-3_RP*xoptsq) then  ! Powell
              if (dnorm*dnorm <= 1.0e-3_RP*inprod(xopt, xopt)) then
                  call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xop&
     &t, xpt)
              end if

! Calculate VLAG and BETA for D.
              call vlagbeta(idz, kopt, bmat, d, xopt, xpt, zmat, beta, v&
     &lag)

! Use the current quadratic model to predict the change in F due
! to the step D.
              call calquad(d, gq, hq, pq, xopt, xpt, vquad)

! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x) ! Set F to NaN. It is necessary.
                  info = NAN_X
                  exit
              end if
              call calfun(x, f)
              nf = int(nf + 1, kind(nf))
              if (iprint >= 3) then
                  call fmssg(iprint, nf, f, x, solver)
              end if
              if (maxfhist >= 1) then
                  khist = mod(nf - 1_IK, maxfhist) + 1_IK
                  fhist(khist) = f
              end if
              if (maxxhist >= 1) then
                  khist = mod(nf - 1_IK, maxxhist) + 1_IK
                  xhist(:, khist) = x
              end if

! FQDIFF is the error of the current model in predicting the change
! in F due to D.
              moderr = f - fsave - vquad

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

! Calculate the reduction ratio and update DELTA accordingly.
              if (is_nan(vquad) .or. vquad >= ZERO) then
                  info = TRSUBP_FAILED
                  exit
              end if
              ratio = (f - fsave)/vquad
              delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ra&
     &tio)
              if (delta <= 1.5_RP*rho) then
                  delta = rho
              end if

! Set KNEW to the index of the interpolation point that will be
! replaced by XNEW. KNEW will ensure that the geometry of XPT
! is "optimal" after the replacement. Note that the information
! of XNEW is included in VLAG and BETA, which are calculated
! according to D = XNEW - XOPT.
! KNEW = 0 means it is not a good idea to replace any current
! interpolation point by XNEW.
              call setremove(idz, kopt, beta, delta, ratio, rho, vlag(1:&
     &npt), xopt, xpt, zmat, knew)

              if (knew > 0) then
! Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation
! point is replaced by XNEW.
                  call updateh(knew, beta, vlag, idz, bmat, zmat)

! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
                  call updateq(idz, knew, bmat(:, knew), moderr, zmat, x&
     &pt(:, knew), gq, hq, pq)

! Include the new interpolation point. This should be done
! after updating the model.
                  fval(knew) = f
                  xpt(:, knew) = xnew

! Update KOPT to KNEW if F < FSAVE (i.e., last FOPT).
                  if (f < fsave) then
                      kopt = knew
                  end if

                  if (delta <= rho) then ! Equivalent to DELTA == RHO.
! Test whether to replace the new quadratic model Q by the
! least Frobenius norm interpolant Q_alt. Perform the
! replacement if certain ceriteria is satisfied. This part
! is optional, but it is crucial for the performance on a
! certain class of problems. See Section 8 of the NEWUOA paper.
! In theory, the FVAL - FSAVE in the following line can be
! replaced by FVAL + C with any constant C. This constant
! will not affect the result in precise arithmetic. Powell
! chose C = - FSAVE.
! Since tryqalt is invoked only when DELTA equals the current
! RHO, why not reset ITEST to 0 when RHO is reduced?
                      call tryqalt(idz, fval - fsave, ratio, bmat(:, 1 :&
     & npt), zmat, itest, gq, hq, pq)
                  end if
              end if

! DNORMSAVE constains the DNORM corresponding to the latest 3
! function evaluations with the current RHO.
              dnormsave = (/ dnorm, dnormsave(1 : size(dnormsave) - 1) /&
     &)
! MODERR is the prediction errors of the latest 3 models.
              moderrsave = (/ moderr, moderrsave(1 : size(moderrsave) - &
     &1) /)
          end if

! The geometry of XPT probably needs improvement if
! 1. The trust region step D is not short but RATIO < TENTH or there
! is no good member in XPT to be replaced by XOPT + D (i.e., KNEW = 0), or
! 2. D is short but the latest three model errors are not small enough
! to render REDUCE_RHO = TRUE.
          if ((.not. shortd .and. (ratio < TENTH .or. knew == 0)) .or. (&
     &shortd .and. .not. reduce_rho)) then
! Find out if the interpolation points are close enough to the
! best point so far, i.e., all the points are within a ball
! centered at XOPT with a radius of 2*DELTA. If not, set KNEW to
! the index of the point that is the farthest away.
              distsq = 4.0_RP*delta*delta
              xdsq = sum((xpt - spread(xopt, dim = 2, ncopies = npt))**2&
     &, dim = 1)
              if (maxval(xdsq) > distsq) then
                  knew = int(maxloc(xdsq, dim = 1), kind(knew))
                  distsq = maxval(xdsq)
              else
                  knew = 0
              end if

! If KNEW is positive (i.e., not all points are close to XOPT),
! then a model step will be taken to ameliorate the geometry of
! the interpolation set and hence improve the model. Otherwise
! (all points are close to XOPT), RHO will be reduced
! (if MAX(DELTA, DNORM) <= RHO and D is "bad") or another
! trust-region step will be taken.
              if (knew > 0) then
! The only possibility that IMPROVE_GEOMETRY = TRUE.
                  improve_geometry = .true.
              else if (max(delta, dnorm) <= rho .and. (ratio <= 0 .or. s&
     &hortd)) then
! The 2nd possibility (out of 2) that REDUCE_RHO = TRUE.
! Even though all points are close to XOPT, a sufficiently
! small trust region does not suggest a good step to improve
! the current iterate. Then we should shrink RHO (i.e., update
! the stadard for defining "closeness" and shortd).
                  reduce_rho = .true.
              end if
          end if

! Before next trust region iteration, we may reduce rho or improve the
! geometry of XPT according to REDUCE_RHO and IMPROVE_GEOMETRY.

          if (reduce_rho) then
! The calculations with the current RHO are complete. Pick the
! next values of RHO and DELTA.
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
                  moderrsave = HUGENUM
                  if (iprint >= 2) then
                      call rhomssg(iprint, nf, fopt, rho, xbase + xopt, &
     &solver)
                 end if
              end if
          end if ! The procedure of reducing RHO ends.

          if (improve_geometry) then
! Save the current FOPT in fsave. It is needed later.
              fsave = fopt

! Set DELBAR, which will be used as the trust region radius for
! the geometry-improving schemes AMELIORGEO. We also need it to
! decide whether to shift XBASE or not.
              delbar = max(min(TENTH*sqrt(distsq), HALF*delta), rho)

! Shift XBASE if XOPT may be too far from XBASE.
              if (delbar*delbar<= 1.0e-3_RP*inprod(xopt, xopt)) then
                  call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xop&
     &t, xpt)
              end if

! Find a step D so that the geometry of XPT will be improved
! when XPT(:, KNEW) is replaced by XOPT + D. The AMELIORGEO
! subroutine will call Powell's BIGLAG and BIGDEN. It will also
! calculate the VLAG and BETA for this D.
              call ameliorgeo(idz, knew, kopt, bmat, delbar, xopt, xpt, &
     &zmat, d, beta, vlag)

! Use the current quadratic model to predict the change in F due
! to the step D.
              call calquad(d, gq, hq, pq, xopt, xpt, vquad)

! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x) ! Set F to NaN. It is necessary.
                  info = NAN_X
                  exit
              end if
              call calfun(x, f)
              nf = int(nf + 1, kind(nf))
              if (iprint >= 3) then
                  call fmssg(iprint, nf, f, x, solver)
              end if
              if (maxfhist >= 1) then
                  khist = mod(nf - 1_IK, maxfhist) + 1_IK
                  fhist(khist) = f
              end if
              if (maxxhist >= 1) then
                  khist = mod(nf - 1_IK, maxxhist) + 1_IK
                  xhist(:, khist) = x
              end if

! FQDIFF is the error of the current model in predicting the
! change in F due to D.
              moderr = f - fsave - vquad

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

! Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation
! point can be moved.
              call updateh(knew, beta, vlag, idz, bmat, zmat)

! Update the quadratic model.
              call updateq(idz, knew, bmat(:, knew), moderr, zmat, xpt(:&
     &, knew), gq, hq, pq)

! Include the new interpolation point. This should be done after
! updating BMAT, ZMAT, and the model.
              fval(knew) = f
              xpt(:, knew) = xnew
              if (f < fsave) then
                  kopt = knew
              end if

! DNORMSAVE constains the DNORM corresponding to the
! latest 3 function evaluations with the current RHO.
!--------------------------------------------------------------!
! Powell's code does not update DNORM. Therefore, DNORM is the
! length of last trust-region trial step, which seems inconsistent
! with what is described in Section 7 (around (7.7)) of the NEWUOA
! paper. Seemingly we should keep DNORM = ||D|| as we do here. The
! value of DNORM will be used when defining REDUCE_RHO.
              dnorm = min(delbar, sqrt(inprod(d, d)))
! In theory, DNORM = DELBAR in this case.
!--------------------------------------------------------------!
              dnormsave = (/ dnorm, dnormsave(1 : size(dnormsave) - 1) /&
     &)
! MODERR is the prediction errors of the latest 3 models.
              moderrsave = (/ moderr, moderrsave(1 : size(moderrsave) - &
     &1) /)
          end if ! The procedure of improving geometry ends.

      end do

! Return from the calculation, after another Newton-Raphson step, if it
! is too short to have been tried before.
! Note that no trust region iteration has been done if MAXTR = 0, and
! hence we should not check whether SHORTD = TRUE but return immediately.
      if (maxtr > 0 .and. shortd .and. nf < maxfun) then
          x = xbase + (xopt + d)
          if (any(is_nan(x))) then
              f = sum(x) ! Set F to NaN. It is necessary.
              info = NAN_X
          else
              call calfun(x, f)
              nf = int(nf + 1, kind(nf))
              if (iprint >= 3) then
                  call fmssg(iprint, nf, f, x, solver)
              end if
              if (maxfhist >= 1) then
                  khist = mod(nf - 1_IK, maxfhist) + 1_IK
                  fhist(khist) = f
              end if
              if (maxxhist >= 1) then
                  khist = mod(nf - 1_IK, maxxhist) + 1_IK
                  xhist(:, khist) = x
              end if
          end if
      end if

! Note that (FOPT .LE. F) is FALSE if F is NaN; When F is NaN, it is also
! necessary to update X and F.
      if (is_nan(f) .or. fopt <= f) then
          x = xbase + xopt
          f = fopt
      end if

! Rearrange FHIST and XHIST so that they are in the chronological order.
      if (maxfhist >= 1 .and. maxfhist < nf) then
          khist = mod(nf - 1_IK, maxfhist) + 1_IK
          fhist = (/ fhist(khist + 1 : maxfhist), fhist(1 : khist) /)
      end if
      if (maxxhist >= 1 .and. maxxhist < nf) then
          khist = mod(nf - 1_IK, maxxhist) + 1_IK
          xhist = reshape((/ xhist(:, khist + 1 : maxxhist), xhist(:, 1 &
     &: khist) /), shape(xhist))
      end if

      if (iprint >= 1) then
          call retmssg(info, iprint, nf, f, x, solver)
      end if

      return
      end subroutine newuob


      end module newuob_mod