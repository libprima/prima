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
! on 08-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! NEWUOB_MOD is a module that performs the major calculations of NEWUOA.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Last Modified: Thursday, July 08, 2021 PM03:05:22

      module newuob_mod

      implicit none
      private
      public :: newuob


      contains


      subroutine newuob(calfun, iprint, maxfun, npt, eta1, eta2, ftarget&
     &, gamma1, gamma2, rhobeg, rhoend, x, nf, f, fhist, xhist, info)
! NEWUOB performs the actual calculations of NEWUOA. The arguments IPRINT, MAXFUN, MAXHIST, NPT,
! ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, X, NF, F, FHIST, XHIST, and INFO are
! identical to the corresponding arguments in subroutine NEWUOA.

! XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
! of the model and Lagrange functions.
! XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
! the least calculated F so far). FOPT = F(XOPT + XBASE).
! D is reserved for trial steps from XOPT.
! XNEW = XOPT + D, corresponding to the vector of variables for the next calculation of F.
! [XPT, FVAL, KOPT] describes the interpolation set:
! XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
! However, there is a delay between the update of XOPT and KOPT. So they are not always consistent
! in the mid of an iteration. See the comment on the update of XOPT for details.
! values of F at the interpolation points; KOPT is the index of XOPT in XPT (XPT(:,KOPT) = XOPT).
! [GQ, HQ, PQ] describes the quadratic model: GQ will hold the gradient of the quadratic model at
! XBASE; HQ will hold the explicit second order derivatives of the quadratic model; PQ will contain
! the parameters of the implicit second order derivatives of the quadratic model.
! [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
! the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
! BMAT will hold the last N ROWs of H; ZMAT will hold the factorization of the leading NPT*NPT
! submatrix of H, this factorization being ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT)=1.
! VLAG will contain the values of the Lagrange functions at a new point X. They are part of a
! product that requires VLAG to be of length NPT + N. Both VLAG and BETA are critical for the
! updating procedure of H, which is detailed formula (4.10)--(4.12) of the NEWUOA paper.
!
! See Section 2 of the NEWUOA paper for more information about these variables.

! Generic modules
      use consts_mod, only : RP, IK, ZERO, TWO, HALF, TENTH, HUGENUM, DE&
     &BUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod

! Solver-specific modules
      use pintrf_mod, only : FUNEVAL
      use initialize_mod, only : initxf, initq, inith
      use trustregion_mod, only : trsapp, trrad
      use geometry_mod, only : setdrop, geostep
      use shiftbase_mod, only : shiftbase
      use vlagbeta_mod, only : vlagbeta
      use update_mod, only : updateh, updateq, tryqalt

      implicit none

! Inputs
      procedure(FUNEVAL) :: calfun
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
      real(RP), intent(inout) :: x(:)
      ! SIZE(X) = N

! Outputs
      integer(IK), intent(out) :: info
      integer(IK), intent(out) :: nf
      real(RP), intent(out) :: f
      real(RP), intent(out) :: fhist(:)
      real(RP), intent(out) :: xhist(:, :)

! Local variables
      integer(IK) :: idz
      integer(IK) :: ij(2, npt)
      integer(IK) :: itest
      integer(IK) :: khist
      integer(IK) :: knew_geo
      integer(IK) :: knew_tr
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
      real(RP) :: dnorm
      real(RP) :: dnormsave(3)
      real(RP) :: fopt
      real(RP) :: fval(npt)
      real(RP) :: gq(size(x))
      real(RP) :: hq(size(x), size(x))
      real(RP) :: moderr
      real(RP) :: moderrsave(size(dnormsave))
      real(RP) :: pq(npt)
      real(RP) :: ratio
      real(RP) :: rho
      real(RP) :: rho_ratio
      real(RP) :: trtol
      real(RP) :: vlag(npt + size(x))
      real(RP) :: vquad
      real(RP) :: xbase(size(x))
      real(RP) :: xdist(npt)
      real(RP) :: xnew(size(x))
      real(RP) :: xopt(size(x))
      real(RP) :: xpt(size(x), npt)
      real(RP) :: zmat(npt, npt - size(x) - 1)
      logical :: improve_geo
      logical :: reduce_rho_1
      logical :: reduce_rho_2
      logical :: shortd
      character(len=6), parameter :: solver = 'NEWUOA'
      character(len=SRNLEN), parameter :: srname = 'NEWUOB'


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
          if (maxfhist * maxxhist > 0 .and. maxfhist /= maxxhist) then
              call errstop(srname, 'FHIST and XHIST are both nonempty bu&
     &t SIZE(FHIST) /= SIZE(XHIST, 2)')
          end if
      end if

      maxtr = maxfun
      ! Maximal number of trust region iterations.

! Initialize FVAL, XBASE, and XPT.
      call initxf(calfun, iprint, ftarget, rhobeg, x, ij, kopt, nf, fhis&
     &t, fval, xbase, xhist, xpt, subinfo)
      xopt = xpt(:, kopt)
      fopt = fval(kopt)
      x = xbase + xopt
      ! Set X.
      f = fopt
      ! Set F.

! Check whether to return after initialization.
      if (subinfo == FTARGET_ACHIEVED .or. subinfo == NAN_X .or. subinfo&
     & == NAN_INF_F) then
! In these cases, pack the data and return immediately.
          info = subinfo
          if (abs(iprint) >= 1) then
              call retmssg(info, iprint, nf, f, x, solver)
          end if
! Rearrange FHIST and XHIST so that they are in the chronological order.
          if (maxfhist >= 1 .and. maxfhist < nf) then
              khist = mod(nf - 1_IK, maxfhist) + 1_IK
              fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
          end if
          if (maxxhist >= 1 .and. maxxhist < nf) then
              khist = mod(nf - 1_IK, maxxhist) + 1_IK
              xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:&
     &khist)], shape(xhist))
! The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
! order of Fortran arrays.
          end if
          return
      end if

! Initialize GQ, HQ, and PQ.
      call initq(ij, fval, xpt, gq, hq, pq, subinfo)

! Initialize BMAT and ZMAT, and IDZ.
      call inith(ij, xpt, idz, bmat, zmat, subinfo)

! After initializing GQ, HQ, PQ, BMAT, ZMAT, one can also choose to return if subinfo = NAN_MODEL
! (NaN occurs in the model). We do not do it here. If such a model is harmful, then it will probably
! lead to other returns (NaN in X, NaN in F, trust region subproblem fails, ...); otherwise, the
! code will continue to run and possibly get rid of the NaN in the model.

! Set some more initial values and parameters.
      rho = rhobeg
      delta = rho
      moderrsave = HUGENUM
      dnormsave = HUGENUM
      itest = 0
      trtol = 1.0E-2_RP
      ! Tolerance used in trsapp.

! Begin the iterative procedure.
! After solving a trust-region subproblem, NEWUOA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust region iteration?
! REDUCE_RHO - Will we reduce rho after the trust region iteration?
! REDUCE_RHO = REDUCE_RHO_1 .OR. REDUCE_RHO_2 (see boxes 14 and 10 of Fig. 1 in the NEWUOA paper).
! NEWUOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
      do tr = 1, maxtr
! Solve the trust region subproblem. In Powell's NEWUOA code, VQUAD is not an output of TRSAPP.
! Here we output it but will NOT use it (for the moment); it will still be calculated later
! by CALQUAD in order to produce the same results as Powell's code.
          call trsapp(delta, gq, hq, pq, trtol, xopt, xpt, crvmin, vquad&
     &, d, subinfo)

! Calculate the length of the trial step D.
          dnorm = min(delta, sqrt(inprod(d, d)))

! SHORTD corresponds to Box 3 of the NEWUOA paper.
          shortd = (dnorm < HALF * rho)
! REDUCE_RHO_1 corresponds to Box 14 of the NEWUOA paper.
          reduce_rho_1 = shortd .and. (maxval(abs(moderrsave)) <= 0.125_&
     &RP * crvmin * rho * rho) .and. (maxval(dnormsave) <= rho)
          if (shortd .and. (.not. reduce_rho_1)) then
! Reduce DELTA. After this, DELTA < DNORM may hold.
              delta = TENTH * delta
              if (delta <= 1.5_RP * rho) then
                  delta = rho
                  ! Set DELTA to RHO when it is close.
              end if
          end if

          if (.not. shortd) then
          ! D is long enough.
! Shift XBASE if XOPT may be too far from XBASE.
!if (inprod(d, d) <= 1.0e-3_RP*inprod(xopt, xopt)) then  ! Powell
              if (dnorm * dnorm <= 1.0E-3_RP * inprod(xopt, xopt)) then
                  call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xop&
     &t, xpt)
              end if

! Calculate VLAG and BETA for D. It makes uses of XOPT, so this is done before updating XOPT.
              call vlagbeta(idz, kopt, bmat, d, xpt, zmat, beta, vlag)

! Use the current quadratic model to predict the change in F due to the step D.
              call calquad(d, gq, hq, pq, xopt, xpt, vquad)

! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)
                  ! Set F to NaN. It is necessary.
                  info = NAN_X
                  exit
              end if
              call calfun(x, f)
              nf = int(nf + 1, kind(nf))
              if (abs(iprint) >= 3) then
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

! DNORMSAVE constains the DNORM of the latest 3 function evaluations with the current RHO.
              dnormsave = [dnorm, dnormsave(1:size(dnormsave) - 1)]

! MODERR is the error of the current model in predicting the change in F due to D.
              moderr = f - fopt - vquad
! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
              moderrsave = [moderr, moderrsave(1:size(moderrsave) - 1)]

! Calculate the reduction ratio and update DELTA accordingly.
              if (is_nan(vquad) .or. vquad >= ZERO) then
                  info = TRSUBP_FAILED
                  exit
              end if
              ratio = (f - fopt) / vquad
! Update DELTA. After this, DELTA < DNORM may hold.
              delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ra&
     &tio)
              if (delta <= 1.5_RP * rho) then
                  delta = rho
              end if

! Update XOPT and FOPT. Before KOPT is updated, XOPT may differ from XPT(:, KOPT), and FOPT
! may differ from FVAL(KOPT). Note that the code may exit before KOPT is updated. See below.
! The updated XOPT is needed by SETDROP.
              if (f < fopt) then
                  xopt = xnew
                  fopt = f
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

! Set KNEW_TR to the index of the interpolation point that will be replaced by XNEW. KNEW_TR
! will ensure that the geometry of XPT is "good enough" after the replacement. Note that the
! information of XNEW is included in VLAG and BETA, which are calculated according to
! D = XNEW - XOPT. KNEW_TR = 0 means it is impossible to obtain a good interpolation set
! by replacing any current interpolation point with XNEW.
              call setdrop(idz, kopt, beta, delta, ratio, rho, vlag(1:np&
     &t), xopt, xpt, zmat, knew_tr)

              if (knew_tr > 0) then
! If KNEW_TR > 0, then update BMAT, ZMAT and IDZ, so that the KNEW_TR-th interpolation
! point is replaced by XNEW. If KNEW_TR = 0, then probably the geometry of XPT needs
! improvement, which will be handled below.
                  call updateh(knew_tr, beta, vlag, idz, bmat, zmat)

! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
                  call updateq(idz, knew_tr, bmat(:, knew_tr), moderr, z&
     &mat, xpt(:, knew_tr), gq, hq, pq)

! Include the new interpolation point.
                  xpt(:, knew_tr) = xnew
                  ! Should be done after UPDATEQ.
                  fval(knew_tr) = f
                  if (fval(knew_tr) < fval(kopt)) then
                      kopt = knew_tr
                  end if
! KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW_TR) = FVAL(KOPT) and
! KNEW_TR < KOPT, then MINLOC(FVAL) = KNEW_TR /= KOPT. Do not change KOPT in this case.
              end if

! Test whether to replace the new quadratic model Q by the least-Frobenius norm interpolant
! Q_alt. Perform the replacement if certain criteria are satisfied.
! N.B.:
! 1. This part is OPTIONAL, but it is crucial for the performance on some problems. See
! Section 8 of the NEWUOA paper.
! 2. TRYQALT is called only after a trust-region step but not after a geometry step, maybe
! because the model is expected to be good after a geometry step.
! 3. If KNEW_TR = 0 after a trust-region step, TRYQALT is not invoked. In this case, the
! interpolation set is unchanged, so it seems reasonable to keep the model unchanged.
! 4. In theory, FVAL - FOPT in the call of TRYQALT can be replaced by FVAL + C with any
! constant C. This constant will not affect the result in precise arithmetic. Powell chose
! C = - FVAL(KOPT_OLD), where KOPT_OLD is the KOPT before the update above (Powell updated
! KOPT after TRYQALT). Here we use C = -FOPT, as it worked slightly better on CUTEst,
! although there is no difference theoretically. Note that FVAL(KOPT_OLD) may not equal
! FOPT_OLD --- it may happen that KNEW_TR = KOPT_OLD so that FVAL(KOPT_OLD) has been revised
! after the last function evaluation.
! 5. Question: Since TRYQALT is invoked only when DELTA equals the current RHO, why not
! reset ITEST to 0 when RHO is reduced?
              if (knew_tr > 0 .and. delta <= rho) then
              ! DELTA = RHO.
                  call tryqalt(idz, fval - fopt, ratio, bmat(:, 1:npt), &
     &zmat, itest, gq, hq, pq)
              end if
          end if
          ! End of if (.not. shortd)

! Before next trust region iteration, we may improve the geometry of XPT or reduce rho
! according to IMPROVE_GEO and REDUCE_RHO. Now we decide these two indicators.

! First define IMPROVE_GEO, which corresponds to Box 8 of the NEWUOA paper.
! The geometry of XPT likely needs improvement if the trust-region step bad --- either too short
! (SHORTD = TRUE) or the reduction ratio is small (RATIO < TENTH). However, if REDUCE_RHO_1 is
! TRUE, meaning that the step is short and the latest model errors have been small, then we do
! not need to improve the geometry; instead, RHO will be reduced.
! To improve the geometry of XPT, we will check whether the interpolation points are all close
! enough to the best point so far, i.e., all the points are within a ball centered at XOPT with
! a radius of 2*DELTA. If not, the farthest point will be replaced with a point selected by
! GEOSTEP, aiming to ameliorate the geometry of the interpolation set.
! N.B.:
! 1. RATIO is set if SHORTD = FALSE. So the expression (SHORTD .OR. RATIO < TENTH) will not
! suffer from unset RATIO.
! 2. If SHORTD = FALSE and KNEW_TR = TRUE, then IMPROVE_GEO = TRUE, because KNEW_TR = TRUE
! necessitates RATIO <= 0 < TENTH. Therefore, IMPROVE_GEO = TRUE if it is impossible to obtain
! a good XPT by replacing a current point with the one suggested by the trust region step. This
! is reasonable.
! 3. If REDUCE_RHO = FALSE and SHORTD = TRUE, then the trust-region step is not tried at all,
! i.e., no function evaluation is invoked at XOPT + D (When REDUCE_RHO = TRUE, the step is not
! tried either, but the same step will be generated again at the next trust-region iteration
! after RHO is reduced and DELTA is updated; see the end of Section 2 of the NEWUOA paper).
! 4. If SHORTD = FALSE and KNEW_TR = 0, then the trust-region step invokes a function evaluation
! at XOPT + D, but [XOPT + D, F(XOPT +D)] is not included into [XPT, FVAL]. In other words, this
! function value is discarded. Note that KNEW_TR = 0 only if RATIO <= 0 (see SETDROP), so that
! a function value that renders a reduction is never discarded.
! 5. If SHORTD = FALSE and KNEW_TR > 0 and RATIO < TENTH, then [XPT, FVAL] is updated so that
! [XPT(KNEW_TR), FVAL(KNEW_TR)] = [XOPT + D, F(XOPT + D)], and the model is updated accordingly,
! but such a model will not be used in the next trust-region iteration, because a geometry step
! will be invoked to improve the geometry of the interpolation set and update the model again.
! 6. DELTA has been updated before arriving here: if REDUCE_RHO = FALSE and SHORTD = TRUE, then
! DELTA was reduced by a factor of 10; if SHORTD = FALSE, then DELTA was updated by TRRAD after
! the trust-region iteration.
! 7. When SHORTD = FALSE and KNEW_TR > 0, then XPT has been updated after the trust-region
! iteration; if RATIO > 0 in addition, then XOPT has been updated as well.
          xdist = sqrt(sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, &
     &dim=1))
          knew_geo = int(maxloc(xdist, dim=1), kind(knew_geo))
          improve_geo = (.not. reduce_rho_1) .and. (shortd .or. ratio < &
     &TENTH) .and. (maxval(xdist) > TWO * delta)
! ------------------------------------------------------------------------------------------!
! Modifying IMPROVE_GEO in the following way seems to make little difference in the performance,
! sometimes worsening, sometimes improving, but never substantially. The advantage of this
! IMPROVE_GEO is that the bound for RATIO is 0, the same as in SETDROP and REDUCE_RHO_2.
!improve_geo = (.not. reduce_rho_1) .and. (shortd .or. ratio <= ZERO) .and. (maxval(xdist) > TWO * delta)
! ------------------------------------------------------------------------------------------!

          if (improve_geo) then
! Set DELBAR, which will be used as the trust region radius for the geometry-improving
! scheme GEOSTEP. We also need it to decide whether to shift XBASE or not.
! Note that DELTA has been updated before arriving here. See the comments above the
! definition of IMPROVE_GEO.
              delbar = max(min(TENTH * maxval(xdist), HALF * delta), rho&
     &)

! Shift XBASE if XOPT may be too far from XBASE.
              if (delbar * delbar <= 1.0E-3_RP * inprod(xopt, xopt)) the&
     &n
                  call shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xop&
     &t, xpt)
              end if

! Find a step D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) is
! replaced by XOPT + D. The GEOSTEP subroutine will call Powell's BIGLAG and BIGDEN. It will
! also calculate the VLAG and BETA for this D.
              call geostep(idz, knew_geo, kopt, bmat, delbar, xpt, zmat,&
     & d, beta, vlag)

! Use the current quadratic model to predict the change in F due to the step D.
              call calquad(d, gq, hq, pq, xopt, xpt, vquad)

! Calculate the next value of the objective function.
              xnew = xopt + d
              x = xbase + xnew
              if (any(is_nan(x))) then
                  f = sum(x)
                  ! Set F to NaN. It is necessary.
                  info = NAN_X
                  exit
              end if
              call calfun(x, f)
              nf = int(nf + 1, kind(nf))
              if (abs(iprint) >= 3) then
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

! DNORMSAVE constains the DNORM of the latest 3 function evaluations with the current RHO.
!------------------------------------------------------------------------------------------!
! Powell's code does not update DNORM. Therefore, DNORM is the length of last trust-region
! trial step, which seems inconsistent with what is described in Section 7 (around (7.7)) of
! the NEWUOA paper. Seemingly we should keep DNORM = ||D|| as we do here. The value of DNORM
! will be used when defining REDUCE_RHO.
              dnorm = min(delbar, sqrt(inprod(d, d)))
! In theory, DNORM = DELBAR in this case.
!------------------------------------------------------------------------------------------!
              dnormsave = [dnorm, dnormsave(1:size(dnormsave) - 1)]

! MODERR is the error of the current model in predicting the change in F due to D.
              moderr = f - fopt - vquad
! MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
              moderrsave = [moderr, moderrsave(1:size(moderrsave) - 1)]

! Update XOPT and FOPT. Before KOPT is updated, XOPT may differ from XPT(:, KOPT), and FOPT
! may differ from FVAL(KOPT). Note that the code may exit before KOPT is updated. See below.
              if (f < fopt) then
                  xopt = xnew
                  fopt = f
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

! Update BMAT, ZMAT and IDZ, so that the KNEW_GEO-th interpolation point can be moved.
              call updateh(knew_geo, beta, vlag, idz, bmat, zmat)

! Update the quadratic model using the updated BMAT, ZMAT, IDZ.
              call updateq(idz, knew_geo, bmat(:, knew_geo), moderr, zma&
     &t, xpt(:, knew_geo), gq, hq, pq)

! Include the new interpolation point.
              xpt(:, knew_geo) = xnew
              ! Should be done after UPDATEQ.
              fval(knew_geo) = f
              if (fval(knew_geo) < fval(kopt)) then
                  kopt = knew_geo
              end if
! KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW_GEO) = FVAL(KOPT) and
! KNEW_GEO < KOPT, then MINLOC(FVAL) = KNEW_GEO /= KOPT. Do not change KOPT in this case.
          end if
          ! The procedure of improving geometry ends.

! If all the interpolation points are close to XOPT and the trust region is small, but the
! trust-region step is "bad" (SHORTD or RATIO <= 0), then we shrink RHO (update the criterion
! for the "closeness" and SHORTD). REDUCE_RHO_2 corresponds to Box 10 of the NEWUOA paper.
! N.B.:
! 1. Even though DNORM gets a new value after the geometry step when IMPROVE_GEO = TRUE, this
! value does not affect REDUCE_RHO_2, because DNORM comes into play only if IMPROVE_GEO = FALSE.
! 2. DELTA < DNORM may hold due to the update of DELTA.
! 3. The following two lines are equivalent.
!reduce_rho_2 = (.not. improve_geo) .and. (max(delta,dnorm)<=rho) .and. (shortd .or. ratio <= 0)
          reduce_rho_2 = (maxval(xdist) <= TWO * delta) .and. (max(delta&
     &, dnorm) <= rho) .and. (shortd .or. ratio <= ZERO)
! ------------------------------------------------------------------------------------------!
! Modifying REDUCE_RHO_2 in the following way seems to make little difference. The advantage
! of this REDUCE_RHO_2 is that RATIO has the same bound as in the definition of IMPROVE_GEO.
!reduce_rho_2 = (maxval(xdist) <= TWO * delta) .and. (max(delta, dnorm) <= rho) .and. (shortd .or. ratio < TENTH)
! ------------------------------------------------------------------------------------------!

          if (reduce_rho_1 .or. reduce_rho_2) then
! The calculations with the current RHO are complete. Pick the next values of RHO and DELTA.
              if (rho <= rhoend) then
                  info = SMALL_TR_RADIUS
                  exit
              else
                  delta = HALF * rho
                  rho_ratio = rho / rhoend
                  if (rho_ratio <= 16.0_RP) then
                      rho = rhoend
                  else if (rho_ratio <= 250.0_RP) then
                      rho = sqrt(rho_ratio) * rhoend
                  else
                      rho = TENTH * rho
                  end if
                  delta = max(delta, rho)
! DNORMSAVE and MODERRSAVE are corresponding to the latest 3 function evaluations with
! the current RHO. Update them after reducing RHO.
                  dnormsave = HUGENUM
                  moderrsave = HUGENUM
                  if (abs(iprint) >= 2) then
                      call rhomssg(iprint, nf, fopt, rho, xbase + xopt, &
     &solver)
                  end if
              end if
          end if
          ! The procedure of reducing RHO ends.

      end do
      ! The iterative procedure ends.

! Return from the calculation, after another Newton-Raphson step, if it is too short to have been
! tried before. Note that no trust region iteration has been done if MAXTR = 0, and hence we
! should not check whether SHORTD = TRUE but return immediately.
      if (maxtr > 0 .and. shortd .and. nf < maxfun) then
          x = xbase + (xopt + d)
          if (any(is_nan(x))) then
              f = sum(x)
              ! Set F to NaN. It is necessary.
              info = NAN_X
          else
              call calfun(x, f)
              nf = int(nf + 1, kind(nf))
              if (abs(iprint) >= 3) then
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

! Note that (FOPT <= F) = FALSE if F is NaN; if F is NaN, it is necessary to update X and F.
      if (is_nan(f) .or. fopt <= f) then
          x = xbase + xopt
          f = fopt
      end if

! Rearrange FHIST and XHIST so that they are in the chronological order.
      if (maxfhist >= 1 .and. maxfhist < nf) then
          khist = mod(nf - 1_IK, maxfhist) + 1_IK
          fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
      end if
      if (maxxhist >= 1 .and. maxxhist < nf) then
          khist = mod(nf - 1_IK, maxxhist) + 1_IK
          xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khis&
     &t)], shape(xhist))
! The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
! order of Fortran arrays.
      end if

      if (abs(iprint) >= 1) then
          call retmssg(info, iprint, nf, f, x, solver)
      end if

      return
      end subroutine newuob


      end module newuob_mod