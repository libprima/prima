!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of cobylb.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 10-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module cobylb_mod


      contains

      subroutine cobylb(m, x, rhobeg, rhoend, iprint, maxfun, con, f, in&
     &fo, ftarget, cstrv)

! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH,&
     & EPS, HUGENUM, DEBUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F, DAMAGING_ROUNDING
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : inprod, matprod, outprod
      use memory_mod, only : cstyle_sizeof

! Solver-specific modules
!use savex_mod, only : savex
      use initialize_mod, only : initialize
      use trustregion_mod, only : trstlp
      use update_mod, only : updatepole, findpole
      use geometry_mod, only : goodgeo, setdrop_geo, setdrop_tr, geostep
      use selectx_mod, only : selectx

      implicit none

! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: m
      integer(IK), intent(in) :: maxfun
      real(RP), intent(in) :: ftarget
      real(RP), intent(in) :: rhobeg
      real(RP), intent(in) :: rhoend

! In-outputs
      real(RP), intent(out) :: con(:)
      ! m+2. Bad name; should be confr
      real(RP), intent(inout) :: x(:)
      ! n

! Outputs
      integer(IK), intent(out) :: info
      real(RP), intent(out) :: f


! Parameters
! NSAVMAX is the maximal number of "dropped X" to save
      integer(IK), parameter :: nsavmax = 2000_IK
! CTOL is the tolerance for constraint violation. A point X is considered to be feasible if its
! constraint violation (CSTRV) is less than CTOL.
      real(RP), parameter :: ctol = EPS

! Local variables

      integer(IK) :: i
      integer(IK) :: tr
      integer(IK) :: maxtr
      integer(IK) :: ifull
      integer(IK) :: j
      integer(IK) :: jdrop
      integer(IK) :: jopt
      integer(IK) :: kopt
      integer(IK) :: n
      integer(IK) :: nf
      integer(IK) :: nsav
      integer(IK) :: subinfo
      real(RP) :: A(size(x), m + 1)
      ! Better name?
! A(:, 1:m) contains the approximate gradient for the constraints, and A(:, m+1) is minus the
! approximate gradient for the objective function.
      real(RP) :: barmu
      real(RP) :: cmax(m)
      real(RP) :: cmin(m)
      real(RP) :: cpen
      ! Penalty parameter for constraint in merit function (PARMU in Powell's code)
      real(RP) :: datmat(m + 2, size(x) + 1)
      ! CONVAL, FVAL, CVVAL
      real(RP) :: datsav(m + 2, max(nsavmax, 0))
      real(RP) :: denom
      real(RP) :: d(size(x))
      real(RP) :: factor_alpha
      real(RP) :: factor_beta
      real(RP) :: factor_delta
      real(RP) :: factor_gamma
      real(RP) :: prerec
      ! Predicted reduction in constraint violation
      real(RP) :: preref
      ! Predicted reduction in objective function
      real(RP) :: prerem
      ! Predicted reduction in merit function
      real(RP) :: cstrv
      real(RP) :: rho
      real(RP) :: sim(size(x), size(x) + 1)
      ! (n, )
      real(RP) :: simi(size(x), size(x))
      ! (n, )
      real(RP) :: simid(size(x))
      real(RP) :: simi_jdrop(size(x))
      real(RP) :: actrem
      real(RP) :: xsav(size(x), max(nsavmax, 0))
      real(RP) :: conopt(size(con))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY
      real(RP), allocatable :: xhist(:, :)
      real(RP), allocatable :: fhist(:)
      real(RP), allocatable :: conhist(:, :)
      real(RP), allocatable :: cstrvhist(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical :: improve_geo
      logical :: good_geo
      logical :: reduce_rho
      logical :: shortd


      character(len=SRNLEN), parameter :: srname = 'COBYLB'

      reduce_rho = .false.

      n = size(x)

! Set the initial values of some parameters. The last column of SIM holds the optimal vertex of the
! current simplex, and the preceding N columns hold the displacements from the optimal vertex to the
! other vertices.  Further, SIMI holds the inverse of the matrix that is contained in the first N
! columns of SIM.
      factor_alpha = QUART
      factor_beta = 2.1E0_RP
      factor_delta = 1.1E0_RP
      factor_gamma = HALF
      rho = rhobeg
      cpen = ZERO

      nsav = 0
      datsav = HUGENUM
      ! This is necessary; otherwise, SELECTX may return an incorrect X.
      datmat = HUGENUM
      ! This is necessary; otherwise, SELECTX may return an incorrect X.

      call initialize(iprint, maxfun, ctol, ftarget, rho, x, nf, datmat,&
     & sim, simi, subinfo)

      if (subinfo == NAN_X .or. subinfo == NAN_INF_F .or. subinfo == FTA&
     &RGET_ACHIEVED .or. subinfo == MAXFUN_REACHED) then
          info = subinfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          sim(:, 1:n) = sim(:, 1:n) + spread(sim(:, n + 1), dim=2, ncopi&
     &es=n)
          !!! TEMPORARY
          xhist = sim
          fhist = datmat(m + 1, :)
          conhist = datmat(1:m, :)
          cstrvhist = datmat(m + 2, :)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          cpen = min(1.0E8_RP, HUGENUM)
! Return the best calculated values of the variables.
          kopt = selectx(cpen, cstrvhist, ctol, fhist)
          x = xhist(:, kopt)
          f = fhist(kopt)
          cstrv = cstrvhist(kopt)
          con = conhist(:, kopt)
          return
      else
          x = sim(:, n + 1)
          f = datmat(m + 1, n + 1)
          cstrv = datmat(m + 2, n + 1)
          con = datmat(:, n + 1)
      end if

      maxtr = huge(tr)
      ! No constraint on the maximal number of trust-region iterations.

! Begin the iterative procedure.
! After solving a trust-region subproblem, COBYLA uses 3 boolean variables to control the work flow.
! SHORTD - Is the trust region trial step too short to invoke a function evaluation?
! IMPROVE_GEO - Will we improve the model after the trust region iteration? If yes, a geometry step
! will be taken, corresponding to the Branch (Delta) in the COBYLA paper.
! REDUCE_RHO - Will we reduce rho after the trust region iteration?
! COBYLA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
      do tr = 1, maxtr

! Before the trust-region step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
          call updatepole(cpen, [(.true., i=1, n + 1)], datmat, sim, sim&
     &i, subinfo)
          if (subinfo == DAMAGING_ROUNDING) then
              info = subinfo
              exit
          end if

! Does the current interpolation set has good geometry? It decides IMPROVE_GEO and REDUCE_RHO.
          good_geo = goodgeo(factor_alpha, factor_beta, rho, sim, simi)

! Calculate the linear approximations to the objective and constraint functions, placing minus
! the objective function gradient after the constraint gradients in the array A.
! N.B.: When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result
! as Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
! finite-precision arithmetic.
! QUESTION: Is it more reasonable to save A transpose instead of A? Better name for A?
          A = transpose(matprod(datmat(1:m + 1, 1:n) - spread(datmat(1:m&
     & + 1, n + 1), dim=2, ncopies=n), simi))
          A(:, m + 1) = -A(:, m + 1)


!!!!!!!!!!!!!!!!! Can this be removed? Is it safe for TRSTLP??????????
          if (any(is_nan(A))) then
              info = -3
              exit
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?????????????????????????????????

! Constraint and objective function values of the optimal vertex.
          conopt = datmat(:, n + 1)

! Calculate the trust-region trial step D.
          call trstlp(n, m, A, -conopt(1:m + 1), rho, d, ifull)
!write (16, *) ifull, d

! Is the trust-region trial step short?
! Is IFULL == 0 necessary ?????????????????????? If no, TRSTLP can be a function.
          shortd = (ifull == 0 .and. inprod(d, d) < QUART * rho * rho)

          if (.not. shortd) then
! Predict the change to F (PREREF) and to the constraint violation (PREREC) due to D.
              preref = inprod(d, A(:, m + 1))
              prerec = datmat(m + 2, n + 1) - maxval([-conopt(1:m) - mat&
     &prod(d, A(:, 1:m)), ZERO])

! Increase CPEN if necessary and branch back if this change alters the optimal vertex.
! Otherwise, PREREM will be set to the predicted reductions in the merit function.
! See the discussions around equation (9) of the COBYLA paper.
              barmu = -preref / prerec
              ! PREREF + BARMU * PREREC = 0
!!!!!!!!!!!!!!! Is it possible that PREREC <= 0????????????? It seems yes.
              if (prerec > ZERO .and. cpen < 1.5E0_RP * barmu) then
                  cpen = min(TWO * barmu, HUGENUM)
                  if (findpole(cpen, [(.true., i=1, n + 1)], datmat) <= &
     &n) then
                      cycle
                  end if
              end if

              prerem = preref + cpen * prerec
              ! Is it positive????

! Set X.
              x = sim(:, n + 1) + d
              if (any(is_nan(x))) then
                  f = sum(x)
                  ! Set F to NaN.
                  con = sum(x)
                  ! Set constraint values and constraint violation to NaN.
                  info = -1
                  exit
              end if

! Evaluate the objective function and constraints at X.
              call calcfc(n, m, x, f, con)
              nf = nf + 1
              cstrv = maxval([-con(1:m), ZERO])
              con(m + 1) = f
              con(m + 2) = cstrv

! Begin the operations that decide whether X should replace one of the vertices of the
! current simplex, the change being mandatory if ACTREM is positive.
              actrem = (datmat(m + 1, n + 1) + cpen * datmat(m + 2, n + &
     &1)) - (f + cpen * cstrv)
              if (cpen <= ZERO .and. abs(f - datmat(m + 1, n + 1)) <= ZE&
     &RO) then
                  prerem = prerec
                  ! Is it positive?????
                  actrem = datmat(m + 2, n + 1) - cstrv
              end if


! Set JDROP to the index of the vertex that is to be replaced by X.
              jdrop = setdrop_tr(actrem, d, factor_alpha, factor_delta, &
     &rho, sim, simi)

! When JDROP=0, the algorithm decides not to include X into the simplex.
              if (jdrop == 0) then
                  call savex(x, con, xsav, datsav, nsav, ctol)
                  !?????
              else
                  call savex(sim(:, n + 1) + sim(:, jdrop), datmat(:, jd&
     &rop), xsav, datsav, nsav, ctol)
! Revise the simplex by updating the elements of SIM, SIMI, and DATMAT.
                  sim(:, jdrop) = d
                  simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d&
     &)
                  simi = simi - outprod(matprod(simi, d), simi_jdrop)
                  simi(jdrop, :) = simi_jdrop
                  datmat(:, jdrop) = con
              end if

              if (is_nan(F) .or. is_posinf(F)) then
                  info = -2
                  exit
              end if
              if (any(is_nan(con(1:m)))) then
                  cstrv = sum(con(1:m))
                  ! Set CSTRV to NaN
                  info = -2
                  exit
              end if
              if (f <= ftarget .and. cstrv <= ctol) then
                  info = 1
                  exit
              end if
              if (nf >= maxfun) then
                  info = MAXFUN_REACHED
                  exit
              end if
          end if

! Should we take a geometry step to improve the geometry of the interpolation set?
          improve_geo = (shortd .or. actrem <= ZERO .or. actrem < TENTH &
     &* prerem) .and. .not. good_geo

! Should we revise RHO (and CPEN)?
          reduce_rho = (shortd .or. actrem <= ZERO .or. actrem < TENTH *&
     & prerem) .and. good_geo

          if (improve_geo) then
! Before the geometry step, call UPDATEPOLE so that SIM(:, N + 1) is the optimal vertex.
              call updatepole(cpen, [(.true., i=1, n + 1)], datmat, sim,&
     & simi, subinfo)
              if (subinfo == DAMAGING_ROUNDING) then
                  info = subinfo
                  exit
              end if

! If the current interpolation set has good geometry, then we skip the geometry step.
! The code has a small difference from the original COBYLA code here: If the current geometry
! is good, then we will continue with a new trust-region iteration; at the beginning of the
! iteration, CPEN may be updated, which may alter the pole point SIM(:, N + 1) by UPDATEPOLE;
! the quality of the interpolation point depends on SIM(:, N + 1), meaning that the same
! interpolation set may have good or bad geometry with respect to different "poles"; if the
! geometry turns out bad with the new pole, the original COBYLA code will take a geometry
! step, but the code here will NOT do it but continue to take a trust region step.
! The argument is this: even if the geometry step is not skipped at the first place, the
! geometry may turn out bad again after the pole is altered due to an update to CPEN; should
! we take another geometry step in that case? If no, why should we do it here? Indeed, this
! distinction makes no practical difference for CUTEst problems with at most 100 variables
! and 5000 constraints, while the algorithm framework is simplified.
              if (.not. goodgeo(factor_alpha, factor_beta, rho, sim, sim&
     &i)) then
! Decide a vertex to drop from the simplex. It will be replaced by SIM(:, N + 1) + D to
! improve acceptability of the simplex. See equations (15) and (16) of the COBYLA paper.
                  jdrop = setdrop_geo(factor_alpha, factor_beta, rho, si&
     &m, simi)

!Calculate the geometry step D.
                  d = geostep(jdrop, cpen, datmat, factor_gamma, rho, si&
     &mi)

! Save the information of the JOPT-th vertex in XSAV and DATSAV.
                  call savex(sim(:, n + 1) + sim(:, jdrop), datmat(:, jd&
     &rop), xsav, datsav, nsav, ctol)

! Update SIM and SIMI.
                  sim(:, jdrop) = d
                  ! Corresponding to the new vertex SIM(:, N + 1) + D
                  simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d&
     &)
                  simi = simi - outprod(matprod(simi, d), simi_jdrop)
                  simi(jdrop, :) = simi_jdrop

! Set X.
                  x = sim(:, n + 1) + d
                  if (any(is_nan(x))) then
                      f = sum(x)
                      ! Set F to NaN.
                      con = sum(x)
                      ! Set constraint values and constraint violation to NaN.
                      info = -1
                      exit
                  end if

! Evaluate the objective function and constraints at X.
                  call calcfc(n, m, x, f, con)
                  nf = nf + 1
                  cstrv = maxval([-con(1:m), ZERO])
                  con(m + 1) = f
                  con(m + 2) = cstrv
                  datmat(:, jdrop) = con

                  if (is_nan(F) .or. is_posinf(F)) then
                      info = -2
                      exit
                  end if
                  if (any(is_nan(con(1:m)))) then
                      cstrv = sum(con(1:m))
                      ! Set CSTRV to NaN
                      info = -2
                      exit
                  end if
                  if (f <= ftarget .and. cstrv <= ctol) then
                      info = 1
                      exit
                  end if
                  if (nf >= maxfun) then
                      info = MAXFUN_REACHED
                      exit
                  end if
              end if
          end if

          if (reduce_rho) then
          ! Update RHO and CPEN.
              if (rho <= rhoend) then
                  info = 0
                  exit
              end if
! See equation (11) in Section 3 of the COBYLA paper for the update of RHO.
              rho = HALF * rho
              if (rho <= 1.5E0_RP * rhoend) then
                  rho = rhoend
              end if
! See equations (12)--(13) in Section 3 of the COBYLA paper for the update of CPEN.
! If the original CPEN = 0, then the updated CPEN is also 0.
              cmin = minval(datmat(1:m, :), dim=2)
              cmax = maxval(datmat(1:m, :), dim=2)
              if (any(cmin < HALF * cmax)) then
                  denom = minval(max(cmax, ZERO) - cmin, mask=(cmin < HA&
     &LF * cmax))
                  cpen = min(cpen, (maxval(datmat(m + 1, :)) - minval(da&
     &tmat(m + 1, :))) / denom)
              else
                  cpen = ZERO
              end if
          end if
      end do

! Return the best calculated values of the variables.
      sim(:, 1:n) = sim(:, 1:n) + spread(sim(:, n + 1), dim=2, ncopies=n&
     &)
      !!! TEMPORARY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make sure that the history includes the last X.
      xhist = reshape([sim, xsav(:, 1:nsav), x], [n, n + nsav + 2])
      fhist = [datmat(m + 1, :), datsav(m + 1, 1:nsav), f]
      conhist = reshape([datmat(1:m, :), datsav(1:m, 1:nsav), con(1:m)],&
     & [m, n + nsav + 2])
      cstrvhist = [datmat(m + 2, :), datsav(m + 2, 1:nsav), cstrv]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cpen = max(cpen, min(1.0E8_RP, HUGENUM))
!write (16, *) fhist
!write (16, *) cstrvhist
      kopt = selectx(cpen, cstrvhist, ctol, fhist)
      x = xhist(:, kopt)
      f = fhist(kopt)
      cstrv = cstrvhist(kopt)
      con = conhist(:, kopt)

!write (16, *) kopt, f, cstrv
!close (16)
      close (11)

      end subroutine cobylb

      subroutine savex(xdrop, datdrop, xsav, datsav, nsav, ctol)
! This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless a vector in XSAV(:, 1:NSAV) is
! better than XDROP. If XDROP is better than some vectors in XSAV(:, 1:NSAV), then these vectors
! will be removed. If XDROP is not better than any of XSAV(:, 1:NSAV) but NSAV=NSAVMAX, then we
! remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NSAV).
!
! When COBYLA calls this subroutine, XDROP is a vector to be "dropped", and  DATDROP contains its
! function/constraint information (constraint value in the first M entries, DATDROP(M+1) = F(XDROP),
! and DATDROP(M+2) = CSTRV(X)). XSAV and DATSAV save at most NSAVMAX vectors "dropped" by COBYLB
! and their function/constraint information. Only XSAV(:, 1:NSAV) and DATSAV(:, 1:NSAV) contains
! such vectors, while XSAV(:, NSAV+1:NSAVMAX) and DATSAV(:, NSAV+1:NSAVMAX) are not initialized yet.
!
! Note: We decide whether X is better than the function/constraint of Y according to the ISBETTER
! function with CPEN = -ONE. Due to the implementation of ISBETTER,
! X is better than Y with CPEN < 0
! ==> X is better than Y with any CPEN >= 0,
! ==> X is better than Y regardless of CPEN.

! Generic modules
      use consts_mod, only : RP, IK, ONE
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod

! Solver-specific modules
      use selectx_mod, only : isbetter

      implicit none

! Inputs
      real(RP), intent(IN) :: ctol
      real(RP), intent(IN) :: datdrop(:)
      ! m+2
      real(RP), intent(IN) :: xdrop(:)
      ! n

! In-outputs
      integer(IK), intent(INOUT) :: nsav
      real(RP), intent(INOUT) :: datsav(:, :)
      ! (M+2, NSAVMAX)
      real(RP), intent(INOUT) :: xsav(:, :)
      ! (N, NSAVMAX)

! Local variables
      integer(IK) :: m
      integer(IK) :: n
      integer(IK) :: nsavmax
      integer(IK) :: i
      real(RP) :: cpen
      logical :: better(nsav)
      logical :: keep(nsav)

      m = size(datdrop) - 2
      n = size(xdrop)
      nsavmax = size(xsav, 2)

      if (nsavmax <= 0) then
          return
          ! Do nothing if NSAVMAX=0
      end if

      cpen = -ONE
      ! See the comments above for why CPEN = -1

! Return immediately if any column of XSAV is better than XDROP.
! BETTER is defined by the array constructor with an implicit do loop.
      better = [(isbetter([datsav(m + 1, i), datsav(m + 2, i)], [datdrop&
     &(m + 1), datdrop(m + 2)], ctol), i=1, nsav)]
      if (any(better)) then
          return
      end if

! Decide which columns of XSAV to keep. We use again the array constructor with an implicit do loop.
      keep = [(.not. isbetter([datdrop(m + 1), datdrop(m + 2)], [datsav(&
     &m + 1, i), datsav(m + 2, i)], ctol), i=1, nsav)]
! If XDROP is not better than any column of XSAV, then we remove the first (oldest) column of XSAV.
      if (count(keep) == nsavmax) then
          keep(1) = .false.
      end if
      xsav(:, 1:count(keep)) = xsav(:, pack([(i, i=1, nsav)], mask=keep)&
     &)
      datsav(:, 1:count(keep)) = datsav(:, pack([(i, i=1, nsav)], mask=k&
     &eep))

! Update NSAV. Note that the update of XSAV and DATSAV used NSAV, so it should be updated afterward.
      nsav = count(keep) + 1

! Save XDROP to XSAV(:, NSAV) and DATDROP to DATSAV(:, NSAV).
      xsav(:, nsav) = xdrop(:)
      datsav(:, nsav) = datdrop(:)

      return

      end subroutine savex


      end module cobylb_mod