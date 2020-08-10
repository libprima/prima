!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of newuoa.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 10-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! NEWUOA_MOD is a module providing a modern Fortran implementation of
! M. J. D. Powell's NEWUOA algorithm described in
!
! M. J. D. Powell, The NEWUOA software for unconstrained optimization
! without derivatives, In Large-Scale Nonlinear Optimization, eds. G. Di
! Pillo and M. Roma, pages 255--297, Springer, New York, US, 2006
!
! NEWUOA seeks the least value of a function of many variables, by a
! trust region method that forms quadratic models by interpolation.
! There can be some freedom in the interpolation conditions, which is
! taken up by minimizing the Frobenius norm of the change to the second
! derivative of the quadratic model, beginning with a zero matrix.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.


      module newuoa_mod

      implicit none
      private
      public :: newuoa


      contains


      subroutine newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfu&
     &n, npt, iprint, eta1, eta2, gamma1, gamma2, xhist, fhist, maxhist,&
     & info)

! Among all the arguments, only CALFUN, X, and F are required. The other
! arguments are OPTIONAL and you can neglect them unless you are farmiliar
! with the algorithm. If you do not specify an optional argument, it will
! be assigned the default value to it as will be explained later. For
! example, it is valid to call NEWUOA by
!
! call newuoa(calfun, x, f)
!
! or
!
! call newuoa(calfun, x, f, rhobeg = 5.0D-1, rhoend = 1.0D-3)
!
! A detailed introduction to the arguments is as follows.
!
! CALFUN
!   Input, subroutine.
!   CALFUN(X, F) should evaluate the objective function at the given real
!   vector X and set the value to the real scalar F. It must be provided
!   by the user.
!
! X
!   Input and outout, real vector.
!   As an input, X should be an N dimensional vector that contains the
!   initial values of the variables, N being the dimension of the problem.
!   As an output, X will be set to an approximate minimizer.
!
! F
!   Output, real scalar.
!   F will be set to the objective function value of the X at exit.
!
! NF
!   Output, integer scalar.
!   NF will be set to the number of function evaluations at exit.
!
! RHOBEG, RHOEND
!   Inputs, real scalars, default: RHOBEG = 1, RHOEND = 10^-6.
!   RHOBEG and RHOEND must be set to the initial and final values of a
!   trust region radius, so both must be positive with RHOEND <= RHOBEG.
!   Typically RHOBEG should be about one tenth of the greatest expected
!   change to a variable, and RHOEND should indicate the accuracy that is
!   required in the final values of the variables.
!
! FTARGET
!   Input, real scalar, default: - Infinity.
!   FTARGET is the target function value. The algorithm will terminate
!   when a point withi a function value <= FTARGET is found.
!
! MAXFUN
!   Input, integer scalar, default: 500N.
!   MAXFUN is the maximal number of function evaluations.
!
! NPT
!   Input, integer scalar, default: 2N + 1.
!   NPT is the number of interpolation conditions for each trust region
!   model. Its value must be in the interval [N+2, (N+1)(N+2)/2].
!
! IPRINT
!   Input, integer scalar, default: 0.
!   The value of IPRINT should be set to 0, 1, 2, 3, or 4, which controls
!   the amount of printing. Specifically, there is no output if IPRINT = 0,
!   and there is output only at the return if IPRINT = 1. Otherwise, each
!   new value of RHO is printed, with the best vector of variables so far
!   and the corresponding value of the objective function. Further, each
!   new value of F with its variables are output if IPRINT = 3. When
!   IPRINT = 4, all the output of IPRINT = 3 will be recorded in a file
!   named NEWUOA.output, which can be costly in terms of time and space;
!   the file will be created if it does not exist; the new output will be
!   appended to the end of this file if it already exists.
!
! ETA1, ETA2, GAMMA1, GAMMA2
!   Input, real scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5,
!   and GAMMA2 = 2.
!   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme
!   of the trust region radius as detailed in the subroutine TRRAD in
!   trustregion.f90. Roughly speaking, the trust region radius is contracted
!   by a factor of GAMMA1 when the reduction ratio is below ETA1, and
!   enlarged by a factor of GAMMA2 when the reduction ratio is above ETA2.
!
! XHIST, FHIST, MAXHIST
!   XHIST: Output, ALLOCATABLE RANK 2 real array;
!   FHIST: Output, ALLOCATABLE RANK 1 real array;
!   MAXHIST: Input, integer scalar, default: equal to MAXFUN.
!   XHIST, if present, will output the history of iterates, while FHIST,
!   if present, will output the histroy function values. MAXHIST should
!   be a nonnegative integer, and XHIST/FHIST will output only the last
!   MAXHIST iterates and/or the corresponding function values. Therefore,
!   MAXHIST = 0 means XHIST/FHIST will output nothing, while setting
!   MAXHIST = MAXFUN ensures that  XHIST/FHIST will output all the history.
!   If XHIST is present, its size at exit will be (N, min(NF, MAXHIST));
!   if FHIST is present, its size at exit will be min(NF, MAXHIST).
!   Note that setting MAXHIST to a large value may be costly in terms of
!   memory. For instance, if N = 1000 and MAXHIST = 100, 000, XHIST will
!   take about 1 GB if we use double precision.
!
! INFO
!   Output, integer scalar.
!   INFO is the exit flag. It can be set to the following values defined
!   in info.F:
!   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
!   FTARGET_ACHIEVED: the target function value is reached;
!   TRSUBP_FAILED: a trust region step failed to reduce the quadratic model;
!   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
!   NAN_X: NaN occurs in x;
!   NAN_INF_F: the objective function returns NaN or nearly infinite value;
!   NAN_MODEL: NaN occurs in the models.


! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, EPS
      use consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT&
     &_DFT, MAXIMAL_HIST
      use infnan_mod, only : is_nan, is_inf
      use memory_mod, only : safealloc

! Solver-specific module
      use prob_mod, only : FUNEVAL
      use newuob_mod, only : newuob

      implicit none

! Dummy variables
      procedure(FUNEVAL) :: calfun
      real(RP), intent(inout) :: x(:)
      real(RP), intent(out) :: f
      integer(IK), intent(out), optional :: nf
      real(RP), intent(in), optional :: rhobeg
      real(RP), intent(in), optional :: rhoend
      real(RP), intent(in), optional :: ftarget
      integer(IK), intent(in), optional :: maxfun
      integer(IK), intent(in), optional :: npt
      integer(IK), intent(in), optional :: iprint
      real(RP), intent(in), optional :: eta1
      real(RP), intent(in), optional :: eta2
      real(RP), intent(in), optional :: gamma1
      real(RP), intent(in), optional :: gamma2
      real(RP), intent(out), optional, allocatable :: fhist(:)
      real(RP), intent(out), optional, allocatable :: xhist(:, :)
      integer(IK), intent(in), optional :: maxhist
      integer(IK), intent(out), optional :: info

! Intermediate variables
      integer(IK) :: iprint_v
      integer(IK) :: maxfun_v
      integer(IK) :: maxfhist
      integer(IK) :: maxhist_v
      integer(IK) :: maxxhist
      integer(IK) :: n
      integer(IK) :: npt_v
      real(RP) :: eta1_v
      real(RP) :: eta2_v
      real(RP), allocatable :: fhist_c(:)
      real(RP) :: ftarget_v
      real(RP) :: gamma1_v
      real(RP) :: gamma2_v
      real(RP) :: rhobeg_v
      real(RP) :: rhoend_v
      real(RP), allocatable :: xhist_c(:, :)

! Get size
      n = int(size(x), kind(n))

! Replace any NaN in X by ZERO.
      where (is_nan(x))
          x = ZERO
      end where

! Verify and possibly revise the inputs. RHOBEG_V is the value of RHOBEG
! after verification. The others are similar.
      rhobeg_v = rhobeg
      rhoend_v = rhoend
      eta1_v = eta1
      eta2_v = eta2
      gamma1_v = gamma1
      gamma2_v = gamma2
      ftarget_v = ftarget
      maxfun_v = maxfun
      npt_v = npt
      iprint_v = iprint

! When the data is passed from the interfaces (e.g., MEX) to the Fortran
! code, RHOBEG, and RHOEND may change a bit. It was oberved in a MATLAB
! test that MEX passed 1 to Fortran as 0.99999999999999978. Therefore,
! if we set RHOEND = RHOBEG in the interfaces, then it may happen that
! RHOEND > RHOBEG, which is considered as an invalid input. To avoid
! this, we force RHOBEG and RHOEND to equal when their difference is tiny.
      if ((rhobeg_v - rhoend_v) < 1.0e2_RP*EPS*max(abs(rhobeg_v), ONE))t&
     &hen
          rhoend_v = rhobeg_v
      end if

      if (rhobeg_v <= 0 .or. is_nan(rhobeg_v) .or. is_inf(rhobeg_v))then
          rhobeg_v = RHOBEG_DFT
      end if
      rhobeg_v = max(EPS, rhobeg_v)

      if (rhoend_v < 0 .or. rhobeg_v < rhoend .or. is_nan(rhoend_v) .or.&
     & is_inf(rhoend_v)) then
          rhoend_v = min(TENTH*rhobeg_v, RHOEND_DFT)
      end if
      rhoend_v = max(EPS, rhoend_v)

      if (eta1_v < 0.0_RP .or. eta1_v > HALF .or. is_nan(eta1_v)) then
          eta1_v = TENTH
      end if

      if (eta2_v < eta1_v .or. eta2_v > 1.0_RP .or. is_nan(eta2_v)) then
         eta2_v = min(1.0_RP, max(eta1_v, 0.7_RP))
      end if

      if (gamma1_v <= 0.0_RP .or. gamma1_v >= 1.0_RP .or. is_nan(gamma1_&
     &v)) then
          gamma1_v = HALF
      end if

      if (gamma2_v <= 1.0_RP .or. is_nan(gamma2_v) .or. is_inf(gamma2_v)&
     &) then
          gamma2_v = TWO
      end if

      if (is_nan(ftarget_v)) then
          ftarget_v = FTARGET_DFT
      end if

      maxfun_v = max(int(n + 3, kind(maxfun_v)), maxfun_v)

      if (npt_v < n + 2 .or. npt > min(maxfun_v - 1, ((n + 2)*(n + 1))/2&
     &)) then
          npt_v = int(min(maxfun_v - 1, 2*n + 1), kind(npt_v))
      end if

      if (iprint_v /= 0 .and. iprint_v /= 1 .and. iprint_v /= 2 .and. ip&
     &rint_v /= 3 .and. iprint_v /= 4) then
          iprint_v = IPRINT_DFT
      end if

      if (present(maxhist)) then
          maxhist_v = max(0_IK, minval((/MAXIMAL_HIST, maxhist, maxfun_v&
     &/)))
      else
          maxhist_v = min(MAXIMAL_HIST, maxfun_v)
      end if

! Allocate memory for the histroy of X. We use XH instead of XHIST,
! which may not be present.
      if (present(xhist)) then
          maxxhist = maxhist_v
      else
          maxxhist = 0
      end if
      call safealloc(xhist_c, n, maxxhist)
! Allocate memory for the histroy of F. We use XH instead of FHIST,
! which may not be present.
      if (present(fhist)) then
          maxfhist = maxhist_v
      else
          maxfhist = 0
      end if
      call safealloc(fhist_c, maxfhist)

      call newuob(calfun, iprint_v, maxfun_v, npt_v, eta1_v, eta2_v, fta&
     &rget_v, gamma1_v, gamma2_v, rhobeg_v, rhoend_v, x, nf, f, fhist_c,&
     & xhist_c, info)

! Copy XH to XHIST and FH to FHIST if needed.
! N. B.: Fortran 2003 supports "automatic (re)allocation of allocatable
! arrays upon intrinsic assignment": if an intrinsic assignment is used,
! an allocatable variable on the left-hand side is automatically
! allocated (if unallocated) or reallocated (if the shape is different).
! In that case, the lines of SAFEALLOC in the following can be removed.
      if (present(xhist)) then
          call safealloc(xhist, n, min(nf, maxxhist))
          xhist = xhist_c(:, 1 : min(nf, maxxhist))
      end if
      deallocate(xhist_c)
      if (present(fhist)) then
          call safealloc(fhist, min(nf, maxfhist))
          fhist = fhist_c(1 : min(nf, maxfhist))
      end if
      deallocate(fhist_c)

      end subroutine newuoa


      end module newuoa_mod