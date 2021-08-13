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
! on 14-Aug-2021.
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
!
! Last Modified: Thursday, August 12, 2021 PM07:55:26

      module newuoa_mod

      implicit none
      private
      public :: newuoa


      contains


      subroutine newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfu&
     &n, npt, iprint, eta1, eta2, gamma1, gamma2, xhist, fhist, maxhist,&
     & info)

! Among all the arguments, only CALFUN, X, and F are obligatory. The others
! are OPTIONAL and you can neglect them unless you are familiar with the
! algorithm. If you do not specify an optional argument, it will be assigned
! the default value  as will be explained later. For example, it is valid to
! call NEWUOA by
!
! call newuoa(calfun, x, f)
!
! or
!
! call newuoa(calfun, x, f, rhobeg = 0.5D0, rhoend = 1.0D-3, maxfun = 100)
!
! See example.f90 for a concrete example.
!
! A detailed introduction to the arguments is as follows.
! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under
! the directory name "common". By default, RP = kind(0.0D0) and IK = kind(0).
! Therefore, REAL(RP) is the double-precision real, and INTEGER(IK) is the
! default integer. For ADVANCED USERS, RP and IK can be defined by specifying
! __REAL_PRECISION__ and __INTEGER_KIND__ in common/ppf.h. Use the default if
! you are unsure.
!
! CALFUN
!   Input, subroutine.
!   CALFUN(X, F) should evaluate the objective function at the given
!   REAL(RP) vector X and set the value to the REAL(RP) scalar F. It
!   must be provided by the user.
!
! X
!   Input and outout, REAL(RP) vector.
!   As an input, X should be an N dimensional vector that contains the
!   initial values of the variables, N being the dimension of the problem.
!   As an output, X will be set to an approximate minimizer.
!
! F
!   Output, REAL(RP) scalar.
!   F will be set to the objective function value of the X at exit.
!
! NF
!   Output, INTEGER(IK) scalar.
!   NF will be set to the number of function evaluations at exit.
!
! RHOBEG, RHOEND
!   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6.
!   RHOBEG and RHOEND must be set to the initial and final values of a
!   trust region radius, so both must be positive with RHOEND <= RHOBEG.
!   Typically RHOBEG should be about one tenth of the greatest expected
!   change to a variable, and RHOEND should indicate the accuracy that is
!   required in the final values of the variables.
!
! FTARGET
!   Input, REAL(RP) scalar, default: - Infinity.
!   FTARGET is the target function value. The algorithm will terminate
!   when a point withi a function value <= FTARGET is found.
!
! MAXFUN
!   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with
!   MAXFUN_DIM_DFT defined in the module CONSTS_MOD (see consts.F90 in
!   the directory named "common").
!   MAXFUN is the maximal number of function evaluations.
!
! NPT
!   Input, INTEGER(IK) scalar, default: 2N + 1.
!   NPT is the number of interpolation conditions for each trust region
!   model. Its value must be in the interval [N+2, (N+1)(N+2)/2].
!
! IPRINT
!   Input, INTEGER(IK) scalar, default: 0.
!   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which
!   controls how much information will be printed during the computation:
!   0: there will be no printing;
!   1: a message will be printed to the screen at the return, showing the
!      best vector of veriables found and its objective function value;
!   2: in addition to 1, each new value of RHO is printed to the screen,
!      with the best vector of variables so far and its objective function
!      value;
!   3: in addition to 2, each function evaluation with its variables will
!      be printed to the screen;
!   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to
!     the screen but to a file named NEWUOA_output.txt; the file will be
!     created if it does not exist; the new output will be appended to
!     the end of this file if it already exists. Note that IPRINT = -3 can
!     be costly in terms of time and space.
!
! ETA1, ETA2, GAMMA1, GAMMA2
!   Input, REAL(RP) scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5,
!   and GAMMA2 = 2.
!   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme
!   of the trust region radius as detailed in the subroutine TRRAD in
!   trustregion.f90. Roughly speaking, the trust region radius is contracted
!   by a factor of GAMMA1 when the reduction ratio is below ETA1, and
!   enlarged by a factor of GAMMA2 when the reduction ratio is above ETA2.
!   It is required that 0 < ETA1 <= ETA2 < 1 and 0 < GAMMA1 < 1 < GAMMA2.
!   Normally, ETA1 <= 0.25. It is NOT recommended to set ETA1 >= 0.5.
!
! XHIST, FHIST, MAXHIST
!   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
!   XHIST, if present, will output the history of iterates, while FHIST,
!   if present, will output the history function values. MAXHIST should
!   be a nonnegative integer, and XHIST/FHIST will output only the last
!   MAXHIST iterates and/or the corresponding function values. Therefore,
!   MAXHIST = 0 means XHIST/FHIST will output nothing, while setting
!   MAXHIST = MAXFUN ensures that  XHIST/FHIST will output all the history.
!   If XHIST is present, its size at exit will be (N, min(NF, MAXHIST));
!   if FHIST is present, its size at exit will be min(NF, MAXHIST).
!
!   Important Notice:
!   Setting MAXHIST to a large value can be costly in terms of memory.
!   For instance, if N = 1000 and MAXHIST = 100, 000, XHIST will take
!   about 1 GB if we use double precision. Therefore, MAXHIST will be
!   reset to a smaller value if the memory needed for XHIST and/or FHIST
!   exceeds MAXMEMORY defined in CONSTS_MOD (see consts.F90 under the
!   directory named "common"; default: 2GB). Use XHIST, FHIST, and MAXHIST
!   with caution!!!
!
! INFO
!   Output, INTEGER(IK) scalar.
!   INFO is the exit flag. It will be set to one of the following values defined
!   in the module INFO_MOD (see info.F90 under the directory named "common"):
!   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
!   FTARGET_ACHIEVED: the target function value is reached;
!   TRSUBP_FAILED: a trust region step failed to reduce the quadratic model;
!   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
!   MAXTR_REACHED: the trust region iteration has been performed MAXTR times,
!       the value of MAXTR being 10*MAXFUN, which is very unlikely to reach;
!   NAN_X: NaN occurs in x;
!   NAN_INF_F: the objective function returns NaN or nearly infinite value;
!   NAN_MODEL: NaN occurs in the models.


! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, E&
     &PS
      use consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT&
     &_DFT, MAXMEMORY, MAXFUN_DIM_DFT
      use infnan_mod, only : is_nan, is_inf, is_finite
      use memory_mod, only : safealloc, cstyle_sizeof

! Solver-specific modules
      use pintrf_mod, only : FUNEVAL
      use preproc_mod, only : preproc
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

! Local variables
      integer :: maximal_hist
      integer(IK) :: info_c
      integer(IK) :: iprint_c
      integer(IK) :: maxfun_c
      integer(IK) :: maxfhist
      integer(IK) :: maxhist_c
      integer(IK) :: maxhist_in
      integer(IK) :: maxxhist
      integer(IK) :: n
      integer(IK) :: nf_c
      integer(IK) :: npt_c
      real(RP) :: eta1_c
      real(RP) :: eta2_c
      real(RP), allocatable :: fhist_c(:)
      real(RP) :: ftarget_c
      real(RP) :: gamma1_c
      real(RP) :: gamma2_c
      real(RP) :: rhobeg_c
      real(RP) :: rhoend_c
      real(RP), allocatable :: xhist_c(:, :)
      character(len=6), parameter :: solver = 'NEWUOA'


! Get size.
      n = int(size(x), kind(n))

! Replace any NaN or Inf in X by ZERO.
      where (is_nan(x) .or. is_inf(x))
          x = ZERO
      end where

! Read the inputs.

! If RHOBEG is present, then RHOBEG_C is a copy of RHOBEG (_C for "copy");
! otherwise, RHOBEG_C takes the default value for RHOBEG, taking the value
! of RHOEND into account. Note that RHOEND is considered only if it is
! present and it is VALID (i.e., finite and positive). The other inputs
! are read in a similar way.
      if (present(rhobeg)) then
          rhobeg_c = rhobeg
      else if (present(rhoend)) then
! Fortran does not take short-circuit evalation of logic expressions.
! Thus it is WRONG to combine the evalation of PRESENT(RHOEND) and the
! evalation of IS_FINITE(RHOEND) as
! "if (present(rhoend) .and. is_finite(rhoend))".
! The compiler may choose the evaluate the IS_FINITE(RHOEND) even if
! PRESENT(RHOEND) is false!
          if (is_finite(rhoend) .and. rhoend > ZERO) then
              rhobeg_c = max(TEN * rhoend, RHOBEG_DFT)
          else
              rhobeg_c = RHOBEG_DFT
          end if
      else
          rhobeg_c = RHOBEG_DFT
      end if

      if (present(rhoend)) then
          rhoend_c = rhoend
      else if (rhobeg_c > 0) then
          rhoend_c = max(EPS, min(TENTH * rhobeg_c, RHOEND_DFT))
      else
          rhoend_c = RHOEND_DFT
      end if

      if (present(ftarget)) then
          ftarget_c = ftarget
      else
          ftarget_c = FTARGET_DFT
      end if

      if (present(maxfun)) then
          maxfun_c = maxfun
      else
          maxfun_c = MAXFUN_DIM_DFT * n
      end if

      if (present(npt)) then
          npt_c = npt
      else if (maxfun_c >= 1) then
          npt_c = int(max(n + 2, min(maxfun_c - 1, 2 * n + 1)), kind(npt&
     &_c))
      else
          npt_c = int(2 * n + 1, kind(npt_c))
      end if

      if (present(iprint)) then
          iprint_c = iprint
      else
          iprint_c = IPRINT_DFT
      end if

      if (present(eta1)) then
          eta1_c = eta1
      else if (present(eta2)) then
          if (eta2 > ZERO .and. eta2 < ONE) then
              eta1_c = max(EPS, eta2 / 7.0_RP)
          end if
      else
          eta1_c = TENTH
      end if

      if (present(eta2)) then
          eta2_c = eta2
      else if (eta1_c > ZERO .and. eta1_c < ONE) then
          eta2_c = (eta1_c + TWO) / 3.0_RP
      else
          eta2_c = 0.7_RP
      end if

      if (present(gamma1)) then
          gamma1_c = gamma1
      else
          gamma1_c = HALF
      end if

      if (present(gamma2)) then
          gamma2_c = gamma2
      else
          gamma2_c = TWO
      end if

      maxhist_in = 0
      ! MAXHIST input by user
      if (present(maxhist)) then
          maxhist_c = maxhist
          maxhist_in = maxhist
      else if (maxfun_c >= n + 3) then
          maxhist_c = maxfun_c
      else
          maxhist_c = MAXFUN_DIM_DFT * n
      end if

! Preprocess the inputs in case some of them are invalid.
      call preproc(n, iprint_c, maxfun_c, maxhist_c, npt_c, eta1_c, eta2&
     &_c, ftarget_c, gamma1_c, gamma2_c, rhobeg_c, rhoend_c)

! Further revise MAXHIST according to MAXMEMORY, i.e., the maximal amount
! of memory allowed for the history.
      if (present(xhist)) then
          maximal_hist = int(MAXMEMORY / ((n + 1) * cstyle_sizeof(0.0_RP&
     &)), kind(maximal_hist))
      else
          maximal_hist = int(MAXMEMORY / (cstyle_sizeof(0.0_RP)), kind(m&
     &aximal_hist))
      end if
      if (maxhist_c > maximal_hist) then
! We cannot simply take MAXHIST_C = MIN(MAXHIST_C, MAXIMAL_HIST)
! becaue they may not be the same kind, and compilers may complain.
! We may convert them to the same kind, but overflow may occur.
          maxhist_c = int(maximal_hist, kind(maxhist_c))
      end if

! Allocate memory for the history of X. We use XHIST_C instead of XHIST,
! which may not be present.
      if (present(xhist)) then
          maxxhist = min(maxhist_c, maxfun_c)
      else
          maxxhist = 0
      end if
      call safealloc(xhist_c, n, maxxhist)

! Allocate memory for the history of F. We use FHIST_C instead of FHIST,
! which may not be present.
      if (present(fhist)) then
          maxfhist = min(maxhist_c, maxfun_c)
      else
          maxfhist = 0
      end if
      call safealloc(fhist_c, maxfhist)

!-------------------- Call NEWUOB, which performs the real calculations. --------------------------!
      call newuob(calfun, iprint_c, maxfun_c, npt_c, eta1_c, eta2_c, fta&
     &rget_c, gamma1_c, gamma2_c, rhobeg_c, rhoend_c, x, nf_c, f, fhist_&
     &c, xhist_c, info_c)
!--------------------------------------------------------------------------------------------------!

! Write the outputs.

      if (present(nf)) then
          nf = nf_c
      end if

      if (present(info)) then
          info = info_c
      end if

! Copy XHIST_C to XHIST if needed.
      if (present(xhist)) then
! The SAFEALLOC line is removable in F2003.
          call safealloc(xhist, n, min(nf_c, maxxhist))
          xhist = xhist_c(:, 1:min(nf_c, maxxhist))
! N.B.:
! 0. Allocate XHIST as long as it is present, even if MAXXHIST = 0;
! otherwise, it will be illeagle to enquire XHIST after exit.
! 1. Even though Fortran 2003 supports automatic (re)allocation of
! allocatable arrays upon intrinsic assignment, we keep the line of
! SAFEALLOC, because some very new compilers (Absoft Fortran 20.0)
! are still not standard-compliant in this respect.
! 2. NF may not be present. Hence we should NOT use NF but NF_C.
! 3. When MAXXHIST > NF_C, which is the normal case in practice,
! XHIST_C contains GARBAGE in XHIST_C(:, NF_C + 1 : MAXXHIST).
! Therefore, we MUST cap XHIST at min(NF_C, MAXXHIST) so that XHIST
! cointains only valid history. For this reason, there is no way to
! avoid allocating two copies of memory for XHIST unless we declare
! it to be a POINTER instead of ALLOCATABLE.
      end if
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet
! we prefer to deallocate them immediately when they finish their jobs.
      deallocate (xhist_c)

! Copy FHIST_C to FHIST if needed.
      if (present(fhist)) then
! The SAFEALLOC line is removable in F2003.
          call safealloc(fhist, min(nf_c, maxfhist))
          fhist = fhist_c(1:min(nf_c, maxfhist))
! The same as XHIST, we must cap FHIST at min(NF_C, MAXFHIST).
      end if
      deallocate (fhist_c)

! If MAXFHIST_IN >= NF_C > MAXFHIST_C, warn that not all history is recorded.
      if ((present(xhist) .or. present(fhist)) .and. maxhist_c < min(nf_&
     &c, maxhist_in)) then
          print '(/1A, I7, 1A)', 'WARNING: '//solver//': due to memory l&
     &imit, MAXHIST is reset to ', maxhist_c, '.'
          print '(1A/)', 'Only the history of the last MAXHIST iteration&
     &s is recoreded.'
      end if

      end subroutine newuoa


      end module newuoa_mod