      module newuoa_mod

      implicit none
      private
      public :: newuoa


      contains

      subroutine newuoa (npt, x, rhobeg, rhoend, iprint, maxfun, f,     &
     & info, ftarget)
      ! NEWUOA seeks the least value of a function of many variables,
      ! by a trust region method that forms quadratic models by
      ! interpolation. There can be some freedom in the interpolation
      ! conditions, which is taken up by minimizing the Frobenius norm
      ! of the change to the second derivative of the quadratic model,
      ! beginning with a zero matrix. The arguments of the subroutine
      ! are as follows.

      ! N must be set to the number of variables.
      ! NPT is the number of interpolation conditions. Its value must be
      ! in the interval [N+2,(N+1)(N+2)/2].
      ! Initial values of the variables must be set in X(1 : N). They
      ! will be changed to the values that give the least calculated F.
      ! RHOBEG and RHOEND must be set to the initial and final values of
      ! a trust region radius, so both must be positive with
      ! RHOEND<=RHOBEG. Typically RHOBEG should be about one tenth of
      ! the greatest expected change to a variable, and RHOEND should
      ! indicate the accuracy that is required in the final values of
      ! the variables.
      ! The value of IPRINT should be set to 0, 1, 2 or 3, which
      ! controls the amount of printing. Specifically, there is no
      ! output if IPRINT = 0 and there is output only at the return if
      ! IPRINT=1. Otherwise, each new value of RHO is printed, with the
      ! best vector of variables so far and the corresponding value of
      ! the objective function. Further, each new value of F with its
      ! variables are output if IPRINT=3.
      ! MAXFUN must be set to the maximal number of calls of CALFUN.
      !
      ! F is the objective function value when the algorithm exit.
      ! INFO is the exit flag, which can be set to:
      ! 0: the lower bound for the trust region radius is reached.
      ! 1: the target function value is reached.
      ! 2: a trust region step has failed to reduce the quadratic model.
      ! 3: the objective function has been evaluated MAXFUN times.
      ! 4: much cancellation in a denominator.
      ! 5: NPT is not in the required interval.
      ! 6: one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
      ! 7: rounding errors are becoming damaging.
      ! 8: rounding errors prevent reasonable changes to X.
      ! 9: the denominator of the updating formule is zero.
      ! 10: N should not be less than 2.
      ! 11: MAXFUN is less than NPT+1.
      ! 12: the gradient of constraint is zero.
      ! -1: NaN occurs in x.
      ! -2: the objective function returns NaN or nearly infinite value.
      ! -3: NaN occurs in the models
      !
      ! Subroutine CALFUN (N,X,F) must be provided by the user. It must
      ! set F to the value of the objective function for the variables
      ! X(1 : N).

      use consts_mod, only : RP, IK,ZERO, ONE, TENTH, EPS
      use consts_mod, only : RHOBEG_DEF, RHOEND_DEF, FTARGET_DEF
      use consts_mod, only : MAXFUN_DIM_DEF, IPRINT_DEF
      use newuob_mod, only : newuob
      use infnan_mod, only : is_nan, is_inf
      implicit none
      
      integer(IK), intent(in) :: npt, iprint, maxfun
      integer(IK), intent(out) :: info
      real(RP), intent(in) :: rhobeg, rhoend, ftarget
      real(RP), intent(out) :: f
      real(RP), intent(inout) :: x(:)

      real(RP) :: rhobeg_v, rhoend_v, ftarget_v
      integer(IK) :: n, nf, maxfun_v, npt_v, iprint_v

      n = size(x)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun, 2020-05-05
      !rhoend = min(rhobeg, rhoend)

      rhobeg_v = rhobeg
      rhoend_v = rhoend
      maxfun_v = maxfun
      npt_v = npt
      iprint_v = iprint
      ftarget_v = ftarget

      where (is_nan(x)) 
          x = ZERO
      end where

      ! When the data is passed from the interfaces (e.g., MEX) to the 
      ! Fortran code, RHOBEG, and RHOEND may change a bit. It was 
      ! oberved in a MATLAB test that MEX passed 1 to Fortran as
      ! 0.99999999999999978. Therefore, if we set RHOEND = RHOBEG in the
      ! interfaces, then it may happen that RHOEND > RHOBEG, which is
      ! considered as an invalid input. To avoid this, we force RHOBEG
      ! and RHOEND to equal when their difference is tiny.
      if ((rhobeg_v-rhoend_v) < 1.0e2_RP*EPS*max(abs(rhobeg_v),ONE))then
          rhoend_v = rhobeg_v
      end if

      if (rhobeg_v <= 0 .or. is_nan(rhobeg_v) .or. is_inf(rhobeg_v))then
          rhobeg_v = RHOBEG_DEF
      end if
      rhobeg_v = max(EPS, rhobeg_v)

      if (rhoend_v < 0 .or.  rhobeg_v < rhoend .or. is_nan(rhoend_v)    &
     & .or. is_inf(rhoend_v)) then
          rhoend_v = min(TENTH*rhobeg_v, RHOEND_DEF)
      end if
      rhoend_v = max(EPS, rhoend_v)

      maxfun_v = max(n + 3, maxfun_v)

      if (npt_v < n + 2 .or. npt > min(maxfun_v-1,((n+2)*(n+1))/2)) then 
          npt_v = min(maxfun_v - 1, 2*n + 1)
      end if

      if (iprint_v /= 0 .and. iprint_v /= 1 .and. iprint_v /= 2         &
     & .and. iprint_v /=3) then
          iprint_v = IPRINT_DEF
      end if

      if (is_nan(ftarget_v)) then
          ftarget_v = FTARGET_DEF
      end if

      call newuob(npt_v, rhobeg_v, rhoend_v, iprint_v, maxfun_v,        &
     & ftarget_v, x, f, nf, info)
      return
      end subroutine newuoa

      end module newuoa_mod
