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

      use consts, only : RP, IK
      use newuob_mod, only : newuob
      implicit none
      
      integer(IK), intent(in) :: npt, iprint, maxfun
      integer(IK), intent(out) :: info
      real(RP), intent(in) :: rhobeg, rhoend, ftarget
      real(RP), intent(out) :: f
      real(RP), intent(inout) :: x(:)

      integer(IK) :: n, nf

      n = size(x)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun, 2020-05-05
      ! When the data is passed from the interfaces to the Fortran code,
      ! RHOBEG, and RHOEND may change a bit (due to rounding ???). It
      ! was oberved in a MATLAB test that MEX passed 1 to Fortran as
      ! 0.99999999999999978. Therefore, if we set RHOEND = RHOBEG in the
      ! interfaces, then it may happen that RHOEND > RHOBEG. That is why
      ! we do the following.
      !rhoend = min(rhobeg, rhoend)

      ! Other checkings should be done as well !!!!
      if (npt < n + 2 .or. npt > ((n + 2)*(n + 1))/2) then
          call calfun(n, x, f)
          info=5
          return
      end if
      
      call newuob(npt, rhobeg, min(rhobeg, rhoend), iprint, maxfun,     &
     & ftarget, x, f, nf, info)
      return
      end
