!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of cobyla.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 27-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module cobyla_mod

      contains
      subroutine cobyla(n, m, x, rhobeg, rhoend, iprint, maxfun, f, info&
     &, ftarget, resmax, con)
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
      use cobylb_mod, only : cobylb

      implicit none

! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: m
      integer(IK), intent(in) :: maxfun
      integer(IK), intent(in) :: n
      real(RP), intent(in) :: ftarget
      real(RP), intent(in) :: resmax
      real(RP), intent(in) :: rhobeg
      real(RP), intent(in) :: rhoend

! In-outputs
      real(RP), intent(inout) :: x(:)

! Outputs
      integer(IK), intent(out) :: info
      real(RP), intent(out) :: con(:)
      real(RP), intent(out) :: f
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      integer(IK) :: ia
      integer(IK) :: icon
      integer(IK) :: idatm
      integer(IK) :: idx
      integer(IK) :: isigb
      integer(IK) :: isim
      integer(IK) :: isimi
      integer(IK) :: iveta
      integer(IK) :: ivsig
      integer(IK) :: iw
      integer(IK) :: k
      integer(IK) :: mpp
      real(RP) :: rhoend_c
      real(RP) :: confr(m + 2)
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This subroutine minimizes an objective function F(X) subject to M
!     inequality constraints on X, where X is a vector of variables that has
!     N components. The algorithm employs linear approximations to the
!     objective and constraint functions, the approximations being formed by
!     linear interpolation at N+1 points in the space of the variables.
!     We regard these interpolation points as vertices of a simplex. The
!     parameter RHO controls the size of the simplex and it is reduced
!     automatically from RHOBEG to RHOEND. For each RHO the subroutine tries
!     to achieve a good vector of variables for the current size, and then
!     RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and
!     RHOEND should be set to reasonable initial changes to and the required
!     accuracy in the variables respectively, but this accuracy should be
!     viewed as a subject for experimentation because it is not guaranteed.
!     The subroutine has an advantage over many of its competitors, however,
!     which is that it treats each constraint individually when calculating
!     a change to the variables, instead of lumping the constraints together
!     into a single penalty function. The name of the subroutine is derived
!     from the phrase Constrained Optimization BY Linear Approximations.
!
!     The user must set the values of N, M, RHOBEG and RHOEND, and must
!     provide an initial vector of variables in X. Further, the value of
!     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
!     printing during the calculation. Specifically, there is no output if
!     IPRINT=0 and there is output only at the end of the calculation if
!     IPRINT=1. Otherwise each new value of RHO and SIGMA is printed.
!     Further, the vector of variables and some function information are
!     given either when RHO is reduced or when each new value of F(X) is
!     computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA
!     is a penalty parameter, it being assumed that a change to X is an
!     improvement if it reduces the merit function
!                F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)),
!     where C1,C2,...,CM denote the constraint functions that should become
!     nonnegative eventually, at least to the precision of RHOEND. In the
!     printed output the displayed term that is multiplied by SIGMA is
!     called MAXCV, which stands for 'MAXimum Constraint Violation'. The
!     argument MAXFUN is an integer variable that must be set by the user to a
!     limit on the number of calls of CALCFC, the purpose of this routine being
!     given below. The value of MAXFUN will be altered to the number of calls
!     of CALCFC that are made. The arguments W and IACT provide real and
!     integer arrays that are used as working space. Their lengths must be at
!     least N*(3*N+2*M+11)+4*M+6 and M+1 respectively.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     F is the objective function value when the algorithm exit.
!     INFO is the exit flag, which can be set to:
!       0: the lower bound for the trust region radius is reached.
!       1: the target function value is reached.
!       2: a trust region step has failed to reduce the quadratic model.
!       3: the objective function has been evaluated MAXFUN times.
!       4: much cancellation in a denominator.
!       5: NPT is not in the required interval.
!       6: one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
!       7: rounding errors are becoming damaging.
!       8: rounding errors prevent reasonable changes to X.
!       9: the denominator of the updating formule is zero.
!       10: N should not be less than 2.
!       11: MAXFUN is less than NPT+1.
!       12: the gradient of constraint is zero.
!       -1: NaN occurs in x.
!       -2: the objective function returns a NaN or nearly infinite
!           value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     In order to define the objective and constraint functions, we require
!     a subroutine that has the name and arguments
!                SUBROUTINE CALCFC (N,M,X,F,CON)
!                DIMENSION X(*),CON(*)  .
!     The values of N and M are fixed and have been defined already, while
!     X is now the current vector of variables. The subroutine should return
!     the objective and constraint functions at X in F and CON(1),CON(2),
!     ...,CON(M). Note that we are trying to adjust X so that F(X) is as
!     small as possible subject to the constraint functions being nonnegative.
!
!     Partition the working space array W to provide the storage that is needed
!     for the main calculation.
!
      mpp = M + 2
      icon = 1
      isim = icon + mpp
      isimi = isim + N * N + N
      idatm = isimi + N * N
      ia = idatm + N * mpp + mpp
      ivsig = ia + M * N + N
      iveta = ivsig + N
      isigb = iveta + N
      idx = isigb + N
      iw = idx + N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! and RHOEND may change a bit (due to rounding ???). It was observed in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOEND = RHOBEG in the interfaces, then it may happen
! that RHOEND > RHOBEG. That is why we do the following.
      rhoend_c = min(rhobeg, rhoend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write (10, *) 'bc'
      call cobylb(n, m, x, rhobeg, rhoend_c, iprint, maxfun, confr, f, i&
     &nfo, ftarget, resmax)
      con(1:m) = confr(1:m)
      write (10, *) 'ac'
      close (10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end subroutine COBYLA

      end module cobyla_mod