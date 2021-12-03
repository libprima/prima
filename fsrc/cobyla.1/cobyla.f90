module cobyla_mod

implicit none
private
public :: cobyla


contains

subroutine cobyla(calcfc, n, m, x, rhobeg, rhoend, iprint, maxfun, f, info, ftarget, cstrv, constr)

! Generic modules
use, non_intrinsic :: pintrf_mod, only : FUNCON
use, non_intrinsic :: consts_mod, only : RP, IK, EPS, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: output_mod, only : retmssg, rhomssg, fmssg
use, non_intrinsic :: memory_mod, only : safealloc

! Solver-specific modules
use, non_intrinsic :: cobylb_mod, only : cobylb

implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: m
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: n
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! Parameters
real(RP) :: ctol

! In-outputs
procedure(FUNCON) :: calcfc
real(RP), intent(inout) :: x(:)

! Outputs
integer(IK), intent(out) :: info
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: cstrv
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
integer(IK) :: maxxhist
integer(IK) :: maxfhist
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: nf_c
real(RP), allocatable :: xhist(:, :)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: conhist(:, :)
real(RP), allocatable :: chist(:)
real(RP), allocatable :: xhist_c(:, :)
real(RP), allocatable :: fhist_c(:)
real(RP), allocatable :: conhist_c(:, :)
real(RP), allocatable :: chist_c(:)

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
integer(IK) :: mpp
real(RP) :: rhoend_c
character(len=*), parameter :: srname = 'COBYLA'
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
!                SUBROUTINE CALCFC (N,M,X,F,CONSTR)
!                DIMENSION X(*),CONSTR(*)  .
!     The values of N and M are fixed and have been defined already, while
!     X is now the current vector of variables. The subroutine should return
!     the objective and constraint functions at X in F and CONSTR(1),CONSTR(2),
!     ...,CONSTR(M). Note that we are trying to adjust X so that F(X) is as
!     small as possible subject to the constraint functions being nonnegative.
!
!     Partition the working space array W to provide the storage that is needed
!     for the main calculation.
!
mpp = m + 2
icon = 1
isim = icon + mpp
isimi = isim + n * n + n
idatm = isimi + n * n
ia = idatm + n * mpp + mpp
ivsig = ia + m * n + n
iveta = ivsig + n
isigb = iveta + n
idx = isigb + n
iw = idx + n
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! and RHOEND may change a bit (due to rounding ???). It was observed in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOEND = RHOBEG in the interfaces, then it may happen
! that RHOEND > RHOBEG. That is why we do the following.
rhoend_c = min(rhobeg, rhoend)
! CTOL is the tolerance for constraint violation. A point X is considered to be feasible if its
! constraint violation (CSTRV) is less than CTOL.
ctol = EPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
maxxhist = maxfun
maxfhist = maxfun
maxchist = maxfun
maxconhist = maxfun
call safealloc(xhist_c, n, maxxhist)
call safealloc(fhist_c, maxfhist)
call safealloc(conhist_c, m, maxconhist)
call safealloc(chist_c, maxchist)
call cobylb(calcfc, iprint, maxfun, ctol, ftarget, rhobeg, rhoend, constr, x, nf_c, chist_c, &
    & conhist_c, cstrv, f, fhist_c, xhist_c, info)
call safealloc(xhist, n, min(nf_c, maxxhist))
xhist = xhist_c(:, 1:min(nf_c, maxxhist))
deallocate (xhist_c)
call safealloc(fhist, min(nf_c, maxfhist))
fhist = fhist_c(1:min(nf_c, maxfhist))
deallocate (fhist_c)
call safealloc(conhist, m, min(nf_c, maxconhist))
conhist = conhist_c(:, 1:min(nf_c, maxconhist))
deallocate (conhist_c)
call safealloc(chist, min(nf_c, maxchist))
chist = chist_c(1:min(nf_c, maxchist))
deallocate (chist_c)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine cobyla

end module cobyla_mod
