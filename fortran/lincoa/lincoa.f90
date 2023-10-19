module lincoa_mod
!--------------------------------------------------------------------------------------------------!
! LINCOA_MOD is a module providing the reference implementation of Powell's LINCOA algorithm.
!
! The algorithm approximately solves
!
!   min F(X) subject to Aineq*X <= Bineq, Aeq*x = Beq, XL <= X <= XU,
!
! where X is a vector of variables that has N components, F is a real-valued objective function,
! Aineq is an Mineq-by-N matrix, Bineq is an Mineq-dimensional real vector, Aeq is an Meq-by-N
! matrix, Beq is an Meq-dimensional real vector, XL is an N-dimensional real vector, and XU is
! an N-dimensional real vector.
!
! It tackles the problem by a trust region method that forms quadratic models by interpolation.
! Usually there is much freedom in each new model after satisfying the interpolation conditions,
! which is taken up by minimizing the Frobenius norm of the change to the second derivative matrix
! of the model. One new function value is calculated on each iteration, usually at a point where
! the current model predicts a reduction in the least value so far of the objective function subject
! to the linear constraints. Alternatively, a new vector of variables may be chosen to replace an
! interpolation point that may be too far away for reliability, and the new point does not have to
! satisfy the constraints.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! and Powell's code, with modernization, bug fixes, and improvements.
!
! N.B.:
! 1. Powell did not publish a paper to introduce the algorithm. The above paper does not describe
! LINCOA but discusses how to solve linearly-constrained trust-region subproblems.
! 2. Powell's code does not accept linear equality constraints or bound constraints.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, October 19, 2023 AM10:34:35
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lincoa


contains


subroutine lincoa(calfun, x, &
    & f, cstrv, &
    & Aineq, bineq, &
    & Aeq, beq, &
    & xl, xu, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, chist, maxhist, maxfilt, info)
!--------------------------------------------------------------------------------------------------!
! Among all the arguments, only CALFUN, and X are obligatory. The others are OPTIONAL and you can
! neglect them unless you are familiar with the algorithm. Any unspecified optional input will take
! the default value detailed below. For instance, we may invoke the solver as follows.
!
! ! First define CALFUN and X, and then do the following.
! call lincoa(calfun, x, f)
!
! or
!
! ! First define CALFUN, X, Aineq, and Bineq, and then do the following.
! call lincoa(calfun, x, f, cstrv, Aineq = Aineq, bineq = bineq, rhobeg = 0.5D0, rhoend = 1.0D-3, maxfun = 100)
!
! See examples/lincoa_exmp.f90 for a concrete example.
!
! A detailed introduction to the arguments is as follows.
! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory named
! "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
! real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
! setting PRIMA_REAL_PRECISION and PRIMA_INTEGER_KIND in common/ppf.h. Use the default if unsure.
!
! CALFUN
!   Input, subroutine.
!   CALFUN(X, F) should evaluate the objective function at the given REAL(RP) vector X and set the
!   value to the REAL(RP) scalar F. It must be provided by the user, and its definition must conform
!   to the following interface:
!   !-------------------------------------------------------------------------!
!    subroutine calfun(x, f)
!    real(RP), intent(in) :: x(:)
!    real(RP), intent(out) :: f
!    end subroutine calfun
!   !-------------------------------------------------------------------------!
!
! X
!   Input and output, REAL(RP) vector.
!   As an input, X should be an N dimensional vector that contains the starting point, N being the
!   dimension of the problem. As an output, X will be set to an approximate minimizer.
!
! F
!   Output, REAL(RP) scalar.
!   F will be set to the objective function value of X at exit.
!
! CSTRV
!   Output, REAL(RP) scalar.
!   CSTRV will be set to the L-infinity constraint violation of X at exit, namely
!   MAXVAL([0, Aineq*X - Bineq, abs(Aeq*X - Beq), XL - X, X - XU])
!   N.B.: We use the original constraints to evaluate CSTRV, even though they may be modified during
!   the computation.
!
! Aineq, Bineq
!   Input, REAL(RP) matrix of size [Mineq, N] and REAL vector of size Mineq unless they are both
!   empty, default: [] and [].
!   Aineq and Bineq represent the linear inequality constraints: Aineq*X <= Bineq.
!
! Aeq, Beq
!   Input, REAL(RP) matrix of size [Meq, N] and REAL vector of size Meq unless they are both
!   empty, default: [] and [].
!   Aeq and Beq represent the linear equality constraints: Aeq*X = Beq.
!
! XL, XU
!   Input, REAL(RP) vectors of size N unless they are both empty, default: [] and [].
!   XL and XU represent the lower and upper bounds of the variables: XL <= X <= XU.
!
! NF
!   Output, INTEGER(IK) scalar.
!   NF will be set to the number of calls of CALFUN at exit.
!
! RHOBEG, RHOEND
!   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
!   the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
!   Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
!   RHOEND should indicate the accuracy that is required in the final values of the variables.
!
! FTARGET
!   Input, REAL(RP) scalar, default: -Inf.
!   FTARGET is the target function value. The algorithm will terminate when a point with a function
!   value <= FTARGET is found.
!
! CTOL
!   Input, REAL(RP) scalar, default: machine epsilon.
!   CTOL is the tolerance of constraint violation. X is considered feasible if CSTRV(X) <= CTOL.
!   N.B.: 1. CTOL is absolute, not relative.
!   2. CTOL is used for choosing the returned X. It does not affect the iterations of the algorithm.
!
! CWEIGHT
!   Input, REAL(RP) scalar, default: CWEIGHT_DFT defined in the module CONSTS_MOD in common/consts.F90.
!   CWEIGHT is the weight that the constraint violation takes in the selection of the returned X.
!
! MAXFUN
!   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
!   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of calls of CALFUN.
!
! NPT
!   Input, INTEGER(IK) scalar, default: 2N + 1.
!   NPT is the number of interpolation conditions for each trust region model. Its value must be in
!   the interval [N+2, (N+1)(N+2)/2]. Typical choices of Powell were NPT=N+6 and NPT=2*N+1. Powell
!   commented that "larger values tend to be highly inefficient when the number of variables is
!   substantial, due to the amount of work and extra difficulty of adjusting more points."
!
! IPRINT
!   Input, INTEGER(IK) scalar, default: 0.
!   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls how much
!   information will be printed during the computation:
!   0: there will be no printing;
!   1: a message will be printed to the screen at the return, showing the best vector of variables
!      found and its objective function value;
!   2: in addition to 1, each new value of RHO is printed to the screen, with the best vector of
!      variables so far and its objective function value;
!   3: in addition to 2, each function evaluation with its variables will be printed to the screen;
!   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
!      named LINCOA_output.txt; the file will be created if it does not exist; the new output will
!      be appended to the end of this file if it already exists.
!   Note that IPRINT = +/-3 can be costly in terms of time and/or space.
!
! ETA1, ETA2, GAMMA1, GAMMA2
!   Input, REAL(RP) scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5, and GAMMA2 = 2.
!   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme of the trust-region radius
!   detailed in the subroutine TRRAD in trustregion.f90. Roughly speaking, the trust-region radius
!   is contracted by a factor of GAMMA1 when the reduction ratio is below ETA1, and enlarged by a
!   factor of GAMMA2 when the reduction ratio is above ETA2. It is required that 0 < ETA1 <= ETA2
!   < 1 and 0 < GAMMA1 < 1 < GAMMA2. Normally, ETA1 <= 0.25. It is NOT advised to set ETA1 >= 0.5.
!
! XHIST, FHIST, CHIST, MAXHIST
!   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   CHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
!   XHIST, if present, will output the history of iterates, while FHIST/CHIST, if present, will output
!   the history function values/constraint violations. MAXHIST should be a nonnegative integer, and
!   XHIST/FHIST/CHIST will output only the history of the last MAXHIST iterations. Therefore,
!   MAXHIST = 0 means XHIST/FHIST/CHIST will output nothing, while setting MAXHIST = MAXFUN requests
!   XHIST/FHIST/CHIST to output all the history.
!   If XHIST is present, its size at exit will be [N, min(NF, MAXHIST)]; if FHIST/CHIST is present,
!   its size at exit will be min(NF, MAXHIST).
!
!   IMPORTANT NOTICE:
!   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
!   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
!   CONSTS_MOD (see consts.F90 under the directory named "common").
!   Use *HIST with caution! (N.B.: the algorithm is NOT designed for large problems).
!
! INFO
!   Output, INTEGER(IK) scalar.
!   INFO is the exit flag. It will be set to one of the following values defined in the module
!   INFOS_MOD (see common/infos.f90):
!   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
!   FTARGET_ACHIEVED: the target function value is reached;
!   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
!   MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
!   NAN_INF_MODEL: NaN or Inf occurs in the model;
!   NAN_INF_X: NaN or Inf occurs in X;
!   DAMAGING_ROUNDING: the rounding error becomes damaging;
!   ZERO_LINEAR_CONSTRAINT: one of the linear constraints has a zero gradient
!   !--------------------------------------------------------------------------!
!   The following case(s) should NEVER occur unless there is a bug.
!   NAN_INF_F: the objective function returns NaN or +Inf;
!   TRSUBP_FAILED: a trust region step failed to reduce the model.
!   !--------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, TEN, TENTH, EPS, REALMAX
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: evaluate_mod, only : moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: preproc_mod, only : preproc
use, non_intrinsic :: selectx_mod, only : isbetter
use, non_intrinsic :: string_mod, only : num2str

! Solver-specific modules
use, non_intrinsic :: lincob_mod, only : lincob

implicit none

! Compulsory arguments
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(inout) :: x(:)  ! X(N)

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfilt
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
integer(IK), intent(in), optional :: npt
real(RP), intent(in), optional :: Aeq(:, :)  ! Aeq(Meq, N)
real(RP), intent(in), optional :: Aineq(:, :)  ! Aineq(Mineq, N)
real(RP), intent(in), optional :: beq(:)  ! Beq(Meq)
real(RP), intent(in), optional :: bineq(:)  ! Bineq(Mineq)
real(RP), intent(in), optional :: ctol
real(RP), intent(in), optional :: cweight
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend
real(RP), intent(in), optional :: xl(:)  ! XL(N)
real(RP), intent(in), optional :: xu(:)  ! XU(N)

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), optional :: cstrv
real(RP), intent(out), optional :: f
real(RP), intent(out), optional, allocatable :: chist(:)  ! CHIST(MAXCHIST)
real(RP), intent(out), optional, allocatable :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out), optional, allocatable :: xhist(:, :)  ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'LINCOA'
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfilt_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: meq
integer(IK) :: mineq
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: npt_loc
real(RP) :: cstrv_loc
real(RP) :: ctol_loc
real(RP) :: cweight_loc
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: f_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP), allocatable :: Aeq_loc(:, :)  ! Aeq_LOC(Meq, N)
real(RP), allocatable :: Aineq_loc(:, :)  ! Aineq_LOC(Mineq, N)
real(RP), allocatable :: amat(:, :)  ! AMAT(N, M); each column corresponds to a constraint
real(RP), allocatable :: beq_loc(:)  ! Beq_LOC(Meq)
real(RP), allocatable :: bineq_loc(:)  ! Bineq_LOC(Mineq)
real(RP), allocatable :: bvec(:)  ! BVEC(M)
real(RP), allocatable :: chist_loc(:)  ! CHIST_LOC(MAXCHIST)
real(RP), allocatable :: fhist_loc(:)  ! FHIST_LOC(MAXFHIST)
real(RP), allocatable :: xhist_loc(:, :)  ! XHIST_LOC(N, MAXXHIST)
real(RP), allocatable :: xl_loc(:)  ! XL_LOC(N)
real(RP), allocatable :: xu_loc(:)  ! XU_LOC(N)

! Sizes
if (present(bineq)) then
    mineq = int(size(bineq), kind(mineq))
else
    mineq = 0
end if
if (present(beq)) then
    meq = int(size(beq), kind(meq))
else
    meq = 0
end if
n = int(size(x), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(mineq >= 0, 'Mineq >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(present(Aineq) .eqv. present(bineq), 'Aineq and Bineq are both present or both absent', srname)
    if (present(Aineq)) then
        call assert((size(Aineq, 1) == mineq .and. size(Aineq, 2) == n) &
            & .or. (size(Aineq, 1) == 0 .and. size(Aineq, 2) == 0 .and. mineq == 0), &
            & 'SIZE(Aineq) == [Mineq, N] unless Aineq and Bineq are both empty', srname)
    end if
    call assert(present(Aeq) .eqv. present(beq), 'Aeq and Beq are both present or both absent', srname)
    if (present(Aeq)) then
        call assert((size(Aeq, 1) == meq .and. size(Aeq, 2) == n) &
            & .or. (size(Aeq, 1) == 0 .and. size(Aeq, 2) == 0 .and. meq == 0), &
            & 'SIZE(Aeq) == [Meq, N] unless Aeq and Beq are both empty', srname)
    end if
    if (present(xl)) then
        call assert(size(xl) == n .or. size(xl) == 0, 'SIZE(XL) == N unless XL is empty', srname)
    end if
    if (present(xu)) then
        call assert(size(xu) == n .or. size(xu) == 0, 'SIZE(XU) == N unless XU is empty', srname)
    end if
end if

! Read the inputs

x = moderatex(x)

call safealloc(Aineq_loc, mineq, n)  ! NOT removable even in F2003, as Aineq may be absent or of size 0-by-0.
if (present(Aineq) .and. mineq > 0) then
    ! We must check Mineq > 0. Otherwise, the size of Aineq_LOC may be changed to 0-by-0 due to
    ! automatic (re)allocation if that is the size of Aineq; we allow Aineq to be 0-by-0, but
    ! Aineq_LOC should be n-by-0.
    Aineq_loc = Aineq
end if

call safealloc(bineq_loc, mineq)  ! NOT removable even in F2003, as Bineq may be absent.
if (present(bineq)) then
    bineq_loc = bineq
end if

call safealloc(Aeq_loc, meq, n)  ! NOT removable even in F2003, as Aeq may be absent or of size 0-by-0.
if (present(Aeq) .and. meq > 0) then
    ! We must check Meq > 0. Otherwise, the size of Aeq_LOC may be changed to 0-by-0 due to
    ! automatic (re)allocation if that is the size of Aeq; we allow Aeq to be 0-by-0, but
    ! Aeq_LOC should be n-by-0.
    Aeq_loc = Aeq
end if

call safealloc(beq_loc, meq)  ! NOT removable even in F2003, as Beq may be absent.
if (present(beq)) then
    beq_loc = beq
end if

call safealloc(xl_loc, n)  ! NOT removable even in F2003, as XL may be absent.
if (present(xl)) then
    xl_loc = xl
else
    xl_loc = -REALMAX
end if

call safealloc(xu_loc, n)  ! NOT removable even in F2003, as XU may be absent.
if (present(xu)) then
    xu_loc = xu
else
    xu_loc = REALMAX
end if

! If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
! value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
! it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
if (present(rhobeg)) then
    rhobeg_loc = rhobeg
elseif (present(rhoend)) then
    ! Fortran does not take short-circuit evaluation of logic expressions. Thus it is WRONG to
    ! combine the evaluation of PRESENT(RHOEND) and the evaluation of IS_FINITE(RHOEND) as
    ! "IF (PRESENT(RHOEND) .AND. IS_FINITE(RHOEND))". The compiler may choose to evaluate the
    ! IS_FINITE(RHOEND) even if PRESENT(RHOEND) is false!
    if (is_finite(rhoend) .and. rhoend > 0) then
        rhobeg_loc = max(TEN * rhoend, RHOBEG_DFT)
    else
        rhobeg_loc = RHOBEG_DFT
    end if
else
    rhobeg_loc = RHOBEG_DFT
end if

if (present(rhoend)) then
    rhoend_loc = rhoend
elseif (rhobeg_loc > 0) then
    rhoend_loc = max(EPS, min(TENTH * rhobeg_loc, RHOEND_DFT))
else
    rhoend_loc = RHOEND_DFT
end if

if (present(ctol)) then
    ctol_loc = ctol
else
    ctol_loc = CTOL_DFT
end if

if (present(cweight)) then
    cweight_loc = cweight
else
    cweight_loc = CWEIGHT_DFT
end if

if (present(ftarget)) then
    ftarget_loc = ftarget
else
    ftarget_loc = FTARGET_DFT
end if

if (present(maxfun)) then
    maxfun_loc = maxfun
else
    maxfun_loc = MAXFUN_DIM_DFT * n
end if

if (present(npt)) then
    npt_loc = npt
elseif (maxfun_loc >= 1) then
    npt_loc = max(n + 2_IK, min(maxfun_loc - 1_IK, 2_IK * n + 1_IK))
else
    npt_loc = 2_IK * n + 1_IK
end if

if (present(iprint)) then
    iprint_loc = iprint
else
    iprint_loc = IPRINT_DFT
end if

if (present(eta1)) then
    eta1_loc = eta1
elseif (present(eta2)) then
    if (eta2 > 0 .and. eta2 < 1) then
        eta1_loc = max(EPS, eta2 / 7.0_RP)
    end if
else
    eta1_loc = TENTH
end if

if (present(eta2)) then
    eta2_loc = eta2
elseif (eta1_loc > 0 .and. eta1_loc < 1) then
    eta2_loc = (eta1_loc + TWO) / 3.0_RP
else
    eta2_loc = 0.7_RP
end if

if (present(gamma1)) then
    gamma1_loc = gamma1
else
    gamma1_loc = HALF
end if

if (present(gamma2)) then
    gamma2_loc = gamma2
else
    gamma2_loc = TWO
end if

if (present(maxhist)) then
    maxhist_loc = maxhist
else
    maxhist_loc = maxval([maxfun_loc, n + 3_IK, MAXFUN_DIM_DFT * n])
end if

if (present(maxfilt)) then
    maxfilt_loc = maxfilt
else
    maxfilt_loc = MAXFILT_DFT
end if

! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & npt=npt_loc, ctol=ctol_loc, cweight=cweight_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, &
    & gamma2=gamma2_loc, maxfilt=maxfilt_loc)

! Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! CHIST = NaN(1, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc, present(chist), chist_loc)

! Wrap the linear and bound constraints into a single constraint: AMAT^T*X <= BVEC.
call get_lincon(Aeq_loc, Aineq_loc, beq_loc, bineq_loc, rhoend_loc, xl_loc, xu_loc, x, amat, bvec)

!-------------------- Call LINCOB, which performs the real calculations. --------------------------!
call lincob(calfun, iprint_loc, maxfilt_loc, maxfun_loc, npt_loc, Aeq_loc, Aineq_loc, amat, &
    & beq_loc, bineq_loc, bvec, ctol_loc, cweight_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, &
    & gamma2_loc, rhobeg_loc, rhoend_loc, xl_loc, xu_loc, x, nf_loc, chist_loc, cstrv_loc, &
    & f_loc, fhist_loc, xhist_loc, info_loc)
!--------------------------------------------------------------------------------------------------!

! Deallocate variables not needed any more. Indeed, automatic allocation will take place at exit.
deallocate (Aineq_loc, Aeq_loc, amat, bineq_loc, beq_loc, bvec, xl_loc, xu_loc)


! Write the outputs.

if (present(f)) then
    f = f_loc
end if

if (present(cstrv)) then
    cstrv = cstrv_loc
end if

if (present(nf)) then
    nf = nf_loc
end if

if (present(info)) then
    info = info_loc
end if

! Copy XHIST_LOC to XHIST if needed.
if (present(xhist)) then
    nhist = min(nf_loc, int(size(xhist_loc, 2), kind(nhist)))
    !----------------------------------------------------!
    call safealloc(xhist, n, nhist)  ! Removable in F2003.
    !----------------------------------------------------!
    xhist = xhist_loc(:, 1:nhist)
    ! N.B.:
    ! 0. Allocate XHIST as long as it is present, even if the size is 0; otherwise, it will be
    ! illegal to enquire XHIST after exit.
    ! 1. Even though Fortran 2003 supports automatic (re)allocation of allocatable arrays upon
    ! intrinsic assignment, we keep the line of SAFEALLOC, because some very new compilers (Absoft
    ! Fortran 21.0) are still not standard-compliant in this respect.
    ! 2. NF may not be present. Hence we should NOT use NF but NF_LOC.
    ! 3. When SIZE(XHIST_LOC, 2) > NF_LOC, which is the normal case in practice, XHIST_LOC contains
    ! GARBAGE in XHIST_LOC(:, NF_LOC + 1 : END). Therefore, we MUST cap XHIST at NF_LOC so that
    ! XHIST contains only valid history. For this reason, there is no way to avoid allocating
    ! two copies of memory for XHIST unless we declare it to be a POINTER instead of ALLOCATABLE.
end if
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer to deallocate
! them immediately when they finish their jobs.
deallocate (xhist_loc)

! Copy FHIST_LOC to FHIST if needed.
if (present(fhist)) then
    nhist = min(nf_loc, int(size(fhist_loc), kind(nhist)))
    !--------------------------------------------------!
    call safealloc(fhist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    fhist = fhist_loc(1:nhist)  ! The same as XHIST, we must cap FHIST at NF_LOC.
end if
deallocate (fhist_loc)

! Copy CHIST_LOC to CHIST if needed.
if (present(chist)) then
    nhist = min(nf_loc, int(size(chist_loc), kind(nhist)))
    !--------------------------------------------------!
    call safealloc(chist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    chist = chist_loc(1:nhist)  ! The same as XHIST, we must cap CHIST at NF_LOC.
end if
deallocate (chist_loc)

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist) .or. present(chist)) .and. maxhist_loc < nf_loc) then
    call warning(solver, 'Only the history of the last '//num2str(maxhist_loc)//' iteration(s) is recorded')
end if

! Postconditions
if (DEBUGGING) then
    call assert(nf_loc <= maxfun_loc, 'NF <= MAXFUN', srname)
    call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
    nhist = min(nf_loc, maxhist_loc)
    if (present(xhist)) then
        call assert(size(xhist, 1) == n .and. size(xhist, 2) == nhist, 'SIZE(XHIST) == [N, NHIST]', srname)
        call assert(.not. any(is_nan(xhist)), 'XHIST does not contain NaN', srname)
    end if
    if (present(fhist)) then
        call assert(size(fhist) == nhist, 'SIZE(FHIST) == NHIST', srname)
        call assert(.not. any(is_nan(fhist) .or. is_posinf(fhist)), 'FHIST does not contain NaN/+Inf', srname)
    end if
    if (present(chist)) then
        call assert(size(chist) == nhist, 'SIZE(CHIST) == NHIST', srname)
        call assert(.not. any(is_nan(chist) .or. is_posinf(chist)), 'CHIST does not contain NaN/+Inf', srname)
    end if
    if (present(fhist) .and. present(chist)) then
        call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f_loc, cstrv_loc, ctol_loc)),&
            & 'No point in the history is better than X', srname)
    end if
end if

end subroutine lincoa


subroutine get_lincon(Aeq, Aineq, beq, bineq, rhoend, xl, xu, x0, amat, bvec)
!--------------------------------------------------------------------------------------------------!
! This subroutine wraps the linear and bound constraints into a single constraint: AMAT^T*X <= BVEC.
! N.B.:
! 1. LINCOA modifies the right hand sides of the constraints to make the starting point feasible if
! it is not. This is not ideal, but Powell's code was implemented in this way. In the
! MATLAB/Python/Julia/R code, we should include a preprocessing subroutine to project the starting
! point to the feasible region if it is infeasible, so that the modification will not occur.
! 2. The linear inequality constraints received by LINCOB is AMAT^T * X <= BVEC. Note that Each
! column of AMAT corresponds to a constraint. This is different from Aineq and Aeq, whose rows
! correspond to constraints. AMAT is defined in this way because it is accessed in columns during
! the computation, and because Fortran saves arrays in the column-major order. In Python/C
! implementations, AMAT should be transposed.
! 3. LINCOA normalizes the linear constraints so that each constraint has a gradient of norm 1. This
! is essential for LINCOA.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, EPS, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: linalg_mod, only : matprod, eye, trueloc
use, non_intrinsic :: memory_mod, only : safealloc

implicit none

! Inputs
real(RP), intent(in) :: Aeq(:, :)
real(RP), intent(in) :: Aineq(:, :)
real(RP), intent(in) :: beq(:)
real(RP), intent(in) :: bineq(:)
real(RP), intent(in) :: rhoend
real(RP), intent(in) :: xl(:)
real(RP), intent(in) :: xu(:)
real(RP), intent(in) :: x0(:)

! Outputs
real(RP), intent(out), allocatable :: amat(:, :)
real(RP), intent(out), allocatable :: bvec(:)

! Local variables
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'GET_LINCON'
integer(IK) :: m
integer(IK) :: meq
integer(IK) :: mineq
integer(IK) :: mxl
integer(IK) :: mxu
integer(IK) :: n
integer(IK), allocatable :: ieq(:)
integer(IK), allocatable :: iineq(:)
integer(IK), allocatable :: ixl(:)
integer(IK), allocatable :: ixu(:)
logical :: constr_modified
real(RP) :: Aeq_norm(size(Aeq, 1))
real(RP) :: Aeqx0(size(Aeq, 1))
real(RP) :: Aineq_norm(size(Aineq, 1))
real(RP) :: Aineqx0(size(Aineq, 1))
real(RP) :: idmat(size(x0), size(x0))
real(RP) :: smallx
real(RP), allocatable :: Anorm(:)

! Sizes
n = int(size(x0), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(Aineq, 1) == size(bineq) .and. size(Aineq, 2) == n, 'SIZE(AINEQ) == [SIZE(BINEQ), N]', srname)
    call assert(size(Aeq, 1) == size(beq) .and. size(Aeq, 2) == n, 'SIZE(AEQ) == [SIZE(BEQ), N]', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == SIZE(XU) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

! Decide the number of nontrivial and valid (gradient is nonzero) constraints.
mxl = int(count(xl > -REALMAX), kind(mxl))
mxu = int(count(xu < REALMAX), kind(mxu))
Aeq_norm = sqrt(sum(Aeq**2, dim=2))
meq = int(count(Aeq_norm > 0), kind(meq))
Aineq_norm = sqrt(sum(Aineq**2, dim=2))
mineq = int(count(Aineq_norm > 0), kind(mineq))
m = mxl + mxu + 2_IK * meq + mineq  ! The final number of linear inequality constraints.

! Print a warning if some constraints are invalid. They will be ignored (Powell's code would stop).
if (meq < size(Aeq, 1) .or. mineq < size(Aineq, 1)) then
    call warning(solver, 'Some linear constraints have zero gradients; they are ignored')
end if

! Allocate memory. Removable in F2003.
call safealloc(ixl, mxl)
call safealloc(ixu, mxu)
call safealloc(ieq, meq)
call safealloc(iineq, mineq)
call safealloc(amat, n, m)
call safealloc(bvec, m)
call safealloc(Anorm, 2_IK * meq + mineq)

! Define the indices of the valid and nontrivial constraints.
ixl = trueloc(xl > -REALMAX)
ixu = trueloc(xu < REALMAX)
ieq = trueloc(Aeq_norm > 0)
iineq = trueloc(Aineq_norm > 0)

! Wrap the linear constraints.
! The bound constraint XL <= X <= XU is handled as two constraints -X <= -XL, X <= XU.
! The equality constraint Aeq*X = Beq is handled as two constraints -Aeq*X <= -Beq, Aeq*X <= Beq.
! N.B.:
! 1. The treatment of the equality constraints is naive. One may choose to eliminate them instead.
! 2. The code below is quite inefficient in terms of memory, but we prefer readability.
idmat = eye(n, n)
amat = reshape(shape=shape(amat), source= &
    & [-idmat(:, ixl), idmat(:, ixu), -transpose(Aeq(ieq, :)), transpose(Aeq(ieq, :)), transpose(Aineq(iineq, :))])
bvec = [-xl(ixl), xu(ixu), -beq(ieq), beq(ieq), bineq(iineq)]
!!MATLAB code:
!!amat = [-idmat(:, ixl), idmat(:, ixu), -Aeq(ieq, :)', Aeq(ieq, :)', Aineq(iineq, :)'];
!!bvec = [-xl(ixl); xu(ixu); -beq(ieq); beq(ieq); bineq(iineq)];

! Modify BVEC if necessary so that the initial point is feasible.
Aeqx0 = matprod(Aeq, x0)
Aineqx0 = matprod(Aineq, x0)
bvec = max(bvec, [-x0(ixl), x0(ixu), -Aeqx0(ieq), Aeqx0(ieq), Aineqx0(iineq)])

! Normalize the linear constraints so that each constraint has a gradient of norm 1.
Anorm = [Aeq_norm(ieq), Aeq_norm(ieq), Aineq_norm(iineq)]
amat(:, mxl + mxu + 1:m) = amat(:, mxl + mxu + 1:m) / spread(Anorm, dim=1, ncopies=n)
bvec(mxl + mxu + 1:m) = bvec(mxl + mxu + 1:m) / Anorm

! Deallocate memory.
deallocate (ixl, ixu, ieq, iineq, Anorm)

! Print a warning if the starting point is sufficiently infeasible and the constraints are modified.
smallx = 1.0E-6_RP * rhoend
constr_modified = (any(x0 + smallx < xl) .or. any(x0 - smallx > xu) .or. &
    & any(abs(Aeqx0 - beq) > smallx * Aeq_norm) .or. any(Aineqx0 - bineq > smallx * Aineq_norm))
if (constr_modified) then
    call warning(solver, 'The starting point is infeasible. '//solver// &
        & ' modified the right-hand sides of the constraints to make it feasible')
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(amat, 1) == size(x0) .and. size(amat, 2) == size(bvec), &
        & 'SIZE(AMAT) == [SIZE(X), SIZE(BVEC)]', srname)
    call assert(all(matprod(x0, amat) - bvec <= max(1.0E-12_RP, 1.0E2 * EPS) * &
        & (ONE + sum(abs(x0)) + sum(abs(bvec)))), 'The starting point is feasible', srname)
end if
end subroutine get_lincon


end module lincoa_mod
