module lincoa_mod
!--------------------------------------------------------------------------------------------------!
! LINCOA_MOD is a module providing the reference implementation of Powell's LINCOA algorithm.
!
! The algorithm approximately solves
!
!   min F(X) subject to A^T*X <= B,
!
! where X is a vector of variables that has N components, F is a real-valued objective function,
! A is an N-by-M matrix, and B is an M-dimensional real vector. It tackles the problem by a trust
! region method that forms quadratic models by interpolation. Usually there is much freedom in each
! new model after satisfying the interpolation conditions, which is taken up by minimizing the
! Frobenius norm of the change to the second derivative matrix of the model. One new function value
! is calculated on each iteration, usually at a point where the current model predicts a reduction
! in the least value so far of the objective function subject to the linear constraints.
! Alternatively, a new vector of variables may be chosen to replace an interpolation point that may
! be too far away for reliability, and then the new point does not have to satisfy the constraints.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! and Powell's code, with modernization, bug fixes, and improvements.
!
! Powell did not publish a paper to introduce the algorithm. The above paper does not describe
! LINCOA but discusses how to solve linearly-constrained trust-region subproblems.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Saturday, March 04, 2023 PM01:41:00
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lincoa


contains


subroutine lincoa(calfun, x, f, &
    & cstrv, &
    & A, b, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, chist, maxhist, maxfilt, info)
!--------------------------------------------------------------------------------------------------!
! Among all the arguments, only CALFUN, X, and F are obligatory. The others are OPTIONAL and you can
! neglect them unless you are familiar with the algorithm. Any unspecified optional input will take
! the default value detailed below. For instance, we may invoke the solver as follows.
!
! ! First define CALFUN and X, and then do the following.
! call lincoa(calfun, x, f)
!
! or
!
! ! First define CALFUN, X, A, and B, and then do the following.
! call lincoa(calfun, x, f, A = A, b = b, rhobeg = 0.5D0, rhoend = 1.0D-3, maxfun = 100)
!
! See examples/lincoa_exmp.f90 for a concrete example.
!
! A detailed introduction to the arguments is as follows.
! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory named
! "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
! real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
! setting __REAL_PRECISION__ and __INTEGER_KIND__ in common/ppf.h. Use the default if unsure.
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
!   CSTRV will be set to the L-infinity constraint violation of X at exit, namely MAXVAL([0, A*X-B])
!   N.B.: Even though A and B are scaled during the computation, we use the original A and B to
!   evaluate CSTRV.
!
! A, B
!   Input, REAL(RP) matrix of size [N, M] and REAL vector of size M unless they are both empty,
!   default: [] and [].
!   A and B represent the linear inequality constraints: A^T*X <= B.
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
!   CTOL is the tolerance of constraint violation. Any X with MAXVAL(A^T*X-B) <= CTOL is
!   considered feasible.
!   N.B.: 1. CTOL is absolute, not relative. 2. CTOL is used only when selecting the returned X.
!   It does not affect the iterations of the algorithm.
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
!   commented in his code that "larger values tend to be highly inefficent when the number of
!   variables is substantial, due to the amount of work and extra difficulty of adjusting more
!   points."
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
!      be appended to the end of this file if it already exists. Note that IPRINT = -3 can be costly
!      in terms of time and space.
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
!   ZERO_LINEAR_CONSTRAINT: one of the linear constraint has a zero gradient
!   !--------------------------------------------------------------------------!
!   The following case(s) should NEVER occur unless there is a bug.
!   NAN_INF_F: the objective function returns NaN or +Inf;
!   TRSUBP_FAILED: a trust region step failed to reduce the model.
!   !--------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, TEN, TENTH, EPS, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: evaluate_mod, only : moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: infos_mod, only : ZERO_LINEAR_CONSTRAINT
use, non_intrinsic :: linalg_mod, only : inprod, norm
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: preproc_mod, only : preproc
use, non_intrinsic :: selectx_mod, only : isbetter

! Solver-specific modules
use, non_intrinsic :: lincob_mod, only : lincob

implicit none

! Compulsory arguments
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(inout) :: x(:)  ! X(N)
real(RP), intent(out) :: f

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfilt
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
integer(IK), intent(in), optional :: npt
real(RP), intent(in), optional :: A(:, :)  ! A(N, M)
real(RP), intent(in), optional :: b(:)  ! B(M)
real(RP), intent(in), optional :: ctol
real(RP), intent(in), optional :: cweight
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), allocatable, optional :: chist(:)  ! CHIST(MAXCHIST)
real(RP), intent(out), allocatable, optional :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out), allocatable, optional :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(out), optional :: cstrv

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'LINCOA'
character(len=MSGLEN) :: wmsg
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: j
integer(IK) :: m
integer(IK) :: maxfilt_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: npt_loc
logical :: constr_modified
real(RP) :: anorm
real(RP) :: ax
real(RP) :: cstrv_loc
real(RP) :: ctol_loc
real(RP) :: cweight_loc
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP) :: smallx
real(RP), allocatable :: A_loc(:, :)  ! A_LOC(N, M)
real(RP), allocatable :: A_normalized(:, :)  ! A_NORMALIZED(N, M)
real(RP), allocatable :: b_loc(:)  ! B_LOC(M)
real(RP), allocatable :: b_normalized(:)  ! B_NORMALIZED(M)
real(RP), allocatable :: chist_loc(:)  ! CHIST_LOC(MAXCHIST)
real(RP), allocatable :: fhist_loc(:)  ! FHIST_LOC(MAXFHIST)
real(RP), allocatable :: xhist_loc(:, :)  ! XHIST_LOC(N, MAXXHIST)

! Sizes
if (present(b)) then
    m = int(size(b), kind(m))
else
    m = 0
end if
n = int(size(x), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(present(A) .eqv. present(b), 'A and B are both present or both absent', srname)
    if (present(A)) then
        call assert((size(A, 1) == n .and. size(A, 2) == m) &
            & .or. (size(A, 1) == 0 .and. size(A, 2) == 0 .and. m == 0), &
            & 'SIZE(A) == [N, M] unless A and B are both empty', srname)
    end if
end if

! Read the inputs

x = moderatex(x)

call safealloc(A_loc, n, m)  ! NOT removable even in F2003, as A may be absent or of size 0-by-0.
if (present(A) .and. m > 0) then
    ! We must check M > 0. Otherwise, the size of A_LOC may be changed to 0-by-0 due to automatic
    ! (re)allocation if that is the size of A; we allow A to be 0-by-0, but A_LOC should be n-by-0.
    A_loc = A
end if

call safealloc(b_loc, m)  ! NOT removable even in F2003, as B may be absent.
if (present(b)) then
    b_loc = b
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


! Normalize the constraints after increasing the right hand sides if necessary so that the starting
! point is feasible.
! N.B.: Yes, LINCOA modifies the right hand sides of the constraints to make the starting point
! feasible if it is not. This is not ideal, but Powell's code was implemented in this way. In the
! MATLAB/Python code, we include a preprocessing subroutine to project the starting point to
! the feasible region if it is infeasible, so that the modification will not occur.
constr_modified = .false.; 
smallx = 1.0E-6_RP * rhoend_loc
call safealloc(A_normalized, n, m)
call safealloc(b_normalized, m)
if (m > 0) then
    do j = 1, m
        ax = inprod(x, A_loc(:, j))
        anorm = norm(A_loc(:, j))
        if (anorm <= 0) then
            if (present(info)) then
                info = ZERO_LINEAR_CONSTRAINT
            end if
            return
        end if
        constr_modified = constr_modified .or. (ax - b_loc(j) > smallx * anorm)
        b_normalized(j) = max(b_loc(j), ax) / anorm
        A_normalized(:, j) = A_loc(:, j) / anorm
    end do
    if (constr_modified) then
        call warning(solver, 'The starting point is infeasible. '//solver// &
            & ' modified the right-hand sides of the constraints to make it feasible')
    end if
end if


!-------------------- Call LINCOB, which performs the real calculations. --------------------------!
call lincob(calfun, iprint_loc, maxfilt_loc, maxfun_loc, npt_loc, A_loc, A_normalized, &
    & b_loc, b_normalized, ctol_loc, cweight_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, &
    & gamma2_loc, rhobeg_loc, rhoend_loc, x, nf_loc, chist_loc, cstrv_loc, f, fhist_loc, xhist_loc, info_loc)
!--------------------------------------------------------------------------------------------------!

! Deallocate variables not needed any more. Indeed, automatic allocation will take place at exit.
deallocate (A_loc)
deallocate (A_normalized)
deallocate (b_loc)
deallocate (b_normalized)


! Write the outputs.

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
    ! XHIST cointains only valid history. For this reason, there is no way to avoid allocating
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
    write (wmsg, ifmt) maxhist_loc
    call warning(solver, 'Only the history of the last '//trim(wmsg)//' iteration(s) is recorded')
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
        call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f, cstrv_loc, ctol_loc)),&
            & 'No point in the history is better than X', srname)
    end if
end if

end subroutine lincoa


end module lincoa_mod
