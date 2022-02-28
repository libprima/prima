module cobyla_mod
!--------------------------------------------------------------------------------------------------!
! COBYLA_MOD is a module providing a modernized and improved Fortran implementation of Powell's
! COBYLA algorithm in
!
! M. J. D. Powell, A direct search optimization method that models the objective and constraint
! functions by linear interpolation, In Advances in Optimization and Numerical Analysis, eds. S.
! Gomez and J. P. Hennart, pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994
!
! COBYLA approximately solves
!
!   min F(X) subject to CONSTR(X) >= 0,
!
! where X is a vector of variables that has N components, F is a real-valued objective function, and
! CONSTR(X) is an M-dimensional vector-valued mapping. The algorithm employs linear approximations
! to the objective and constraint functions, the approximations being formed by linear interpolation
! at N + 1 points in the space of the variables. We regard these interpolation points as vertices of
! a simplex. The parameter RHO controls the size of the simplex and it is reduced automatically from
! RHOBEG to RHOEND. For each RHO the subroutine tries to achieve a good vector of variables for the
! current size, and then RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and
! RHOEND should be set to reasonable initial changes to and the required accuracy in the variables
! respectively, but this accuracy should be viewed as a subject for experimentation because it is
! not guaranteed.  The subroutine has an advantage over many of its competitors, however, which is
! that it treats each constraint individually when calculating a change to the variables, instead of
! lumping the constraints together into a single penalty function. The name of the subroutine is
! derived from the phrase Constrained Optimization BY Linear Approximations.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Monday, February 28, 2022 PM08:49:38
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobyla


contains


subroutine cobyla(calcfc, m, x, f, &
    & cstrv, constr, &
    & f0, constr0, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, chist, conhist, maxhist, maxfilt, info)
!--------------------------------------------------------------------------------------------------!
! Among all the arguments, only CALCFC, X, and F are obligatory. The others are OPTIONAL and you can
! neglect them unless you are familiar with the algorithm. If you do not specify an optional input,
! it will be assigned the default value detailed below. For instance, we may invoke the solver by
!
! call cobyla(calcfc, x, f)
!
! or
!
! call cobyla(calcfc, x, f, rhobeg = 0.5D0, rhoend = 1.0D-3, maxfun = 100)
!
! See examples/cobyla_exmp.f90 for a concrete example.
!
! A detailed introduction to the arguments is as follows.
! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory name
! "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
! real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
! setting __REAL_PRECISION__ and __INTEGER_KIND__ in common/ppf.h. Use the default if unsure.
!
! CALCFC
!   Input, subroutine.
!   CALCFC(X, F, CONSTR) should evaluate the objective function and constraints at the given
!   REAL(RP) vector X; it should set the objective function value to the REAL(RP) scalar F and the
!   constraint value to the REAL(RP) vector CONSTR. It must be provided by the user, and its
!   definition must conform to the following interface:
!   !-------------------------------------------------------------------------!
!    subroutine calcfc(x, f, constr)
!    real(RP), intent(in) :: x(:)
!    real(RP), intent(out) :: f
!    real(RP), intent(out) :: constr(:)
!    end subroutine calcfc
!   !-------------------------------------------------------------------------!
!   Besides, the size of CONSTR must be M, which is the second compulsory argument (see below).
!
! M
!   Input, INTEGER(IK) scalar.
!   M must be set to the number of constraints, namely the size (length) of CONSTR(X).
!
! X
!   Input and output, REAL(RP) vector.
!   As an input, X should be an N-dimensional vector that contains the starting point, N being the
!   dimension of the problem. As an output, X will be set to an approximate minimizer.
!
! F
!   Output, REAL(RP) scalar.
!   F will be set to the objective function value of X at exit.
!
! CSTRV
!   Output, REAL(RP) scalar.
!   CSTRV will be set to the constraint violation of X at exit, i.e., MAXVAL([-CONSTR(X), 0]).
!
! CONSTR
!   Output, ALLOCATABLE REAL(RP) vector.
!   CONSTR will be set to the constraint value of X at exit.
!
! F0
!   Input, REAL(RP) scalar.
!   F0, if present, should be set to the objective function value of the starting X.
!
! CONSTR0
!   Input, REAL(RP) vector.
!   CONSTR0, if present, should be set to the constraint value of the starting X; in addition,
!   SIZE(CONSTR0) must be M, or the solver will abort.
!
! NF
!   Output, INTEGER(IK) scalar.
!   NF will be set to the number of calls of CALCFC at exit.
!
! RHOBEG, RHOEND
!   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
!   the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
!   Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
!   RHOEND should indicate the accuracy that is required in the final values of the variables.
!
! FTARGET
!   Input, REAL(RP) scalar, default: -Inf.
!   FTARGET is the target function value. The algorithm will terminate when a feasible point with a
!   function value <= FTARGET is found.
!
! CTOL
!   Input, REAL(RP) scalar, default: machine epsilon.
!   CTOL is the tolerance of constraint violation. Any X with MAXVAL(-CONSTR(X)) <= CTOL is
!   considered as feasible.
!
! CWEIGHT
!   Input, REAL(RP) scalar, default: CWEIGHT_DFT defined in the module CONSTS_MOD in common/consts.F90.
!   CWEIGHT is the weight that the constraint violation takes in the selection of the returned X.
!
! MAXFUN
!   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
!   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of calls of CALCFC.
!
! NPT
!   Input, INTEGER(IK) scalar, default: 2N + 1.
!   NPT is the number of interpolation conditions for each trust region model. Its value must be in
!   the interval [N+2, (N+1)(N+2)/2].
!
! IPRINT
!   Input, INTEGER(IK) scalar, default: 0.
!   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls how much
!   information will be printed during the computation:
!   0: there will be no printing;
!   1: a message will be printed to the screen at the return, showing the best vector of variables
!      found and its objective function value;
!   2: in addition to 1, each new value of RHO is printed to the screen, with the best vector of
!      variables so far and its objective function value; each new value of CPEN is also printed;
!   3: in addition to 2, each function evaluation with its variables will be printed to the screen;
!   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
!      named COBYLA_output.txt; the file will be created if it does not exist; the new output will
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
! XHIST, FHIST, CHIST, CONHIST, MAXHIST
!   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   CHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   CONHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
!   XHIST, if present, will output the history of iterates; FHIST, if present, will output the
!   history function values; CHIST, if present, will output the history of constraint violations;
!   CONHIST, if present, will output the history of constraint values; MAXHIST should be a
!   nonnegative integer, and XHIST/FHIST/CHIST/CONHIST will output only the history of the last
!   MAXHIST iterations. Therefore, MAXHIST= 0 means XHIST/FHIST/CONHIST/CHIST will output nothing,
!   while setting MAXHIST = MAXFUN requests XHIST/FHIST/CHIST/CONHIST to output all the history.
!   If XHIST is present, its size at exit will be (N, min(NF, MAXHIST)); if FHIST is present, its
!   size at exit will be min(NF, MAXHIST); if CHIST is present, its size at exit will be
!   min(NF, MAXHIST); if CONHIST is present, its size at exit will be (M, min(NF, MAXHIST)).
!
!   Important Notice:
!   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
!   For instance, if N = 1000 and MAXHIST = 100, 000, XHIST will take up to 1 GB if we use double
!   precision. MAXHIST will be reset to a smaller value if the memory needed exceeds MAXMEMORY
!   defined in CONSTS_MOD (see consts.F90 under the directory named "common"; default: 2GB).
!   Use *HIST with caution!!! (N.B.: the algorithm is NOT designed for large problems).
!
! MAXFILT
!   Input, INTEGER(IK) scalar.
!   MAXFILT is a nonnegative integer indicating the maximal length of the filter used for selecting
!   the returned solution; default: MAXFILT_DFT (a value lower than MIN_MAXFILT is not recommended);
!   see common/consts.F90 for the definitions of MAXFILT_DFT and MIN_MAXFILT.
!
! INFO
!   Output, INTEGER(IK) scalar.
!   INFO is the exit flag. It will be set to one of the following values defined in the module
!   INFO_MOD (see common/info.F90):
!   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
!   FTARGET_ACHIEVED: the target function value is reached;
!   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
!   MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
!   NAN_INF_X: NaN or Inf occurs in X;
!   DAMAGING_ROUNDING: rounding errors are becoming damaging.
!   !--------------------------------------------------------------------------!
!   The following case(s) should NEVER occur unless there is a bug.
!   NAN_INF_F: the objective function returns NaN or +Inf;
!   NAN_MODEL: NaN occurs in the model;
!   TRSUBP_FAILED: a trust region step failed to reduce the model
!   !--------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, errstop, warning
use, non_intrinsic :: evaluate_mod, only : evaluate, moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_neginf, is_posinf
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJCON
use, non_intrinsic :: selectx_mod, only : isbetter
use, non_intrinsic :: preproc_mod, only : preproc

! Solver-specific modules
use, non_intrinsic :: cobylb_mod, only : cobylb

implicit none

! Compulsory arguments
procedure(OBJCON) :: calcfc
real(RP), intent(inout) :: x(:)
real(RP), intent(out) :: f
integer(IK), intent(in) :: m

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfilt
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
real(RP), intent(in), optional :: constr0(:)
real(RP), intent(in), optional :: ctol
real(RP), intent(in), optional :: cweight
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: f0
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), allocatable, optional :: chist(:)
real(RP), intent(out), allocatable, optional :: conhist(:, :)
real(RP), intent(out), allocatable, optional :: constr(:)
real(RP), intent(out), allocatable, optional :: fhist(:)
real(RP), intent(out), allocatable, optional :: xhist(:, :)
real(RP), intent(out), optional :: cstrv

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
character(len=*), parameter :: solver = 'COBYLA'
character(len=*), parameter :: srname = 'COBYLA'
character(len=MSGLEN) :: wmsg
integer(IK) :: i
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfilt_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
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
real(RP), allocatable :: chist_loc(:)
real(RP), allocatable :: conhist_loc(:, :)
real(RP), allocatable :: constr_loc(:)
real(RP), allocatable :: fhist_loc(:)
real(RP), allocatable :: xhist_loc(:, :)

! Preconditions
if (DEBUGGING) then
    call assert(present(f0) .eqv. present(constr0), 'F0 and CONSTR0 are both present or both absent', srname)
end if

! Sizes
n = int(size(x), kind(n))

! Exit if the size of CONSTR0 is inconsistent with M.
if (present(constr0)) then
    if (size(constr0) /= m) then
        if (DEBUGGING) then
            call errstop(srname, 'SIZE(CONSTR0) /= M. Exiting')
        else
            call warning(srname, 'SIZE(CONSTR0) /= M. Exiting')
            return
        end if
    end if
end if

! Allocate memory for CONSTR_LOC, since M is now available.
call safealloc(constr_loc, m)  ! NOT removable even in F2003!

! If the user provides the function & constraint value at X0, then set up [F, CONSTR_LOC] to them.
! Otherwise, set [F, CONSTR_LOC] = [F(X0), CONSTR(X0)], so that COBYLB only needs a single interface.
if (present(f0) .and. present(constr0) .and. all(is_finite(x))) then
    f = f0
    constr_loc = constr0
else
    ! Replace any NaN in X by ZERO and Inf/-Inf in X by HUGENUM/-HUGENUM.
    x = moderatex(x)
    call evaluate(calcfc, x, f, constr_loc, cstrv_loc) ! Indeed, CSTRV_LOC needs not to be evaluated.
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
    if (is_finite(rhoend) .and. rhoend > ZERO) then
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

if (present(iprint)) then
    iprint_loc = iprint
else
    iprint_loc = IPRINT_DFT
end if

if (present(eta1)) then
    eta1_loc = eta1
elseif (present(eta2)) then
    if (eta2 > ZERO .and. eta2 < ONE) then
        eta1_loc = max(EPS, eta2 / 7.0_RP)
    end if
else
    eta1_loc = TENTH
end if

if (present(eta2)) then
    eta2_loc = eta2
elseif (eta1_loc > ZERO .and. eta1_loc < ONE) then
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
    maxhist_loc = maxval([maxfun_loc, n + 2_IK, MAXFUN_DIM_DFT * n])
end if

if (present(maxfilt)) then
    maxfilt_loc = maxfilt
else
    maxfilt_loc = MAXFILT_DFT
end if

! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & m=m, ctol=ctol_loc, cweight=cweight_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, &
    & gamma2=gamma2_loc, maxfilt=maxfilt_loc)

! Further revise MAXHIST_LOC according to MAXMEMORY, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! CHIST = NaN(1, MAXFUN), CONHIST = NaN(M, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc, &
    & present(chist), chist_loc, m, present(conhist), conhist_loc)


!-------------------- Call COBYLB, which performs the real calculations. --------------------------!
call cobylb(calcfc, iprint_loc, maxfilt_loc, maxfun_loc, ctol_loc, cweight_loc, eta1_loc, eta2_loc, &
    & ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, constr_loc, f, x, nf_loc, &
    & chist_loc, conhist_loc, cstrv_loc, fhist_loc, xhist_loc, info_loc)
!--------------------------------------------------------------------------------------------------!


! Write the outputs.

! Copy CONSTR_LOC to CONSTR if needed.
if (present(constr)) then
    !--------------------------------------------------!
    call safealloc(constr, m)  ! Removable in F2003.
    !--------------------------------------------------!
    constr = constr_loc
end if
deallocate (constr_loc)

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
    nhist = min(nf_loc, int(size(xhist_loc, 2), IK))
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
    nhist = min(nf_loc, int(size(fhist_loc), IK))
    !--------------------------------------------------!
    call safealloc(fhist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    fhist = fhist_loc(1:nhist)  ! The same as XHIST, we must cap FHIST at NF_LOC.
end if
deallocate (fhist_loc)

! Copy CHIST_LOC to CHIST if needed.
if (present(chist)) then
    nhist = min(nf_loc, int(size(chist_loc), IK))
    !--------------------------------------------------!
    call safealloc(chist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    chist = chist_loc(1:nhist)  ! The same as XHIST, we must cap CHIST at NF_LOC.
end if
deallocate (chist_loc)

! Copy CONHIST_LOC to CONHIST if needed.
if (present(conhist)) then
    nhist = min(nf_loc, int(size(conhist_loc, 2), IK))
    !----------------------------------------------------------!
    call safealloc(conhist, m, nhist)  ! Removable in F2003.
    !----------------------------------------------------------!
    conhist = conhist_loc(:, 1:nhist)  ! The same as XHIST, we must cap CONHIST at NF_LOC.
end if
deallocate (conhist_loc)

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist) .or. present(chist) .or. present(conhist)) .and. maxhist_loc < nf_loc) then
    write (wmsg, ifmt) maxhist_loc
    call warning(solver, 'Only the history of the last '//trim(wmsg)//' iteration(s) is recoreded')
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
    if (present(conhist)) then
        call assert(size(conhist, 1) == m .and. size(conhist, 2) == nhist, 'SIZE(CONHIST) == [M, NHIST]', srname)
        call assert(.not. any(is_nan(conhist) .or. is_neginf(conhist)), 'CONHIST does not contain NaN/-Inf', srname)
    end if
    if (present(fhist) .and. present(chist)) then
        call assert(.not. any([(isbetter([fhist(i), chist(i)], [f, cstrv_loc], ctol_loc), i=1, nhist)]),&
            & 'No point in the history is better than X', srname)
    end if
end if

end subroutine cobyla


end module cobyla_mod
