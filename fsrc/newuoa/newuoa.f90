module newuoa_mod
!--------------------------------------------------------------------------------------------------!
! NEWUOA_MOD is a module providing a modern Fortran implementation of Powell's NEWUOA algorithm in
!
! M. J. D. Powell, The NEWUOA software for unconstrained optimization without derivatives, In Large-
! Scale Nonlinear Optimization, eds. G. Di Pillo and M. Roma, 255--297, Springer, New York, 2006
!
! NEWUOA seeks the least value of a function of many variables, by a trust region method that forms
! quadratic models by interpolation. There can be some freedom in the interpolation conditions,
! which is taken up by minimizing the Frobenius norm of the change to the second derivative of the
! quadratic model, beginning with a zero matrix.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Friday, December 17, 2021 PM03:34:51
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: newuoa


contains


subroutine newuoa(calfun, x, f, &
    & nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, &
    & eta1, eta2, gamma1, gamma2, xhist, fhist, maxhist, info)
!--------------------------------------------------------------------------------------------------!
! Among all the arguments, only CALFUN, X, and F are obligatory. The others are OPTIONAL and you can
! neglect them unless you are familiar with the algorithm. If you do not specify an optional
! argument, it will be assigned the default value  as will be explained later. For example, it is
! valid to call NEWUOA by
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
! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory name
! "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
! real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
! setting __REAL_PRECISION__ and __INTEGER_KIND__ in common/ppf.h. Use the default if unsure.
!
! CALFUN
!   Input, subroutine.
!   CALFUN(X, F) should evaluate the objective function at the given REAL(RP) vector X and set the
!   value to the REAL(RP) scalar F. It must be provided by the user.
!
! X
!   Input and outout, REAL(RP) vector.
!   As an input, X should be an N dimensional vector that contains the starting point, N being the
!   dimension of the problem. As an output, X will be set to an approximate minimizer.
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
!   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
!   the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
!   Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
!   RHOEND should indicate the accuracy that is required in the final values of the variables.
!
! FTARGET
!   Input, REAL(RP) scalar, default: - Inf.
!   FTARGET is the target function value. The algorithm will terminate when a point with a function
!   value <= FTARGET is found.
!
! MAXFUN
!   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
!   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of function evaluations.
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
!      variables so far and its objective function value;
!   3: in addition to 2, each function evaluation with its variables will be printed to the screen;
!   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
!      named NEWUOA_output.txt; the file will be created if it does not exist; the new output will
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
! XHIST, FHIST, MAXHIST
!   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
!   XHIST, if present, will output the history of iterates, while FHIST, if present, will output the
!   history function values. MAXHIST should be a nonnegative integer, and XHIST/FHIST will output
!   only the last MAXHIST iterates and/or the corresponding function values. Therefore, MAXHIST
!   = 0 means XHIST/FHIST will output nothing, while setting MAXHIST = MAXFUN ensures that
!   XHIST/FHIST will output all the history. If XHIST is present, its size at exit will be (N,
!   min(NF, MAXHIST)); if FHIST is present, its size at exit will be min(NF, MAXHIST).
!
!   Important Notice:
!   Setting MAXHIST to a large value can be costly in terms of memory for problems with a large N.
!   For instance, if N = 1000 and MAXHIST = 100, 000, XHIST will take up to 1 GB if we use double
!   precision. MAXHIST will be  reset to a smaller value if the memory needed for XHIST and/or FHIST
!   exceeds MAXMEMORY defined in CONSTS_MOD (see consts.F90 under the directory named "common";
!   default: 2GB). Use XHIST, FHIST, and MAXHIST with caution!!! (Fortunately, the algorithm is NOT
!   designed for large problems).
!
! INFO
!   Output, INTEGER(IK) scalar.
!   INFO is the exit flag. It will be set to one of the following values defined in the module
!   INFO_MOD (see common/info.F90):
!   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
!   FTARGET_ACHIEVED: the target function value is reached;
!   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
!   MAXTR_REACHED: the trust region iteration has been performed MAXTR times,
!       the value of MAXTR being 10*MAXFUN, which is UNLIKELY to reach;
!   NAN_INF_X: NaN or Inf occurs in x;
!   !--------------------------------------------------------------------------!
!   The following cases should NEVER occur unless there is a bug, because the
!   code will try to continue in the corresponding scenarios.
!   NAN_INF_F: the objective function returns NaN or +Inf
!   NAN_MODEL: NaN occurs in the model;
!   TRSUBP_FAILED: a trust region step failed to reduce the quadratic model
!   !--------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXMEMORY, MAXFUN_DIM_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_inf, is_finite
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof
use, non_intrinsic :: pintrf_mod, only : FUN
use, non_intrinsic :: preproc_mod, only : preproc

! Solver-specific modules
use, non_intrinsic :: newuob_mod, only : newuob

implicit none

! Arguments
procedure(FUN) :: calfun
! N.B.: The INTENT attribute cannot be specified for a dummy procedure without the POINTER attribute
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
character(len=*), parameter :: solver = 'NEWUOA'
character(len=*), parameter :: srname = 'NEWUOA'
integer :: maximal_hist
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfhist
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: maxhist_in
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: npt_loc
integer(IK) :: unit_memo
logical :: output_hist
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP), allocatable :: fhist_loc(:)
real(RP), allocatable :: xhist_loc(:, :)

! Sizes
n = int(size(x), kind(n))

! Replace any NaN or Inf in X by ZERO.
where (is_nan(x) .or. is_inf(x))
    x = ZERO
end where

! Read the inputs.

! If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
! value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
! it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
if (present(rhobeg)) then
    rhobeg_loc = rhobeg
else if (present(rhoend)) then
    ! Fortran does not take short-circuit evaluation of logic expressions. Thus it is WRONG to
    ! combine the evaluation of PRESENT(RHOEND) and the evaluation of IS_FINITE(RHOEND) as
    ! "if (present(rhoend) .and. is_finite(rhoend))". The compiler may choose the evaluate the
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
else if (rhobeg_loc > 0) then
    rhoend_loc = max(EPS, min(TENTH * rhobeg_loc, RHOEND_DFT))
else
    rhoend_loc = RHOEND_DFT
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
else if (maxfun_loc >= 1) then
    npt_loc = int(max(n + 2, min(maxfun_loc - 1, 2 * n + 1)), kind(npt_loc))
else
    npt_loc = int(2 * n + 1, kind(npt_loc))
end if

if (present(iprint)) then
    iprint_loc = iprint
else
    iprint_loc = IPRINT_DFT
end if

if (present(eta1)) then
    eta1_loc = eta1
else if (present(eta2)) then
    if (eta2 > ZERO .and. eta2 < ONE) then
        eta1_loc = max(EPS, eta2 / 7.0_RP)
    end if
else
    eta1_loc = TENTH
end if

if (present(eta2)) then
    eta2_loc = eta2
else if (eta1_loc > ZERO .and. eta1_loc < ONE) then
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

maxhist_in = 0_IK  ! MAXHIST input by user
if (present(maxhist)) then
    maxhist_in = maxhist
    maxhist_loc = maxhist
else if (maxfun_loc >= n + 3) then
    maxhist_loc = maxfun_loc
else
    maxhist_loc = MAXFUN_DIM_DFT * n
end if
output_hist = (present(xhist) .or. present(fhist))
if (.not. output_hist) then
    maxhist_in = 0
    maxhist_loc = 0
end if

! Preprocess the inputs in case some of them are invalid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, &
    & gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc)

! Further revise MAXHIST according to MAXMEMORY, i.e., the maximal memory allowed for the history.
unit_memo = 0_IK
if (present(xhist)) then
    unit_memo = unit_memo + n
end if
if (present(fhist)) then
    unit_memo = unit_memo + 1_IK
end if
maximal_hist = int(MAXMEMORY / max(1_IK, unit_memo * cstyle_sizeof(0.0_RP)), kind(maximal_hist))
if (maxhist_loc > maximal_hist) then
    ! We cannot simply take MAXHIST_LOC = MIN(MAXHIST_LOC, MAXIMAL_HIST), as they may not have the
    ! same kind, and compilers may complain. We may convert them to the same, but overflow may occur
    maxhist_loc = int(maximal_hist, kind(maxhist_loc))
end if

! Allocate memory for the history of X. We use XHIST_LOC instead of XHIST, which may not be present.
if (present(xhist)) then
    maxxhist = maxhist_loc
else
    maxxhist = 0_IK
end if
call safealloc(xhist_loc, n, maxxhist)

! Allocate memory for the history of F. We use FHIST_LOC instead of FHIST, which may not be present.
if (present(fhist)) then
    maxfhist = maxhist_loc
else
    maxfhist = 0_IK
end if
call safealloc(fhist_loc, maxfhist)

!-------------------- Call NEWUOB, which performs the real calculations. --------------------------!
call newuob(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, &
    & gamma2_loc, rhobeg_loc, rhoend_loc, x, nf_loc, f, fhist_loc, xhist_loc, info_loc)
!--------------------------------------------------------------------------------------------------!

! Write the outputs.

if (present(nf)) then
    nf = nf_loc
end if

if (present(info)) then
    info = info_loc
end if

! Copy XHIST_LOC to XHIST if needed.
if (present(xhist)) then
    !--------------------------------------------------!
    !---- The SAFEALLOC line is removable in F2003. ---!
    call safealloc(xhist, n, min(nf_loc, maxxhist))
    !--------------------------------------------------!
    xhist = xhist_loc(:, 1:min(nf_loc, maxxhist))
    ! N.B.:
    ! 0. Allocate XHIST as long as it is present, even if MAXXHIST = 0; otherwise, it will be
    ! illegal to enquire XHIST after exit.
    ! 1. Even though Fortran 2003 supports automatic (re)allocation of allocatable arrays upon
    ! intrinsic assignment, we keep the line of SAFEALLOC, because some very new compilers (Absoft
    ! Fortran 20.0) are still not standard-compliant in this respect.
    ! 2. NF may not be present. Hence we should NOT use NF but NF_LOC.
    ! 3. When MAXXHIST > NF_LOC, which is the normal case in practice, XHIST_LOC contains GARBAGE in
    ! XHIST_LOC(:, NF_LOC + 1 : MAXXHIST). Therefore, we MUST cap XHIST at min(NF_LOC, MAXXHIST) so
    ! that XHIST cointains only valid history. For this reason, there is no way to avoid allocating
    ! two copies of memory for XHIST unless we declare it to be a POINTER instead of ALLOCATABLE.
end if
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer to deallocate
! them immediately when they finish their jobs.
deallocate (xhist_loc)

! Copy FHIST_LOC to FHIST if needed.
if (present(fhist)) then
    !--------------------------------------------------!
    !---- The SAFEALLOC line is removable in F2003. ---!
    call safealloc(fhist, min(nf_loc, maxfhist))
    !--------------------------------------------------!
    fhist = fhist_loc(1:min(nf_loc, maxfhist))
    ! The same as XHIST, we must cap FHIST at min(NF_LOC, MAXFHIST).
end if
deallocate (fhist_loc)

! If MAXFHIST_IN >= NF_LOC > MAXFHIST_LOC, warn that not all history is recorded.
if (output_hist .and. maxhist_loc < min(nf_loc, maxhist_in)) then
    print '(/1A, I7, 1A)', 'WARNING: '//solver//': due to memory limit, MAXHIST is reset to ', maxhist_loc, '.'
    print '(1A/)', 'Only the history of the last MAXHIST iterations is recoreded.'
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
        call assert(.not. any(fhist < f), 'F is the smallest in FHIST', srname)
    end if
end if

end subroutine newuoa


end module newuoa_mod
