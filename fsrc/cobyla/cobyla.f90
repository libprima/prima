module cobyla_mod
!--------------------------------------------------------------------------------------------------!
! COBYLA_MOD is a module providing a modern Fortran implementation of Powell's COBYLA algorithm in
!
! M. J. D. Powell, A direct search optimization method that models the objective and constraint
! functions by linear interpolation, In Advances in Optimization and Numerical Analysis, eds. S.
! Gomez and J. P. Hennart, pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994
!
! COBYLA minimizes an objective function F(X) subject to M inequality constraints on X, where X is
! a vector of variables that has N components. The algorithm employs linear approximations to the
! objective and constraint functions, the approximations being formed by linear interpolation at N+1
! points in the space of the variables. We regard these interpolation points as vertices of
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
! Started: July 2021
!
! Last Modified: Friday, December 17, 2021 AM11:36:47
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobyla


contains


subroutine cobyla(calcfc, x, f, &
    & cstrv, constr, &
    & m, f0, constr0, &
    & nf, rhobeg, rhoend, ftarget, ctol, maxfun, iprint, &
    & xhist, fhist, conhist, chist, maxhist, info)
!& eta1, eta2, gamma1, gamma2, xhist, fhist, conhist, chist, maxhist, info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXMEMORY, MAXFUN_DIM_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, FTARGET_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TEN, TENTH, EPS
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_inf, is_finite, is_neginf, is_posinf
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof
use, non_intrinsic :: pintrf_mod, only : FUNCON
use, non_intrinsic :: selectx_mod, only : isbetter
!use, non_intrinsic :: preproc_mod, only : preproc

! Solver-specific modules
use, non_intrinsic :: cobylb_mod, only : cobylb

implicit none

! Compulsory arguments
procedure(FUNCON) :: calcfc
real(RP), intent(inout) :: x(:)
real(RP), intent(out) :: f

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: m
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
real(RP), intent(in), optional :: constr0(:)
real(RP), intent(in), optional :: ctol
real(RP), intent(in), optional :: f0
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), allocatable, optional :: chist(:)
real(RP), intent(out), allocatable, optional :: conhist(:, :)
real(RP), intent(out), allocatable, optional :: fhist(:)
real(RP), intent(out), allocatable, optional :: xhist(:, :)
real(RP), intent(out), optional :: constr(:)
real(RP), intent(out), optional :: cstrv

! Local variables
character(len=*), parameter :: solver = 'COBYLA'
character(len=*), parameter :: srname = 'COBYLA'
integer(IK) :: i
integer(IK) :: iprint_loc
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: maxhist_in
integer(IK) :: maxxhist
integer(IK) :: maximal_hist
integer(IK) :: m_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: info_loc
real(RP) :: cstrv_loc
real(RP) :: ctol_loc
real(RP) :: ftarget_loc
real(RP), allocatable :: chist_loc(:)
real(RP), allocatable :: conhist_loc(:, :)
!real(RP), allocatable :: constr0_loc(:)
real(RP), allocatable :: constr_loc(:)
real(RP), allocatable :: fhist_loc(:)
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP), allocatable :: xhist_loc(:, :)

! Sizes
n = int(size(x), kind(n))
if (present(m) .and. present(constr0)) then
    call assert(m == size(constr0), 'M == SIZE(CONSTR0)', srname)
end if
if (present(m)) then
    m_loc = m
elseif (present(constr0)) then
    m_loc = int(size(constr0), kind(m_loc))
else
    m_loc = 0_IK
end if

! Allocate memory for CONSTR_LOC
call safealloc(constr_loc, m_loc)

! Replace any NaN or Inf in X by ZERO.
where (is_nan(x) .or. is_inf(x))
    x = ZERO
end where

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

if (present(ctol)) then
    ctol_loc = ctol
else
    ctol_loc = CTOL_DFT
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

maxhist_in = 0_IK  ! MAXHIST input by user
if (present(maxhist)) then
    maxhist_loc = maxhist
    maxhist_in = maxhist
else if (maxfun_loc >= n + 2) then
    maxhist_loc = maxfun_loc
else
    maxhist_loc = MAXFUN_DIM_DFT * n
end if

! Preprocess the inputs in case some of them are invalid.
!call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc)

! Further revise MAXHIST according to MAXMEMORY, i.e., the maximal memory allowed for the history.
if (present(xhist) .and. present(conhist)) then
    maximal_hist = int(MAXMEMORY / ((n + m_loc + 2) * cstyle_sizeof(0.0_RP)), kind(maximal_hist))
elseif (present(xhist)) then
    maximal_hist = int(MAXMEMORY / ((n + 2) * cstyle_sizeof(0.0_RP)), kind(maximal_hist))
elseif (present(conhist)) then
    maximal_hist = int(MAXMEMORY / ((m_loc + 2) * cstyle_sizeof(0.0_RP)), kind(maximal_hist))
else
    maximal_hist = int(MAXMEMORY / (2 * cstyle_sizeof(0.0_RP)), kind(maximal_hist))
end if
if (maxhist_loc > maximal_hist) then
    ! We cannot simply take MAXHIST_LOC = MIN(MAXHIST_LOC, MAXIMAL_HIST), as they may not have the
    ! same kind, and compilers may complain. We may convert them to the same, but overflow may occur
    maxhist_loc = int(maximal_hist, kind(maxhist_loc))
end if

! Allocate memory for the history of X. We use XHIST_LOC instead of XHIST, which be absent.
if (present(xhist)) then
    maxxhist = maxhist_loc
else
    maxxhist = 0_IK
end if
call safealloc(xhist_loc, n, maxxhist)

! Allocate memory for the history of F. We use FHIST_LOC instead of FHIST, which be absent.
if (present(fhist)) then
    maxfhist = maxhist_loc
else
    maxfhist = 0_IK
end if
call safealloc(fhist_loc, maxfhist)

! Allocate memory for the history of CONSTR. Use CONHIST_LOC instead of CONHIST, which may be absent.
if (present(conhist)) then
    maxconhist = maxhist_loc
else
    maxconhist = 0_IK
end if
call safealloc(conhist_loc, m_loc, maxconhist)

! Allocate memory for the history of CSTRV. We use CHIST_LOC instead of CHIST, which may be absent.
if (present(chist)) then
    maxchist = maxhist_loc
else
    maxchist = 0_IK
end if
call safealloc(chist_loc, maxchist)

!-------------------- Call COBYLB, which performs the real calculations. --------------------------!
call cobylb(calcfc, iprint_loc, maxfun_loc, ctol_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & constr_loc, x, nf_loc, chist_loc, conhist_loc, cstrv_loc, f, fhist_loc, xhist_loc, info_loc)
!--------------------------------------------------------------------------------------------------!

! Write the outputs.

! Copy CONSTR_LOC to CONSTR if needed.
if (present(constr)) then
    !--------------------------------------------------!
    !---- The SAFEALLOC line is removable in F2003. ---!
    !!call safealloc(constr, m_loc)
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

! Copy CONHIST_LOC to CONHIST if needed.
if (present(conhist)) then
    !--------------------------------------------------!
    !---- The SAFEALLOC line is removable in F2003. ---!
    call safealloc(conhist, m_loc, min(nf_loc, maxconhist))
    !--------------------------------------------------!
    conhist = conhist_loc(:, 1:min(nf_loc, maxconhist))
    ! The same as XHIST, we must cap CONHIST at min(NF_LOC, MAXCONHIST).
end if
deallocate (conhist_loc)

! Copy CHIST_LOC to CHIST if needed.
if (present(chist)) then
    !--------------------------------------------------!
    !---- The SAFEALLOC line is removable in F2003. ---!
    call safealloc(chist, min(nf_loc, maxchist))
    !--------------------------------------------------!
    chist = chist_loc(1:min(nf_loc, maxchist))
    ! The same as XHIST, we must cap CHIST at min(NF_LOC, MAXCHIST).
end if
deallocate (chist_loc)

! If MAXFHIST_IN >= NF_LOC > MAXFHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist)) .and. maxhist_loc < min(nf_loc, maxhist_in)) then
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
        call assert(.not. any(is_nan(fhist) .or. is_posinf(fhist)), 'FHIST does not contain NaN/+Inf', srname)
    end if
    if (present(conhist)) then
        call assert(size(conhist, 1) == m_loc .and. size(conhist, 2) == nhist, 'SIZE(CONHIST) == [M, NHIST]', srname)
        call assert(.not. any(is_nan(conhist) .or. is_neginf(conhist)), 'CONHIST does not contain NaN/-Inf', srname)
    end if
    if (present(chist)) then
        call assert(size(chist) == nhist, 'SIZE(CHIST) == NHIST', srname)
        call assert(.not. any(is_nan(chist) .or. is_posinf(chist)), 'CHIST does not contain NaN/+Inf', srname)
    end if
    if (present(fhist) .and. present(chist)) then
        call assert(.not. any([(isbetter([fhist(i), chist(i)], [f, cstrv], ctol), i=1, nhist)]), &
            & 'No point in the filter is better than X', srname)
    end if
end if

end subroutine cobyla


end module cobyla_mod
