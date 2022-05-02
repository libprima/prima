module bobyqa_mod
!--------------------------------------------------------------------------------------------------!
! BOBYQA_MOD is a module providing a modernized and improved Fortran implementation of Powell's
! BOBYQA algorithm in
!
! M. J. D. Powell, The BOBYQA algorithm for bound constrained optimization without derivatives,
! Technical Report DAMTP 2009/NA06, Department of Applied Mathematics and Theoretical Physics,
! Cambridge University, Cambridge, UK, 2009
!
! BOBYQA approximately solves
!
!   min F(X) subject to XL <= X <= XU,
!
! where X is a vector of variables that has N components, and F is a real-valued objective function.
! XL and XU are a pair of N-dimensional vectors indicating the lower and upper bounds of X. The
! algorithm assumes that XL < XU entrywise. It tackles the problem by applying a trust region method
! that forms quadratic models by interpolation. There is usually some freedom in the interpolation
! conditions, which is taken up by minimizing the Frobenius norm of the change to the second
! derivative of the model, beginning with the ZERO matrix. The values of the variables are
! constrained by upper and lower bounds. The arguments of the subroutine are as follows.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, May 02, 2022 AM09:45:19
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: bobyqa


contains


subroutine bobyqa(calfun, x, f, &
    & xl, xu, &
    & nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, maxhist, honour_x0, info)
!--------------------------------------------------------------------------------------------------!
!     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
!     F to the value of the objective function for the current values of the
!     variables X(1),X(2),...,X(N), which are generated automatically in a
!     way that satisfies the bounds given in XL and XU.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in
!       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
!       recommended.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!       9: the denominator of the updating formula is ZERO.
!       10: N should not be less than 2.
!       11: MAXFUN is less than NPT+1.
!       12: the gradient of constraint is ZERO.
!       -1: NaN occurs in x.
!       -2: the objective function returns a NaN or nearly infinite value.
!       -3: NaN occurs in BMAT or ZMAT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH
use, non_intrinsic :: consts_mod, only : EPS, HUGEBOUND, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: evaluate_mod, only : moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: linalg_mod, only : trueloc
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: preproc_mod, only : preproc

! Solver-specific modules
use, non_intrinsic :: bobyqb_mod, only : bobyqb

implicit none

! Compulsory arguments
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(inout) :: x(:)  ! X(N)
real(RP), intent(out) :: f

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
integer(IK), intent(in), optional :: npt
logical, intent(in), optional :: honour_x0
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
real(RP), intent(out), allocatable, optional :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(out), allocatable, optional :: xhist(:, :)  ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
character(len=*), parameter :: solver = 'BOBYQA'
character(len=*), parameter :: srname = 'BOBYQA'
character(len=MSGLEN) :: wmsg
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: npt_loc
logical :: has_rhobeg
logical :: honour_x0_loc
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP) :: sl(size(x))
real(RP) :: su(size(x))
real(RP) :: xl_loc(size(x))
real(RP) :: xu_loc(size(x))
real(RP), allocatable :: fhist_loc(:)  ! FHIST_LOC(MAXFHIST)
real(RP), allocatable :: xhist_loc(:, :)  ! XHIST_LOC(N, MAXXHIST)

! Sizes
n = int(size(x), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    if (present(xl)) then
        call assert(size(xl) == n .or. size(xl) == 0, 'SIZE(XL) == N unless XL is empty', srname)
    end if
    if (present(xu)) then
        call assert(size(xu) == n .or. size(xu) == 0, 'SIZE(XU) == N unless XU is empty', srname)
    end if
end if

! Read the inputs

if (present(xl)) then
    xl_loc = xl
else
    xl_loc = -HUGEBOUND
end if
xl_loc(trueloc(is_nan(xl_loc) .or. xl_loc < -HUGEBOUND)) = -HUGEBOUND

if (present(xu)) then
    xu_loc = xu
else
    xu_loc = HUGEBOUND
end if
xu_loc(trueloc(is_nan(xu_loc) .or. xu_loc > HUGEBOUND)) = HUGEBOUND

x = max(xl_loc, min(xu_loc, moderatex(x)))

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
    maxhist_loc = maxval([maxfun_loc, n + 3_IK, MAXFUN_DIM_DFT * n])
end if

has_rhobeg = present(rhobeg)
if (present(honour_x0)) then
    honour_x0_loc = honour_x0
else
    honour_x0_loc = (.not. has_rhobeg)
end if


!write (16, *) xl_loc
!write (16, *) xu_loc
!write (16, *) x
!write (16, *) rhobeg_loc
!write (16, *) '---'

! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & npt=npt_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, gamma2=gamma2_loc, &
    & has_rhobeg=has_rhobeg, honour_x0=honour_x0_loc, xl=xl_loc, xu=xu_loc, x0=x)

! Further revise MAXHIST_LOC according to MAXMEMORY, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc)

! The lower and upper bounds on moves from the updated X are set now, in the ISL and ISU partitions
! of W, in order to provide useful and exact information about components of X that become within
! distance RHOBEG from their bounds.
sl = xl_loc - x
su = xu_loc - x
! After PREPROC, the nonzero entries of SL should be less than -RHOBEG_LOC, while those of SU should
! be larger than RHOBEG_LOC. However, this may not be true due to rounding. The following lines
! revise SL and SU to ensure it.
where (sl < 0)
    sl = min(sl, -rhobeg_loc)
elsewhere
    x = xl_loc
    sl = ZERO
    su = xu_loc - xl_loc
end where
where (su > 0)
    su = max(su, rhobeg_loc)
elsewhere
    x = xu_loc
    sl = xl_loc - xu_loc
    su = ZERO
end where
!!MATLAB code for revising X, SL, and SU:
!!sl(sl < 0) = min(sl(sl < 0), -rhobeg_loc);
!!x(sl >= 0) = xl_loc(sl >= 0);
!!sl(sl >= 0) = 0;
!!su(sl >= 0) = xu_loc(sl >= 0) - xl_loc(sl >= 0);
!!su(su > 0) = max(su(su > 0), rhobeg_loc);
!!x(su <= 0) = xu_loc(su <= 0);
!!sl(su <= 0) = xl_loc(su <= 0) - xu_loc(su <= 0);
!!su(su <= 0) = 0;


!write (16, *) rhobeg_loc
!write (16, *) x
!write (16, *) sl
!write (16, *) su

!-------------------- Call BOBYQB, which performs the real calculations. --------------------------!
call bobyqb(calfun, iprint_loc, maxfun_loc, npt_loc, eta1_loc, eta2_loc, ftarget_loc, &
    & gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, sl, su, xl_loc, xu_loc, x, nf_loc, f, &
    & fhist_loc, xhist_loc, info_loc)
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

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist)) .and. maxhist_loc < nf_loc) then
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
        call assert(.not. any(fhist < f), 'F is the smallest in FHIST', srname)
    end if
end if

end subroutine bobyqa


end module bobyqa_mod
