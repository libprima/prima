module cobyla_mod
!--------------------------------------------------------------------------------------------------!
! Classical mode. Not maintained. Not recommended. Please use the modernized version instead.
!
! The usage is the same as the modernized version.
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobyla


contains


subroutine cobyla(calcfc, m, x, f, &
    & cstrv, constr, &
    & f0, constr0, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, &
    & eta1, eta2, gamma1, gamma2, xhist, fhist, chist, conhist, maxhist, maxfilt, info)

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, errstop, warning
use, non_intrinsic :: evaluate_mod, only : evaluate, moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJCON
!use, non_intrinsic :: selectx_mod, only : isbetter
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
real(RP), allocatable :: constr_loc(:)
real(RP), allocatable :: chist_loc(:)
real(RP), allocatable :: conhist_loc(:, :)
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

!! If the user provides the function & constraint value at X0, then set up F_X0 and CONSTR_X0.
!if (present(f0) .and. present(constr0)) then
!    fc_x0_provided = .true.
!    !--------------------------------------------------!
!    call safealloc(x0, n)  ! Removable in F2003.
!    call safealloc(constr_x0, m)  ! Removable in F2003.
!    !--------------------------------------------------!
!    x0 = x
!    f_x0 = f0
!    constr_x0 = constr0
!else
!    call evaluate(calcfc, x, f_x0, constr_x0)
!end if

if (present(f0) .and. present(constr0) .and. all(is_finite(x))) then
    f = f0
    constr_loc = constr0
    cstrv_loc = maxval([ZERO, -constr_loc])
else
    ! Replace any NaN in X by ZERO and Inf/-Inf in X by HUGENUM/-HUGENUM.
    x = moderatex(x)
    call evaluate(calcfc, x, f, constr_loc, cstrv_loc)
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

!--------------------------------------------------------------------------------------------------!
!-------------------- Call COBYLB, which performs the real calculations. --------------------------!
!!!! ETA1, ETA2, GAMAA1, GAMMA2, MAXFILT, CTOL, CWEIGHT are not used in the classical mode. !!!!
call cobylb(calcfc, iprint_loc, maxfun_loc, rhobeg_loc, rhoend_loc, constr_loc, x, cstrv_loc, f, info_loc, &
    & nf_loc, xhist_loc, fhist_loc, chist_loc, conhist_loc)
!--------------------------------------------------------------------------------------------------!
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
    call warning(solver, 'Only the history of the last '//trim(wmsg)//' iteration(s) is recorded')
end if

end subroutine cobyla


end module cobyla_mod
