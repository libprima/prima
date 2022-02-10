module newuoa_mod
!--------------------------------------------------------------------------------------------------!
! Classical mode. Not maintained. Not recommended. Please use the modernized version instead.
!
! The usage is the same as the modernized version.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code.
!
! Started: January 2022
!
! Last Modified: Thursday, February 10, 2022 PM05:19:33
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: newuoa


contains


subroutine newuoa(calfun, x, f, &
    & nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, &
    & eta1, eta2, gamma1, gamma2, xhist, fhist, maxhist, info)

! Generic modules
!use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: evaluate_mod, only : moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_inf, is_finite
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: preproc_mod, only : preproc

! Solver-specific modules
use, non_intrinsic :: newuob_mod, only : newuob

implicit none

! Arguments
procedure(OBJ) :: calfun
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
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
character(len=*), parameter :: solver = 'NEWUOA'
!character(len=*), parameter :: srname = 'NEWUOA'
character(len=MSGLEN) :: wmsg
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK) :: npt_loc
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

! Replace any NaN in X by ZERO and Inf/-Inf in X by HUGENUM/-HUGENUM.
x = moderatex(x)

! Read the inputs.

! If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
! value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
! it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
if (present(rhobeg)) then
    rhobeg_loc = rhobeg
elseif (present(rhoend)) then
    ! Fortran does not take short-circuit evaluation of logic expressions. Thus it is WRONG to
    ! combine the evaluation of PRESENT(RHOEND) and the evaluation of IS_FINITE(RHOEND) as
    ! "if (present(rhoend) .and. is_finite(rhoend))". The compiler may choose to evaluate the
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

! Preprocess the inputs in case some of them are invalid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & npt=npt_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, gamma2=gamma2_loc)

! Further revise MAXHIST_LOC according to MAXMEMORY, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN) if they are requested; replace MAXFUN with 0 for
! the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc)

!--------------------------------------------------------------------------------------------------!
!-------------------- Call NEWUOB, which performs the real calculations. --------------------------!
call newuob(calfun, iprint_loc, maxfun_loc, npt_loc, ftarget_loc, rhobeg_loc, rhoend_loc, x, f, info_loc, &
    & nf_loc, xhist_loc, fhist_loc)
!--------------------------------------------------------------------------------------------------!
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

! If MAXFHIST_IN >= NF_LOC > MAXFHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist)) .and. maxhist_loc < nf_loc) then
    write (wmsg, ifmt) maxhist_loc
    call warning(solver, 'Only the history of the last '//trim(wmsg)//' iteration(s) is recoreded')
end if

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! The classical mode cannot pass the tests, of course.
!! Postconditions
!if (DEBUGGING) then
!    call assert(nf_loc <= maxfun_loc, 'NF <= MAXFUN', srname)
!    call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
!    nhist = min(nf_loc, maxhist_loc)
!    if (present(xhist)) then
!        call assert(size(xhist, 1) == n .and. size(xhist, 2) == nhist, 'SIZE(XHIST) == [N, NHIST]', srname)
!        call assert(.not. any(is_nan(xhist)), 'XHIST does not contain NaN', srname)
!    end if
!    if (present(fhist)) then
!        call assert(size(fhist) == nhist, 'SIZE(FHIST) == NHIST', srname)
!        call assert(.not. any(fhist < f), 'F is the smallest in FHIST', srname)
!    end if
!end if
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


end subroutine newuoa


end module newuoa_mod
