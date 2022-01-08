module preproc_mod
!--------------------------------------------------------------------------------------------------!
! PREPROC_MOD is a module that preprocesses the inputs.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and papers.
!
! Started: July 2020
!
! Last Modified: Saturday, January 01, 2022 PM12:56:55
!--------------------------------------------------------------------------------------------------!

! N.B.: If all the inputs are valid, then PREPROC should do nothing.

implicit none
private
public :: preproc


contains


subroutine preproc(solver, n, iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, m, npt, maxfilt, &
        & ctol, eta1, eta2, gamma1, gamma2, is_constrained)
!--------------------------------------------------------------------------------------------------!
! This subroutine preprocesses the inputs. It does nothing to the inputs that are valid.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, TENTH, EPS, MAXMEMORY, MSGLEN, DEBUGGING
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, ETA1_DFT, ETA2_DFT, GAMMA1_DFT, GAMMA2_DFT
use, non_intrinsic :: consts_mod, only : CTOL_DFT, FTARGET_DFT, IPRINT_DFT, MIN_MAXFILT, MAXFILT_DFT
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: infnan_mod, only : is_nan, is_inf, is_finite
use, non_intrinsic :: memory_mod, only : cstyle_sizeof
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Compulsory inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: n

! Optional inputs
integer(IK), intent(in), optional :: m

! Compulsory in-outputs
integer(IK), intent(inout) :: iprint
integer(IK), intent(inout) :: maxfun
integer(IK), intent(inout) :: maxhist
real(RP), intent(inout) :: ftarget
real(RP), intent(inout) :: rhobeg
real(RP), intent(inout) :: rhoend

! Optional in-outputs
integer(IK), intent(inout), optional :: npt
integer(IK), intent(inout), optional :: maxfilt
logical, intent(in), optional :: is_constrained
real(RP), intent(inout), optional :: ctol
real(RP), intent(inout), optional :: eta1
real(RP), intent(inout), optional :: eta2
real(RP), intent(inout), optional :: gamma1
real(RP), intent(inout), optional :: gamma2

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! Format of integers; use the minimum number of digits
character(len=*), parameter :: rfmt = '(1PD15.6)'  ! Format of reals
character(len=*), parameter :: srname = 'PREPROC'
character(len=MSGLEN) :: min_maxfun_str
character(len=MSGLEN) :: wmsg
integer(IK) :: m_loc
integer(IK) :: maxfilt_in
integer(IK) :: min_maxfun
integer(IK) :: unit_memo
logical :: is_constrained_loc
real(RP) :: eta1_loc
real(RP) :: eta2_loc

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    if (present(m)) then
        call assert(m >= 0, 'M >= 0', srname)
        call assert(m == 0 .or. lower(solver) == 'cobyla', 'M == 0 unless the solver is COBYLA', srname)
    end if
    if (lower(solver) == 'cobyla' .and. present(m) .and. present(is_constrained)) then
        call assert(m == 0 .or. is_constrained, 'For COBYLA, M == 0 unless the problem is constrained', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! Read M, if necessary
if (lower(solver) == 'cobyla' .and. present(m)) then
    m_loc = m
else
    m_loc = 0_IK
end if

! Decide whether the problem is truly constrained
if (present(is_constrained)) then
    is_constrained_loc = is_constrained
else
    is_constrained_loc = (m_loc > 0)
end if


! Validate IPRINT
if (abs(iprint) > 3) then
    iprint = IPRINT_DFT
    write (wmsg, ifmt) iprint
    call warning(solver, 'Invalid IPRINT; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to '//trimstr(wmsg))
end if

! Validate MAXFUN
select case (lower(solver))
case ('uobyqa')
    min_maxfun = (n + 1_IK) * (n + 2_IK) / 2_IK + 1_IK
    min_maxfun_str = '(N+1)(N+2)/2 + 1'
case ('cobyla')
    min_maxfun = n + 2_IK
    min_maxfun_str = 'N + 2'
case default  ! CASE ('NEWUOA', 'BOBYQA', 'LINCOA')
    min_maxfun = n + 3_IK
    min_maxfun_str = 'N + 3'
end select
if (maxfun < min_maxfun) then
    maxfun = min_maxfun
    write (wmsg, ifmt) maxfun
    call warning(solver, 'Invalid MAXFUN; it should be at least '//trimstr(min_maxfun_str)//'; it is set to '//trimstr(wmsg))
end if

! Validate MAXHIST
if (maxhist < 0) then
    maxhist = maxfun
    write (wmsg, ifmt) maxhist
    call warning(solver, 'Invalid MAXHIST; it should be a nonnegative integer; it is set to '//trimstr(wmsg))
end if
maxhist = min(maxhist, maxfun)  ! MAXHIST > MAXFUN is never needed.

! Validate FTARGET
if (is_nan(ftarget)) then
    ftarget = FTARGET_DFT
    write (wmsg, rfmt) ftarget
    call warning(solver, 'Invalid FTARGET; it should be a real number; it is set to '//trimstr(wmsg))
end if

! Validate NPT
if ((lower(solver) == 'newuoa' .or. lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa') &
    & .and. present(npt)) then
    if (npt < n + 2 .or. npt > min(maxfun - 1, ((n + 2) * (n + 1)) / 2)) then
        npt = int(min(maxfun - 1, 2 * n + 1), kind(npt))
        write (wmsg, ifmt) npt
        call warning(solver, 'Invalid NPT; it should be an integer in the interval [N+2, (N+1)(N+2)/2], '// &
            & 'and it should be less than MAXFUN; it is set to '//trimstr(wmsg))
    end if
end if

! Validate MAXFILT
if (present(maxfilt) .and. (lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla')) then
    maxfilt_in = maxfilt
    if (maxfilt < 1) then
        maxfilt = MAXFILT_DFT  ! The inputted MAXFILT is obviously wrong.
    else
        maxfilt = max(MIN_MAXFILT, maxfilt)  ! The inputted MAXFILT is too small.
    end if
    ! Further revise MAXFILT according to MAXMEMORY.
    select case (lower(solver))
    case ('lincoa')
        unit_memo = (n + 2_IK) * cstyle_sizeof(0.0_RP)
    case ('cobyla')
        unit_memo = (m_loc + n + 2_IK) * cstyle_sizeof(0.0_RP)
    case default
        unit_memo = 1_IK
    end select
    ! We cannot simply set MAXFILT = MIN(MAXFILT, MAXMEMORY/...), as they may not have
    ! the same kind, and compilers may complain. We may convert them, but overflow may occur.
    if (maxfilt > MAXMEMORY / unit_memo) then
        maxfilt = int(MAXMEMORY / unit_memo, kind(maxfilt))
    end if
    maxfilt = min(maxfun, max(MIN_MAXFILT, maxfilt))
    if (is_constrained_loc) then
        write (wmsg, ifmt) maxfilt
        if (maxfilt_in < 1) then
            call warning(solver, 'Invalid MAXFILT; it should be a positive integer; it is set to ' &
                & //trimstr(wmsg))
        elseif (maxfilt_in < min(maxfun, MIN_MAXFILT)) then
            call warning(solver, 'MAXFILT is too small; it is set to '//trimstr(wmsg))
        elseif (maxfilt < min(maxfilt_in, maxfun)) then
            call warning(solver, 'MAXFILT is set to '//trimstr(wmsg)//' due to memory limit.')
        end if
    end if
end if

! Validate ETA1 and ETA2
if (present(eta1)) then
    eta1_loc = eta1
else
    eta1_loc = ETA1_DFT
end if
if (present(eta2)) then
    eta2_loc = eta2
else
    eta2_loc = ETA2_DFT
end if

! When the difference between ETA1 and ETA2 is tiny, we force them to equal.
! See the explanation around RHOBEG and RHOEND for the reason.
if (present(eta1) .and. present(eta2)) then
    if (abs(eta1 - eta2) < 1.0E2_RP * EPS * max(abs(eta1), ONE)) then
        eta2 = eta1
    end if
end if

if (present(eta1)) then
    if (is_nan(eta1)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        eta1 = ETA1_DFT
    elseif (eta1 < 0 .or. eta1 >= 1) then
        ! Take ETA2 into account if it has a valid value.
        if (present(eta2) .and. eta2_loc > 0 .and. eta2_loc <= 1) then
            eta1 = max(EPS, eta2 / 7.0_RP)
        else
            eta1 = ETA1_DFT
        end if
        write (wmsg, rfmt) eta1
        call warning(solver, 'Invalid ETA1; it should be in the interval [0, 1) and not more than ETA2;'// &
            & ' it is set to '//trimstr(wmsg))
    end if
end if

if (present(eta2)) then
    if (is_nan(eta2)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        eta2 = ETA2_DFT
    elseif (present(eta1) .and. (eta2 < eta1_loc .or. eta2 > 1)) then
        eta2 = (eta1 + TWO) / 3.0_RP
        write (wmsg, rfmt) eta2
        call warning(solver, 'Invalid ETA2; it should be in the interval [0, 1) and not less than ETA1;'// &
            & ' it is set to '//trimstr(wmsg))
    end if
end if

! Validate GAMMA1 and GAMMA2
if (present(gamma1)) then
    if (is_nan(gamma1)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        gamma1 = GAMMA1_DFT
    elseif (gamma1 <= 0 .or. gamma1 >= 1) then
        gamma1 = GAMMA1_DFT
        write (wmsg, rfmt) gamma1
        call warning(solver, 'Invalid GAMMA1; it should in the interval (0, 1); it is set to '//trimstr(wmsg))
    end if
end if

if (present(gamma2)) then
    if (is_nan(gamma2)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        gamma2 = GAMMA2_DFT
    elseif (gamma2 < 1 .or. is_inf(gamma2)) then
        gamma2 = GAMMA2_DFT
        write (wmsg, rfmt) gamma2
        call warning(solver, 'Invalid GAMMA2; it should be a real number not less than 1; it is set to '//trimstr(wmsg))
    end if
end if

! Validate RHOBEG and RHOEND
if (abs(rhobeg - rhoend) < 1.0E2_RP * EPS * max(abs(rhobeg), ONE)) then
    ! When the data is passed from the interfaces (e.g., MEX) to the Fortran code, RHOBEG, and RHOEND
    ! may change a bit. It was observed in a MATLAB test that MEX passed 1 to Fortran as
    ! 0.99999999999999978. Therefore, if we set RHOEND = RHOBEG in the interfaces, then it may happen
    ! that RHOEND > RHOBEG, which is considered as an invalid input. To avoid this situation, we
    ! force RHOBEG and RHOEND to equal when the difference is tiny.
    rhoend = rhobeg
end if

if (rhobeg <= 0 .or. is_nan(rhobeg) .or. is_inf(rhobeg)) then
    ! Take RHOEND into account if it has a valid value.
    if (is_finite(rhoend) .and. rhoend > 0) then
        rhobeg = max(TEN * rhoend, RHOBEG_DFT)
    else
        rhobeg = RHOBEG_DFT
    end if
    write (wmsg, rfmt) rhobeg
    call warning(solver, 'Invalid RHOBEG; it should be a positive number; it is set to '//trimstr(wmsg))
end if

if (rhoend <= 0 .or. rhobeg < rhoend .or. is_nan(rhoend) .or. is_inf(rhoend)) then
    rhoend = max(EPS, min(TENTH * rhobeg, RHOEND_DFT))
    write (wmsg, rfmt) rhoend
    call warning(solver, 'Invalid RHOEND; it should be a positive number and RHOEND <= RHOBEG; '// &
        & 'it is set to '//trimstr(wmsg))
end if

! Validate CTOL
if (present(ctol)) then
    if (is_nan(ctol) .or. ctol < 0) then
        ctol = CTOL_DFT
        if (is_constrained_loc) then
            write (wmsg, rfmt) ctol
            call warning(solver, 'Invalid CTOL; it should be a positive number; it is set to '//trimstr(wmsg))
        end if
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', solver)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', solver)
    if (present(npt)) then
        call assert(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', solver)
        call assert(npt >= 3, 'NPT >= 3', solver)
    end if
    if (present(maxfilt)) then
        call assert(maxfilt >= min(MIN_MAXFILT, maxfun) .and. maxfilt <= maxfun, &
            & 'MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN', solver)
    end if
    if (present(eta1) .and. present(eta2)) then
        call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', solver)
    end if
    if (present(gamma1) .and. present(gamma2)) then
        call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', solver)
    end if
    call assert(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', solver)
    if (present(ctol)) then
        call assert(ctol >= 0, 'CTOL >= 0', solver)
    end if
end if

end subroutine preproc


end module preproc_mod
