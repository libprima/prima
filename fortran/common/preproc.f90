module preproc_mod
!--------------------------------------------------------------------------------------------------!
! PREPROC_MOD is a module that preprocesses the inputs.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
!
! Started: July 2020
!
! Last Modified: Thursday, October 12, 2023 PM09:03:20
!--------------------------------------------------------------------------------------------------!

! N.B.:
! 1. If all the inputs are valid, then PREPROC should do nothing.
! 2. In PREPROC, we use VALIDATE instead of ASSERT, so that the parameters are validated even if we
! are not in debug mode.

implicit none
private
public :: preproc


contains


subroutine preproc(solver, n, iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, m, npt, maxfilt, &
        & ctol, cweight, eta1, eta2, gamma1, gamma2, is_constrained, has_rhobeg, honour_x0, xl, xu, x0)
!--------------------------------------------------------------------------------------------------!
! This subroutine preprocesses the inputs. It does nothing to the inputs that are valid.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, TENTH, HALF, EPS, MAXFUN_DIM_DFT, MAXHISTMEM, DEBUGGING
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, ETA1_DFT, ETA2_DFT, GAMMA1_DFT, GAMMA2_DFT
use, non_intrinsic :: consts_mod, only : CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT, IPRINT_DFT, MIN_MAXFILT, MAXFILT_DFT
use, non_intrinsic :: debug_mod, only : validate, warning
use, non_intrinsic :: infnan_mod, only : is_nan, is_inf, is_finite
use, non_intrinsic :: linalg_mod, only : trueloc, falseloc
use, non_intrinsic :: memory_mod, only : cstyle_sizeof
use, non_intrinsic :: string_mod, only : lower, num2str
implicit none

! Compulsory inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: n

! Optional inputs
integer(IK), intent(in), optional :: m
logical, intent(in), optional :: has_rhobeg
logical, intent(in), optional :: honour_x0
logical, intent(in), optional :: is_constrained
real(RP), intent(in), optional :: xl(:)
real(RP), intent(in), optional :: xu(:)

! Compulsory in-outputs
integer(IK), intent(inout) :: iprint
integer(IK), intent(inout) :: maxfun
integer(IK), intent(inout) :: maxhist
real(RP), intent(inout) :: eta1
real(RP), intent(inout) :: eta2
real(RP), intent(inout) :: ftarget
real(RP), intent(inout) :: gamma1
real(RP), intent(inout) :: gamma2
real(RP), intent(inout) :: rhobeg
real(RP), intent(inout) :: rhoend

! Optional in-outputs
integer(IK), intent(inout), optional :: npt
integer(IK), intent(inout), optional :: maxfilt
real(RP), intent(inout), optional :: ctol
real(RP), intent(inout), optional :: cweight
real(RP), intent(inout), optional :: x0(:)

! Local variables
character(len=*), parameter :: srname = 'PREPROC'
character(len=:), allocatable :: min_maxfun_str
integer(IK) :: m_loc
integer(IK) :: maxfilt_in
integer(IK) :: min_maxfun
integer(IK) :: unit_memo
logical :: is_constrained_loc
logical :: lbx(n)
logical :: ubx(n)
real(RP) :: rhobeg_default
real(RP) :: rhobeg_old
real(RP) :: rhoend_default
real(RP) :: x0_old(n)

! Preconditions
if (DEBUGGING) then
    call validate(n >= 1, 'N >= 1', srname)
    if (present(m)) then
        call validate(m >= 0, 'M >= 0', srname)
        call validate(m == 0 .or. lower(solver) == 'cobyla', 'M == 0 unless the solver is COBYLA', srname)
    end if
    if (lower(solver) == 'cobyla' .and. present(m) .and. present(is_constrained)) then
        call validate(m == 0 .or. is_constrained, 'For COBYLA, M == 0 unless the problem is constrained', srname)
    end if
    if (lower(solver) == 'bobyqa') then
        call validate(present(xl) .and. present(xu), 'XL and XU are present if the solver is BOBYQA', srname)
        call validate(all(xu - xl >= TWO * EPS), 'MINVAL(XU-XL) > 2*EPS', srname)
    end if
    if (present(honour_x0)) then
        call validate(lower(solver) == 'bobyqa' .and. present(has_rhobeg) .and. present(xl) .and. present(xu) &
            & .and. present(x0), 'If HONOUR_X0 is present, then so are XL, XU, and X0, and the solver is BOBYQA', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! Read M, if necessary
if (lower(solver) == 'cobyla' .and. present(m)) then
    m_loc = m
else
    m_loc = 0
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
    call warning(solver, 'Invalid IPRINT; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to '//num2str(iprint))
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
if (maxfun <= max(0_IK, min_maxfun)) then
    if (maxfun > 0) then
        maxfun = min_maxfun
    else  ! We assume that nonpositive values of MAXFUN are produced by overflow.
        if (MAXFUN_DIM_DFT >= huge(maxfun) / n) then  ! Avoid overflow when N is large.
            maxfun = huge(maxfun)
        else
            maxfun = MAXFUN_DIM_DFT * n
        end if
    end if
    call warning(solver, 'Invalid MAXFUN; it should be at least '//min_maxfun_str//'; it is set to '//num2str(maxfun))
end if

! Validate MAXHIST
if (maxhist < 0) then
    maxhist = maxfun
    call warning(solver, 'Invalid MAXHIST; it should be a nonnegative integer; it is set to '//num2str(maxhist))
end if
maxhist = min(maxhist, maxfun)  ! MAXHIST > MAXFUN is never needed.

! Validate FTARGET. We do not raise a warning if FTARGET is NaN.
if (is_nan(ftarget)) then
    ftarget = FTARGET_DFT
end if

! Validate NPT
if ((lower(solver) == 'newuoa' .or. lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa') &
    & .and. present(npt)) then
    if (npt < n + 2 .or. npt > min(maxfun - 1, ((n + 2) * (n + 1)) / 2)) then
        npt = int(min(maxfun - 1, 2 * n + 1), kind(npt))
        call warning(solver, 'Invalid NPT; it should be an integer in the interval [N+2, (N+1)(N+2)/2]'// &
            & ' and less than MAXFUN; it is set to '//num2str(npt))
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
    ! Further revise MAXFILT according to MAXHISTMEM.
    select case (lower(solver))
    case ('lincoa')
        unit_memo = (n + 2_IK) * cstyle_sizeof(0.0_RP)
    case ('cobyla')
        unit_memo = (m_loc + n + 2_IK) * cstyle_sizeof(0.0_RP)
    case default
        unit_memo = 1
    end select
    ! We cannot simply set MAXFILT = MIN(MAXFILT, MAXHISTMEM/...), as they may not have
    ! the same kind, and compilers may complain. We may convert them, but overflow may occur.
    if (maxfilt > MAXHISTMEM / unit_memo) then
        maxfilt = int(MAXHISTMEM / unit_memo, kind(maxfilt))  ! Integer division.
    end if
    maxfilt = min(maxfun, max(MIN_MAXFILT, maxfilt))
    if (is_constrained_loc) then
        if (maxfilt_in < 1) then
            call warning(solver, 'Invalid MAXFILT; it should be a positive integer; it is set to ' &
                & //num2str(maxfilt))
        elseif (maxfilt_in < min(maxfun, MIN_MAXFILT)) then
            call warning(solver, 'MAXFILT is too small; it is set to '//num2str(maxfilt))
        elseif (maxfilt < min(maxfilt_in, maxfun)) then
            call warning(solver, 'MAXFILT is set to '//num2str(maxfilt)//' due to memory limit')
        end if
    end if
end if

! Validate ETA1 and ETA2
if (is_nan(eta1) .or. eta1 < 0 .or. eta1 >= 1) then
    ! Take ETA2 into account if it has a valid value.
    if (eta2 >= 0 .and. eta2 < 1) then
        eta1 = eta2 / 7.0_RP
    else
        eta1 = ETA1_DFT
    end if
    call warning(solver, 'Invalid ETA1; it should be in the interval [0, 1) and not more than ETA2;'// &
        & ' it is set to '//num2str(eta1))
end if

if (is_nan(eta2) .or. eta2 < eta1 .or. eta2 < 0 .or. eta2 >= 1) then
    ! Take ETA1 into account if it has a valid value.
    if (eta1 >= 0 .and. eta1 < 1) then
        eta2 = (eta1 + TWO) / 3.0_RP
    else
        eta2 = ETA2_DFT
    end if
    call warning(solver, 'Invalid ETA2; it should be in the interval [0, 1) and not less than ETA1;'// &
        & ' it is set to '//num2str(eta2))
end if

! When the difference between ETA1 and ETA2 is tiny, we force them to equal.
! See the explanation around RHOBEG and RHOEND for the reason.
if (abs(eta1 - eta2) < 1.0E2_RP * EPS * max(abs(eta1), ONE)) then
    eta2 = eta1
end if

! Validate GAMMA1 and GAMMA2
if (is_nan(gamma1)) then
    ! In this case, we take the value hard coded in Powell's original code
    ! without any warning. It is useful when interfacing with MATLAB/Python.
    gamma1 = GAMMA1_DFT
elseif (gamma1 <= 0 .or. gamma1 >= 1) then
    gamma1 = GAMMA1_DFT
    call warning(solver, 'Invalid GAMMA1; it should in the interval (0, 1); it is set to '//num2str(gamma1))
end if

if (is_nan(gamma2)) then
    ! In this case, we take the value hard coded in Powell's original code
    ! without any warning. It is useful when interfacing with MATLAB/Python.
    gamma2 = GAMMA2_DFT
elseif (gamma2 < 1 .or. is_inf(gamma2)) then
    gamma2 = GAMMA2_DFT
    call warning(solver, 'Invalid GAMMA2; it should be a real number not less than 1; it is set to '//num2str(gamma2))
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

! Revise the default values for RHOBEG/RHOEND according to the solver.
if (lower(solver) == 'bobyqa') then
    rhobeg_default = max(EPS, min(RHOBEG_DFT, minval(xu - xl) / 4.0_RP))
    rhoend_default = max(EPS, min(TENTH * rhobeg_default, RHOEND_DFT))
else
    rhobeg_default = RHOBEG_DFT
    rhoend_default = RHOEND_DFT
end if

if (lower(solver) == 'bobyqa') then
    ! Do NOT merge the IF below into the ELSEIF above! Otherwise, XU and XL may be accessed even if
    ! the solver is not BOBYQA, because the logical evaluation is not short-circuit.
    if (rhobeg > minval(xu - xl) / TWO) then
        ! Do NOT make this revision if RHOBEG not positive or not finite, because otherwise RHOBEG
        ! will get a huge value when XU or XL contains huge values that indicate unbounded variables.
        rhobeg = minval(xu - xl) / 4.0_RP  ! Here, we do not take RHOBEG_DEFAULT.
        call warning(solver, 'Invalid RHOBEG; '//solver//' requires 0 < RHOBEG <= MINVAL(XU-XL)/2;' &
            & //' it is set to MINVAL(XU-XL)/4')
    end if
end if
if (rhobeg <= 0 .or. is_nan(rhobeg) .or. is_inf(rhobeg)) then
    ! Take RHOEND into account if it has a valid value. We do not do this if the solver is BOBYQA,
    ! which requires that RHOBEG <= (XU-XL)/2.
    if (is_finite(rhoend) .and. rhoend > 0 .and. lower(solver) /= 'bobyqa') then
        rhobeg = max(TEN * rhoend, rhobeg_default)
    else
        rhobeg = rhobeg_default
    end if
    call warning(solver, 'Invalid RHOBEG; it should be a positive number; it is set to '//num2str(rhobeg))
end if

if (rhoend <= 0 .or. rhobeg < rhoend .or. is_nan(rhoend) .or. is_inf(rhoend)) then
    rhoend = max(EPS, min(TENTH * rhobeg, rhoend_default))
    call warning(solver, 'Invalid RHOEND; it should be a positive number and RHOEND <= RHOBEG; '// &
        & 'it is set to '//num2str(rhoend))
end if

! For BOBYQA, revise X0 or RHOBEG so that the distance between X0 and the inactive bounds is at
! least RHOBEG. If HONOUR_X0 == TRUE, revise RHOBEG if needed; otherwise, revise HONOUR_X0 if needed.
if (present(honour_x0)) then
    if (honour_x0) then
        rhobeg_old = rhobeg; 
        lbx = (is_finite(xl) .and. x0 - xl <= EPS * max(ONE, abs(xl))) ! X0 essentially equals XL
        ubx = (is_finite(xu) .and. x0 - xu >= -EPS * max(ONE, abs(xu))) ! X0 essentially equals XU
        x0(trueloc(lbx)) = xl(trueloc(lbx))
        x0(trueloc(ubx)) = xu(trueloc(ubx))
        rhobeg = max(EPS, minval([rhobeg, x0(falseloc(lbx)) - xl(falseloc(lbx)), xu(falseloc(ubx)) - x0(falseloc(ubx))]))
        if (rhobeg_old - rhobeg > EPS * max(ONE, rhobeg_old)) then
            rhoend = max(EPS, min(TENTH * rhobeg, rhoend)) ! We do not revise RHOEND unless RHOBEG is truly revised.
            if (has_rhobeg) then
                call warning(solver, 'RHOBEG is revised to '//num2str(rhobeg)//' and RHOEND to at most 0.1*RHOBEG'// &
                    & ' so that the distance between X0 and the inactive bounds is at least RHOBEG')
            end if
        else
            rhoend = min(rhoend, rhobeg)  ! This may update RHOEND slightly.
        end if
    else
        x0_old = x0  ! Recorded to see whether X0 is really revised.
        ! N.B.: The following revision is valid only if XL <= X0 <= XU and RHOBEG <= MINVAL(XU-XL)/2,
        ! which should hold at this point due to the revision of RHOBEG and moderation of X0.
        ! The cases below are mutually exclusive in precise arithmetic as MINVAL(XU-XL) >= 2*RHOBEG.
        where (x0 <= xl + HALF * rhobeg)
            x0 = xl
        elsewhere(x0 < xl + rhobeg)
            x0 = xl + rhobeg
        end where
        where (x0 >= xu - HALF * rhobeg)
            x0 = xu
        elsewhere(x0 > xu - rhobeg)
            x0 = xu - rhobeg
        end where
        !!MATLAB code:
        !!lbx = (x0 <= xl + 0.5 * rhobeg);
        !!lbx_plus = (x0 > xl + 0.5 * rhobeg .and. x0 < xl + rhobeg);
        !!ubx = (x0 >= xu - 0.5 * rhobeg);
        !!ubx_minus = (x0 < xu - 0.5 * rhobeg .and. x0 > xu - rhobeg);
        !!x0(lbx) = xl(lbx);
        !!x0(lbx_plus) = xl(lbx_plus) + rhobeg;
        !!x0(ubx) = xu(ubx);
        !!x0(ubx_minus) = xu(ubx_minus) - rhobeg;

        if (any(abs(x0_old - x0) > 0)) then
            call warning(solver, 'X0 is revised so that the distance between X0 and the inactive bounds is at least RHOBEG; '// &
                  & 'set HONOUR_X0 to .TRUE. if you prefer to keep X0 unchanged')
        end if
    end if
end if

! Validate CTOL (it can be 0)
if (present(ctol)) then
    if (is_nan(ctol) .or. ctol < 0) then
        ctol = CTOL_DFT
        if (is_constrained_loc) then
            call warning(solver, 'Invalid CTOL; it should be a nonnegative number; it is set to '//num2str(ctol))
        end if
    end if
end if

! Validate CWEIGHT (it can be +Inf)
if (present(cweight)) then
    if (is_nan(cweight) .or. cweight < 0) then
        cweight = CWEIGHT_DFT
        if (is_constrained_loc) then
            call warning(solver, 'Invalid CWEIGHT; it should be a nonnegative number; it is set to '//num2str(cweight))
        end if
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call validate(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', solver)
    call validate(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', solver)
    if (present(npt)) then
        call validate(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', solver)
        call validate(npt >= 3, 'NPT >= 3', solver)
    end if
    if (present(maxfilt)) then
        call validate(maxfilt >= min(MIN_MAXFILT, maxfun) .and. maxfilt <= maxfun, &
            & 'MIN(MIN_MAXFILT, MAXFUN) <= MAXFILT <= MAXFUN', solver)
    end if
    call validate(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', solver)
    call validate(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', solver)
    call validate(rhobeg >= rhoend .and. rhoend > 0, 'RHOBEG >= RHOEND > 0', solver)
    if (lower(solver) == 'bobyqa') then
        call validate(all(rhobeg <= (xu - xl) / TWO), 'RHOBEG <= MINVAL(XU-XL)/2', solver)
        call validate(all(is_finite(x0)), 'X0 is finite', solver)
        call validate(all(x0 >= xl .and. (x0 <= xl .or. x0 >= xl + rhobeg)), 'X0 == XL or X0 >= XL + RHOBEG', solver)
        call validate(all(x0 <= xu .and. (x0 >= xu .or. x0 <= xu - rhobeg)), 'X0 == XU or X0 >= XU - RHOBEG', solver)
    end if
    if (present(ctol)) then
        call validate(ctol >= 0, 'CTOL >= 0', solver)
    end if
end if

end subroutine preproc


end module preproc_mod
