module preproc_mod
!--------------------------------------------------------------------------------------------------!
! PREPROC_MOD is a module that preprocesses the inputs.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
!
! Started: July 2020
!
! Last Modified: Wed 10 Sep 2025 02:16:58 AM CST
!--------------------------------------------------------------------------------------------------!

! N.B.:
! 1. If all the inputs are valid, then PREPROC should do nothing.
! 2. In PREPROC, we use VALIDATE instead of ASSERT, so that the parameters are validated even if we
!    are not in debug mode.

implicit none
private
public :: preproc


contains


subroutine preproc(solver, n, iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, m, npt, maxfilt, &
        & ctol, cweight, eta1, eta2, gamma1, gamma2, is_constrained, has_rhobeg, honour_x0, xl, xu, x0)
!--------------------------------------------------------------------------------------------------!
! This subroutine preprocesses the inputs. It does nothing to the inputs that are valid.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, HALF, EPS, MAXHISTMEM, DEBUGGING
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, ETA1_DFT, ETA2_DFT, GAMMA1_DFT, GAMMA2_DFT
use, non_intrinsic :: consts_mod, only : CTOL_DFT, CWEIGHT_DFT, IPRINT_DFT, MIN_MAXFILT, MAXFILT_DFT
use, non_intrinsic :: debug_mod, only : validate, warning
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
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
integer :: min_maxfun  ! INTEGER(IK) may overflow if IK corresponds to the 16-bit integer.
integer :: unit_memo  ! INTEGER(IK) may overflow if IK corresponds to the 16-bit integer.
integer(IK) :: iprint_in
integer(IK) :: m_loc
integer(IK) :: maxfilt_in
integer(IK) :: maxfun_in
integer(IK) :: maxhist_in
integer(IK) :: npt_in
logical :: is_constrained_loc
logical :: lbx(n)
logical :: ubx(n)
real(RP) :: ctol_in
real(RP) :: cweight_in
real(RP) :: eta1_in
real(RP) :: eta2_in
real(RP) :: gamma1_in
real(RP) :: gamma2_in
real(RP) :: rhobeg_default
real(RP) :: rhobeg_in
real(RP) :: rhoend_default
real(RP) :: rhoend_in
real(RP) :: x0_in(n)

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
    call validate(present(maxfilt) .eqv. (lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla'), &
        & 'MAXFILT is present if and only if the solver is LINCOA or COBYLA', srname)
    if (lower(solver) == 'bobyqa') then
        call validate(present(xl) .and. present(xu), 'XL and XU are present if the solver is BOBYQA', srname)
        call validate(all(xu - xl >= TWO * EPS), 'MINVAL(XU-XL) > 2*EPS', srname)
    end if
    call validate((present(honour_x0) .eqv. present(x0)) .and. (present(honour_x0) .eqv. present(has_rhobeg)), &
        & 'HONOUR_X0, X0, and HAS_RHOBEG are present or absent simultaneously', srname)
    call validate(present(honour_x0) .eqv. lower(solver) == 'bobyqa', &
        & 'HONOUR_X0 is present if and only if the solver is BOBYQA', srname)
    ! N.B.: LINCOA and COBYLA will have HONOUR_X0 as well if we intend to make them respect bounds.
    ! !call validate(present(honour_x0) .eqv. &
    ! !    & (lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla'), &
    ! !    & 'HONOUR_X0 is present if and only if the solver is BOBYQA, LINCOA, or COBYLA', srname)
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
    iprint_in =  iprint
    iprint = IPRINT_DFT
    call warning(solver, 'Invalid IPRINT: '//num2str(iprint_in)// &
        & '; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to '// num2str(iprint))
end if

! Validate MAXFUN
! N.B.: The INT(N), INT(N+1), and INT(N+2) below convert integers to the default integer kind,
! which is the kind of MIN_MAXFUN. Fortran compilers may complain without the conversion. It is
! not needed in Python/MATLAB/Julia/R.
select case (lower(solver))
case ('uobyqa')
    min_maxfun = (int(n + 1) * int(n + 2)) / 2 + 1  ! INT(*) avoids overflow when IK is 16-bit.
    min_maxfun_str = '(N+1)(N+2)/2 + 1'
case ('cobyla')
    min_maxfun = int(n) + 2
    min_maxfun_str = 'N + 2'
case default  ! CASE ('NEWUOA', 'BOBYQA', 'LINCOA')
    min_maxfun = int(n) + 3
    min_maxfun_str = 'N + 3'
end select
if (maxfun <= max(0, min_maxfun - 1)) then
    maxfun_in = maxfun
    if (maxfun > 0) then
        maxfun = int(min_maxfun, kind(maxfun))
    else  ! We assume that non-positive values of MAXFUN are produced by overflow.
        maxfun = int(max(min_maxfun, 10**min(range(maxfun), 5)), kind(maxfun))  !!MATLAB: maxfun =  max(min_maxfun, 10^5);
        ! N.B.: Do NOT set MAXFUN to HUGE(MAXFUN), as it may cause overflow and infinite cycling
        ! when used as the upper bound of DO loops. This occurred on 20240225 with gfortran 13. See
        ! https://fortran-lang.discourse.group/t/loop-variable-reaching-integer-huge-causes-infinite-loop
        ! https://fortran-lang.discourse.group/t/loops-dont-behave-like-they-should
    end if
    call warning(solver, 'Invalid MAXFUN: '//num2str(maxfun_in)// &
        & '; it should be at least '//min_maxfun_str//' with N = '//num2str(n)//'; it is set to '//num2str(maxfun))
end if

! Validate MAXHIST
if (maxhist <= 0) then
    maxhist_in = maxhist
    maxhist = maxfun
    call warning(solver, 'Invalid MAXHIST: '//num2str(maxhist_in)// &
        & '; it should be a positive integer; it is set to '//num2str(maxhist))
end if
maxhist = min(maxhist, maxfun)  ! MAXHIST > MAXFUN is never needed.

! Validate FTARGET
if (is_nan(ftarget)) then  ! No warning if FTARGET is NaN, which is interpreted as no target function value is provided.
    ftarget = -huge(ftarget)
end if

! Validate NPT
if ((lower(solver) == 'newuoa' .or. lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa') &
    & .and. present(npt)) then
    if (npt < n + 2 .or. npt >= maxfun .or. 2 * npt > int(n + 2) * int(n + 1)) then  ! INT(*) avoids overflow when IK is 16-bit.
        npt_in = npt
        npt = int(min(maxfun - 1, 2 * n + 1), kind(npt))
        call warning(solver, 'Invalid NPT: '//num2str(npt_in)// &
            & '; it should be an integer in the interval [N+2, (N+1)(N+2)/2] with N = '//num2str(n)// &
            & ' and less than MAXFUN = '//num2str(maxfun)//'; it is set to '//num2str(npt))
    end if
end if

! Validate MAXFILT
if (present(maxfilt) .and. (lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla')) then
    maxfilt_in = maxfilt
    if (maxfilt <= 0) then
        maxfilt = MAXFILT_DFT
    else
        maxfilt = max(MIN_MAXFILT, maxfilt)  ! The inputted MAXFILT is too small.
    end if
    ! Further revise MAXFILT according to MAXHISTMEM.
    select case (lower(solver))
    case ('lincoa')
        unit_memo = int(n + 2) * int(cstyle_sizeof(0.0_RP))  ! INT(*) avoids overflow when IK is 16-bit.
    case ('cobyla')
        unit_memo = int(m_loc + n + 2) * int(cstyle_sizeof(0.0_RP))  ! INT(*) avoids overflow when IK is 16-bit.
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
        if (maxfilt_in <= 0) then
            call warning(solver, 'Invalid MAXFILT: '//num2str(maxfilt_in)// &
                & '; it should be a positive integer; it is set to '//num2str(maxfilt))
        elseif (maxfilt_in < min(maxfun, MIN_MAXFILT)) then
            call warning(solver, 'MAXFILT = '//num2str(maxfilt_in)//' is too small; it is set to '//num2str(maxfilt))
        elseif (maxfilt < min(maxfilt_in, maxfun)) then
            call warning(solver, 'MAXFILT is reduced from '//num2str(maxfilt_in)//' to '//num2str(maxfilt)//' due to memory limit')
        end if
    end if
end if

! Validate ETA1 and ETA2
if (.not. (eta1 >= 0 .and. eta1 < 1)) then  ! ETA1 = NaN falls into this case.
    eta1_in = eta1
    ! Take ETA2 into account if it has a valid value.
    if (eta2 >= 0 .and. eta2 < 1) then
        eta1 = eta2 / 7.0_RP
    else
        eta1 = ETA1_DFT
    end if
    call warning(solver, 'Invalid ETA1: '//num2str(eta1_in)// &
        & '; it should be in the interval [0, 1) and not more than ETA2 = '//num2str(eta2)//'; it is set to '//num2str(eta1))
end if

if (.not. (eta2 >= eta1 .and. eta2 < 1)) then  ! ETA2 = NaN falls into this case.
    eta2_in = eta2
    ! Take ETA1 into account if it has a valid value.
    if (eta1 >= 0 .and. eta1 < 1) then
        eta2 = (eta1 + TWO) / 3.0_RP
    else
        eta2 = ETA2_DFT
    end if
    call warning(solver, 'Invalid ETA2: '//num2str(eta2_in)// &
        & '; it should be in the interval [0, 1) and not less than ETA1 = '//num2str(eta1)//'; it is set to '//num2str(eta2))
end if

! The following revision may update ETA1 slightly. It prevents ETA1 > ETA2 due to rounding
! errors, which would not be accepted by the solvers.
eta1 = min(eta1, eta2)

! Validate GAMMA1 and GAMMA2
if (.not. (gamma1 > 0 .and. gamma1 < 1)) then  ! GAMMA1 = NaN falls into this case.
    gamma1_in = gamma1
    gamma1 = GAMMA1_DFT
    call warning(solver, 'Invalid GAMMA1: '//num2str(gamma1_in)// &
        & '; it should in the interval (0, 1); it is set to '//num2str(gamma1))
end if

if (.not. (is_finite(gamma2) .and. gamma2 >= 1)) then  ! GAMMA2 = NaN falls into this case.
    gamma2_in = gamma2
    gamma2 = GAMMA2_DFT
    call warning(solver, 'Invalid GAMMA2: '//num2str(gamma2_in)// &
        & '; it should be a real number not less than 1; it is set to '//num2str(gamma2))
end if

! Validate RHOBEG and RHOEND

rhobeg_in = rhobeg
rhoend_in = rhoend

! Revise the default values for RHOBEG/RHOEND according to the solver.
if (lower(solver) == 'bobyqa') then
    rhobeg_default = max(EPS, min(RHOBEG_DFT, minval(xu - xl) / 4.0_RP))
    rhoend_default = max(EPS, min((RHOEND_DFT / RHOBEG_DFT) * rhobeg_default, RHOEND_DFT))
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
        call warning(solver, 'Invalid RHOBEG: '//num2str(rhobeg_in)// &
            & '; '//solver//' requires 0 < RHOBEG <= MINVAL(XU-XL)/2 = '//num2str(minval(xu - xl) / 2.0_RP)// &
            & '; it is set to '//num2str(rhobeg))
    end if
end if

if (.not. (is_finite(rhobeg) .and. rhobeg > 0)) then  ! RHOBEG = NaN falls into this case.
    ! Take RHOEND into account if it has a valid value. We do not do this if the solver is BOBYQA,
    ! which requires that RHOBEG <= (XU-XL)/2.
    if (is_finite(rhoend) .and. rhoend > 0 .and. lower(solver) /= 'bobyqa') then
        rhobeg = max(TEN * rhoend, rhobeg_default)
    else
        rhobeg = rhobeg_default
    end if
    call warning(solver, 'Invalid RHOBEG: '//num2str(rhobeg_in)// &
        & '; it should be a positive number; it is set to '//num2str(rhobeg))
end if

if (.not. (is_finite(rhoend) .and. rhoend >= 0 .and. rhoend <= rhobeg)) then  ! RHOEND = NaN falls into this case.
    rhoend = max(EPS, min((RHOEND_DFT / RHOBEG_DFT) * rhobeg, rhoend_default))
    call warning(solver, 'Invalid RHOEND: '//num2str(rhoend_in)// &
        & '; we should have '//num2str(rhobeg)//' = RHOBEG >= RHOEND >= 0; it is set to '//num2str(rhoend))
end if

! For BOBYQA, revise X0 or RHOBEG so that the distance between X0 and the inactive bounds is at
! least RHOBEG. If HONOUR_X0 == FALSE, revise X0 if needed; then revise RHOBEG if needed.
! N.B.: We should do the same for LINCOA and COBYLA if we make them respect the bounds in the future.
! !if (lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa' .or. lower(solver) == 'cobyla') then
if (lower(solver) == 'bobyqa') then
    ! Revise X0 if allowed and needed.
    if (.not. honour_x0) then
        x0_in = x0  ! Recorded to see whether X0 is really revised.
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

        if (any(abs(x0_in - x0) > 0)) then
            call warning(solver, 'X0 is revised so that the distance between X0 and the inactive bounds is at least RHOBEG = '// &
                & num2str(rhobeg)//'; revise RHOBEG or set HONOUR_X0 to .TRUE. if you prefer to keep X0 unchanged')
        end if
    end if

    ! Revise RHOBEG if needed.
    ! N.B.: If X0 has been revised above (i.e., HONOUR_X0 is FALSE), then the following revision
    ! is unnecessary in precise arithmetic. However, it may still be needed due to rounding errors.
    lbx = (is_finite(xl) .and. x0 - xl <= EPS * max(ONE, abs(xl))) ! X0 essentially equals XL
    ubx = (is_finite(xu) .and. x0 - xu >= -EPS * max(ONE, abs(xu))) ! X0 essentially equals XU
    x0(trueloc(lbx)) = xl(trueloc(lbx))
    x0(trueloc(ubx)) = xu(trueloc(ubx))
    rhobeg = max(EPS, minval([rhobeg, x0(falseloc(lbx)) - xl(falseloc(lbx)), xu(falseloc(ubx)) - x0(falseloc(ubx))]))
    if (rhobeg_in - rhobeg > EPS * max(ONE, rhobeg_in)) then
        rhoend = max(EPS, min((rhoend / rhobeg_in) * rhobeg, rhoend)) ! We do not revise RHOEND unless RHOBEG is truly revised.
        if (has_rhobeg) then
            call warning(solver, 'RHOBEG is revised from '//num2str(rhobeg_in)//' to '//num2str(rhobeg)// &
                & ' and RHOEND from '//num2str(rhoend_in)//' to '//num2str(rhoend)// &
                & ' so that the distance between X0 and the inactive bounds is at least RHOBEG')
        end if
    end if
end if

! The following revision may update RHOBEG and RHOEND slightly. It particularly prevents
! RHOEND > RHOBEG due to rounding errors, which would not be accepted by the solvers.
rhobeg = max(rhobeg, EPS)
rhoend = min(max(rhoend, EPS), rhobeg)

! Validate CTOL (it can be 0)
if (present(ctol)) then
    if (.not. (ctol >= 0)) then  ! CTOL = NaN falls into this case.
        ctol_in = ctol
        ctol = CTOL_DFT
        if (is_constrained_loc) then
            call warning(solver, 'Invalid CTOL: '//num2str(ctol_in)// &
                & '; it should be a nonnegative number; it is set to '//num2str(ctol))
        end if
    end if
end if

! Validate CWEIGHT (it can be +Inf)
if (present(cweight)) then
    if (.not. (cweight >= 0)) then  ! CWEIGHT = NaN falls into this case.
        cweight_in = cweight
        cweight = CWEIGHT_DFT
        if (is_constrained_loc) then
            call warning(solver, 'Invalid CWEIGHT: '//num2str(cweight_in)// &
                & '; it should be a nonnegative number; it is set to '//num2str(cweight))
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
    call validate(maxfun >= min_maxfun, 'MAXFUN >= MIN_MAXFUN', solver)
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
        call validate(all(x0 >= xl .and. (x0 <= xl .or. x0 - xl >= rhobeg)), 'X0 == XL or X0 - XL >= RHOBEG', solver)
        call validate(all(x0 <= xu .and. (x0 >= xu .or. xu - x0 >= rhobeg)), 'X0 == XU or XU - X0 >= RHOBEG', solver)
    end if
    if (present(ctol)) then
        call validate(ctol >= 0, 'CTOL >= 0', solver)
    end if
end if

end subroutine preproc


end module preproc_mod
