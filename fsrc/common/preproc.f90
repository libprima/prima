module preproc_mod
!--------------------------------------------------------------------------------------------------!
! PREPROC_MOD is a module that preprocesses the inputs.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and papers.
!
! Started: July 2020
!
! Last Modified: Saturday, December 18, 2021 AM12:00:20
!--------------------------------------------------------------------------------------------------!

! N.B.: If all the inputs are valid, then PREPROC should do nothing.

implicit none
private
public :: preproc


contains


subroutine preproc(solver, n, iprint, maxfun, maxhist, ftarget, rhobeg, rhoend, npt, ctol, eta1, eta2, gamma1, gamma2)

use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, TENTH, EPS, DEBUGGING
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, ETA1_DFT, ETA2_DFT, GAMMA1_DFT, GAMMA2_DFT
use, non_intrinsic :: consts_mod, only : CTOL_DFT, FTARGET_DFT, IPRINT_DFT
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_inf, is_finite
use, non_intrinsic :: string_mod, only : lower

implicit none

! Inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: n

! Compulsory in-outputs
integer(IK), intent(inout) :: iprint
integer(IK), intent(inout) :: maxfun
integer(IK), intent(inout) :: maxhist
real(RP), intent(inout) :: ftarget
real(RP), intent(inout) :: rhobeg
real(RP), intent(inout) :: rhoend

! Optional in-outputs
integer(IK), intent(inout), optional :: npt
real(RP), intent(inout), optional :: ctol
real(RP), intent(inout), optional :: eta1
real(RP), intent(inout), optional :: eta2
real(RP), intent(inout), optional :: gamma1
real(RP), intent(inout), optional :: gamma2

! Local variables
integer(IK) :: min_maxfun
real(RP) :: eta1_loc
real(RP) :: eta2_loc


if (iprint /= 0 .and. abs(iprint) /= 1 .and. abs(iprint) /= 2 .and. abs(iprint) /= 3) then
    iprint = IPRINT_DFT
    print '(/1A, I2, 1A)', solver//': invalid IPRINT; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to ', iprint, '.'
end if

if (lower(solver) == 'newuoa' .or. lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa') then
    min_maxfun = n + 3_IK
else if (lower(solver) == 'uobyqa') then
    min_maxfun = (n + 1_IK) * (n + 2_IK) / 2_IK + 1_IK
else
    min_maxfun = n + 2_IK
end if
if (maxfun < min_maxfun) then
    maxfun = min_maxfun
    print '(/1A, I8, 1A)', solver//': invalid MAXFUN; it should an integer at least N + 3 ; it is set to ', maxfun, '.'
end if

if (maxhist < 0) then
    maxhist = maxfun
    print '(/1A, I8, 1A)', solver//': invalid MAXHIST; it should be a nonnegative integer; it is set to ', maxhist, '.'
end if
maxhist = min(maxhist, maxfun)  ! MAXHIST > MAXFUN is never needed.

if (is_nan(ftarget)) then
    ftarget = FTARGET_DFT
    print '(/1A, 1PD15.6, 1A)', solver//': invalid FTARGET; it should a real number; it is set to ', ftarget, '.'
end if

if ((lower(solver) == 'newuoa' .or. lower(solver) == 'bobyqa' .or. lower(solver) == 'lincoa') &
    & .and. present(npt)) then
    if (npt < n + 2 .or. npt > min(maxfun - 1, ((n + 2) * (n + 1)) / 2)) then
        npt = int(min(maxfun - 1, 2 * n + 1), kind(npt))
        print '(/1A, I6, 1A)', solver//': invalid NPT; it should an integer in the interval [N+2, (N+1)(N+2)/2], '// &
            & 'and it should be less than MAXFUN; it is set to ', npt, '.'
    end if
end if

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
    else if (eta1 < 0 .or. eta1 >= 1) then
        ! Take ETA2 into account if it has a valid value.
        if (present(eta2) .and. eta2_loc > 0 .and. eta2_loc <= 1) then
            eta1 = max(EPS, eta2 / 7.0_RP)
        else
            eta1 = ETA1_DFT
        end if
        print '(/1A, 1PD15.6, 1A)', solver//': invalid ETA1; it should be in the interval [0, 1) and not more than ETA2;'// &
            & ' it is set to ', eta1, '.'
    end if
end if

if (present(eta2)) then
    if (is_nan(eta2)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        eta2 = ETA2_DFT
    else if (present(eta1) .and. (eta2 < eta1_loc .or. eta2 > 1)) then
        eta2 = (eta1 + TWO) / 3.0_RP
        print '(/1A, 1PD15.6, 1A)', solver//': invalid ETA2; it should be in the interval [0, 1] and not less than ETA1;'// &
            & ' it is set to ', eta2, '.'
    end if
end if

if (present(gamma1)) then
    if (is_nan(gamma1)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        gamma1 = GAMMA1_DFT
    else if (gamma1 <= 0 .or. gamma1 >= 1) then
        gamma1 = GAMMA1_DFT
        print '(/1A, 1PD15.6, 1A)', solver//': invalid GAMMA1; it should in the interval (0, 1); it is set to ', gamma1, '.'
    end if
end if

if (present(gamma2)) then
    if (is_nan(gamma2)) then
        ! In this case, we take the value hard coded in Powell's orginal code
        ! without any warning. It is useful when interfacing with MATLAB/Python.
        gamma2 = GAMMA2_DFT
    else if (gamma2 < 1 .or. is_inf(gamma2)) then
        gamma2 = GAMMA2_DFT
        print '(/1A, 1PD15.6, 1A)', solver//': invalid GAMMA2; it should a real number not less than 1; it is set to ', &
            & gamma2, '.'
    end if
end if

if (abs(rhobeg - rhoend) < 1.0E2_RP * EPS * max(abs(rhobeg), ONE)) then
! When the data is passed from the interfaces (e.g., MEX) to the Fortran
! code, RHOBEG, and RHOEND may change a bit. It was observed in a MATLAB
! test that MEX passed 1 to Fortran as 0.99999999999999978. Therefore,
! if we set RHOEND = RHOBEG in the interfaces, then it may happen that
! RHOEND > RHOBEG, which is considered as an invalid input. To avoid this
! situation, we force RHOBEG and RHOEND to equal when the difference is tiny.
    rhoend = rhobeg
end if

if (rhobeg <= 0 .or. is_nan(rhobeg) .or. is_inf(rhobeg)) then
    ! Take RHOEND into account if it has a valid value.
    if (is_finite(rhoend) .and. rhoend > 0) then
        rhobeg = max(TEN * rhoend, RHOBEG_DFT)
    else
        rhobeg = RHOBEG_DFT
    end if
    print '(/1A, 1PD15.6, 1A)', solver//': invalid RHOBEG; it should be a positive number; it is set to ', rhobeg, '.'
end if

if (rhoend <= 0 .or. rhobeg < rhoend .or. is_nan(rhoend) .or. is_inf(rhoend)) then
    rhoend = max(EPS, min(TENTH * rhobeg, RHOEND_DFT))
    print '(/1A, 1PD15.6, 1A)', solver//': invalid RHOEND; it should be a positive number and RHOEND <= RHOBEG; '// &
        & 'it is set to ', rhoend, '.'
end if

if (present(ctol)) then
    if (is_nan(ctol) .or. ctol < 0) then
        ctol = CTOL_DFT
        print '(/1A, 1PD15.6, 1A)', solver//': invalid CTOL; it should be a positive number; '// &
            & 'it is set to ', ctol, '.'
    end if
end if

! Postconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', solver)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', solver)
    if (present(npt)) then
        call assert(maxfun >= npt + 1, 'MAXFUN >= NPT + 1', solver)
        call assert(npt >= 3, 'NPT >= 3', solver)
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
