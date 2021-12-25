module output_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some subroutines concerning output to terminal/files. Note that these output
! operations are sequential in nature. In case parallelism is desirable (especially during
! initializaton), the subroutines may have to be modified or disabled.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and papers.
!
! Started: July 2020
!
! Last Modified: Saturday, December 25, 2021 PM05:14:04
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: retmssg, rhomssg, fmssg

!--------------------------------------------------------------------------------------------------!
! Formats.
! Format for F at return: 16 digits for base, 4 digits for exponent.
character(len=*), parameter :: ffmt = '1PE25.16E4'
! Format for intermediate F: 10 digits for base, 4 digits for exponent.
character(len=*), parameter :: ffmt_intermediate = '1PE19.10E4'
! Format for X: 8 digits for base, 4 digits for exponent, 4 components in each line.
character(len=*), parameter :: xfmt = '/(1P, 4E19.8E4)'
! Format for RHO: 4 digits for base, 2 digits for exponent.
character(len=*), parameter :: rfmt = '1PE11.4E2'
! Format for integers: 10 digits.
character(len=*), parameter :: ifmt = 'I0'  ! Use the minimum number of digits needed to print integers
! Separating spaces: 3 spaces.
character(len=*), parameter :: spaces = '3X'
! Format for NF at return.
character(len=*), parameter :: nf_fmt = '(1A, '//spaces//', 1A, '//ifmt//')'
! Format for F and X at return and when RHO is updated.
character(len=*), parameter :: fx_fmt = '(1A, '//ffmt//', '//spaces//', 1A, '//xfmt//')'
character(len=*), parameter :: fcx_fmt = '(1A, '//ffmt//', '//spaces//', 1A, '//ffmt//', '//spaces//', 1A, '//xfmt//')'
! Format for constraint violation and constraint value at return and when RHO is updated.
character(len=*), parameter :: cc_fmt = fx_fmt
character(len=*), parameter :: cv_fmt = '(1A, '//ffmt//')'
! Format for RHO and NF when RHO is updated.
character(len=*), parameter :: rnf_fmt = '(/1A, '//rfmt//', '//spaces//', 1A, '//ifmt//')'
character(len=*), parameter :: rpnf_fmt = '(/1A, '//rfmt//', '//spaces//', 1A, '//rfmt//', '//spaces//', 1A, '//ifmt//')'
! Format for NF, F, and X during iterations.
character(len=*), parameter :: nffx_fmt = &
    & '(/1A, '//ifmt//', '//spaces//', 1A, '//ffmt_intermediate//', '//spaces//', 1A, '//xfmt//')'
!--------------------------------------------------------------------------------------------------!


contains


subroutine retmssg(solver, info, iprint, nf, f, x, cstrv, constr)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages at return.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, MSSGLEN, FNAMELEN, OUTUNIT, STDOUT, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED
use, non_intrinsic :: info_mod, only : SMALL_TR_RADIUS, TRSUBP_FAILED
use, non_intrinsic :: info_mod, only : NAN_INF_X, NAN_INF_F, NAN_MODEL
implicit none

! Compulsory inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: info
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)

! Optional inputs
real(RP), intent(in), optional :: cstrv
real(RP), intent(in), optional :: constr(:)

! Local variables
character(len=*), parameter :: srname = 'RETMSSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
character(len=MSSGLEN) :: mssg
integer :: ios  ! IO status of the writing. Should be an integer of default kind.
integer :: wunit ! Logical unit for the writing. Should be an integer of default kind.
integer(IK), parameter :: valid_exit_flags(7) = [FTARGET_ACHIEVED, MAXFUN_REACHED, SMALL_TR_RADIUS, &
    & TRSUBP_FAILED, NAN_INF_F, NAN_INF_X, NAN_MODEL]
logical :: fexist
logical :: is_constrained
real(RP) :: cstrv_loc

! Preconditions
if (DEBUGGING) then
    call assert(any(info == valid_exit_flags), 'The exit flag is valid', srname)
end if

!====================!
! Calculation starts !
!====================!

if (abs(iprint) < 1) then
    return  ! No printing
elseif (iprint > 0) then
    wunit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    wunit = OUTUNIT
    fout = trim(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=ios, action='write')
    if (ios /= 0) then
        call warning(srname, 'Failed to open file '//fout)
        return
    end if
end if

! Decide whether the problem is truly constrained.
if (present(constr)) then
    is_constrained = (size(constr) > 0)
else
    is_constrained = present(cstrv)
end if

! Decide the constraint violation.
if (present(cstrv)) then
    cstrv_loc = cstrv
elseif (present(constr)) then
    cstrv_loc = maxval([ZERO, -constr])
else
    cstrv_loc = ZERO
end if

! Decide the exit message.
select case (info)
case (FTARGET_ACHIEVED)
    mssg = 'the target function value is achieved.'
case (MAXFUN_REACHED)
    mssg = 'the objective function has been evaluated MAXFUN times.'
case (SMALL_TR_RADIUS)
    mssg = 'the trust region radius reaches its lower bound.'
case (TRSUBP_FAILED)
    mssg = 'a trust region step has failed to reduce the quadratic model.'
case (NAN_INF_X)
    mssg = 'NaN or Inf occurs in x.'
case (NAN_INF_F)
    mssg = 'the objective function returns NaN/+Inf.'
case (NAN_MODEL)
    mssg = 'NaN occurs in the models.'
case default
    mssg = 'UNKNOWN EXIT FLAG'
end select

! Print the message.
if (abs(iprint) >= 3) then
    write (wunit, '(1X)')
end if
write (wunit, '(/1A)') 'Return from '//solver//' because '//trim(mssg)
write (wunit, nf_fmt) 'At the return from '//solver, 'Number of function evaluations = ', nf
write (wunit, fx_fmt) 'Least function value = ', f, 'The corresponding X is:', x
if (is_constrained) then
    if (present(constr)) then
        write (wunit, cc_fmt) 'Constraint violation = ', cstrv_loc, 'The constraint value is:', constr
    else
        write (wunit, cv_fmt) 'Constraint violation = ', cstrv_loc
    end if
end if
write (wunit, '(1X)')

if (iprint < 0) then
    close (wunit)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine retmssg


subroutine rhomssg(solver, iprint, nf, f, rho, x, cstrv, constr, cpen)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages when RHO is updated.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, OUTUNIT, STDOUT, FNAMELEN
use, non_intrinsic :: debug_mod, only : warning
implicit none

! Compulsory inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: rho
real(RP), intent(in) :: x(:)

! Optional inputs
real(RP), intent(in), optional :: cstrv
real(RP), intent(in), optional :: constr(:)
real(RP), intent(in), optional :: cpen

! Local variables
character(len=*), parameter :: srname = 'RHOMSSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
integer :: ios  ! IO status of the writing. Should be an integer of default kind.
integer :: wunit ! Logical unit for the writing. Should be an integer of default kind.
logical :: fexist
logical :: is_constrained
real(RP) :: cstrv_loc

!====================!
! Calculation starts !
!====================!

if (abs(iprint) < 2) then
    return  ! No printing
elseif (iprint > 0) then
    wunit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    wunit = OUTUNIT
    fout = trim(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=ios, action='write')
    if (ios /= 0) then
        call warning(srname, 'Failed to open file '//fout)
        return
    end if
end if

! Decide whether the problem is truly constrained.
if (present(constr)) then
    is_constrained = (size(constr) > 0)
else
    is_constrained = present(cstrv)
end if

! Decide the constraint violation.
if (present(cstrv)) then
    cstrv_loc = cstrv
elseif (present(constr)) then
    cstrv_loc = maxval([ZERO, -constr])
else
    cstrv_loc = ZERO
end if

! Print the message.
if (abs(iprint) >= 3) then
    write (wunit, '(1X)')
end if
if (present(cpen)) then
    write (wunit, rpnf_fmt) 'New RHO = ', rho, '  CPEN = ', cpen, 'Number of function evaluations = ', nf
else
    write (wunit, rnf_fmt) 'New RHO = ', rho, 'Number of function evaluations = ', nf
end if
if (is_constrained) then
    write (wunit, fcx_fmt) 'Least function value = ', f, 'Constraint violation = ', cstrv_loc, &
        & 'The corresponding X is:', x
else
    write (wunit, fx_fmt) 'Least function value = ', f, 'The corresponding X is:', x
end if

if (iprint < 0) then
    close (wunit)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine rhomssg


subroutine fmssg(solver, iprint, nf, f, x, cstrv, constr)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages at each iteration.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, OUTUNIT, STDOUT, FNAMELEN
use, non_intrinsic :: debug_mod, only : warning
implicit none

! Inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)

! Optional inputs
real(RP), intent(in), optional :: cstrv
real(RP), intent(in), optional :: constr(:)

! Local variables
character(len=*), parameter :: srname = 'FMSSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
integer :: ios  ! IO status of the writing. Should be an integer of default kind.
integer :: wunit ! Logical unit for the writing. Should be an integer of default kind.
logical :: fexist
logical :: is_constrained
real(RP) :: cstrv_loc

!====================!
! Calculation starts !
!====================!

if (abs(iprint) < 3) then
    return  ! No printing
elseif (iprint > 0) then
    wunit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    wunit = OUTUNIT
    fout = trim(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=ios, action='write')
    if (ios /= 0) then
        call warning(srname, 'Failed to open file '//fout)
        return
    end if
end if

! Decide whether the problem is truly constrained.
if (present(constr)) then
    is_constrained = (size(constr) > 0)
else
    is_constrained = present(cstrv)
end if

! Decide the constraint violation.
if (present(cstrv)) then
    cstrv_loc = cstrv
elseif (present(constr)) then
    cstrv_loc = maxval([ZERO, -constr])
else
    cstrv_loc = ZERO
end if

! Print the message.
write (wunit, nffx_fmt) 'Function number ', nf, 'F = ', f, 'The corresponding X is:', x
if (is_constrained) then
    if (present(constr)) then
        write (wunit, cc_fmt) 'Constraint violation = ', cstrv_loc, 'The constraint value is:', constr
    else
        write (wunit, cv_fmt) 'Constraint violation = ', cstrv_loc
    end if
end if

if (iprint < 0) then
    close (OUTUNIT)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine fmssg


end module output_mod
