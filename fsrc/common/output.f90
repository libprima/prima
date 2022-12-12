module output_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some subroutines concerning output to terminal/files. Note that these output
! operations are sequential in nature. In case parallelism is desirable (especially during
! initializaton), the subroutines may have to be modified or disabled.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
!
! Started: July 2020
!
! Last Modified: Monday, December 12, 2022 PM03:18:06
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: retmsg, rhomsg, fmsg, cpenmsg

!--------------------------------------------------------------------------------------------------!
! Formats.
! Separating spaces: 3 spaces.
character(len=*), parameter :: spaces = '3X'
! Format for F: 16 digits for base, 4 digits for exponent.
character(len=*), parameter :: ffmt = '1PE25.16E4'
character(len=*), parameter :: f_fmt = '(1A, '//ffmt//')'
! Format for CSTRV: the same as F.
character(len=*), parameter :: cstrv_fmt = f_fmt
! Format for X: 8 digits for base, 4 digits for exponent, 4 components in each line.
character(len=*), parameter :: xfmt = '/(1P, 4E19.8E4)'
character(len=*), parameter :: x_fmt = '(1A, '//xfmt//')'
! Format for CONSTR: the same as X.
character(len=*), parameter :: constr_fmt = x_fmt
! Format for integers: Use the minimum number of digits needed to print integers.
character(len=*), parameter :: ifmt = 'I0'
! Format for NF during the iterations.
character(len=*), parameter :: nf_fmt = '(/1A, '//ifmt//')'
! Format for NF at return.
character(len=*), parameter :: retnf_fmt = '(1A, '//spaces//', 1A, '//ifmt//')'
! Format for RHO: 8 digits for base, 4 digits for exponent.
character(len=*), parameter :: rfmt = '1PE17.8E4'
! Format for RHO and NF when RHO is updated, with or without CPEN.
character(len=*), parameter :: rnf_fmt = '(/1A, '//rfmt//', '//spaces//', /1A, '//ifmt//')'
character(len=*), parameter :: rpnf_fmt = '(/1A, '//rfmt//', '//spaces//', 1A, '//rfmt//', '//spaces//', /1A, '//ifmt//')'
! Format for the constraint penalty parameter: the same as RHO.
character(len=*), parameter :: pfmt = rfmt
! Format for the constraint penalty parameter when it is updated.
character(len=*), parameter :: p_fmt = '(/1A, '//pfmt//')'
!--------------------------------------------------------------------------------------------------!


contains


subroutine retmsg(solver, info, iprint, nf, f, x, cstrv, constr)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages at return.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, MSGLEN, FNAMELEN, OUTUNIT, STDOUT, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: infos_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, MAXTR_REACHED, &
    & SMALL_TR_RADIUS, TRSUBP_FAILED, NAN_INF_X, NAN_INF_F, NAN_INF_MODEL, DAMAGING_ROUNDING, &
    & NO_SPACE_BETWEEN_BOUNDS, ZERO_LINEAR_CONSTRAINT
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
character(len=*), parameter :: srname = 'RETMSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
character(len=MSGLEN) :: msg
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
integer :: wunit ! Logical unit for the writing. Should be an integer of default kind.
integer(IK), parameter :: valid_exit_flags(11) = [FTARGET_ACHIEVED, MAXFUN_REACHED, MAXTR_REACHED, &
    & SMALL_TR_RADIUS, TRSUBP_FAILED, NAN_INF_F, NAN_INF_X, NAN_INF_MODEL, DAMAGING_ROUNDING, &
    & NO_SPACE_BETWEEN_BOUNDS, ZERO_LINEAR_CONSTRAINT]
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
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
    if (iostat /= 0) then
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
    cstrv_loc = maxval([ZERO, -constr])  ! N.B.: We assume that the constraint is CONSTR >= 0.
else
    cstrv_loc = ZERO
end if

! Decide the exit message.
select case (info)
case (FTARGET_ACHIEVED)
    msg = 'the target function value is achieved.'
case (MAXFUN_REACHED)
    msg = 'the objective function has been evaluated MAXFUN times.'
case (MAXTR_REACHED)
    msg = 'the maximal number of trust region iterations has been reached.'
case (SMALL_TR_RADIUS)
    msg = 'the trust region radius reaches its lower bound.'
case (TRSUBP_FAILED)
    msg = 'a trust region step has failed to reduce the quadratic model.'
case (NAN_INF_X)
    msg = 'NaN or Inf occurs in x.'
case (NAN_INF_F)
    msg = 'the objective function returns NaN/+Inf.'
case (NAN_INF_MODEL)
    msg = 'NaN or Inf occurs in the models.'
case (DAMAGING_ROUNDING)
    msg = 'rounding errors are becoming damaging.'
case (NO_SPACE_BETWEEN_BOUNDS)
    msg = 'there is no space between the lower and upper bounds of variable.'
case (ZERO_LINEAR_CONSTRAINT)
    msg = 'one of the linear constraints has a zero gradient'
case default
    msg = 'UNKNOWN EXIT FLAG'
end select

! Print the message.
if (abs(iprint) >= 3) then
    write (wunit, '(1X)')
end if
write (wunit, '(/1A)') 'Return from '//solver//' because '//trim(msg)
write (wunit, retnf_fmt) 'At the return from '//solver, 'Number of function evaluations = ', nf
write (wunit, f_fmt) 'Least function value = ', f
if (is_constrained) then
    write (wunit, cstrv_fmt) 'Constraint violation = ', cstrv_loc
end if
write (wunit, x_fmt) 'The corresponding X is:', x
if (is_constrained .and. present(constr)) then
    write (wunit, constr_fmt) 'The constraint value is:', constr
end if
write (wunit, '(1X)')

if (iprint < 0) then
    close (wunit)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine retmsg


subroutine rhomsg(solver, iprint, nf, f, rho, x, cstrv, constr, cpen)
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
character(len=*), parameter :: srname = 'RHOMSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
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
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
    if (iostat /= 0) then
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
    cstrv_loc = maxval([ZERO, -constr])  ! N.B.: We assume that the constraint is CONSTR >= 0.
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
write (wunit, f_fmt) 'Least function value = ', f
if (is_constrained) then
    write (wunit, cstrv_fmt) 'Constraint violation = ', cstrv_loc
end if
write (wunit, x_fmt) 'The corresponding X is:', x
if (is_constrained .and. present(constr)) then
    write (wunit, constr_fmt) 'The constraint value is:', constr
end if

if (iprint < 0) then
    close (wunit)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine rhomsg


subroutine cpenmsg(solver, iprint, cpen)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints a message when CPEN is updated.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, FNAMELEN, OUTUNIT, STDOUT
use, non_intrinsic :: debug_mod, only : warning
implicit none

! Compulsory inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: iprint

! Optional inputs
real(RP), intent(in), optional :: cpen

! Local variables
character(len=*), parameter :: srname = 'CPENMSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
integer :: wunit ! Logical unit for the writing. Should be an integer of default kind.
logical :: fexist

if (abs(iprint) < 2) then
    return  ! No printing
elseif (iprint > 0) then
    wunit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    wunit = OUTUNIT
    fout = trim(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
    if (iostat /= 0) then
        call warning(srname, 'Failed to open file '//fout)
        return
    end if
end if
! Print the message.
write (wunit, p_fmt) 'Set CPEN to ', cpen
if (iprint < 0) then
    close (wunit)
end if
end subroutine cpenmsg


subroutine fmsg(solver, iprint, nf, f, x, cstrv, constr)
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
character(len=*), parameter :: srname = 'FMSG'
character(len=3) :: fstat  ! 'OLD' or 'NEW'
character(len=FNAMELEN) :: fout
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
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
    open (unit=wunit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
    if (iostat /= 0) then
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
    cstrv_loc = maxval([ZERO, -constr])  ! N.B.: We assume that the constraint is CONSTR >= 0.
else
    cstrv_loc = ZERO
end if

! Print the message.
write (wunit, nf_fmt) 'Function number ', nf
write (wunit, f_fmt) 'Function value = ', f
if (is_constrained) then
    write (wunit, cstrv_fmt) 'Constraint violation = ', cstrv_loc
end if
write (wunit, x_fmt) 'The corresponding X is:', x
if (is_constrained .and. present(constr)) then
    write (wunit, constr_fmt) 'The constraint value is:', constr
end if

if (iprint < 0) then
    close (OUTUNIT)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine fmsg


end module output_mod
