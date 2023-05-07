#include "fintrf.h"


module output_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some subroutines concerning output to terminal/files. Note that these output
! operations are sequential in nature. In case parallelism is desirable (especially during
! initialization), the subroutines may have to be modified or disabled.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
!
! Started: July 2020
!
! Last Modified: Sunday, May 07, 2023 PM12:45:33
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
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, OUTUNIT, STDOUT, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: infos_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, MAXTR_REACHED, &
    & SMALL_TR_RADIUS, TRSUBP_FAILED, NAN_INF_X, NAN_INF_F, NAN_INF_MODEL, DAMAGING_ROUNDING, &
    & NO_SPACE_BETWEEN_BOUNDS, ZERO_LINEAR_CONSTRAINT
use, non_intrinsic :: string_mod, only : strip, num2str
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
character(len=:), allocatable :: fstat  ! 'OLD' or 'NEW'
character(len=:), allocatable :: fout
character(len=:), allocatable :: msg
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
integer :: funit ! File storage unit for the writing. Should be an integer of default kind.
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
    funit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    funit = OUTUNIT
    fout = strip(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=funit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
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
    write (funit, '(1X)')
end if
write (funit, '(/1A)') 'Return from '//solver//' because '//strip(msg)
write (funit, '(1A)') 'At the return from '//solver//'   Number of function evaluations = '//num2str(nf)

write (funit, '(1A)') 'Least function value = '//num2str(f)
if (is_constrained) then
    write (funit, '(1A)') 'Constraint violation = '//num2str(cstrv_loc)
end if
write (funit, '(1A)') 'The corresponding X is:'//new_line('')//num2str(x)
if (is_constrained .and. present(constr)) then
    write (funit, '(1A)') 'The constraint value is:'//new_line('')//num2str(constr)
end if
write (funit, '(1X)')

if (iprint < 0) then
    close (funit)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine retmsg


subroutine rhomsg(solver, iprint, nf, f, rho, x, cstrv, constr, cpen)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages when RHO is updated.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, OUTUNIT, STDOUT
use, non_intrinsic :: debug_mod, only : warning
use, non_intrinsic :: string_mod, only : strip
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
character(len=:), allocatable :: fstat  ! 'OLD' or 'NEW'
character(len=:), allocatable :: fout
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
integer :: funit ! Logical unit for the writing. Should be an integer of default kind.
logical :: fexist
logical :: is_constrained
real(RP) :: cstrv_loc

!====================!
! Calculation starts !
!====================!

if (abs(iprint) < 2) then
    return  ! No printing
elseif (iprint > 0) then
    funit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    funit = OUTUNIT
    fout = strip(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=funit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
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
    write (funit, '(1X)')
end if
if (present(cpen)) then
    write (funit, rpnf_fmt) 'New RHO = ', rho, '  CPEN = ', cpen, 'Number of function evaluations = ', nf
else
    write (funit, rnf_fmt) 'New RHO = ', rho, 'Number of function evaluations = ', nf
end if
write (funit, f_fmt) 'Least function value = ', f
if (is_constrained) then
    write (funit, cstrv_fmt) 'Constraint violation = ', cstrv_loc
end if
write (funit, x_fmt) 'The corresponding X is:', x
if (is_constrained .and. present(constr)) then
    write (funit, constr_fmt) 'The constraint value is:', constr
end if

if (iprint < 0) then
    close (funit)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine rhomsg


subroutine cpenmsg(solver, iprint, cpen)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints a message when CPEN is updated.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, OUTUNIT, STDOUT
use, non_intrinsic :: debug_mod, only : warning
use, non_intrinsic :: string_mod, only : strip
implicit none

! Compulsory inputs
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: iprint

! Optional inputs
real(RP), intent(in), optional :: cpen

! Local variables
character(len=*), parameter :: srname = 'CPENMSG'
character(len=:), allocatable :: fstat  ! 'OLD' or 'NEW'
character(len=:), allocatable :: fout
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
integer :: funit ! Logical unit for the writing. Should be an integer of default kind.
logical :: fexist

if (abs(iprint) < 2) then
    return  ! No printing
elseif (iprint > 0) then
    funit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    funit = OUTUNIT
    fout = strip(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=funit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
    if (iostat /= 0) then
        call warning(srname, 'Failed to open file '//fout)
        return
    end if
end if
! Print the message.
write (funit, p_fmt) 'Set CPEN to ', cpen
if (iprint < 0) then
    close (funit)
end if
end subroutine cpenmsg


subroutine fmsg(solver, iprint, nf, f, x, cstrv, constr)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages at each iteration.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, OUTUNIT, STDOUT
use, non_intrinsic :: debug_mod, only : warning
use, non_intrinsic :: string_mod, only : strip
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
character(len=:), allocatable :: fstat  ! 'OLD' or 'NEW'
character(len=:), allocatable :: fout
integer :: iostat  ! IO status of the writing. Should be an integer of default kind.
integer :: funit ! Logical unit for the writing. Should be an integer of default kind.
logical :: fexist
logical :: is_constrained
real(RP) :: cstrv_loc

!====================!
! Calculation starts !
!====================!

if (abs(iprint) < 3) then
    return  ! No printing
elseif (iprint > 0) then
    funit = STDOUT  ! Print the message to the standard out.
else  ! Print the message to a file named FOUT with the writing unit being OUTUNIT.
    funit = OUTUNIT
    fout = strip(solver)//'_output.txt'
    inquire (file=fout, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    open (unit=funit, file=fout, status=fstat, position='append', iostat=iostat, action='write')
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
write (funit, nf_fmt) 'Function number ', nf
write (funit, f_fmt) 'Function value = ', f
if (is_constrained) then
    write (funit, cstrv_fmt) 'Constraint violation = ', cstrv_loc
end if
write (funit, x_fmt) 'The corresponding X is:', x
if (is_constrained .and. present(constr)) then
    write (funit, constr_fmt) 'The constraint value is:', constr
end if

if (iprint < 0) then
    close (OUTUNIT)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine fmsg


!subroutine fprint(string, funit, fname, faction)
!use, non_intrinsic :: consts_mod, only : IK, OUTUNIT, STDIN, STDOUT, STDERR, DEBUGGING
!use, non_intrinsic :: debug_mod, only : assert, warning
!use, non_intrinsic :: string_mod, only : strip, num2str
!implicit none

!! Inputs
!character(len=*), intent(in) :: string
!integer, intent(in), optional :: funit
!character(len=*), intent(in), optional :: fname
!character(len=*), intent(in), optional :: faction

!! Local variables
!character(len=*), parameter :: srname = 'FPRINT'
!character(len=:), allocatable :: fname_loc
!character(len=:), allocatable :: fstat
!character(len=:), allocatable :: position
!integer :: funit_loc
!integer :: iostat
!logical :: fexist

!! Preconditions
!if (DEBUGGING) then
!    call assert(OUTUNIT > 0 .and. all(OUTUNIT /= [STDIN, STDOUT, STDERR]), &
!        & 'OUTUNIT is positive and not STDIN, STDOUT, or STDERR', srname)
!    if (present(funit)) then
!        call assert(funit /= STDIN, 'The file unit is not STDIN', srname)
!        if (present(fname)) then
!            call assert(funit /= STDOUT .and. funit /= STDERR, &
!                & 'If the file name is present, then the file unit is neither STDOUT nor STDERR', srname)
!        end if
!    end if
!    if (present(fname)) then
!        call assert(len(strip(fname)) > 0, 'The file name is nonempty and does not contain only spaces', srname)
!    end if
!    call assert(present(fname) .or. .not. present(faction), 'FACTION is present only if FNAME is specified', srname)
!    if (present(faction)) then
!        call assert(faction == 'write' .or. faction == 'w' .or. faction == 'append' .or. faction == 'a', &
!            & 'FACTION is either "write (w)" or "append (a)"', srname)
!    end if
!end if

!!====================!
!! Calculation starts !
!!====================!

!! Decide the file storage unit.
!if (present(funit)) then
!    funit_loc = funit
!else
!    if (present(fname)) then
!        funit_loc = OUTUNIT  ! Print the message to the writing unit OUTUNIT.
!    else
!        funit_loc = STDOUT  ! Print the message to the standard out.
!    end if
!end if

!! Decide the file name.
!if (present(fname)) then
!    fname_loc = fname
!elseif (funit_loc /= STDOUT .and. funit_loc /= STDERR) then
!    fname_loc = 'fort.'//num2str(int(funit_loc, IK))
!else
!    fname_loc = ''
!end if

!if (DEBUGGING) then
!    call assert(len(fname_loc) > 0 .eqv. (funit_loc /= STDOUT .and. funit_loc /= STDERR), &
!        & 'The file name is nonempty if and only if the file unit is neither STDOUT nor STDERR', srname)
!end if

!! Decide the position of printing.
!position = 'append'
!if (present(faction)) then
!    select case (faction)
!    case ('write', 'w')
!        position = 'rewind'
!    case ('append', 'a')
!        position = 'append'
!    case default
!        call warning(srname, 'Unknown file action "'//faction//'"')
!    end select
!end if

!! Open the file if necessary.
!iostat = 0
!if (len(fname_loc) > 0) then
!    inquire (file=fname_loc, exist=fexist)
!    fstat = merge(tsource='old', fsource='new', mask=fexist)
!    open (unit=funit_loc, file=fname_loc, status=fstat, position=position, iostat=iostat, action='write')
!    if (iostat /= 0) then
!        call warning(srname, 'Failed to open file '//fname_loc)
!        return
!    end if
!end if

!! Print the string.
!write (funit_loc, '(1A)') string

!! Close the file if necessary
!if (len(fname_loc) > 0 .and. iostat == 0) then
!    close (funit_loc)
!end if

!!====================!
!!  Calculation ends  !
!!====================!
!end subroutine fprint


end module output_mod
