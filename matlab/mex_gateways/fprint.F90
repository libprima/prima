#include "fintrf.h"

module fprint_mod
!--------------------------------------------------------------------------------------------------!
! This module provides a subroutine that prints a string to STDOUT, STDERR, or a normal file.
!
! N.B.: When interfacing the code with MATLAB, we need to use the MEX function mexPrintf instead of
! WRITE. This is because the Fortran WRITE cannot write to the
! STDOUT when the code is interfaced with MATLAB, since the STDOUT is hijacked by MEX. See
! https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
! https://www.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
!
! Coded by Zaikun ZHANG (www.zhangzk.net)
!
! Started: July 2020
!
! Last Modified: Monday, May 08, 2023 PM06:15:37
!--------------------------------------------------------------------------------------------------!

! N.B.: INT32_MEX is indeed INT32, i.e., the kind of INTEGER*4. We name it INT32_MEX instead of
! INT32 so that it is easily locatable.
use, non_intrinsic :: consts_mod, only : INT32_MEX => INT32
implicit none
private
public :: fprint


! Specify the interfaces of mexPrintf and mexString.
! We may use those specified in fmxapi_mod. However, that would make fprint.F90 depend on fmxapi.F90,
! and make it impossible to use fprint_mod in fmxapi.F90, which may become necessary in the future.
interface
    function mexPrintf(message)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mexPrintf
    character*(*), intent(in) :: message
    end function mexPrintf

    function mexEvalString(command)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mexEvalString
    character*(*), intent(in) :: command
    end function mexEvalString
end interface


contains


subroutine fprint(string, funit, fname, faction)
use, non_intrinsic :: consts_mod, only : IK, OUTUNIT, STDIN, STDOUT, STDERR, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: string_mod, only : num2str
implicit none

! Inputs
character(len=*), intent(in) :: string
integer, intent(in), optional :: funit
character(len=*), intent(in), optional :: fname
character(len=*), intent(in), optional :: faction

! Local variables
character(len=*), parameter :: srname = 'FPRINT'
character(len=:), allocatable :: fname_loc
character(len=:), allocatable :: fstat
character(len=:), allocatable :: position
integer :: funit_loc
integer :: iostat
integer(INT32_MEX) :: k
logical :: fexist

! Preconditions
if (DEBUGGING) then
    call assert(OUTUNIT > 0 .and. all(OUTUNIT /= [STDIN, STDOUT, STDERR]), &
        & 'OUTUNIT is positive and not STDIN, STDOUT, or STDERR', srname)
    if (present(funit)) then
        call assert(funit /= STDIN, 'The file unit is not STDIN', srname)
        if (present(fname)) then
            call assert(len(fname) == 0 .eqv. (funit == STDOUT .or. funit == STDERR), &
                & 'The file name is empty if and only if the file unit is either STDOUT or STDERR', srname)
        end if
    end if
    if (present(faction)) then
        call assert(faction == 'write' .or. faction == 'w' .or. faction == 'append' .or. faction == 'a', &
            & 'FACTION is either "write (w)" or "append (a)"', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! Decide the file storage unit.
if (present(funit)) then
    funit_loc = funit
else
    if (present(fname)) then
        funit_loc = OUTUNIT  ! Print the message to the writing unit OUTUNIT.
    else
        funit_loc = STDOUT  ! Print the message to the standard out.
    end if
end if

! Decide the file name.
if (present(fname)) then
    fname_loc = fname
elseif (funit_loc /= STDOUT .and. funit_loc /= STDERR) then
    fname_loc = 'fort.'//num2str(int(funit_loc, IK))
else
    fname_loc = ''
end if

if (DEBUGGING) then
    call assert(len(fname_loc) == 0 .eqv. (funit_loc == STDOUT .or. funit_loc == STDERR), &
        & 'The file name is empty if and only if the file unit is either STDOUT or STDERR', srname)
end if

! Print the string.
if (len(fname_loc) == 0) then
    ! N.B.: We append a trailing new line to the string to be printed. This is because mexPrintf
    ! does not print strings to the standard output on new lines as of MATLAB R2023a. Ideally, we
    ! should add a new line to the beginning of the string. However, this may lead to the phenomenon
    ! that messages printed by other functions are appended to the end of strings by FPRINT. Note
    ! that all the strings received by FPRINT from MESSAGE have a leading new line.
    k = mexPrintf(string//new_line(''))
    if (k /= len(string) + 1) then
        call warning(srname, 'mexPrintf failed to print a string to the standard output')
        return
    end if
    k = mexEvalString('drawnow;')  ! Without this, the string printed above may not show immediately.
    if (k /= 0) then
        call warning(srname, 'Failed to call mexEvalString(''drawnow;'')')
    end if
else
    ! Decide the position for OPEN. This is the only place where FACTION is used.
    position = 'append'
    if (present(faction)) then
        select case (faction)
        case ('write', 'w')
            position = 'rewind'
        case ('append', 'a')
            position = 'append'
        case default
            call warning(srname, 'Unknown file action "'//faction//'"')
        end select
    end if
    ! Check whether the file is already existing.
    inquire (file=fname_loc, exist=fexist)
    fstat = merge(tsource='old', fsource='new', mask=fexist)
    ! Open the file.
    open (unit=funit_loc, file=fname_loc, status=fstat, position=position, iostat=iostat, action='write')
    if (iostat /= 0) then
        call warning(srname, 'Failed to open file '//fname_loc)
        return
    end if
    ! Print the string.
    write (funit_loc, '(1A)') string
    ! Close the file.
    close (funit_loc)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine fprint


end module fprint_mod
