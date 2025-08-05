module fprint_mod
!--------------------------------------------------------------------------------------------------!
! This module provides a subroutine that prints a string to STDOUT, STDERR, or a normal file.
!
! N.B.: When interfacing the code with MATLAB, this module needs to be revised to use the MATLAB
! MEX function mexPrintf instead of WRITE. This is because the Fortran WRITE cannot write to the
! STDOUT when the code is interfaced with MATLAB, since the STDOUT is hijacked by MEX. See
! https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
! https://www.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and papers.
!
! Started: July 2020
!
! Last Modified: Sunday, May 21, 2023 AM01:29:25
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: fprint


contains


subroutine fprint(string, funit, fname, faction)
use, non_intrinsic :: consts_mod, only : IK, STDIN, STDOUT, STDERR, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: string_mod, only : num2str
implicit none

! Inputs
character(len=*), intent(in) :: string
integer, intent(in), optional :: funit
character(len=*), intent(in), optional :: fname
character(len=*), intent(in), optional :: faction

! Local variables
character(len=*), parameter :: newline = new_line('')
character(len=*), parameter :: srname = 'FPRINT'
character(len=:), allocatable :: fname_loc
character(len=:), allocatable :: fstat
character(len=:), allocatable :: position
integer :: funit_loc
integer :: i
integer :: iostat
integer :: j
integer :: slen
logical :: fexist

! Preconditions
if (DEBUGGING) then
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
        funit_loc = -1  ! This value will not be used.
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

! Open the file if necessary.
iostat = 0
if (len(fname_loc) > 0) then
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
    if (present(funit)) then
        open (unit=funit_loc, file=fname_loc, status=fstat, position=position, iostat=iostat, action='write')
    else
        open (newunit=funit_loc, file=fname_loc, status=fstat, position=position, iostat=iostat, action='write')
    end if
    if (iostat /= 0) then
        call warning(srname, 'Failed to open file '//fname_loc)
        return
    end if
end if

! Print the string.
! N.B.: `WRITE (FUNIT_LOC, '(A)') STRING` would do what we want, but it causes "Buffer overflow on
! output" if string is long. This did occur with NAG Fortran Compiler R7.1(Hanzomon) Build 7122.
! To avoid this problem, we print the string line by line, separated by newlines.
i = 1
j = index(string, newline)  ! Index of the first newline in the string.
slen = len(string)
do while (j >= i)  ! J < I: No more newline in the string.
    write (funit_loc, '(A)') string(i:j - 1)  ! Print the string before the current newline.
    i = j + 1  ! Index of the character after the current newline.
    j = i + index(string(i:slen), newline) - 1  ! Index of the next newline.
end do
if (string(i:slen) /= '') then  ! Print the string after the last newline.
    write (funit_loc, '(A)') string(i:slen)
end if

! Close the file if necessary.
if (len(fname_loc) > 0 .and. iostat == 0) then
    close (funit_loc)
end if

!====================!
!  Calculation ends  !
!====================!
end subroutine fprint


end module fprint_mod
