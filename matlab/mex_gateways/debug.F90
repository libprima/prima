#include "fintrf.h"
#include "ppf.h"

module debug_mod
!--------------------------------------------------------------------------------------------------!
! This is a module defining some procedures concerning debugging, errors, and warnings.
! When interfacing Fortran code with MATLAB, we use this module to replace the one defined in
! fortran/common/debug.F90, in order to generate errors or warnings that are native to MATLAB.
!
! Coded by Zaikun ZHANG (www.zhangzk.net)
!
! Started in July 2020
!
! Last Modified: Sunday, May 21, 2023 PM05:46:15
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: assert, validate, wassert, backtr, warning, errstop


! Specify the interfaces of mexWarnMsgIdAndTxt and mexErrMsgIdAndTxt.
! We cannot use those specified in fmxapi_mod. This is because fmxapi.F90 depends on debug.F90.
interface
    subroutine mexErrMsgIdAndTxt(eid, emsg)
    implicit none
    character*(*), intent(in) :: eid, emsg
    end subroutine mexErrMsgIdAndTxt
    subroutine mexWarnMsgIdAndTxt(wid, wmsg)
    implicit none
    character*(*), intent(in) :: wid, wmsg
    end subroutine mexWarnMsgIdAndTxt
end interface


contains


subroutine assert(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether ASSERTION is true.
! If no but DEBUGGING is true, print the following message and then stop the program:
! SRNAME // 'Assertion failed: ' // DESCRIPTION
! MATLAB analogue: assert(condition, sprintf('%s: Assertion failed: %s', srname, description))
! Python analogue: assert condition, srname + ': Assertion failed: ' + description
! C analogue: assert(condition)  /* An error message will be produced by the compiler */
!--------------------------------------------------------------------------------------------------!
! N.B.: As in C, we design ASSERT to operate only in the debug mode, i.e., when PRIMA_DEBUGGING == 1;
! when PRIMA_DEBUGGING == 0, ASSERT does nothing. For the checking that should take effect in both
! the debug and release modes, use VALIDATE (see below) instead. In the optimized mode of Python
! (python -O), the Python `assert` will also be ignored. MATLAB does not behave in this way.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DEBUGGING
implicit none
logical, intent(in) :: condition  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the condition in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure
if (DEBUGGING .and. .not. condition) then
    call errstop(trim(adjustl(srname)), 'Assertion failed: '//trim(adjustl(description)))
end if
end subroutine assert


subroutine validate(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether CONDITION is true.
! If no, print the following message to and then stop the program:
! SRNAME // 'Validation failed: ' // DESCRIPTION
! MATLAB analogue: assert(condition, sprintf('%s: Validation failed: %s', srname, description))
! In Python or C, VALIDATE can be implemented following the Fortran implementation below.
! N.B.: ASSERT checks the condition only when debugging, but VALIDATE does it always.
!--------------------------------------------------------------------------------------------------!
implicit none
logical, intent(in) :: condition  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the condition in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure
if (.not. condition) then
    call errstop(trim(adjustl(srname)), 'Validation failed: '//trim(adjustl(description)))
end if
end subroutine validate


subroutine wassert(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether CONDITION is true.
! If no but DEBUGGING is true, print the following message to STDERR (but do not stop the program):
! SRNAME // 'Assertion failed: ' // DESCRIPTION
! MATLAB analogue:
! !if ~condition
! !    warning(sprintf('%s: Assertion failed: %s', srname, description))
! !end
! In Python or C, WASSERT can be implemented following the Fortran implementation below.
! N.B.: When DEBUGGING is true, ASSERT stops the program with an error if the condition is false,
! but WASSERT only raises a warning.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DEBUGGING
implicit none
logical, intent(in) :: condition  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the condition in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure
if (DEBUGGING .and. .not. condition) then
    call backtr()
    call warning(trim(adjustl(srname)), 'Assertion failed: '//trim(adjustl(description)))
end if
end subroutine wassert


subroutine errstop(srname, msg, code)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints 'ERROR: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.', then stop.
! It also calls BACKTR to print the backtrace.
!--------------------------------------------------------------------------------------------------!
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: msg
integer, intent(in), optional :: code

character(len=128) :: code_str
character(len=:), allocatable :: eid
character(len=:), allocatable :: emsg

call backtr()
eid = 'FMXAPI:'//trim(adjustl(srname))
if (present(code)) then
    write (code_str, '(I0)') code
    emsg = trim(adjustl(srname))//': '//trim(adjustl(msg))//'. The error code is '//trim(adjustl(code_str))//'.'
else
    emsg = trim(adjustl(srname))//': '//trim(adjustl(msg))//'.'
end if
call mexErrMsgIdAndTxt(eid, emsg)
end subroutine errstop


subroutine backtr
!--------------------------------------------------------------------------------------------------!
! This subroutine calls a compiler-dependent intrinsic to show a backtrace if we are in the
! debugging mode, i.e., PRIMA_DEBUGGING == 1.
! N.B.:
! 1. The intrinsic is compiler-dependent and does not exist in all compilers. Indeed, it is not
! standard-conforming. Therefore, compilers may warn that a non-standard intrinsic is in use.
! 2. More seriously, if the compiler is instructed to conform to the standards (e.g., gfortran with
! the option -std=f2018) while PRIMA_DEBUGGING is set to 1, then the compilation may FAIL when
! linking, complaining that a subroutine cannot be found (e.g., backtrace for gfortran). In that
! case, we must set PRIMA_DEBUGGING to 0 in ppf.h. This is also why in this subroutine we do not use
! the constant DEBUGGING defined in the consts_mod module but use the macro PRIMA_DEBUGGING in ppf.h.
! 3. Zaikun 2020412/20230521: Even with the `-g` option, MEX does not always print the line numbers
! in the backtrace when calling `backtrace`. No idea why.
!--------------------------------------------------------------------------------------------------!
#if PRIMA_DEBUGGING == 1

#if defined __GFORTRAN__
implicit none
call backtrace
#elif defined __INTEL_COMPILER
use, non_intrinsic :: ifcore, only : tracebackqq
implicit none
call tracebackqq(user_exit_code=-1)
! According to "Intel Fortran Compiler 19.1 Developer Guide and Reference", item "TRACEBACKQQ":
! By specifying a user exit code of -1, control returns to the calling program. Specifying a user
! exit code with a positive value requests that specified value be returned to the operating system.
! The default value is 0, which causes the application to abort execution.
#endif

#endif
end subroutine backtr


subroutine warning(srname, msg)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints 'Warning: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.'
!--------------------------------------------------------------------------------------------------!
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: msg

character(len=:), allocatable :: wid
character(len=:), allocatable :: wmsg

wid = 'FMXAPI:'//trim(adjustl(srname))
wmsg = trim(adjustl(srname))//': '//trim(adjustl(msg))//'.'
call mexWarnMsgIdAndTxt(wid, wmsg)
end subroutine warning


end module debug_mod
