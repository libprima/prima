#include "fintrf.h"
#include "ppf.h"

module debug_mod
!--------------------------------------------------------------------------------------------------!
! This is a module defining some procedures concerning debugging, errors, and warnings.
! When interfacing Fortran code with MATLAB, we use this module to replace the one defined in
! fsrc/common/debug.F90, in order to generate errors or warnings that are native to MATLAB.
!
! Authors:
!   Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!   and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!   Department of Applied Mathematics,
!   The Hong Kong Polytechnic University
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015)
!
! Started in March 2020
!
! Last Modified: Thursday, June 23, 2022 PM04:11:18
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: assert, validate, wassert, backtr, warning, errstop

! Specify the interfaces of mexWarnMsgIdAndTxt and mexErrMsgIdAndTxt.
! We may use those specified in fmxapi_mod. However, that would make debug.F90 depend on fmxapi.F90,
! and make it impossible to use debug_mod in fmxapi.F90.
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
!! N.B.: As in C, we design ASSERT to operate only in the debug mode, i.e., when __DEBUGGING__ == 1;
!! when __DEBUGGING__ == 0, ASSERT does nothing. For the checking that should take effect in both
!! the debug and release modes, use VALIDATE (see below) instead. In the optimized mode of Python
!! (python -O), the Python `assert` will also be ignored. MATLAB does not behave in this way.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DEBUGGING
implicit none
logical, intent(in) :: condition  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the condition in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure
if (DEBUGGING .and. .not. condition) then
    call errstop(trim(srname), 'Assertion failed: '//trim(description))
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
    call errstop(trim(srname), 'Validation failed: '//trim(description))
end if
end subroutine validate


subroutine wassert(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether CONDITION is true.
! If no but DEBUGGING is true, print the following message to STDERR (but do not stop the program):
! SRNAME // 'Assertion failed: ' // DESCRIPTION
! MATLAB analogue:
!!if ~condition
!!    warning(sprintf('%s: Assertion failed: %s', srname, description))
!!end
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
    call warning(trim(srname), 'Assertion failed: '//trim(description))
end if
end subroutine wassert


subroutine errstop(srname, msg)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints 'ERROR: '//TRIM(SRNAME)//': '//TRIM(MSG)//'.', then stop.
! It also calls BACKTR to print the backtrace.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : MSGLEN
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: msg

character(len=MSGLEN) :: eid
character(len=MSGLEN) :: emsg

call backtr()
eid = 'FMXAPI:'//trim(srname)
emsg = trim(srname)//': '//trim(msg)//'.'
call mexErrMsgIdAndTxt(trim(eid), trim(emsg))
end subroutine errstop


subroutine backtr
!--------------------------------------------------------------------------------------------------!
! This subroutine calls a compiler-dependent intrinsic to show a backtrace if we are in the
! debugging mode, i.e., __DEBUGGING__ == 1.
! N.B.:
! 1. The intrinsic is compiler-dependent and does not exist in all compilers. Indeed, it is not
! standard-conforming. Therefore, compilers may warn that a non-standard intrinsic is in use.
! 2. More seriously, if the compiler is instructed to conform to the standards (e.g., gfortran with
! the option -std=f2003) while __DEBUGGING__ is set to 1, then the compilation may FAIL when linking,
! complaining that a subroutine cannot be found (e.g., backtrace for gfortran). In that case, we
! must set __DEBUGGING__ to 0 in ppf.h.
!
! Zaikun 20220412: MEX does not print the line numbers in the backtrace even with the '-g' option.
!--------------------------------------------------------------------------------------------------!
#if __DEBUGGING__ == 1

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
! This subroutine prints 'Warning: '//TRIM(SRNAME)//': '//TRIM(MSG)//'.'
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : MSGLEN
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: msg

character(len=MSGLEN) :: wid
character(len=MSGLEN) :: wmsg

wid = 'FMXAPI:'//trim(srname)
wmsg = trim(srname)//': '//trim(msg)//'.'
call mexWarnMsgIdAndTxt(trim(wid), trim(wmsg))
end subroutine warning


end module debug_mod
