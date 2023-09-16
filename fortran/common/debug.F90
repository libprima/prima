#include "ppf.h"

module debug_mod
!--------------------------------------------------------------------------------------------------!
! This is a module defining some procedures concerning debugging, errors, and warnings.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Saturday, September 16, 2023 AM09:14:37
!--------------------------------------------------------------------------------------------------!
implicit none
private
public :: assert, validate, wassert, backtr, warning, errstop


contains


subroutine assert(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether ASSERTION is true.
! If no but DEBUGGING is true, print the following message to STDERR and then stop the program:
! 'ERROR: ' // SRNAME // 'Assertion fails: ' // DESCRIPTION
! MATLAB analogue: assert(condition, sprintf('%s: Assertion fails: %s', srname, description))
! Python analogue: assert condition, srname + ': Assertion fails: ' + description
! C analogue: assert(condition)  /* An error message will be produced by the compiler */
! N.B.: As in C, we design ASSERT to operate only in the debug mode, i.e., when PRIMA_DEBUGGING == 1;
! when PRIMA_DEBUGGING == 0, ASSERT does nothing. For the checking that should take effect in both
! the debug and release modes, use VALIDATE (see below) instead. In the optimized mode of Python
! (python -O), the Python `assert` will also be ignored. MATLAB does not behave in this way.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: infos_mod, only : ASSERTION_FAILS
implicit none
logical, intent(in) :: condition  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the condition in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure
if (DEBUGGING .and. .not. condition) then
    call errstop(trim(adjustl(srname)), 'Assertion fails: '//trim(adjustl(description)), ASSERTION_FAILS)
end if
end subroutine assert


subroutine validate(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether CONDITION is true.
! If no, print the following message to STDERR and then stop the program:
! 'ERROR: ' // SRNAME // 'Validation fails: ' // DESCRIPTION
! MATLAB analogue: assert(condition, sprintf('%s: Validation fails: %s', srname, description))
! In Python or C, VALIDATE can be implemented following the Fortran implementation below.
! N.B.: ASSERT checks the condition only when debugging, but VALIDATE does it always.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: infos_mod, only : VALIDATION_FAILS
implicit none
logical, intent(in) :: condition  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the condition in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure
if (.not. condition) then
    call errstop(trim(adjustl(srname)), 'Validation fails: '//trim(adjustl(description)), VALIDATION_FAILS)
end if
end subroutine validate


subroutine wassert(condition, description, srname)
!--------------------------------------------------------------------------------------------------!
! This subroutine checks whether CONDITION is true.
! If no but DEBUGGING is true, print the following message to STDERR (but do not stop the program):
! 'Warning: ' // SRNAME // 'Assertion fails: ' // DESCRIPTION
! MATLAB analogue:
! !if ~condition
! !    warning(sprintf('%s: Assertion fails: %s', srname, description))
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
    call warning(trim(adjustl(srname)), 'Assertion fails: '//trim(adjustl(description)))
end if
end subroutine wassert


subroutine errstop(srname, msg, code)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints 'ERROR: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.' to STDERR, then stop.
! It also calls BACKTR to print the backtrace.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : STDERR
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: msg
integer, intent(in), optional :: code

! `backtr` prints a backtrace. As of gfortran 12.1.0, even without calling `backtrace`, a backtrace
! is printed when the program is stopped by an error stop. Therefore, here, we do not call `backtr`
! if the compiler is gfortran.
#if !defined __GFORTRAN__
call backtr()
#endif

write (STDERR, '(/A/)') 'ERROR: '//trim(adjustl(srname))//': '//trim(adjustl(msg))//'.'
if (present(code)) then
    ! N.B.: In Fortran 2008, stop code must be a scalar default character or integer CONSTANT
    ! expression, but Fortran 2018 lifts the requirement on constancy. gfortran is strict in this
    ! aspect. Consequently, for gfortran, compile with either `-std=f2018` or no `-std` at all.
    error stop code
else
    error stop
end if
! N.B.
! 1. ERROR STOP means to stop the whole program.
! 2. (Zaikun 20230410): We prefer ERROR STOP to STOP, as the former has been allowed in PURE
! procedures since F2018. Later, when F2018 is better supported, we should take advantage of this
! feature to make our subroutines PURE whenever possible.
end subroutine errstop


subroutine backtr()
!--------------------------------------------------------------------------------------------------!
! This subroutine calls a compiler-dependent intrinsic to show a backtrace if we are in the
! debugging mode, i.e., PRIMA_DEBUGGING == 1.
! N.B.:
! 1. The intrinsic is compiler-dependent and does not exist in all compilers. Indeed, it is not
! standard-conforming. Therefore, compilers may warn that a non-standard intrinsic is in use.
! 2. More seriously, if the compiler is instructed to conform to the standards (e.g., gfortran with
! the option -std=f2018) while PRIMA_DEBUGGING is set to 1, then the compilation may FAIL when
! linking, complaining that a subroutine cannot be found (e.g., `backtrace` for gfortran). In that
! case, we must either use the `-fall-intrinsics` option of `gfortran`, or set PRIMA_DEBUGGING to 0
! in ppf.h. This is also why in this subroutine we do not use the constant DEBUGGING defined in the
! consts_mod module but use the macro PRIMA_DEBUGGING in ppf.h.
! 3. As of gfortran 12.1.0, even without calling `backtrace`, a backtrace is printed when the
! program is stopped by an error stop. Therefore, in `errstop`, we do not call `backtr` if the
! compiler is gfortran. However, we cannot remove `backtrace` in `backtr`, because `backtr` is
! also invoked in `wassert`, where `backtrace` is still needed as error stop is not involved.
!--------------------------------------------------------------------------------------------------!
#if PRIMA_DEBUGGING == 1

#if defined __GFORTRAN__
implicit none
call backtrace  ! gfortran: if `-std=f20xy` is imposed, then `-fall-intrinsics` is needed.
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
! This subroutine prints 'Warning: '//STRIP(SRNAME)//': '//STRIP(MSG)//'.' to STDERR.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : STDERR
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: msg

write (STDERR, '(/A/)') 'Warning: '//trim(adjustl(srname))//': '//trim(adjustl(msg))//'.'
end subroutine warning


end module debug_mod
