! DEBUG_MOD is a module defining some procedures concerning debugging,
! errors, and warnings.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Friday, December 17, 2021 PM04:33:25


#include "ppf.h"


module debug_mod

implicit none
private
public :: assert, backtr, warning, errstop, verisize

interface verisize
    module procedure verisize_real_1, verisize_real_2
    module procedure verisize_int_1, verisize_int_2
    module procedure verisize_logical_1, verisize_logical_2
end interface verisize


contains


subroutine assert(assertion, description, srname)
! This subroutine checks whether ASSERTION is true.
! If no, print SRNAME // 'Assertion failed: ' // DESCRIPTION
! and then stop the program by calling ERRSTOP.
! MATLAB analogue: assert(assertion, sprintf('%s: Assertion failed: %s', srname, description))
! Python analogue: assert assertion, srname + ': Assertion failed: ' + description
! C analogue: assert(assertion)  /* An error message will be produced by the compiler */
implicit none
logical, intent(in) :: assertion  ! A condition that is expected to be true
character(len=*), intent(in) :: description  ! Description of the assertion in human language
character(len=*), intent(in) :: srname  ! Name of the subroutine that calls this procedure

if (.not. assertion) then
    call errstop(trim(srname), 'Assertion failed: '//trim(description))
end if
end subroutine assert


subroutine warning(srname, mssg)
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: mssg

call backtr()
print '(/1A/)', 'WARNING: '//trim(srname)//': '//trim(mssg)//'.'
end subroutine warning


subroutine errstop(srname, mssg)
implicit none
character(len=*), intent(in) :: srname
character(len=*), intent(in) :: mssg

call backtr()
print '(/1A/)', 'ERROR: '//trim(srname)//': '//trim(mssg)//'.'
stop  ! This means to stop the whole program.
end subroutine errstop


subroutine backtr
! BACKTR calls a compiler-dependent intrinsic to show a backtrace if we
! are in the debuge mode, i.e., __DEBUGGING__ == 1.
! N.B.:
! 1. The intrinsic is compiler-dependent and does not exist in all
! compilers. Indeed, it is not standard-conforming. Therefore, compilers
! may warn that a non-standard intrinsic is in use.
! 2. More seriously, if the compiler is instructed to conform to the
! standards (e.g., gfortran with the option -std=f2003) while __DEBUGGING__
! is set to 1, then the compilation may FAIL at the linking stage,
! complaining that a subroutine cannot be found (e.g., backtrace for
! gfortran). In that case, we must set __DEBUGGING__ to 0 in ppf.h.
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


subroutine verisize_real_1(x, n)
! VERISIZE_REAL_1 verifies whether SIZE(X) = N.
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
real(RP), intent(in) :: x(:)
integer(IK), intent(in) :: n

character(len=*), parameter :: srname = 'VERISIZE_REAL_1'

if (int(size(x), IK) /= n) then
    call errstop(srname, 'SIZE(X) /= N')
end if
end subroutine verisize_real_1


subroutine verisize_real_2(x, m, n)
! VERISIZE_REAL_2 verifies whether SIZE(X, 1) = M, SIZE(X, 2) = N.
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
real(RP), intent(in) :: x(:, :)
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n

character(len=*), parameter :: srname = 'VERISIZE_REAL_2'

if (int(size(x, 1), IK) /= m) then
    call errstop(srname, 'SIZE(X, 1) /= M')
end if
if (int(size(x, 2), IK) /= n) then
    call errstop(srname, 'SIZE(X, 2) /= N')
end if
end subroutine verisize_real_2


subroutine verisize_int_1(x, n)
! VERISIZE_INT_1 verifies whether SIZE(X) = N.
use, non_intrinsic :: consts_mod, only : IK
implicit none
integer(IK), intent(in) :: x(:)
integer(IK), intent(in) :: n

character(len=*), parameter :: srname = 'VERISIZE_INT_1'

if (int(size(x), IK) /= n) then
    call errstop(srname, 'SIZE(X) /= N')
end if
end subroutine verisize_int_1


subroutine verisize_int_2(x, m, n)
! VERISIZE_INT_2 verifies whether SIZE(X, 1) = M, SIZE(X, 2) = N.
use, non_intrinsic :: consts_mod, only : IK
implicit none
integer(IK), intent(in) :: x(:, :)
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n

character(len=*), parameter :: srname = 'VERISIZE_INT_2'

if (int(size(x, 1), IK) /= m) then
    call errstop(srname, 'SIZE(X, 1) /= M')
end if
if (int(size(x, 2), IK) /= n) then
    call errstop(srname, 'SIZE(X, 2) /= N')
end if
end subroutine verisize_int_2


subroutine verisize_logical_1(x, n)
! VERISIZE_LOGICAL_1 verifies whether SIZE(X) = N.
use, non_intrinsic :: consts_mod, only : IK
implicit none
logical, intent(in) :: x(:)
integer(IK), intent(in) :: n

character(len=*), parameter :: srname = 'VERISIZE_LOGICAL_1'

if (int(size(x), IK) /= n) then
    call errstop(srname, 'SIZE(X) /= N')
end if
end subroutine verisize_logical_1


subroutine verisize_logical_2(x, m, n)
! VERISIZE_LOGICAL_2 verifies whether SIZE(X, 1) = M, SIZE(X, 2) = N.
use, non_intrinsic :: consts_mod, only : IK
implicit none
logical, intent(in) :: x(:, :)
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n

character(len=*), parameter :: srname = 'VERISIZE_LOGICAL_2'

if (int(size(x, 1), IK) /= m) then
    call errstop(srname, 'SIZE(X, 1) /= M')
end if
if (int(size(x, 2), IK) /= n) then
    call errstop(srname, 'SIZE(X, 2) /= N')
end if
end subroutine verisize_logical_2


end module debug_mod
