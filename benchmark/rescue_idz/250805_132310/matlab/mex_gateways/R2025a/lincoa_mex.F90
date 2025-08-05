!--------------------------------------------------------------------------------------------------!
! The MEX gateway for LINCOA
! This is an adapted MEX gateway to circumvent the MATLAB R2025a bug that it segfaults on Linux
! if the Fortran MEX function contains an internal procedure that is passed as an actual argument.
! This MEX uses a module variable FUN_PTR to store the function handle, which is essentially a
! global variable and is not thread-safe or recursion-safe.
! See MathWorks Technical Support Case 07931486 and
! https://www.mathworks.com/matlabcentral/answers/2178414-bug-matlab-2025a-segfaults-on-ubuntu-when-handling-fortran-mex-files-with-internal-subroutines
! https://stackoverflow.com/questions/79699706/matlab-2025a-vs-fortran-mex-files-with-internal-subroutines
! https://fortran-lang.discourse.group/t/implementation-of-a-parametrized-objective-function-without-using-module-variables-or-internal-subroutines
! https://stackoverflow.com/questions/79705107/fortran-implementating-a-parametrized-objective-function-without-using-module-v
!
! Authors:
!   Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!   and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!   Department of Applied Mathematics, 
!   The Hong Kong Polytechnic University
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015)
!
! Started in July 2020
!
! Last Modified: Sun 20 Jul 2025 08:31:41 AM PDT
!--------------------------------------------------------------------------------------------------!

#include "fintrf.h"


module calfun_mod
implicit none
private
public :: fun_ptr, calfun

mwPointer :: fun_ptr  ! Pointer to the objective function handle

contains

subroutine calfun(x, f)
use, non_intrinsic :: cbfun_mod, only : evalcb
use, non_intrinsic :: consts_mod, only : RP
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
call evalcb(fun_ptr, x, f)
end subroutine calfun

end module calfun_mod


subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! This is the entry point to the Fortran MEX function. If the compiled MEX file is named as
! FUNCTION_NAME.mex*** (extension depends on the platform), then in MATLAB we can call:
! [x, f, cstrv, info, nf, xhist, fhist, chist] = ...
!   FUNCTION_NAME(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, rhobeg, rhoend, eta1, eta2, ...
!   gamma1, gamma2, ftarget, ctol, cweight, maxfun, npt, iprint, maxhist, output_xhist, maxfilt)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: memory_mod, only : safealloc

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific modules
use, non_intrinsic :: calfun_mod, only : fun_ptr, calfun
use, non_intrinsic :: lincoa_mod, only : lincoa

implicit none

! mexFunction arguments nargout and nargin are of type INTEGER in MATLAB 2019a documents.
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Local variables
integer(IK) :: info
integer(IK) :: iprint
integer(IK) :: maxfilt
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: nf
integer(IK) :: npt
logical :: output_xhist
real(RP) :: cstrv
real(RP) :: ctol
real(RP) :: cweight
real(RP) :: eta1
real(RP) :: eta2
real(RP) :: f
real(RP) :: ftarget
real(RP) :: gamma1
real(RP) :: gamma2
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: Aeq(:, :)
real(RP), allocatable :: Aineq(:, :)
real(RP), allocatable :: beq(:)
real(RP), allocatable :: bineq(:)
real(RP), allocatable :: chist(:)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: lb(:)
real(RP), allocatable :: ub(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

! Validate the number of arguments
call fmxVerifyNArgin(nargin, 23)
call fmxVerifyNArgout(nargout, 8)

! Verify that input 1 is a function handle; the other inputs will be verified when read.
call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read the inputs
fun_ptr = pinput(1)  ! FUN_PTR is a pointer to the function handle
call fmxReadMPtr(pinput(2), x)
call fmxReadMPtr(pinput(3), Aineq)
call fmxReadMPtr(pinput(4), bineq)
call fmxReadMPtr(pinput(5), Aeq)
call fmxReadMPtr(pinput(6), beq)
call fmxReadMPtr(pinput(7), lb)
call fmxReadMPtr(pinput(8), ub)
call fmxReadMPtr(pinput(9), rhobeg)
call fmxReadMPtr(pinput(10), rhoend)
call fmxReadMPtr(pinput(11), eta1)
call fmxReadMPtr(pinput(12), eta2)
call fmxReadMPtr(pinput(13), gamma1)
call fmxReadMPtr(pinput(14), gamma2)
call fmxReadMPtr(pinput(15), ftarget)
call fmxReadMPtr(pinput(16), ctol)
call fmxReadMPtr(pinput(17), cweight)
call fmxReadMPtr(pinput(18), maxfun)
call fmxReadMPtr(pinput(19), npt)
call fmxReadMPtr(pinput(20), iprint)
call fmxReadMPtr(pinput(21), maxhist)
call fmxReadMPtr(pinput(22), output_xhist)
call fmxReadMPtr(pinput(23), maxfilt)

! Call the Fortran code
! There are different cases because XHIST may or may not be passed to the Fortran code.
if (output_xhist) then
    call lincoa(calfun, x, f, cstrv, Aineq, bineq, Aeq, beq, lb, ub, nf, rhobeg, rhoend, ftarget, &
        & ctol, cweight, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, xhist = xhist, fhist = fhist, &
        & chist = chist, maxhist = maxhist, maxfilt = maxfilt, info = info)
else
    call lincoa(calfun, x, f, cstrv, Aineq, bineq, Aeq, beq, lb, ub, nf, rhobeg, rhoend, ftarget, &
        & ctol, cweight, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, fhist = fhist, chist = chist, &
        & maxhist = maxhist, maxfilt = maxfilt, info = info)
end if

! After the Fortran code, XHIST may not be allocated, because it may not have been passed to the
! Fortran code. We allocate it here. Otherwise, fmxWriteMPtr will fail.
if (.not. allocated(xhist)) then
    call safealloc(xhist, int(size(x), IK), 0_IK)
end if

! Write the outputs
call fmxWriteMPtr(x, poutput(1))
call fmxWriteMPtr(f, poutput(2))
call fmxWriteMPtr(cstrv, poutput(3))
call fmxWriteMPtr(info, poutput(4))
call fmxWriteMPtr(nf, poutput(5))
call fmxWriteMPtr(xhist(:, 1:min(nf, int(size(xhist, 2), IK))), poutput(6))
call fmxWriteMPtr(fhist(1:min(nf, int(size(fhist), IK))), poutput(7), 'row')
call fmxWriteMPtr(chist(1:min(nf, int(size(chist), IK))), poutput(8), 'row')
! N.B.:
! It can happen that 0 < SIZE(XHIST, 2) < MAXHIST or 0 < SIZE(FHIST) < MAXHIST due to the memory
! limit in the Fortran code. Similar for CHIST.

! Free memory. We prefer explicit deallocation to the automatic one.
deallocate (x)  ! Allocated by fmxReadMPtr.
deallocate (Aineq)  ! Allocated by fmxReadMPtr.
deallocate (bineq)  ! Allocated by fmxReadMPtr.
deallocate (Aeq)  ! Allocated by fmxReadMPtr.
deallocate (beq)  ! Allocated by fmxReadMPtr.
deallocate (lb)  ! Allocated by fmxReadMPtr.
deallocate (ub)  ! Allocated by fmxReadMPtr.
deallocate (xhist)  ! Allocated by the solver
deallocate (fhist)  ! Allocated by the solver
deallocate (chist)  ! Allocated by the solver

end subroutine mexFunction
