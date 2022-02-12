!--------------------------------------------------------------------------------------------------!
! The MEX gateway for UOBYQA
!
! **********************************************************************
!   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!               Department of Applied Mathematics,
!               The Hong Kong Polytechnic University
! **********************************************************************
!
! Started in March 2020
!
! Last Modified: Saturday, February 12, 2022 AM11:28:38
!--------------------------------------------------------------------------------------------------!

#include "fintrf.h"

subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! This is the entry point to the Fortran MEX function. If the compiled MEX file is named as
! FUNCTION_NAME.mex*** (extension depends on the platform), then in MATLAB we can call:
! [x, f, info, nf, xhist, fhist] = ...
!   FUNCTION_NAME(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, ...
!   iprint, maxhist, output_xhist)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: memory_mod, only : safealloc

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific modules
use, non_intrinsic :: uobyqa_mod, only : uobyqa

implicit none

! mexFunction arguments nargout and nargin are of type INTEGER in MATLAB 2019a documents.
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Local variables
integer(IK) :: info
integer(IK) :: iprint
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: nf
logical :: output_xhist
mwPointer :: fun_ptr
real(RP) :: eta1
real(RP) :: eta2
real(RP) :: f
real(RP) :: ftarget
real(RP) :: gamma1
real(RP) :: gamma2
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

! Validate the number of arguments
call fmxVerifyNArgin(nargin, 13)
call fmxVerifyNArgout(nargout, 6)

! Verify that input 1 is a function handle; the other inputs will be verified when read.
call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read the inputs
fun_ptr = pinput(1)  ! FUN_PTR is a pointer to the function handle
call fmxReadMPtr(pinput(2), x)
call fmxReadMPtr(pinput(3), rhobeg)
call fmxReadMPtr(pinput(4), rhoend)
call fmxReadMPtr(pinput(5), eta1)
call fmxReadMPtr(pinput(6), eta2)
call fmxReadMPtr(pinput(7), gamma1)
call fmxReadMPtr(pinput(8), gamma2)
call fmxReadMPtr(pinput(9), ftarget)
call fmxReadMPtr(pinput(10), maxfun)
call fmxReadMPtr(pinput(11), iprint)
call fmxReadMPtr(pinput(12), maxhist)
call fmxReadMPtr(pinput(13), output_xhist)

! Call the Fortran code
! There are different cases because XHIST may or may not be passed to the Fortran code.
if (output_xhist) then
    call uobyqa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint, &
        & eta1, eta2, gamma1, gamma2, xhist=xhist, fhist=fhist, maxhist=maxhist, info=info)
else
    call uobyqa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint, &
        & eta1, eta2, gamma1, gamma2, fhist=fhist, maxhist=maxhist, info=info)
end if

! After the Fortran code, XHIST may not be allocated, because it may not have been passed to the
! Fortran code. We allocate it here. Otherwise, fmxWriteMPtr will fail.
if (.not. allocated(xhist)) then
    call safealloc(xhist, int(size(x), IK), 0_IK)
end if

! Write the outputs
call fmxWriteMPtr(x, poutput(1))
call fmxWriteMPtr(f, poutput(2))
call fmxWriteMPtr(info, poutput(3))
call fmxWriteMPtr(nf, poutput(4))
call fmxWriteMPtr(xhist(:, 1:min(int(nf), size(xhist, 2))), poutput(5))
call fmxWriteMPtr(fhist(1:min(int(nf), size(fhist))), poutput(6), 'row')
! N.B.:
! 1. INT(NF) converts NF to the default integer type; if not, MIN may complain.
! 2. It can happen that 0 < SIZE(XHIST, 2) < MAXHIST or 0 < SIZE(FHIST) < MAXHIST due to the memory
! limit in the Fortran code.

! Free memory. Indeed, automatic deallocation would take place.
deallocate (x) ! Allocated by fmxReadMPtr.
deallocate (xhist)  ! Allocated by the solver
deallocate (fhist)  ! Allocated by the solver

!--------------------------------------------------------------------!
contains

subroutine calfun(x_sub, f_sub)
! This is an internal procedure that defines CALFUN. Since F2008, we
! can pass internal procedures as actual arguments. See Note 12.18
! on page 290 of WD 1539-1 J3/10-007r1 (F2008 Working Document).
! We implement CALFUN internally so that FUN_PTR is visible to it.
! Do NOT pass FUN_PTR by a module variable, which is thread-unsafe.
use, non_intrinsic :: cbfun_mod, only : evalcb
implicit none
real(RP), intent(in) :: x_sub(:)
real(RP), intent(out) :: f_sub
call evalcb(fun_ptr, x_sub, f_sub)
end subroutine calfun
!--------------------------------------------------------------------!
end subroutine mexFunction
