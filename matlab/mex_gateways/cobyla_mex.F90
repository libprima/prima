!--------------------------------------------------------------------------------------------------!
! The MEX gateway for COBYLA
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
! Last Modified: Tuesday, January 18, 2022 PM09:09:13
!--------------------------------------------------------------------------------------------------!

#include "fintrf.h"

subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! This is the entry point to the Fortran MEX function. If the compiled MEX file is named as
! FUNCTION_NAME.mex*** (extension depends on the platform), then in MATLAB we can call:
! [x, f, cstrv, constr, info, nf, xhist, fhist, chist, conhist] = ...
!   FUNCTION_NAME(funcon, x0, f0, constr0, rhobeg, rhoend, ftarget, ctol, maxfun, iprint, ...
!   maxhist, output_xhist, output_conhist, maxfilt)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK
!use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific modules
use, non_intrinsic :: cobyla_mod, only : cobyla
use, non_intrinsic :: prob_mod, only : funcon_ptr, calcfc

implicit none

! mexFunction arguments nargout and nargin are of type INTEGER in MATLAB 2019a documents.
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Local variables
!real(RP) :: eta1
!real(RP) :: eta2
!real(RP) :: gamma1
!real(RP) :: gamma2
integer(IK) :: info
integer(IK) :: iprint
integer(IK) :: m
integer(IK) :: maxfilt
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: nf
logical :: output_conhist
logical :: output_xhist
real(RP) :: cstrv
real(RP) :: ctol
real(RP) :: f
real(RP) :: f0
real(RP) :: ftarget
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: chist(:)
real(RP), allocatable :: conhist(:, :)
real(RP), allocatable :: constr(:)
real(RP), allocatable :: constr0(:)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

! Validate the number of arguments
call fmxVerifyNArgin(nargin, 14)
call fmxVerifyNArgout(nargout, 10)

! Verify that input 1 is a function handle; the other inputs will be verified when read.
call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read the inputs
funcon_ptr = pinput(1)  ! FUNCON_PTR is a pointer to the function handle
call fmxReadMPtr(pinput(2), x)
call fmxReadMPtr(pinput(3), f0)
call fmxReadMPtr(pinput(4), constr0)
call fmxReadMPtr(pinput(5), rhobeg)
call fmxReadMPtr(pinput(6), rhoend)
call fmxReadMPtr(pinput(7), ftarget)
call fmxReadMPtr(pinput(8), ctol)
call fmxReadMPtr(pinput(9), maxfun)
call fmxReadMPtr(pinput(10), iprint)
call fmxReadMPtr(pinput(11), maxhist)
call fmxReadMPtr(pinput(12), output_xhist)
call fmxReadMPtr(pinput(13), output_conhist)
call fmxReadMPtr(pinput(14), maxfilt)

! Get the sizes
m = int(size(constr0), kind(m))  ! M is a compulsory input of the Fortran code.

! Call the Fortran code
! There are different cases because XHIST/CONHIST may or may not be passed to the Fortran code.
if (output_xhist .and. output_conhist) then
    call cobyla(calcfc, m, x, f, cstrv, constr, f0, constr0, nf, rhobeg, rhoend, ftarget, ctol, &
        & maxfun, iprint, xhist=xhist, fhist=fhist, chist=chist, conhist=conhist, maxhist=maxhist, maxfilt=maxfilt, info=info)
elseif (output_xhist) then
    call cobyla(calcfc, m, x, f, cstrv, constr, f0, constr0, nf, rhobeg, rhoend, ftarget, ctol, &
        & maxfun, iprint, xhist=xhist, fhist=fhist, chist=chist, maxhist=maxhist, maxfilt=maxfilt, info=info)
elseif (output_conhist) then
    call cobyla(calcfc, m, x, f, cstrv, constr, f0, constr0, nf, rhobeg, rhoend, ftarget, ctol, &
        & maxfun, iprint, fhist=fhist, chist=chist, conhist=conhist, maxhist=maxhist, maxfilt=maxfilt, info=info)
else
    call cobyla(calcfc, m, x, f, cstrv, constr, f0, constr0, nf, rhobeg, rhoend, ftarget, ctol, &
        & maxfun, iprint, fhist=fhist, chist=chist, maxhist=maxhist, maxfilt=maxfilt, info=info)
end if

! After the Fortran code, XHIST or CONHIST may not be allocated, because it may not have been passed
! to the Fortran code. We allocate it here. Otherwise, fmxWriteMPtr will fail.
if (.not. allocated(xhist)) then
    call safealloc(xhist, int(size(x), IK), 0_IK)
end if
if (.not. allocated(conhist)) then
    call safealloc(conhist, int(size(constr), IK), 0_IK)
end if

! Write the outputs
call fmxWriteMPtr(x, poutput(1))
call fmxWriteMPtr(f, poutput(2))
call fmxWriteMPtr(cstrv, poutput(3))
call fmxWriteMPtr(constr, poutput(4))
call fmxWriteMPtr(info, poutput(5))
call fmxWriteMPtr(nf, poutput(6))
call fmxWriteMPtr(xhist(:, 1:min(int(nf), size(xhist, 2))), poutput(7))
call fmxWriteMPtr(fhist(1:min(int(nf), size(fhist))), poutput(8), 'row')
call fmxWriteMPtr(chist(1:min(int(nf), size(chist))), poutput(9), 'row')
call fmxWriteMPtr(conhist(:, 1:min(int(nf), size(conhist, 2))), poutput(10))
! N.B.:
! 1. INT(NF) converts NF to the default integer type; if not, MIN may complain.
! 2. It can happen that 0 < SIZE(XHIST, 2) < MAXHIST or 0 < SIZE(FHIST) < MAXHIST due to the memory
! limit in the Fortran code. Similar for CHIST and CONHIST.

! Free memory. Indeed, automatic deallocation would take place.
deallocate (x) ! Allocated by fmxReadMPtr.
deallocate (constr)
deallocate (constr0)
deallocate (xhist)
deallocate (fhist)
deallocate (chist)
deallocate (conhist)
end subroutine mexFunction
