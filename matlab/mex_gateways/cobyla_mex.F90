!--------------------------------------------------------------------------------------------------!
! The MEX gateway for COBYLA
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
! Last Modified: Friday, August 04, 2023 AM01:12:50
!--------------------------------------------------------------------------------------------------!

#include "fintrf.h"

subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! This is the entry point to the Fortran MEX function. If the compiled MEX file is named as
! FUNCTION_NAME.mex*** (extension depends on the platform), then in MATLAB we can call:
! [x, f, cstrv, nlconstr, info, nf, xhist, fhist, chist, nlchist] = ...
!   FUNCTION_NAME(funcon, x0, f0, nlconstr0, Aineq, bineq, Aeq, beq, lb, ub, rhobeg, rhoend, ...
!   eta1, eta2, gamma1, gamma2, ftarget, ctol, cweight, maxfun, iprint, maxhist, output_xhist, ...
!   output_nlchist, maxfilt)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: memory_mod, only : safealloc

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific modules
use, non_intrinsic :: cobyla_mod, only : cobyla

implicit none

! mexFunction arguments nargout and nargin are of type INTEGER in MATLAB 2019a documents.
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Local variables
integer(IK) :: info
integer(IK) :: iprint
integer(IK) :: m_nlcon
integer(IK) :: maxfilt
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: nf
logical :: output_nlchist
logical :: output_xhist
mwPointer :: funcon_ptr
real(RP) :: cstrv
real(RP) :: ctol
real(RP) :: cweight
real(RP) :: eta1
real(RP) :: eta2
real(RP) :: f
real(RP) :: f0
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
real(RP), allocatable :: nlchist(:, :)
real(RP), allocatable :: nlconstr(:)
real(RP), allocatable :: nlconstr0(:)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: lb(:)
real(RP), allocatable :: ub(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

! Validate the number of arguments
call fmxVerifyNArgin(nargin, 25)
call fmxVerifyNArgout(nargout, 10)

! Verify that input 1 is a function handle; the other inputs will be verified when read.
call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read the inputs
funcon_ptr = pinput(1)  ! FUNCON_PTR is a pointer to the function handle
call fmxReadMPtr(pinput(2), x)
call fmxReadMPtr(pinput(3), f0)
call fmxReadMPtr(pinput(4), nlconstr0)
call fmxReadMPtr(pinput(5), Aineq)
call fmxReadMPtr(pinput(6), bineq)
call fmxReadMPtr(pinput(7), Aeq)
call fmxReadMPtr(pinput(8), beq)
call fmxReadMPtr(pinput(9), lb)
call fmxReadMPtr(pinput(10), ub)
call fmxReadMPtr(pinput(11), rhobeg)
call fmxReadMPtr(pinput(12), rhoend)
call fmxReadMPtr(pinput(13), eta1)
call fmxReadMPtr(pinput(14), eta2)
call fmxReadMPtr(pinput(15), gamma1)
call fmxReadMPtr(pinput(16), gamma2)
call fmxReadMPtr(pinput(17), ftarget)
call fmxReadMPtr(pinput(18), ctol)
call fmxReadMPtr(pinput(19), cweight)
call fmxReadMPtr(pinput(20), maxfun)
call fmxReadMPtr(pinput(21), iprint)
call fmxReadMPtr(pinput(22), maxhist)
call fmxReadMPtr(pinput(23), output_xhist)
call fmxReadMPtr(pinput(24), output_nlchist)
call fmxReadMPtr(pinput(25), maxfilt)

! Get the sizes
m_nlcon = int(size(nlconstr0), kind(m_nlcon))  ! M_NLCON is a compulsory input of the Fortran code.

! Call the Fortran code
! There are different cases because XHIST/CONHIST may or may not be passed to the Fortran code.
if (output_xhist .and. output_nlchist) then
    call cobyla(calcfc, m_nlcon, x, f, cstrv, nlconstr, Aineq, bineq, Aeq, beq, lb, ub, &
        & f0, nlconstr0, nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, eta1, eta2, &
        & gamma1, gamma2, xhist=xhist, fhist=fhist, chist=chist, nlchist=nlchist, maxhist=maxhist, &
        & maxfilt=maxfilt, info=info)
elseif (output_xhist) then
    call cobyla(calcfc, m_nlcon, x, f, cstrv, nlconstr, Aineq, bineq, Aeq, beq, lb, ub, &
        & f0, nlconstr0, nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, eta1, eta2, &
        & gamma1, gamma2, xhist=xhist, fhist=fhist, chist=chist, maxhist=maxhist, maxfilt=maxfilt, &
        & info=info)
elseif (output_nlchist) then
    call cobyla(calcfc, m_nlcon, x, f, cstrv, nlconstr, Aineq, bineq, Aeq, beq, lb, ub, &
        & f0, nlconstr0, nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, eta1, eta2, &
        & gamma1, gamma2, fhist=fhist, chist=chist, nlchist=nlchist, maxhist=maxhist, &
        & maxfilt=maxfilt, info=info)
else
    call cobyla(calcfc, m_nlcon, x, f, cstrv, nlconstr, Aineq, bineq, Aeq, beq, lb, ub, &
        & f0, nlconstr0, nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, eta1, eta2, &
        & gamma1, gamma2, fhist=fhist, chist=chist, maxhist=maxhist, maxfilt=maxfilt, info=info)
end if

! After the Fortran code, XHIST or CONHIST may not be allocated, because it may not have been passed
! to the Fortran code. We allocate it here. Otherwise, fmxWriteMPtr will fail.
if (.not. allocated(xhist)) then
    call safealloc(xhist, int(size(x), IK), 0_IK)
end if
if (.not. allocated(nlchist)) then
    call safealloc(nlchist, int(size(nlconstr), IK), 0_IK)
end if

! Write the outputs
call fmxWriteMPtr(x, poutput(1))
call fmxWriteMPtr(f, poutput(2))
call fmxWriteMPtr(cstrv, poutput(3))
call fmxWriteMPtr(nlconstr, poutput(4))
call fmxWriteMPtr(info, poutput(5))
call fmxWriteMPtr(nf, poutput(6))
call fmxWriteMPtr(xhist(:, 1:min(int(nf), size(xhist, 2))), poutput(7))
call fmxWriteMPtr(fhist(1:min(int(nf), size(fhist))), poutput(8), 'row')
call fmxWriteMPtr(chist(1:min(int(nf), size(chist))), poutput(9), 'row')
call fmxWriteMPtr(nlchist(:, 1:min(int(nf), size(nlchist, 2))), poutput(10))
! N.B.:
! 1. INT(NF) converts NF to the default integer type; if not, MIN may complain.
! 2. It can happen that 0 < SIZE(XHIST, 2) < MAXHIST or 0 < SIZE(FHIST) < MAXHIST due to the memory
! limit in the Fortran code. Similar for CHIST and CONHIST.

! Free memory. Indeed, automatic deallocation would take place.
deallocate (x) ! Allocated by fmxReadMPtr.
deallocate (nlconstr0)  ! Allocated by fmxReadMPtr.
deallocate (Aineq) ! Allocated by fmxReadMPtr.
deallocate (bineq) ! Allocated by fmxReadMPtr.
deallocate (Aeq) ! Allocated by fmxReadMPtr.
deallocate (beq) ! Allocated by fmxReadMPtr.
deallocate (lb) ! Allocated by fmxReadMPtr.
deallocate (ub) ! Allocated by fmxReadMPtr.
deallocate (nlconstr)  ! Allocated by the solver
deallocate (xhist)  ! Allocated by the solver
deallocate (fhist)  ! Allocated by the solver
deallocate (chist)  ! Allocated by the solver
deallocate (nlchist)  ! Allocated by the solver

!--------------------------------------------------------------------!
contains

subroutine calcfc(x_sub, f_sub, nlconstr_sub)
! This is an internal procedure that defines CALCFC. Since F2008, we
! can pass internal procedures as actual arguments. See Note 12.18
! on page 290 of WD 1539-1 J3/10-007r1 (F2008 Working Document).
! We implement CALCFC internally so that FUNCON_PTR is visible to it.
! Do NOT pass FUNCON_PTR by a module variable, which is thread-unsafe.
use, non_intrinsic :: cbfun_mod, only : evalcb
implicit none
real(RP), intent(in) :: x_sub(:)
real(RP), intent(out) :: f_sub
real(RP), intent(out) :: nlconstr_sub(:)
call evalcb(funcon_ptr, x_sub, f_sub, nlconstr_sub)
end subroutine calcfc
!--------------------------------------------------------------------!
end subroutine mexFunction
