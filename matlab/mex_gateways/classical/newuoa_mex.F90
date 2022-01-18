!--------------------------------------------------------------------------------------------------!
! The MEX gateway for NEWUOA (classical version)
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
! Last Modified: Tuesday, January 18, 2022 AM12:13:00
!--------------------------------------------------------------------------------------------------!

#include "fintrf.h"

module newuoacl_mod

use, non_intrinsic :: fmxcl_mod, only : RP_CL, IK_CL
implicit none
private
public :: fun_ptr, nf, fhist, xhist
public :: newuoa
public :: solver

! Pointer to objective function
mwPointer :: fun_ptr
! Number of function evaluations
integer(IK_CL) :: nf
! History of evaluations
real(RP_CL), allocatable :: xhist(:, :), fhist(:)
! Solver name
character(len=*), parameter :: solver = 'NEWUOA'

interface
    subroutine newuoa(n, npt, x, rhobeg, rhoend, iprint, maxfun, f, info, ftarget)
    use fmxcl_mod, only : RP_CL, IK_CL
    implicit none
    integer(IK_CL), intent(in) :: n, npt, iprint, maxfun
    integer(IK_CL), intent(out) :: info
    real(RP_CL), intent(in) :: rhobeg, rhoend, ftarget
    real(RP_CL), intent(out) :: f
    real(RP_CL), intent(inout) :: x(n)  ! x(:) will not work !!!
    end subroutine newuoa
end interface

end module newuoacl_mod


subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! This is the entry point to the Fortran MEX function. If the compiled MEX file is named as
! FUNCTION_NAME.mex*** (extension depends on the platform), then in MATLAB we can call:
! [x, f, info, nf, xhist, fhist] = ...
!   FUNCTION_NAME(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, npt, ...
!   iprint, maxhist, output_xhist)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape
use, non_intrinsic :: fmxcl_mod, only : IK_CL, RP_CL, MAXMEMORY_CL
use, non_intrinsic :: fmxcl_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific modules
use, non_intrinsic :: newuoacl_mod, only : fun_ptr, nf, xhist, fhist, newuoa, solver

implicit none

! mexFunction arguments (dummy variables)
! nargout and nargin are of type INTEGER in MATLAB 2019a documents
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Local variables
integer :: maximal_hist
integer(IK_CL) :: info
integer(IK_CL) :: iprint
integer(IK_CL) :: khist
integer(IK_CL) :: maxfhist
integer(IK_CL) :: maxfun
integer(IK_CL) :: maxhist
integer(IK_CL) :: maxxhist
integer(IK_CL) :: n
integer(IK_CL) :: npt
logical :: output_xhist
real(RP_CL) :: f
real(RP_CL) :: ftarget
real(RP_CL) :: rhobeg
real(RP_CL) :: rhoend
real(RP_CL), allocatable :: x(:)

! Validate number of arguments
call fmxVerifyNArgin(nargin, 10)
call fmxVerifyNArgout(nargout, 6)

! Verify that input 1 is a function handle; the other inputs will be verified when read.
call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read the inputs
fun_ptr = pinput(1)  ! Pointer to the function handle
call fmxReadMPtr(pinput(2), x)
call fmxReadMPtr(pinput(3), rhobeg)
call fmxReadMPtr(pinput(4), rhoend)
call fmxReadMPtr(pinput(5), ftarget)
call fmxReadMPtr(pinput(6), maxfun)
call fmxReadMPtr(pinput(7), npt)
call fmxReadMPtr(pinput(8), iprint)
call fmxReadMPtr(pinput(9), maxhist)
call fmxReadMPtr(pinput(10), output_xhist)

! Get the size
n = int(size(x), kind(n))

! Decide the maximal length of history according to MEXMEMORY in CONSTS_MOD
if (output_xhist) then
    maximal_hist = int(MAXMEMORY_CL / ((n + 1) * cstyle_sizeof(0.0_RP_CL)), kind(maximal_hist))
    maxxhist = max(0_IK_CL, min(maxfun, maxhist))
    ! We cannot simply take MAXXHIST = MIN(MAXXHIST, MAXIMAL_HIST),
    ! as they may not be the same kind, and compilers may complain.
    ! We may convert them to the same kind, but overflow may occur.
    if (maxxhist > maximal_hist) then
        maxxhist = int(maximal_hist, kind(maxxhist))
    end if
else
    maximal_hist = int(MAXMEMORY_CL / (cstyle_sizeof(0.0_RP_CL)), kind(maximal_hist))
    maxxhist = 0
end if
maxfhist = max(0_IK_CL, min(maxfun, maxhist))
if (maxfhist > maximal_hist) then
    maxfhist = int(maximal_hist, kind(maxfhist))
end if

! Initialize NF and the history
nf = 0
call safealloc(xhist, int(n, IK), int(maxxhist, IK)) ! Not removable
call safealloc(fhist, int(maxfhist, IK)) ! Not removable

! Call the Fortran code
call newuoa(n, npt, x, rhobeg, rhoend, iprint, maxfun, f, info, ftarget)

! If necessary, rearrange XHIST and FHIST so that they are in the chronological order.
if (maxxhist >= 1 .and. maxxhist < nf) then
    khist = modulo(nf - 1_IK_CL, maxxhist) + 1_IK_CL
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
end if
if (maxfhist >= 1 .and. maxfhist < nf) then
    khist = modulo(nf - 1_IK_CL, maxfhist) + 1_IK_CL
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if

! If MAXFHIST_IN >= NF_C > MAXFHIST_C, warn that not all history is recorded.
if (maxfhist < min(nf, maxhist)) then
    print '(/1A, I7, 1A)', 'WARNING: '//solver//': due to memory limit, MAXHIST is reset to ', maxfhist, '.'
    print '(1A/)', 'Only the history of the last MAXHIST iterations is recoreded.'
end if

! Write the outputs
call fmxWriteMPtr(x, poutput(1))
call fmxWriteMPtr(f, poutput(2))
call fmxWriteMPtr(info, poutput(3))
call fmxWriteMPtr(nf, poutput(4))
call fmxWriteMPtr(xhist(:, 1:min(nf, maxxhist)), poutput(5))
call fmxWriteMPtr(fhist(1:min(nf, maxfhist)), poutput(6), 'row')

! Free memory. Indeed, automatic deallocation would take place.
deallocate (x)
deallocate (xhist)
deallocate (fhist)
end subroutine mexFunction


subroutine calfun(n, x, f)
!--------------------------------------------------------------------------------------------------!
! This is the Fortran subroutine that evaluates the objective function.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: debug_mod, only : validate

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
use, non_intrinsic :: fmxapi_mod, only : fmxCallMATLAB, fmxIsDoubleScalar
use, non_intrinsic :: fmxcl_mod, only : RP_CL, IK_CL
use, non_intrinsic :: fmxcl_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
use, non_intrinsic :: newuoacl_mod, only : fun_ptr, nf, fhist, xhist, solver

implicit none

! Inputs
integer(IK_CL), intent(in) :: n
real(RP_CL), intent(in) :: x(n)

! Outputs
real(RP_CL), intent(out) :: f

! Local variables
integer(IK_CL) :: maxfhist, maxxhist, khist
mwPointer :: pinput(1), poutput(1)

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Verify the class and shape of outputs (even if not debugging)
call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function returns a scalar', solver)

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)

! Destroy the matrix created by fmxWriteMPtr for X. This must be done.
call mxDestroyArray(pinput(1))
! Destroy POUTPUT(:).
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
call mxDestroyArray(poutput(1))  ! Destroy it even though it is a scalar

! Update global variables
nf = nf + int(1, kind(nf))

maxxhist = int(size(xhist, 2), kind(maxxhist))
if (maxxhist >= 1) then
    khist = modulo(nf - 1_IK_CL, maxxhist) + 1_IK_CL
    xhist(:, khist) = x
end if
maxfhist = int(size(fhist), kind(maxfhist))
if (maxfhist >= 1) then
    khist = modulo(nf - 1_IK_CL, maxfhist) + 1_IK_CL
    fhist(khist) = f
end if

end subroutine calfun
