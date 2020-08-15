! The mex gateway for NEWUOA (classical version)
!
! Coded by Zaikun Zhang in July 2020.

#include "fintrf.h"

module newuoacl_mod

use fmxcl_mod, only : RP_CL, IK_CL
implicit none
private
public :: fun_ptr, nf, fhist, xhist
public :: newuoa
public :: solver

! Some global veriables
! Pointer to bjective function
mwPointer :: fun_ptr 
! Number of function evaluations
integer(IK_CL) :: nf 
! History of evaluations
real(RP_CL), allocatable :: xhist(:, :), fhist(:)
! Solver name 
character(len = 6), parameter :: solver = 'NEWUOA'

interface
    subroutine newuoa(n, npt, x, rhobeg, rhoend, iprint, maxfun, w, f, info, ftarget)
    use fmxcl_mod, only : RP_CL, IK_CL
    implicit none
    integer(IK_CL), intent(in) :: n, npt, iprint, maxfun
    integer(IK_CL), intent(out) :: info
    real(RP_CL), intent(in) :: rhobeg, rhoend, ftarget
    real(RP_CL), intent(out) :: f
    real(RP_CL), intent(inout) :: x(n), w(*)  
    ! x(:) or w(:) does not  work !!!
    end subroutine newuoa
end interface

end module newuoacl_mod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Entry point to Fortran MEX function
subroutine mexFunction(nargout, poutput, nargin, pinput)
! If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
! extension depends on the platform), then the following function is
! callable in matlab:
! [xopt, fopt, info, nf, xhist, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, ftarget, maxfun, npt, iprint, maxhist, output_xhist)

! Generic modules
use consts_mod, only : MSSGLEN
use memory_mod, only : cstyle_sizeof
use fmxapi_mod, only : mexErrMsgIdAndTxt
use fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use fmxapi_mod, only : fmxVerifyClassShape
use fmxcl_mod, only : IK_CL, RP_CL, MAXMEMORY_CL
use fmxcl_mod, only : fmxAllocate, fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
use newuoacl_mod, only : fun_ptr, nf, xhist, fhist, newuoa, solver

implicit none

! mexFunction arguments (dummy variables)
! nargout and nargin are of type INTEGER in MATLAB 2019a documents
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Intermediate variables
integer :: maximal_hist
integer :: n_int
integer :: npt_int
integer :: nw
integer(IK_CL) :: info
integer(IK_CL) :: iprint
integer(IK_CL) :: khist
integer(IK_CL) :: maxfun
integer(IK_CL) :: maxfhist
integer(IK_CL) :: maxhist
integer(IK_CL) :: maxxhist
integer(IK_CL) :: n
integer(IK_CL) :: npt
integer(IK_CL) :: output_xhist
real(RP_CL) :: f
real(RP_CL) :: ftarget
real(RP_CL) :: rhobeg
real(RP_CL) :: rhoend
real(RP_CL), allocatable :: w(:)
real(RP_CL), allocatable :: x(:)
character(len = MSSGLEN) :: eid, mssg

! Validate number of arguments
call fmxVerifyNArgin(nargin, 10)
call fmxVerifyNArgout(nargout, 6)

! Verify that input 1 is a function handle; 
! the other inputs will be verified when read.
call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read inputs (there are 10)
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

! Get size
n = int(size(x), kind(n))

! Allocate workspace
n_int = int(n, kind(n_int))
npt_int = int(npt, kind(npt_int))
nw = (npt_int+13)*(npt_int+n_int)+3*n_int*(n_int+3)/2 + 1
if (nw > MAXMEMORY_CL/cstyle_sizeof(0.0_RP_CL)) then
    ! Without this checking, W may take too much memory, 
    ! or, more seriously, NW may overflow and cause a Segmentation Falt!
    eid = solver // ':WorkspaceTooLarge'
    mssg = solver // ': Workspace exceeds the largest memory allowed.'
    call mexErrMsgIdAndTxt(eid, mssg)
end if
call fmxAllocate(w, int(nw, IK_CL))

! Decide the maximal length of history according to MEXMEMORY in CONSTS_MOD
if (output_xhist > 0) then
    maximal_hist = MAXMEMORY_CL/((n+1)*cstyle_sizeof(0.0_RP_CL))
    maxxhist = max(0_IK_CL, min(maxfun, maxhist))
    ! We cannot simply take MAXXHIST = MIN(MAXXHIST, MAXIMAL_HIST),
    ! becaue they may not be the same kind, and compilers may complain.
    ! We may convert them to the same kind, but overflow may occur.
    if (maxxhist > maximal_hist) then
        maxxhist = int(maximal_hist, kind(maxxhist))
    end if
else
    maximal_hist = MAXMEMORY_CL/(cstyle_sizeof(0.0_RP_CL))
    maxxhist = 0
end if
maxfhist = max(0_IK_CL, min(maxfun, maxhist))
if (maxfhist > maximal_hist) then
    maxfhist = int(maximal_hist, kind(maxfhist))
end if

! Initialize global variables
nf = 0
call fmxAllocate(xhist, n, maxxhist)
call fmxAllocate(fhist, maxfhist)

! Call NEWUOA
call newuoa(n, npt, x, rhobeg, rhoend, iprint, maxfun, w, f, info, ftarget)
! If necessary, rearrange XHIST and FHIST so that they are in the 
! chronological order.
if (maxxhist >= 1 .and. maxxhist < nf) then
    khist = mod(nf - 1_IK_CL, maxxhist) + 1_IK_CL
    xhist = reshape((/ xhist(:, khist + 1 : maxxhist), xhist(:, 1 : khist) /), shape(xhist))
end if
if (maxfhist >= 1 .and. maxfhist < nf) then
    khist = mod(nf - 1_IK_CL, maxfhist) + 1_IK_CL
    fhist = (/ fhist(khist + 1 : maxfhist), fhist(1 : khist) /)
end if

! If MAXFHIST_IN >= NF_C > MAXFHIST_C, warn that not all history is recorced.
if (maxfhist < min(nf, maxhist)) then
    print '(/1A, I7, 1A)', 'WARNING: ' // solver // ': due to memory limit, MAXHIST is reset to ', maxfhist, '.'
    print '(1A/)', 'Only the history of the last MAXHIST iterations is recoreded.' 
end if

! Write outputs
call fmxWriteMPtr(x, poutput(1))
call fmxWriteMPtr(f, poutput(2))
call fmxWriteMPtr(info, poutput(3))
call fmxWriteMPtr(nf, poutput(4))
call fmxWriteMPtr(xhist(:, 1 : min(nf, maxxhist)), poutput(5))
call fmxWriteMPtr(fhist(1 : min(nf, maxfhist)), poutput(6), 'row')


! Free memory
deallocate (x)
deallocate (w)
deallocate (xhist)
deallocate (fhist)

return
end subroutine mexFunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The Fortran subroutine that evaluates the objective function
subroutine calfun(n, x, funval)

! Generic modules
use consts_mod, only : INT32_MEX, MSSGLEN
use fmxapi_mod, only : mxGetM, mxGetN, mxIsDouble
use fmxapi_mod, only : mxDestroyArray
use fmxapi_mod, only : mexErrMsgIdAndTxt
use fmxapi_mod, only : fmxCallMATLAB
use fmxcl_mod, only : RP_CL, IK_CL
use fmxcl_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
use newuoacl_mod, only : fun_ptr, nf, fhist, xhist, solver

implicit none

! Inputs
integer(IK_CL), intent(in) :: n
real(RP_CL), intent(in) :: x(n)

! Output
real(RP_CL), intent(out) :: funval

! Intermediate variables
mwPointer :: pinput(1), poutput(1) 
mwSize :: row, col
integer(IK_CL) :: maxfhist, maxxhist, khist
integer(INT32_MEX) :: isdble
character(len = MSSGLEN) :: eid, mssg

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Verify the class and shape of outputs. 
row = mxGetM(poutput(1)) 
col = mxGetN(poutput(1))
isdble = mxIsDouble(poutput(1))
if (row*col /= 1 .or. isdble /= 1) then
    eid = solver // ':ObjectiveNotScalar'
    mssg = solver // ': Objective function does not return a scalar.'
    call mexErrMsgIdAndTxt(eid, mssg)
end if

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), funval)

! Destroy the matrix created by fmxWriteMPtr for X. This must be done. 
call mxDestroyArray(pinput(1))  

! Update global variables
nf = nf + int(1, kind(nf))

maxxhist = int(size(xhist, 2), kind(maxxhist))
if (maxxhist >= 1) then
    khist = mod(nf - 1_IK_CL, maxxhist) + 1_IK_CL
    xhist(:, khist) = x 
end if

maxfhist = int(size(fhist), kind(maxfhist))
if (maxfhist >= 1) then
    khist = mod(nf - 1_IK_CL, maxfhist) + 1_IK_CL
    fhist(khist) = funval
end if

end subroutine calfun
