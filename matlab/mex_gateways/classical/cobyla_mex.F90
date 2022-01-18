!--------------------------------------------------------------------------------------------------!
! The MEX gateway for COBYLA (classical version)
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
! Last Modified: Tuesday, January 18, 2022 PM10:42:31
!--------------------------------------------------------------------------------------------------!

#include "fintrf.h"

module cobylacl_mod

use, non_intrinsic :: fmxcl_mod, only : RP_CL, IK_CL
implicit none
private
public :: funcon_ptr, nf, fhist, xhist, chist, conhist
public :: cobyla
public :: solver

! Pointer to objective function
mwPointer :: funcon_ptr
! Number of function evaluations
integer(IK_CL) :: nf
! History of evaluations
real(RP_CL), allocatable :: xhist(:, :), fhist(:), chist(:), conhist(:, :)
! Solver name
character(len=*), parameter :: solver = 'COBYLA'

interface
    subroutine cobyla(n, m, x, rhobeg, rhoend, iprint, maxfun, f, info, ftarget, cstrv, constr)
    use fmxcl_mod, only : RP_CL, IK_CL
    implicit none
    integer(IK_CL), intent(in) :: n, m, iprint, maxfun
    integer(IK_CL), intent(out) :: info
    real(RP_CL), intent(in) :: rhobeg, rhoend, ftarget
    real(RP_CL), intent(out) :: f, cstrv, constr(m)
    real(RP_CL), intent(inout) :: x(n)  ! x(:) will not work !!!
    end subroutine cobyla
end interface

end module cobylacl_mod


subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! This is the entry point to the Fortran MEX function. If the compiled MEX file is named as
! FUNCTION_NAME.mex*** (extension depends on the platform), then in MATLAB we can call:
! [x, f, cstrv, constr, info, nf, xhist, fhist, chist, conhist] = ...
!   FUNCTION_NAME(funcon, x0, f0, constr0, rhobeg, rhoend, ftarget, maxfun, iprint, ...
!   maxhist, output_xhist, output_conhist)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape
use, non_intrinsic :: fmxcl_mod, only : IK_CL, RP_CL, MAXMEMORY_CL
use, non_intrinsic :: fmxcl_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
use, non_intrinsic :: cobylacl_mod, only : funcon_ptr, nf, xhist, fhist, chist, conhist, cobyla, solver

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
integer(IK_CL) :: m
integer(IK_CL) :: maxchist
integer(IK_CL) :: maxconhist
integer(IK_CL) :: maxfhist
integer(IK_CL) :: maxfun
integer(IK_CL) :: maxhist
integer(IK_CL) :: maxhist_in
integer(IK_CL) :: maxxhist
integer(IK_CL) :: n
integer(IK_CL) :: unit_memo
logical :: output_conhist
logical :: output_xhist
real(RP_CL) :: cstrv
real(RP_CL) :: f
real(RP_CL) :: f0
real(RP_CL) :: ftarget
real(RP_CL) :: rhobeg
real(RP_CL) :: rhoend
real(RP_CL), allocatable :: x(:)
real(RP_CL), allocatable :: constr(:)
real(RP_CL), allocatable :: constr0(:)

! Validate the number of arguments
call fmxVerifyNArgin(nargin, 12)
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
call fmxReadMPtr(pinput(8), maxfun)
call fmxReadMPtr(pinput(9), iprint)
call fmxReadMPtr(pinput(10), maxhist)
call fmxReadMPtr(pinput(11), output_xhist)
call fmxReadMPtr(pinput(12), output_conhist)

! Get the sizes
n = int(size(x), kind(n))
m = int(size(constr0), kind(m))

! Decide the maximal amount of history to record
! MERGE(TSOURCE, FSOURCE, MASK) = TSOURCE if MASK is .TRUE., or FSOURCE if MASK is .FALSE.
unit_memo = int((n * merge(1, 0, output_xhist) + m * merge(1, 0, output_conhist) + 2) * cstyle_sizeof(0.0_RP_CL), kind(unit_memo))
maximal_hist = int(MAXMEMORY_CL / unit_memo, kind(maximal_hist))
maxhist_in = maxhist
maxhist = max(0_IK_CL, min(maxfun, maxhist))
if (maxhist > maximal_hist) then
    ! We cannot simply take MAXXHIST = MIN(MAXXHIST, MAXIMAL_HIST), as they may not be of the same
    ! kind, and compilers may complain. We may convert them to the same kind, but overflow may occur.
    maxhist = int(maximal_hist, kind(maxhist))
end if
maxfhist = maxhist
maxxhist = merge(maxhist, 0_IK_CL, output_xhist)
maxchist = maxhist
maxconhist = merge(maxhist, 0_IK_CL, output_conhist)

! Initialize NF and the history
nf = 0
call safealloc(xhist, int(n, IK), int(maxxhist, IK)) ! Not removable
call safealloc(fhist, int(maxfhist, IK)) ! Not removable
call safealloc(chist, int(maxchist, IK)) ! Not removable
call safealloc(conhist, int(m, IK), int(maxconhist, IK)) ! Not removable

! Allocate memory for CONSTR (X has been allocated by fmxReadMPtr)
call safealloc(constr, int(m, IK))

! Call the Fortran code
call cobyla(n, m, x, rhobeg, rhoend, iprint, maxfun, f, info, ftarget, cstrv, constr)

! If necessary, rearrange XHIST and FHIST so that they are in the chronological order.
if (maxxhist >= 1 .and. maxxhist < nf) then
    khist = modulo(nf - 1_IK_CL, maxxhist) + 1_IK_CL
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
end if
if (maxfhist >= 1 .and. maxfhist < nf) then
    khist = modulo(nf - 1_IK_CL, maxfhist) + 1_IK_CL
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if
if (maxchist >= 1 .and. maxchist < nf) then
    khist = modulo(nf - 1_IK_CL, maxchist) + 1_IK_CL
    chist = [chist(khist + 1:maxchist), chist(1:khist)]
end if
if (maxconhist >= 1 .and. maxconhist < nf) then
    khist = modulo(nf - 1_IK_CL, maxconhist) + 1_IK_CL
    conhist = reshape([conhist(:, khist + 1:maxconhist), conhist(:, 1:khist)], shape(conhist))
end if

! If MAXHIST_IN >= NF > MAXHIST, warn that not all history is recorded.
if (maxhist < min(nf, maxhist_in)) then
    print '(/1A, I7, 1A)', 'WARNING: '//solver//': due to memory limit, MAXHIST is reset to ', maxhist, '.'
    print '(1A/)', 'Only the history of the last MAXHIST iterations is recoreded.'
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

! Free memory. Indeed, automatic deallocation would take place.
deallocate (x) ! Allocated by fmxReadMPtr.
deallocate (constr)
deallocate (constr0)
deallocate (xhist)
deallocate (fhist)
deallocate (chist)
deallocate (conhist)
end subroutine mexFunction


subroutine calcfc(n, m, x, f, constr)
!--------------------------------------------------------------------------------------------------!
! The Fortran subroutine that evaluates the objective&constraint functions in COBYLA.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : HUGEFUN, HUGECON
use, non_intrinsic :: debug_mod, only : validate
use, non_intrinsic :: infnan_mod, only : is_nan

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
use, non_intrinsic :: fmxapi_mod, only : fmxCallMATLAB, fmxIsDoubleScalar, fmxIsDoubleVector
use, non_intrinsic :: fmxcl_mod, only : RP_CL, IK_CL
use, non_intrinsic :: fmxcl_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific modules
use, non_intrinsic :: cobylacl_mod, only : funcon_ptr, nf, fhist, xhist, chist, conhist, solver

implicit none

! Inputs
integer(IK_CL), intent(in) :: n
integer(IK_CL), intent(in) :: m
real(RP_CL), intent(in) :: x(n)

! Outputs
real(RP_CL), intent(out) :: f
real(RP_CL), intent(out) :: constr(m)

! Local variables
character(len=*), parameter :: srname = 'CALCFC'
integer(IK_CL) :: maxchist, maxconhist, maxfhist, maxxhist, khist
mwPointer :: pinput(1), poutput(2)
real(RP_CL) :: cstrv
real(RP_CL), allocatable :: constr_loc(:)

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
call fmxCallMATLAB(funcon_ptr, pinput, poutput)

! Verify the class and shape of outputs (even if not debugging).
call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function returns a real scalar', solver)
call validate(fmxIsDoubleVector(poutput(2)), 'Constriant function returns a real vector', solver)

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)
call fmxReadMPtr(poutput(2), constr_loc)
! Check that the size of CONSTR_LOC is correct (even if not debugging)
call validate(size(constr_loc) == size(constr), 'SIZE(CONSTR_LOC) == SIZE(CONSTR)', srname)
constr = constr_loc

! Destroy the matrix created by fmxWriteMPtr for X. This must be done.
call mxDestroyArray(pinput(1))
! Destroy POUTPUT(:).
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
call mxDestroyArray(poutput(1))  ! Destroy it even though it is a scalar
call mxDestroyArray(poutput(2))

! Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
! finite value. This is naive, and better approaches surely exist.
if (f > HUGEFUN .or. is_nan(f)) then
    f = HUGEFUN
end if
where (constr < -HUGECON .or. is_nan(constr))
    ! The constraint is CONSTR(X) >= 0, so NaN should be replaced with a large negative value.
    constr = -HUGECON  ! MATLAB code: constr(constr < -HUGECON | isnan(constr)) = -HUGECON
end where

! Moderate huge positive values of CONSTR, or they may lead to Inf/NaN in subsequent calculations.
! This is NOT an extreme barrier.
!constr = min(HUGECON, constr)
where (constr > HUGECON)
    constr = HUGECON
end where
!! We may moderate F similarly, but we decide not to.
!!f = max(-HUGEFUN, f)

! Evaluate the constraint violation for constraints CONSTR(X) >= 0.
!cstrv = maxval([-constr, ZERO])
cstrv = maxval([-constr, 0.0_RP_CL])

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
maxchist = int(size(chist), kind(maxchist))
if (maxchist >= 1) then
    khist = modulo(nf - 1_IK_CL, maxchist) + 1_IK_CL
    chist(khist) = cstrv
end if
maxconhist = int(size(conhist, 2), kind(maxconhist))
if (maxconhist >= 1) then
    khist = modulo(nf - 1_IK_CL, maxconhist) + 1_IK_CL
    conhist(:, khist) = constr
end if

end subroutine calcfc
