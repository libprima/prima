#include "fintrf.h"

module prob_mod
!--------------------------------------------------------------------------------------------------!
! PROB_MOD is a module defining the optimization problem. In particular, it implements CALFUN and
! CALCFC.  CALFUN evaluates the objective function for unconstrained, bound constrained, and
! linearly constrained problems; CALCFC evaluates the objective function and constraint for
! nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Tuesday, January 18, 2022 PM11:29:35
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : DP
implicit none
private
public :: fun_ptr, calfun
public :: funcon_ptr, calcfc
public :: xhist, fhist, chist, conhist
public :: nf

mwPointer :: fun_ptr ! Pointer to objective function, used by UOBYQA, NEWUOA, BOBYQA, and LINCOA
mwPointer :: funcon_ptr ! Pointer to objective&constraint functions, used by COBYLA
real(DP), allocatable :: xhist(:, :), fhist(:), chist(:), conhist(:, :)
integer :: nf


contains


subroutine calfun(n, x, f)
!--------------------------------------------------------------------------------------------------!
! The Fortran subroutine that evaluates the objective function in UOBYQA, NEWUOA, BOBYQA, and LINCOA
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


end module prob_mod
