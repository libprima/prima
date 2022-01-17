! PROB_MOD is a module defining the optimization problem. In particular, it implements CALFUN and
! CALCFC.  CALFUN evaluates the objective function for unconstrained, bound constrained, and
! linearly constrained problems; CALCFC evaluates the objective function and constraint for
! nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Tuesday, January 18, 2022 AM12:14:54


#include "fintrf.h"

module prob_mod

implicit none
private
public :: fun_ptr, calfun
public :: funcon_ptr, calcfc

mwPointer :: fun_ptr ! Pointer to objective function, used by UOBYQA, NEWUOA, BOBYQA, and LINCOA
mwPointer :: funcon_ptr ! Pointer to objective&constraint functions, used by COBYLA


contains


! The Fortran subroutine that evaluates the objective function in UOBYQA, NEWUOA, BOBYQA, and LINCOA
subroutine calfun(x, f)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : validate

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
use, non_intrinsic :: fmxapi_mod, only : fmxIsDoubleScalar
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr, fmxCallMATLAB

implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
character(len=*), parameter :: srname = 'CALFUN'
mwPointer :: pinput(1), poutput(1)

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Verify the class and shape of outputs. Indeed, fmxReadMPtr does also the verification. We do it
! here in order to print a more informative error message when the verification fails.
call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function returns a scalar', srname)

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)

! Destroy the matrix created by fmxWriteMPtr for X. This must be done. Otherwise, the matrix will be
! destroyed only when the MEX function terminates. However, this subroutine will be called maybe
! thousands of times before that.
call mxDestroyArray(pinput(1))

! Destroy POUTPUT(:).
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
call mxDestroyArray(poutput(1))  ! Destroy it even though it is a scalar

end subroutine calfun


! The Fortran subroutine that evaluates the objective&constraint functions in COBYLA.
subroutine calcfc(x, f, constr)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : validate
use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
use, non_intrinsic :: fmxapi_mod, only : fmxIsDoubleScalar, fmxIsDoubleVector
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr, fmxCallMATLAB

implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)

! Local variables
mwPointer :: pinput(1), poutput(2)
character(len=*), parameter :: srname = 'CALCFC'
real(RP), allocatable :: constr_loc(:)

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective&constraint functions
call fmxCallMATLAB(funcon_ptr, pinput, poutput)

! Verify the class and shape of outputs (even not debugging). Indeed, fmxReadMPtr does also the
! verification. We do it here in order to print a more informative error message in case of failure.
call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function returns a real scalar', srname)
call validate(fmxIsDoubleVector(poutput(2)), 'Constriant function returns a real vector', srname)

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)
call fmxReadMPtr(poutput(2), constr_loc)
! Check that the size of CONSTR_LOC is correct (even if not debugging)
call validate(size(constr_loc) == size(constr), 'SIZE(CONSTR_LOC) == SIZE(CONSTR)', srname)
constr = constr_loc

! Destroy the matrix created by fmxWriteMPtr for X. This must be done. Otherwise, the matrix will be
! destroyed only when the MEX function terminates. However, this subroutine will be called maybe
! thousands of times before that.
call mxDestroyArray(pinput(1))

! Destroy POUTPUT(:).
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
call mxDestroyArray(poutput(1))  ! Destroy it even though it is a scalar
call mxDestroyArray(poutput(2))

end subroutine calcfc


end module prob_mod
