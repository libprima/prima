! PROB_MOD is a module defining the optimization problem. In particular, it implements CALFUN and
! CALCFC.  CALFUN evaluates the objective function for unconstrained, bound constrained, and
! linearly constrained problems; CALCFC evaluates the objective function and constraint for
! nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Wednesday, January 12, 2022 PM08:35:53


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
use, non_intrinsic :: consts_mod, only : RP, MSGLEN
use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
use, non_intrinsic :: fmxapi_mod, only : mexErrMsgIdAndTxt
use, non_intrinsic :: fmxapi_mod, only : fmxIsDoubleScalar
use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr, fmxCallMATLAB

implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
mwPointer :: pinput(1), poutput(1)
character(len=MSGLEN) :: eid, msg

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Verify the class and shape of outputs. Indeed, fmxReadMPtr does also the verification. We do it
! here in order to print a more informative error message when the verification fails.
if (.not. fmxIsDoubleScalar(poutput(1))) then
    eid = 'PROBLEM:ObjectiveNotScalar'
    msg = 'PROBLEM: Objective value is not a scalar.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)

! Destroy the matrix created by fmxWriteMPtr for X. This must be done. Otherwise, the matrix will be
! destroyed only when the MEX function terminates. However, this subroutine will be called maybe
! thousands of times before that.
call mxDestroyArray(pinput(1))

end subroutine calfun


! The Fortran subroutine that evaluates the objective&constraint functions in COBYLA.
subroutine calcfc(x, f, constr)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, MSGLEN
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
use, non_intrinsic :: fmxapi_mod, only : mexErrMsgIdAndTxt
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
character(len=MSGLEN) :: eid, msg
real(RP), allocatable :: constr_loc(:)

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective&constraint functions
call fmxCallMATLAB(funcon_ptr, pinput, poutput)

! Verify the class and shape of outputs. Indeed, fmxReadMPtr does also the verification. We do it
! here in order to print a more informative error message when the verification fails.
if (.not. fmxIsDoubleScalar(poutput(1))) then
    eid = 'PROBLEM:ObjectiveNotScalar'
    msg = 'PROBLEM: Objective value is not a scalar.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if
if (.not. fmxIsDoubleVector(poutput(2))) then
    eid = 'PROBLEM:ConstraintNotVector'
    msg = 'PROBLEM: Constraint value is not a vector.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)
call fmxReadMPtr(poutput(2), constr_loc)
call assert(size(constr_loc) == size(constr), 'SIZE(CONSTR_LOC) == SIZE(CONSTR)', srname)
constr = constr_loc

! Destroy the matrix created by fmxWriteMPtr for X. This must be done. Otherwise, the matrix will be
! destroyed only when the MEX function terminates. However, this subroutine will be called maybe
! thousands of times before that.
call mxDestroyArray(pinput(1))

! Should we destroy POUTPUT(2)?
call mxDestroyArray(poutput(2))
!MATLAB allocates dynamic memory to store the arrays in plhs for mexCallMATLAB. MATLAB automatically
!deallocates the dynamic memory when you exit the MEX file. However, if heap space is at a premium,
!call mxDestroyArray when you are finished with the arrays in plhs.
!https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html

end subroutine calcfc


end module prob_mod
