! PROB_MOD is a module defining the optimization problem. In particular, it implements CALFUN and
! CALCFC.  CALFUN evaluates the objective function for unconstrained, bound constrained, and
! linearly constrained problems; CALCFC evaluates the objective function and constraint for
! nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Tuesday, August 31, 2021 AM12:06:40


#include "fintrf.h"

module prob_mod

implicit none
private
public :: fun_ptr, calfun

mwPointer :: fun_ptr ! Pointer to objective function


contains


! The Fortran subroutine that evaluates the objective function.
subroutine calfun(x, f)

! Generic modules
use consts_mod, only : RP, HUGEFUN, MSSGLEN
use infnan_mod, only : is_nan
use fmxapi_mod, only : mxDestroyArray
use fmxapi_mod, only : mexErrMsgIdAndTxt
use fmxapi_mod, only : fmxIsDoubleScalar
use fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr, fmxCallMATLAB

implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Output
real(RP), intent(out) :: f

! Local variables
mwPointer :: pinput(1), poutput(1)
character(len=MSSGLEN) :: eid, mssg

! Associate X with INPUT(1)
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Verify the class and shape of outputs. Indeed, fmxReadMPtr does also the verification. We do it
! here in order to print a more informative error message when the verification fails.
if (.not. fmxIsDoubleScalar(poutput(1))) then
    eid = 'PROBLEM:ObjectiveNotScalar'
    mssg = 'PROBLEM: Objective function does not return a scalar.'
    call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
end if

! Read the data in OUTPUT
call fmxReadMPtr(poutput(1), f)

! Destroy the matrix created by fmxWriteMPtr for X. This must be done. Otherwise, the matrix will be
! destroyed only when the MEX function terminates. However, this subroutine will be called maybe
! thousands of times before that.
call mxDestroyArray(pinput(1))

end subroutine calfun


end module prob_mod
