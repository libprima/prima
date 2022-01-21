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
! Last Modified: Friday, January 21, 2022 AM01:58:09
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: fun_ptr, calfun
public :: funcon_ptr, calcfc

mwPointer :: fun_ptr ! Pointer to objective function, used by UOBYQA, NEWUOA, BOBYQA, and LINCOA
mwPointer :: funcon_ptr ! Pointer to objective&constraint functions, used by COBYLA


contains


subroutine calfun(x, f)
!--------------------------------------------------------------------------------------------------!
! The Fortran subroutine that evaluates the objective function in UOBYQA, NEWUOA, BOBYQA, and LINCOA
!--------------------------------------------------------------------------------------------------!

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
integer :: i
mwPointer :: pinput(1), poutput(1)

! Associate the input with INPUT.
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function.
call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Destroy the arrays in PINPUT(:).
! This must be done. Otherwise, the array created for X by fmxWriteMPtr will be destroyed only when
! the MEX function terminates, but this subroutine will be called maybe thousands of times before that.
do i = 1, size(pinput)
    call mxDestroyArray(pinput(i))
end do

! Read the data in POUTPUT.
! First, verify the class & shape of outputs (even not debugging). Indeed, fmxReadMPtr does also the
! verification. We do it here in order to print a more informative error message in case of failure.
call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function returns a scalar', srname)
! Second, copy the data.
call fmxReadMPtr(poutput(1), f)
! Third, destroy the arrays in POUTPUT.
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
do i = 1, size(poutput)
    call mxDestroyArray(poutput(i))
end do

end subroutine calfun


subroutine calcfc(x, f, constr)
!--------------------------------------------------------------------------------------------------!
! The Fortran subroutine that evaluates the objective&constraint functions in COBYLA.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : validate

! Fortran MEX API modules
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
character(len=*), parameter :: srname = 'CALCFC'
integer :: i
mwPointer :: pinput(1), poutput(2)
real(RP), allocatable :: constr_loc(:)

! Associate the input with PINPUT.
call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function.
call fmxCallMATLAB(funcon_ptr, pinput, poutput)

! Destroy the arrays in PINPUT.
! This must be done. Otherwise, the array created for X by fmxWriteMPtr will be destroyed only when
! the MEX function terminates, but this subroutine will be called maybe thousands of times before that.
do i = 1, size(pinput)
    call mxDestroyArray(pinput(i))
end do

! Read the data in POUTPUT.
! First, verify the class & shape of outputs (even not debugging). Indeed, fmxReadMPtr does also the
! verification. We do it here in order to print a more informative error message in case of failure.
call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function returns a real scalar', srname)
call validate(fmxIsDoubleVector(poutput(2)), 'Constriant function returns a real vector', srname)
! Second, copy the data.
call fmxReadMPtr(poutput(1), f)
call fmxReadMPtr(poutput(2), constr_loc)
! Third, destroy the arrays in POUTPUT.
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
do i = 1, size(poutput)
    call mxDestroyArray(poutput(i))
end do

! Copy CONSTR_LOC to CONSTR.
! Before copying, check that the size of CONSTR_LOC is correct (even if not debugging).
call validate(size(constr_loc) == size(constr), 'SIZE(CONSTR_LOC) == SIZE(CONSTR)', srname)
constr = constr_loc
! Deallocate CONSTR_LOC, allocated by fmxReadMPtr. Indeed, it would be deallocated automatically.
deallocate (constr_loc)

end subroutine calcfc


end module prob_mod
