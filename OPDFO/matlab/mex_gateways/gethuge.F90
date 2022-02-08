!--------------------------------------------------------------------------------------------------!
! GETHUGE subroutine
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
! Last Modified: Saturday, January 22, 2022 PM08:16:50
!--------------------------------------------------------------------------------------------------!

! N.B.:
!
! 1. Be careful with the "kind" and storage size for integer-type (integer, mwSize, mwIndex)
! variables/functions. Some of them may be 32bit, while the others may be 64bit, depending on the
! machine, the version of matlab, and the compilation option of mex. Do NOT assume any two of them
! to be the same. If ever a Segmentation Fault occurs, check these variables first.
!
! 2. Be careful with the line width limit. After preprocessing (macro substitution), some lines may
! become too long and hence get truncated.

#include "fintrf.h"

subroutine mexFunction(nargout, poutput, nargin, pinput)
!--------------------------------------------------------------------------------------------------!
! Usage: data_huge = gethuge(data_type)
! This function returns the largest value of data_type on the current platform. The possible values
! of data_type are 'integer',  'int', 'float', 'real', 'single', 'double', 'mwSI', 'fun', 'function',
! 'con', 'constraint'.
!
! In previous versions, GETHUGE accepted 'mwSize' or 'mwIndex' as inputs, but it is not the case
! anymore. This is because mwSize and mwIndex are macros defined in fintrf.h, and they will be
! replaced by other strings after preprocessing. Instead of 'mwSize' and 'mwIndex', we now use
! 'mwSI' to get the smaller value between huge(msZero) and huge(miZero).
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, SP, DP, RP, HUGEFUN, HUGECON, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, validate
use, non_intrinsic :: string_mod, only : lower

! Fortran MEX API modules
use, non_intrinsic :: fmxapi_mod, only : mxGetN, mxGetString
use, non_intrinsic :: fmxapi_mod, only : mxCreateDoubleScalar
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use, non_intrinsic :: fmxapi_mod, only : fmxVerifyClassShape

implicit none

! mexFunction arguments:
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Local variables
mwSize :: cols  ! Size of the input
! The largest length of the input string.
integer, parameter :: maxlen = 50
! The input string, which specifies the data type
character(len=maxlen) :: data_type
! Integer zero used by the Fortran code
integer(IK), parameter :: intZero = 0
! Integer zero used by MEX for sizes
mwSize, parameter :: msZero = 0
! Integer zero used by MEX for indices
mwIndex, parameter :: miZero = 0
! The huge value that will be returned
real(DP) :: hugeValue
! The success indicator of mxGetString
logical :: success
! The validity indicator of the input
logical :: valid_input
! Name of the current subroutine.
character(len=*), parameter :: srname = 'GETHUGE'

! Check inputs
call fmxVerifyNArgin(nargin, 1)
call fmxVerifyNArgout(nargout, 1)
call fmxVerifyClassShape(pinput(1), 'char', 'row')

! Get the input string. Perform validation even if not debugging.
cols = int(mxGetN(pinput(1)), kind(cols))
call validate(cols <= maxlen, 'COLS <= MAXLEN', srname)
success = (mxGetString(pinput(1), data_type, cols) == 0)
call validate(success, 'Get the input successfully', srname)

! Define hugeValue.
! Note that the REAL values passed to MATLAB via MEX can only be doubles. Therefore, we may need to
! cap the huge values by taking min to ensure that they will not overflow when cast to doubles.
hugeValue = -real(1.0, kind(hugeValue))
valid_input = .true.
select case (lower(data_type))
case ('float')  ! No overflow in this case
    hugeValue = real(huge(0.0), DP)
case ('single')  ! No overflow in this case
    hugeValue = real(huge(0.0_SP), DP)
case ('double')  ! No overflow in this case
    hugeValue = huge(0.0_DP)
case ('real')  ! HUGE(0.0_RP) may be bigger than HUGE(0.0_DP) in this case
    if (RP == DP) then
        hugeValue = huge(0.0_DP)
    else
        hugeValue = 10.0_DP**(min(real(log10(huge(0.0_RP)), DP), log10(huge(0.0_DP))) - 1.0_DP)
    end if
case ('integer', 'int')  ! No overflow in this case
    hugeValue = real(huge(intZero), DP)
case ('mwsi')  ! No overflow in this case
    hugeValue = min(real(huge(msZero), DP), real(huge(miZero), DP))
case ('fun', 'function')  ! HUGEFUN < huge(0.0_DP) according to the definition in CONSTS_MOD.
    hugeValue = real(HUGEFUN, DP)
case ('con', 'constraint')  ! HUGECON < huge(0.0_DP) according to the definition in CONSTS_MOD.
    hugeValue = real(HUGECON, DP)
case default
    valid_input = .false.
end select

! Check that the input is valid (even if not debugging).
call validate(valid_input, 'Input is valid', srname)

! Write output.
! Do NOT use fmxWriteMPtr; when DP /= RP, there is no version of fmxWriteMPtr available in FMXAPI_MOD.
poutput(1) = mxCreateDoubleScalar(hugeValue)

! Postconditions
if (DEBUGGING) then
    call assert(hugeValue > 0, 'HUGEVALUE > 0', srname)
end if
end subroutine mexFunction
