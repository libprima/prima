#include "fintrf.h"

module fmxapi_mod
!--------------------------------------------------------------------------------------------------!
! FMXAPI_MOD is a module that does the following.
! 1. Define some constants to be used in MEX gateways.
! 2. Declare the interfaces of some MEX API subroutine/functions provided by MathWorks.
! 3. Define some user-friendly subroutines for interfacing Fortran with MATLAB. Note that we suppose
! that the REAL type used in the Fortran code is REAL(RP), and the INTEGER type is INTEGER(IK).
!
! Authors:
!   Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!   and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!   Department of Applied Mathematics,
!   The Hong Kong Polytechnic University
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015)
!
! Started in March 2020
!
! Last Modified: Saturday, February 12, 2022 PM02:36:44
!--------------------------------------------------------------------------------------------------!

! N.B.:
! 1. MathWorks may change its APIs in the future!!!
! 2. Make sure that everything is identical to the description in the official documentation of
! MathWorks. Otherwise, failure or unexpected behavior may occur!!!
! 3. Be careful with the "kind" and storage size for integer-type (integer, mwSize, mwIndex)
! variables/functions. Some of them may be 32bit, while the others may be 64bit, depending on the
! machine, the version of MATLAB, and the compilation option of mex. Do NOT assume any two of them
! to be the same. If ever a Segmentation Fault occurs, check these variables first.
! 4. Note that MEX generally use double precision for real values. It is not necessarily the case
! in the Fortran code. Therefore, explicit type conversion is necessary whenever real values are
! exchanged between Fortran and MATLAB. Type mismatch will lead to errors like Segmentation Fault.
! 5. Be careful with the line width limit. After preprocessing (macro expansion), some lines may
! become too long and hence get truncated. For the same reason, do NOT have any continued line
! involving macros, because the & may not appear at the correct position after macro expansion. This
! is why, for example, we define EID and MSG in the subroutines to avoid line continuation involving
! mexErrMsgIdAndTxt.

use, non_intrinsic :: consts_mod, only : DP, RP
implicit none
private

public :: notComplex, mwOne, intOne, intTwo, cvsnTol

! MEX API subroutines
public :: mexErrMsgIdAndTxt
public :: mxCopyPtrToReal8
public :: mxCopyReal8ToPtr
public :: mxDestroyArray

! MEX API functions
public :: mexCallMATLAB
public :: mxCreateDoubleMatrix
public :: mxCreateDoubleScalar
public :: mxGetDoubles
public :: mxGetM
public :: mxGetN
public :: mxGetPr
public :: mxGetString
public :: mxIsChar
public :: mxIsClass
public :: mxIsDouble

! MEX API subroutines/functions defined in this module
public :: fmxGetDble
public :: fmxVerifyNArgin
public :: fmxVerifyNArgout
public :: fmxVerifyClassShape
public :: fmxReadMPtr
public :: fmxWriteMPtr
public :: fmxCallMATLAB
public :: fmxIsDoubleScalar
public :: fmxIsDoubleVector

! Comments on INT32_MEX:
! 1. INT32_MEX is indeed INT32, i.e., the kind of INTEGER*4. It is needed when using some MEX API
! subroutines of MathWorks, e.g., mexCallMATLAB. We do not define it elsewhere (e.g., in consts_mod)
! since it is needed only here. We name it INT32_MEX instead of INT32 so that it is easily locatable.
! 2. For gfortran, SELECTED_REAL_KIND(K) returns INT32 with K = 5--9.
! 3. In F2008, INT32 can be obtained by:
!!use, intrinsic :: iso_fortran_env, only : INT32_MEX => INT32
! 4. Do not write `use consts_mod, only : INT32_MEX => INT32`, as INT32 may not be defined.
integer, parameter :: INT32_MEX = selected_int_kind(7)

! notComplex is used in mxCreateDoubleMatrix
integer(INT32_MEX), parameter :: notComplex = 0
! intOne and intTwo may be used when calling mexCallMATLAB
integer(INT32_MEX), parameter :: intOne = 1, intTwo = 2
! mwOne may be used in mxCreateDoubleMatrix and mxCopyPtrToReal8
mwSize, parameter :: mwOne = 1 ! Integer 1 with type mwSize
! cvsnTol is the tolerance of difference due to conversion between REAL(RP) and REAL(DP)
real(DP), parameter :: cvsnTol = 1.0E1_DP * max(epsilon(0.0_DP), real(epsilon(0.0_RP), DP))

interface fmxReadMPtr
    ! fmxReadMPtr reads the numeric data associated with an mwPointer. It verifies the class and
    ! shape of the data and converts it to REAL(RP) or INTEGER(IK).
    module procedure read_rscalar, read_rvector, read_rmatrix
    module procedure read_iscalar
    module procedure read_lscalar
end interface fmxReadMPtr

interface fmxWriteMPtr
    ! fmxWriteMPtr associates numeric data with an mwPointer. It converts the data to REAL(DP), and
    ! allocates space if the data is a vector or matrix. Therefore, it is necessary to call
    ! mxDestroyArray when the usage of the vector/matrix terminates.
    module procedure write_rscalar, write_rvector, write_rmatrix
    module procedure write_iscalar
end interface fmxWriteMPtr


interface
! Here we declare the interfaces of MEX API subroutines/functions provided by MathWorks. MathWorks
! may change the interfaces in the future!!! Make sure that the interfaces are identical to those
! described in the official documentation of MathWorks!!!
! In particular, pay attention to the following.
! 1. What is the type of an array? Is it automatic (like y(n)), assumed shape (like y(:)), or
! assumed size (like y(*))?
! 2. What is the kind of an integer argument? Is it INT32, INT64, or default INTEGER?
! 3. What is the kind of a real argument? Is it REAL32, REAL64, or default REAL?
! 4. The return values of IsClass, IsChar, and IsDouble, etc., are INTEGER*4 (here we use INT32_MEX
! to represent it). MathWorks may change them in the future to, e.g., logical or default INTEGER.
! 5. Very weirdly, according to MATLAB 2020a documentation, the signature of mexFunction (entry
! point to Fortran MEX function) is
!
!    !---------------------------------------------!
!    subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!    integer nlhs, nrhs
!    mwPointer plhs(*), prhs(*)
!    !---------------------------------------------!
!
! while that of mexCallMATLAB is
!
!    !------------------------------------------------------------!
!    integer*4 mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)
!    integer*4 nlhs, nrhs
!    mwPointer plhs(*), prhs(*)
!    character*(*) functionName
!    !------------------------------------------------------------!
!
! Note that the NLHS/NRHS in the two signatures DO NOT have the same type (INTEGER v.s. INTEGER*4).
! This does not cause any problem, but very bizarre! MathWorks may well modify this later ---
! for example, change all the INTEGER*4 to INTEGER. In that case, we would have to replace all the
! INTEGER(INT32_MEX) by INTEGER.

! MEX subroutines
    subroutine mexErrMsgIdAndTxt(eid, emsg)
    implicit none
    character*(*), intent(in) :: eid, emsg
    end subroutine mexErrMsgIdAndTxt

    subroutine mexWarnMsgIdAndTxt(wid, wmsg)
    implicit none
    character*(*), intent(in) :: wid, wmsg
    end subroutine mexWarnMsgIdAndTxt

    subroutine mxCopyPtrToReal8(px, y, n)
    use, non_intrinsic :: consts_mod, only : DP
    implicit none
    mwPointer, intent(in) :: px
    mwSize, intent(in) :: n
    real(DP), intent(out) :: y(n)
    end subroutine mxCopyPtrToReal8

    subroutine mxCopyReal8ToPtr(y, px, n)
    use, non_intrinsic :: consts_mod, only : DP
    implicit none
    mwPointer, intent(in) :: px
    mwSize, intent(in) :: n
    real(DP), intent(in) :: y(n)
    end subroutine mxCopyReal8ToPtr

    subroutine mxDestroyArray(pm)
    implicit none
    mwPointer, intent(in) :: pm
    end subroutine mxDestroyArray


! MEX functions
    function mexCallMATLAB(nout, pout, nin, pin, f)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mexCallMATLAB
    integer(INT32_MEX), intent(in) :: nout, nin
    ! N.B.: Segmentation Fault will occur if we write POUT(:) or PIN(:)
    mwPointer, intent(in) :: pin(*)
    mwPointer, intent(out) :: pout(*)
    character*(*), intent(in) :: f
    end function mexCallMATLAB

    function mxCreateDoubleMatrix(m, n, ComplexFlag)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    mwPointer :: mxCreateDoubleMatrix
    mwSize, intent(in) :: m, n
    integer(INT32_MEX), intent(in) :: ComplexFlag
    end function mxCreateDoubleMatrix

    function mxCreateDoubleScalar(x)
    use, non_intrinsic :: consts_mod, only : DP
    implicit none
    mwPointer :: mxCreateDoubleScalar
    real(DP), intent(in) :: x
    end function mxCreateDoubleScalar

    function mxGetDoubles(pa)
    implicit none
    mwPointer :: mxGetDoubles
    mwPointer, intent(in) :: pa
    end function mxGetDoubles

    function mxGetM(pm)
    implicit none
    ! The type of mxGetM/N is mwPointer by MATLAB R2020a documentation. Shouldn't it be mwSize?
    mwPointer :: mxGetM
    mwPointer, intent(in) :: pm
    end function mxGetM

    function mxGetN(pm)
    implicit none
    mwPointer :: mxGetN
    mwPointer, intent(in) :: pm
    end function mxGetN

    function mxGetPr(pa)
    implicit none
    mwPointer :: mxGetPr
    mwPointer, intent(in) :: pa
    end function mxGetPr

    function mxGetString(pm, str, strlen)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mxGetString
    mwPointer, intent(in) :: pm
    character*(*), intent(out) :: str
    mwSize, intent(in) :: strlen
    end function mxGetString

    function mxIsClass(pm, classname)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mxIsClass
    mwPointer, intent(in) :: pm
    character*(*), intent(in) :: classname
    end function mxIsClass

    function mxIsChar(pm)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mxIsChar
    mwPointer, intent(in) :: pm
    end function mxIsChar

    function mxIsDouble(pm)
    import :: INT32_MEX  ! Without IMPORT, INT32_MEX will not be available in this interface.
    implicit none
    integer(INT32_MEX) :: mxIsDouble
    mwPointer, intent(in) :: pm
    end function mxIsDouble

end interface


contains


! Here we define some API subroutines/functions for interfacing Fortran code with MATLAB.

function fmxGetDble(pa)
implicit none
mwPointer :: fmxGetDble
mwPointer, intent(in) :: pa
! fmxGetDble gets the pointer pointing to a real array. It is nothing but a wrapper of the
! mxGetDoubles or mxGetPr subroutine defined in fintrf.h. We use mxGetDoubles instead of mxGetPr
! if possible, the former being available since MATLAB R2018b. The macros below should be put after
! fintrf.h is included, because mxGetDoubles is defined in it.
#if defined mxGetDoubles
fmxGetDble = mxGetDoubles(pa)
#else
fmxGetDble = mxGetPr(pa)
#endif
end function


subroutine fmxVerifyNArgin(nin, expected_nin)
! fmxVerifyNArgin verifies that nin = expected_nin.
use, non_intrinsic :: consts_mod, only : MSGLEN
implicit none
integer, intent(in) :: nin  ! NARGIN is of type INTEGER
integer, intent(in) :: expected_nin

character(len=MSGLEN) :: eid, msg

if (nin /= expected_nin) then
    eid = 'FMXAPI:nInput'
    msg = 'fmxVerifyNArgin: Incorrect number of input arguments.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if
end subroutine fmxVerifyNArgin

subroutine fmxVerifyNArgout(nout, expected_nout)
! fmxVerifyNArgout verifies that nout <= expected_nout.
use, non_intrinsic :: consts_mod, only : MSGLEN
implicit none
integer, intent(in) :: nout  ! NARGOUT is of type INTEGER
integer, intent(in) :: expected_nout

character(len=MSGLEN) :: eid, msg

if (nout > expected_nout) then
    eid = 'FMXAPI:nOutput'
    msg = 'fmxVerifyNArgout: Too many output arguments.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if
end subroutine fmxVerifyNArgout


subroutine fmxVerifyClassShape(px, class_name, shape_type)
! fmxVerifyClassShape verifies the class and shape of the data associated with mwPointer px.
use, non_intrinsic :: consts_mod, only : MSGLEN
use, non_intrinsic :: string_mod, only : lower
implicit none
mwPointer, intent(in) :: px
character(len=*), intent(in) :: class_name
character(len=*), intent(in) :: shape_type

mwSize :: m, n
character(len=MSGLEN) :: eid, msg

if (px == 0) then
    eid = 'FMXAPI:NULLPointer'
    msg = 'fmxVerifyClassShape: NULL pointer received.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if

if (mxIsClass(px, class_name) /= 1) then
    eid = 'FMXAPI:WrongInput'
    msg = 'fmxVerifyClassShape: A variable of invalid class received when an argument of class "'//class_name//'" is expected.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if

! Check fmxGetDble(px) if px is associated with a double
if (lower(class_name) == 'double' .and. (lower(shape_type) == 'scalar' .or. lower(shape_type) == 'rank0')) then
    if (fmxGetDble(px) == 0) then  ! This can happen if px is associated with an empty array
        eid = 'FMXAPI:NULLPointer'
        msg = 'fmxVerifyClassShape: NULL pointer returned by fmxGetDble.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end if
end if

m = mxGetM(px)
n = mxGetN(px)

select case (lower(shape_type))
case ('rank0', 'scalar')
    if (m /= 1 .or. n /= 1) then
        eid = 'FMXAPI:WrongInput'
        msg = 'fmxVerifyClassShape: A variable of invalid shape received when an array of rank 0 (scalar) is expected.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end if
case ('rank1', 'vector')
    ! We accept a 0-by-0 array as a vector.
    if ((m /= 1 .or. n < 0) .and. (m < 0 .or. n /= 1) .and. (m /= 0 .or. n /= 0)) then
        eid = 'FMXAPI:WrongInput'
        msg = 'fmxVerifyClassShape: A variable of invalid shape received when an array of rank 1 (vector) is expected.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end if
case ('rank2', 'matrix')
    ! We accept a 0-by-0 array as a matrix.
    if (m < 0 .or. n < 0) then
        eid = 'FMXAPI:WrongInput'
        msg = 'fmxVerifyClassShape: A variable of invalid shape received when an array of rank 2 (matrix) is expected.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end if
case ('column', 'col')
    ! We accept a 0-by-0 array as a column.
    if ((m < 0 .or. n /= 1) .and. (m /= 0 .or. n /= 0)) then
        eid = 'FMXAPI:WrongInput'
        msg = 'fmxVerifyClassShape: A variable of invalid shape received when a column vector is expected.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end if
case ('row')
    ! We accept a 0-by-0 array as a row.
    if ((m /= 1 .or. n < 0) .and. (m /= 0 .or. n /= 0)) then
        eid = 'FMXAPI:WrongInput'
        msg = 'fmxVerifyClassShape: A variable of invalid shape received when a row vector is expected.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end if
case default
    eid = 'FMXAPI:WrongShapeType'
    msg = 'fmxVerifyClassShape: An invalid shape type "'//shape_type//'" received.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end select
end subroutine fmxVerifyClassShape


subroutine read_rscalar(px, x)
! READ_RSCALAR reads the double scalar associated with an mwPointer PX and saves the data in X,
! which is a REAL(RP) scalar.
use, non_intrinsic :: consts_mod, only : RP, DP, ONE, MSGLEN
implicit none

! Input
mwPointer, intent(in) :: px

! Output
real(RP), intent(out) :: x

! Local variables
real(DP) :: x_dp(1)
character(len=MSGLEN) :: wid, msg

! Check input type and size
call fmxVerifyClassShape(px, 'double', 'scalar')

! Read the input
call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

! Convert the input to the type expected by the Fortran code
x = real(x_dp(1), kind(x))
! Check whether the type conversion is proper
if (kind(x) /= kind(x_dp)) then
    if (abs(x - x_dp(1)) > cvsnTol * max(abs(x), ONE)) then
        wid = 'FMXAPI:LargeConversionError'
        msg = 'READ_RSCALAR: Large error occurs when converting REAL(DP) to REAL(RP) (maybe due to overflow).'
        call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
    end if
end if
end subroutine read_rscalar


subroutine read_rvector(px, x)
! READ_RVECTOR reads the double vector associated with an mwPointer PX and saves the data in X,
! which is a REAL(RP) allocatable vector and should have size mxGetM(PX)*mxGetN(PX) at return.
use, non_intrinsic :: consts_mod, only : RP, DP, IK, ONE, MSGLEN
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Input
mwPointer, intent(in) :: px

! Output
real(RP), allocatable, intent(out) :: x(:)

! Local variables
real(DP), allocatable :: x_dp(:)
integer(IK) :: n
mwSize :: n_mw
character(len=MSGLEN) :: wid, msg

! Check input type and size
call fmxVerifyClassShape(px, 'double', 'vector')

! Get size
n_mw = int(mxGetM(px) * mxGetN(px), kind(n_mw))
n = int(n_mw, kind(n))

! Copy input to X_DP
call safealloc(x_dp, n) ! NOT removable
call mxCopyPtrToReal8(fmxGetDble(px), x_dp, n_mw)

! Convert X_DP to the type expected by the Fortran code
call safealloc(x, n) ! Removable in F2003
x = real(x_dp, kind(x))
! Check whether the type conversion is proper
if (kind(x) /= kind(x_dp)) then
    if (maxval(abs(x - x_dp)) > cvsnTol * max(maxval(abs(x)), ONE)) then
        wid = 'FMXAPI:LargeConversionError'
        msg = 'READ_RVECTOR: Large error occurs when converting REAL(DP) to REAL(RP) (maybe due to overflow).'
        call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
    end if
end if

! Deallocate X_DP. Indeed, automatic deallocation would take place.
deallocate (x_dp)
end subroutine read_rvector


subroutine read_rmatrix(px, x)
! READ_RMATRIX reads the double matrix associated with an mwPointer PX and saves the data in X,
! which is a REAL(RP) allocatable matrix and should have size [mxGetM(PX), mxGetN(PX)] at return.
use, non_intrinsic :: consts_mod, only : RP, DP, IK, ONE, MSGLEN
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Input
mwPointer, intent(in) :: px

! Output
real(RP), allocatable, intent(out) :: x(:, :)

! Local variables
real(DP), allocatable :: x_dp(:, :)
integer(IK) :: m, n
mwSize :: xsize
character(len=MSGLEN) :: wid, msg

! Check input type and size
call fmxVerifyClassShape(px, 'double', 'matrix')

! Get size
m = int(mxGetM(px), kind(m))
n = int(mxGetN(px), kind(n))
xsize = int(m * n, kind(xsize))

! Copy input to X_DP
call safealloc(x_dp, m, n) ! NOT removable
call mxCopyPtrToReal8(fmxGetDble(px), x_dp, xsize)

! Convert X_DP to the type expected by the Fortran code
call safealloc(x, m, n) ! Removable in F2003
x = real(x_dp, kind(x))
! Check whether the type conversion is proper
if (kind(x) /= kind(x_dp)) then
    if (maxval(abs(x - x_dp)) > cvsnTol * max(maxval(abs(x)), ONE)) then
        wid = 'FMXAPI:LargeConversionError'
        msg = 'READ_RMATRIX: Large error occurs when converting REAL(DP) to REAL(RP) (maybe due to overflow).'
        call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
    end if
end if

! Deallocate X_DP. Indeed, automatic deallocation would take place.
deallocate (x_dp)
end subroutine read_rmatrix


subroutine read_iscalar(px, x)
! READ_ISCALAR reads a MEX input X that is a double scalar with an integer value. Such a value will
! be passed to the Fortran code as an integer. In MEX, data is passed by pointers, but there are
! only very limited functions that can read an integer value from a pointer or write an integer
! value to a pointer (mxCopyPtrToInteger1, mxCopyInteger1ToPtr, mxCopyPtrToInteger2,
! mxCopyInteger2ToPtr, mxCopyPtrToInteger4, mxCopyInteger4ToPtr; no function for INTEGER*8). This
! makes it impossible to pass integer data properly unless we know the kind of the integer.
! Therefore, in general, it is recommended to pass integers as double variables and then cast them
! back to integers before using them in the Fortran code. Indeed, in MATLAB, even if we define
! X = 1000, the class of X is double! To get an integer X, we would have to define convert it to an
! integer explicitly!
use, non_intrinsic :: consts_mod, only : DP, IK, MSGLEN
implicit none

! Input
mwPointer, intent(in) :: px

! Output
integer(IK), intent(out) :: x

! Local variables
real(DP) :: x_dp(1)
character(len=MSGLEN) :: wid, msg

! Check input type and size
call fmxVerifyClassShape(px, 'double', 'scalar')

! Read the input
call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

! Convert the input to the type expected by the Fortran code
x = int(x_dp(1), kind(x))

! Check whether the type conversion is proper
if (abs(x - x_dp(1)) > epsilon(x_dp) * max(abs(x), 1_IK)) then
    wid = 'FMXAPI:LargeConversionError'
    msg = 'READ_ISCALAR: Large error occurs when converting REAL(DP) to INTEGER(IK) ' &
        & //'(maybe due to overflow, or the input is not an integer).'
    call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
end if
end subroutine read_iscalar


subroutine read_lscalar(px, x)
! READ_LSCALAR reads a MEX input X that is a double scalar with a boolean value. Such a value will
! be passed to the Fortran code as a logical. In MEX, data is passed by pointers, but there is no
! functions that can read a boolean value from a pointer. Therefore, in general, it is recommended
! to pass logicals as double variables and then cast them back to logicals before using them in the
! Fortran code.
use, non_intrinsic :: consts_mod, only : DP, MSGLEN
implicit none

! Input
mwPointer, intent(in) :: px

! Output
logical, intent(out) :: x

! Local variables
character(len=MSGLEN) :: wid, msg
integer :: x_int
real(DP) :: x_dp(1)

! Check input type and size
call fmxVerifyClassShape(px, 'double', 'scalar')

! Read the input
call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

! Convert the input to the type expected by the Fortran code
x_int = int(x_dp(1))

! Check whether the type conversion is proper
if (abs(x_int - x_dp(1)) > epsilon(x_dp) * max(abs(x_int), 1)) then
    wid = 'FMXAPI:LargeConversionError'
    msg = 'READ_LSCALAR: Large error occurs when converting REAL(DP) to INTEGER ' &
        & //'(maybe due to overflow, or the input is not an integer).'
    call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
end if
if (x_int /= 0 .and. x_int /= 1) then
    wid = 'FMXAPI:InputNotBoolean'
    msg = 'READ_LSCALAR: The input should be boolean, either 0 or 1.'
    call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
end if

x = (x_int /= 0)

end subroutine read_lscalar


subroutine write_rscalar(x, px)
! WRITE_RSCALAR associates a REAL(RP) scalar X with an mwPointer PX, after which X can be passed to
! MATLAB either as an output of mexFunction or an input of mexCallMATLAB.
use, non_intrinsic :: consts_mod, only : RP, DP, ONE, MSGLEN
implicit none

! Input
real(RP), intent(in) :: x

! Output
mwPointer, intent(out) :: px

! Local variables
real(DP) :: x_dp
character(len=MSGLEN) :: wid, msg

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
x_dp = real(x, kind(x_dp))
! Check whether the type conversion is proper
if (kind(x_dp) /= kind(x)) then
    if (abs(x - x_dp) > cvsnTol * max(abs(x), ONE)) then
        wid = 'FMXAPI:LargeConversionError'
        msg = 'WRITE_RSCALAR: Large error occurs when converting REAL(RP) to REAL(DP) (maybe due to overflow).'
        call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
    end if
end if

px = mxCreateDoubleScalar(x_dp)

end subroutine write_rscalar


subroutine write_rvector(x, px, shape_type)
! WRITE_RVECTOR associates a REAL(RP) vector X with an mwPointer PX, after which X can be passed to
! MATLAB either as an output of mexFunction or an input of mexCallMATLAB. If ROWCOL = 'row', then
! the vector is passed as a row vector, otherwise, it will be a column vector.
use, non_intrinsic :: consts_mod, only : DP, RP, IK, ONE, MSGLEN
use, non_intrinsic :: string_mod, only : lower
implicit none

! Input
real(RP), intent(in) :: x(:)
character(len=*), intent(in), optional :: shape_type

! Output
mwPointer, intent(out) :: px

! Local variables
real(DP) :: x_dp(size(x))
integer(IK) :: n
mwSize :: n_mw
logical :: row
character(len=MSGLEN) :: wid, msg

! Get size of X
n_mw = int(size(x), kind(n_mw))
n = int(n_mw, kind(n))

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
x_dp = real(x, kind(x_dp))
! Check whether the type conversion is proper
if (kind(x) /= kind(x_dp)) then
    if (maxval(abs(x - x_dp)) > cvsnTol * max(maxval(abs(x)), ONE)) then
        wid = 'FMXAPI:LargeConversionError'
        msg = 'WRITE_RVECTOR: Large error occurs when converting REAL(RP) to REAL(DP) (maybe due to overflow).'
        call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
    end if
end if

row = .false.
if (present(shape_type)) then
    row = (lower(shape_type) == 'row')
end if
! Create a MATLAB matrix using the data in X_DP
if (row) then
    px = mxCreateDoubleMatrix(mwOne, n_mw, notComplex)
else
    px = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
end if
call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), n_mw)

end subroutine write_rvector


subroutine write_rmatrix(x, px)
! WRITE_RMATRIX associates a REAL(RP) matrix X with an mwPointer PX, after which X can be passed to
! MATLAB either as an output of mexFunction or an input of mexCallMATLAB.
use, non_intrinsic :: consts_mod, only : DP, RP, IK, ONE, MSGLEN
implicit none

! Input
real(RP), intent(in) :: x(:, :)

! Output
mwPointer, intent(out) :: px

! Local variables
real(DP) :: x_dp(size(x, 1), size(x, 2))
integer(IK) :: m, n
mwSize :: m_mw, n_mw
character(len=MSGLEN) :: wid, msg

! Get size of X
m = int(size(x, 1), kind(m))
n = int(size(x, 2), kind(n))
m_mw = int(m, kind(m_mw))
n_mw = int(n, kind(n_mw))

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
x_dp = real(x, kind(x_dp))
! Check whether the type conversion is proper
if (kind(x) /= kind(x_dp)) then
    if (maxval(abs(x - x_dp)) > cvsnTol * max(maxval(abs(x)), ONE)) then
        wid = 'FMXAPI:LargeConversionError'
        msg = 'WRITE_RMATRIX: Large error occurs when converting REAL(RP) to REAL(DP) (maybe due to overflow).'
        call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
    end if
end if

! Create a MATLAB matrix using the data in X_DP
px = mxCreateDoubleMatrix(m_mw, n_mw, notComplex)
call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), m_mw * n_mw)

end subroutine write_rmatrix


subroutine write_iscalar(x, px)
! WRITE_RSCALAR associates an INTEGER(IK) scalar X with an mwPointer PX, after which X can be passed
! to MATLAB either as an output of mexFunction or an input of mexCallMATLAB.
use, non_intrinsic :: consts_mod, only : DP, IK, MSGLEN
implicit none

! Input
integer(IK), intent(in) :: x

! Output
mwPointer, intent(out) :: px

! Local variables
real(DP) :: x_dp
character(len=MSGLEN) :: wid, msg

! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
x_dp = real(x, kind(x_dp))
if (abs(x - x_dp) > epsilon(x_dp) * max(abs(x), 1_IK)) then
    wid = 'FMXAPI:LargeConversionError'
    msg = 'WRITE_ISCALAR: Large error occurs when converting INTEGER(IK) to REAL(DP) (maybe due to overflow).'
    call mexWarnMsgIdAndTxt(trim(wid), trim(msg))
end if

px = mxCreateDoubleScalar(x_dp)

end subroutine write_iscalar


subroutine fmxCallMATLAB(fun_ptr, pin, pout)
! fmxCallMATLAB executes matlab command
! output = feval(fun, input),
! where fun_ptr is an mwPointer pointing to the function handle of fun, while pin/pout are mwPointer
! arrays associated with the inputs/outputs
use, non_intrinsic :: consts_mod, only : MSGLEN
implicit none

mwPointer, intent(in) :: fun_ptr
mwPointer, intent(in) :: pin(:)
mwPointer, intent(out) :: pout(:)

mwPointer :: aug_pin(size(pin) + 1)
integer(INT32_MEX) :: nin
integer(INT32_MEX) :: nout
character(5), parameter :: FEVAL = 'feval'

integer(INT32_MEX) :: r
character(len=MSGLEN) :: eid, msg

! Augment the input to include FUN_PTR
aug_pin = [fun_ptr, pin]

! Get number of inputs and number of outputs
nin = int(size(aug_pin), kind(nin))
nout = int(size(pout), kind(nout))

! If mexCallMATLAB returns 0, the execution is successful.
r = mexCallMATLAB(nout, pout, nin, aug_pin, FEVAL)
if (r /= 0) then
    eid = 'FMXAPI:UnsuccessfulCall'
    msg = 'fmxCallMATLAB: MEX fails to call a MATLAB function.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if

if (any(pout == 0)) then
    eid = 'FMXAPI:NULLPointer'
    msg = 'fmxCallMATLAB: NULL pointer returned when MEX calls a MATLAB function.'
    call mexErrMsgIdAndTxt(trim(eid), trim(msg))
end if

end subroutine fmxCallMATLAB


function fmxIsDoubleScalar(px) result(y)
implicit none
mwPointer, intent(in) :: px
logical :: y
y = (mxIsDouble(px) == 1 .and. mxGetM(px) * mxGetN(px) == 1)
end function fmxIsDoubleScalar


function fmxIsDoubleVector(px, shape_type) result(y)
use, non_intrinsic :: consts_mod, only : MSGLEN, IK
use, non_intrinsic :: string_mod, only : lower
implicit none
mwPointer, intent(in) :: px
character(len=*), intent(in), optional :: shape_type
logical :: y
character(len=MSGLEN) :: eid, msg
integer(IK) :: m, n

m = int(mxGetM(px), kind(m))
n = int(mxGetN(px), kind(n))
y = ((m == 1 .and. n >= 0) .or. (m >= 0 .and. n == 1))

if (present(shape_type)) then
    select case (lower(shape_type))
    case ('row')
        y = (y .and. (m == 1))
    case ('column', 'col')
        y = (y .and. (n == 1))
    case default
        eid = 'FMXAPI:WrongShapeType'
        msg = 'fmxIsDoubleVector: An invalid shape type "'//shape_type//'" received.'
        call mexErrMsgIdAndTxt(trim(eid), trim(msg))
    end select
end if

! We accept a 0-by-0 array as a vector/row/column because the MATLAB code may pass an empty vector
! in such a way.
y = (y .or. (m == 0 .and. n == 0))
y = (y .and. (mxIsDouble(px) == 1))
end function fmxIsDoubleVector


end module fmxapi_mod

! Remark: What is the distinction between mx and mex prefixes?
! 1. According to matlab.izmiran.ru/help/techdoc/matlab_external/ch03cre5.html, "Routines in the API
! that are prefixed with mx allow you to create, access, manipulate, and destroy mxArrays. Routines
! prefixed with mex perform operations back in the MATLAB environment."
! 2. We use "fmx" prefix for all subroutines defined by us, e.g., fmxReadMPtr. Here "fmx" indicates
! "Fortran" and "MEX".
