! MEXAPI_MOD is a module that does the following. 
! 1. Define some constants to be used in MEX gateways.
! 2. Declare the interfaces of some MEX API subroutine/functions 
!    provided by MathWorks.
! 3. Define some user-friendly subroutines for interfacing Fortran with
!    MATLAB. We suppose that the REAL type used in the Fortran code is
!    REAL(RP), and the INTEGER type is INTEGER(IK).
!
! N.B.: 
! 1. MathWorks may change its APIs in the future!!!
! 2. Make sure that everything is identical to the description in the
!    official documentation of MathWorks. Otherwise, failure or 
!    unexpected behavior may occur!!! 
! 3. Be careful with the "kind" and storage size for integer-type
!    (integer, mwSize, mwIndex) variables/functions. Some of them may be
!    32bit, while they others may be 64bit, depending on the machine,
!    the version of matlab, and the compilation option of mex. Do NOT
!    assume any two of them to be the same. If ever a Segmentation Fault
!    occurs, check these variables first.
! 4. Note that MEX generally use double precesion for real values. It is
!    not necessarily the case in the Fortran code. Therefore, explicit
!    type conversion is necessary whenever real values are exchanged
!    between Fortran and MATLAB. Type mismatch will lead to errors like
!    Segmentation Fault.
! 5. Be careful with the line width limit. After preprocessing (macro
!    expansion), some lines may become too long and hence get
!    truncated. For the same reason, do NOT have any continued line
!    involving macros, because the & may not appear at the correct
!    position after macro expansion. This is why, for example, we define
!    EID and MSSG in the following subroutines to avoid line
!    continuatuion involving mexErrMsgIdAndTxt.

! Coded by Zaikun ZHANG, July 2020


! Do we intend to use quad precision in the Fortran code (1) or not (0)? 
! 1. It is rarely a good idea to use REAL128 as the working precision,
!     which is probably inefficient and unnecessary. 
! 2. Do NOT change __USE_QP__ to 1 unless you are sure that you intend
!    to use quad precision and such a precision is available on the 
!    current platform.
! 3. To use quad precision, you also need to modify ppf.h to set 
!    __QP_AVAILABLE__ to 1 and __REAL_PRECISION__ to 128. Otherwise, it
!    will not work.  
! 4. __USE_QP__ affects only subroutines alloc_rvector_qp and
!    alloc_rmatrix_qp. 
#if defined __USE_QP__
#undef __USE_QP__
#endif
#define __USE_QP__ 0  


#include "fintrf.h"


      module mexapi_mod

      use consts_mod, only : INT32, DP, RP
      implicit none
      private

      public :: notComplex, mwOne, intOne, intTwo, convTol

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
      public :: fmxAllocate
      public :: fmxReadMPtr
      public :: fmxWriteMPtr
      public :: fmxCallMATLAB

      
      ! notComplex is used in mxCreateDoubleMatrix
      integer(INT32), parameter :: notComplex = 0
      ! intOne and intTwo may be used when calling mexCallMATLAB 
      integer(INT32), parameter :: intOne=1, intTwo=2 
      ! mwOne may be used in mxCreateDoubleMatrix and mxCopyPtrToReal8
      mwSize, parameter :: mwOne = 1 ! Integer 1 with type mwSize
      ! convTol is the tolerance of difference due to conversion between 
      ! REAL(RP) and REAL(DP)
      real(DP), parameter :: convTol = 1.0e1_DP*max(epsilon(0.0_DP),    &
     & real(epsilon(0.0_RP), DP))

      interface fmxAllocate
          ! fmxAllocate allocates the space for a vector/matrix
          module procedure alloc_rvector_sp, alloc_rmatrix_sp
          module procedure alloc_rvector_dp, alloc_rmatrix_dp
#if __USE_QP__
          module procedure alloc_rvector_qp, alloc_rmatrix_qp
#endif
      end interface fmxAllocate 

      interface fmxReadMPtr
          ! fmxReadMPtr reads the numeric data associated with an
          ! mwPointer. It verifies the class and shape of the data and
          ! converts it to REAL(RP).
          module procedure read_rscalar, read_rvector, read_rmatrix
          module procedure read_iscalar
      end interface fmxReadMPtr

      interface fmxWriteMPtr
          ! fmxWriteMPtr associates numeric data with an mwPointer.
          ! It converts the data to REAL(DP), and allocates space if the
          ! data is a vector or matrix. Therefore, it is necessary to
          ! call mxDestroyArray when the usage of the vector/matrix
          ! terminates.
          module procedure write_rscalar, write_rmatrix
          module procedure write_rvector, write_rvector_rowcol
          module procedure write_iscalar
      end interface fmxWriteMPtr
      

      interface 
      ! Here we declare the interfaces of MEX API subroutines/functions
      ! provided by MathWorks. MathWorks may change the interfaces in 
      ! the future!!! 
      ! Make sure that the interfaces are identical to those described
      ! in the official documentation of MathWorks!!!
      ! In particular, pay attention to the following.
      ! 1. What is the type of an array? Is it automatic (like y(n)),
      !    asumed shape (like y(:)), or assumed size (like y(*))?
      ! 2. What is the kind of an integer argument? Is it INT32, INT64,
      !    or default INTEGER?
      ! 3. What is the kind of a real argument? Is it REAL32, REAL64, or
      !    default REAL?
      ! 4. The return values of IsClass, IsChar, and IsDouble, etc., 
      !    are INT32. MathWorks may change them in the future to, e.g.,
      !    logical or default INTEGER.

      ! MEX subroutines
          subroutine mexErrMsgIdAndTxt(errorid, errormsg)
          implicit none
          character*(*), intent(in) ::  errorid, errormsg
          end subroutine mexErrMsgIdAndTxt

          subroutine  mxCopyPtrToReal8(px, y, n)
          use consts_mod, only : DP
          implicit none
          mwPointer, intent(in) :: px
          mwSize, intent(in) :: n
          real(DP), intent(out) :: y(n)
          end subroutine mxCopyPtrToReal8

          subroutine mxCopyReal8ToPtr(y, px, n)
          use consts_mod, only : DP
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
          use consts_mod, only : INT32
          implicit none
          integer(INT32) :: mexCallMATLAB
          integer(INT32), intent(in) :: nout, nin
          ! N.B.: 
          ! Segmentation Fault will occur if we write pout(:)
          ! or pin(:) 
          mwPointer, intent(in) :: pin(*)
          mwPointer, intent(out) ::  pout(*) 
          character*(*), intent(in) :: f
          end function mexCallMATLAB

          function mxCreateDoubleMatrix(m, n, ComplexFlag)
          use consts_mod, only : INT32
          implicit none
          mwPointer :: mxCreateDoubleMatrix
          mwSize, intent(in) :: m, n
          integer(INT32), intent(in) :: ComplexFlag
          end function mxCreateDoubleMatrix

          function mxCreateDoubleScalar(x)
          use consts_mod, only : DP
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
          ! The type of mxGetM/N is mwPointer according to MATLAB
          ! R2020a documentation. Shouldn't it be mwSize?
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
          use consts_mod, only : INT32
          implicit none
          integer(INT32) :: mxGetString
          mwPointer, intent(in) :: pm
          character*(*), intent(out) :: str
          mwSize, intent(in) :: strlen
          end function mxGetString

          function mxIsClass(pm, classname)
          use consts_mod, only : INT32
          implicit none
          integer(INT32) :: mxIsClass
          mwPointer, intent(in) :: pm
          character*(*), intent(in) :: classname
          end function mxIsClass

          function mxIsChar(pm)
          use consts_mod, only : INT32
          implicit none
          integer(INT32) :: mxIsChar
          mwPointer, intent(in) :: pm
          end function mxIsChar

          function mxIsDouble(pm)
          use consts_mod, only : INT32
          implicit none
          integer(INT32) :: mxIsDouble
          mwPointer, intent(in) :: pm
          end function mxIsDouble
    
      end interface


      contains
      ! Here we define some API subroutines/functions for interfacing
      ! Fortran code with MATLAB.


      function fmxGetDble(pa)    
      implicit none
      mwPointer :: fmxGetDble
      mwPointer, intent(in) :: pa
      ! fmxGetDble gets the pointer pointing to a real array. It is
      ! nothing but a wrapper of the mxGetDoubles or mxGetPr subroutine
      ! defined in fintrf.h. We use mxGetDoubles instead of mxGetPr if
      ! possible, the former being available since MATLAB R2018b. The 
      ! macros below should be put after fintrf.h is included, because 
      ! mxGetDoubles is defined in it. 
#if defined mxGetDoubles
      fmxGetDble = mxGetDoubles(pa)
#else
      fmxGetDble = mxGetPr(pa)
#endif
      end function 


      subroutine fmxVerifyNArgin(nin, expected_nin)
      ! fmxVerifyNArgin verifies that nin = expected_nin.
      use consts_mod, only : MSSGLEN
      implicit none
      integer, intent(in) :: nin  ! NARGIN is of type INTEGER
      integer, intent(in) :: expected_nin

      character(len = MSSGLEN) :: eid, mssg

      if (nin /= expected_nin) then
          eid = 'MEXAPI:nInput'
          mssg = 'fmxVerifyNArgin: Incorrect number of input arguments.'
          call mexErrMsgIdAndTxt (trim(eid), trim(mssg))
      end if
      end subroutine fmxVerifyNArgin

      subroutine fmxVerifyNArgout(nout, expected_nout)
      ! fmxVerifyNArgout verifies that nout <= expected_nout.
      use consts_mod, only : MSSGLEN
      implicit none
      integer, intent(in) :: nout  ! NARGOUT is of type INTEGER
      integer, intent(in) :: expected_nout

      character(len = MSSGLEN) :: eid, mssg

      if (nout > expected_nout) then
          eid = 'MEXAPI:nOutput'
          mssg = 'fmxVerifyNArgout: Too many output arguments.'
          call mexErrMsgIdAndTxt (trim(eid), trim(mssg))
      end if
      end subroutine fmxVerifyNArgout

      
      subroutine fmxVerifyClassShape(px, class_name, shape_type)
      ! fmxVerifyClassShape verifies the class and shape of the data
      ! associated with mwPointer px. 
      use consts_mod, only : MSSGLEN
      implicit none
      mwPointer, intent(in) :: px
      character(len = *), intent(in) :: class_name
      character(len = *), intent(in) :: shape_type

      mwSize :: m, n
      character(len = MSSGLEN) :: eid, mssg

      if (px == 0) then
          eid = 'MEXAPI:NULLPointer'
          mssg = 'fmxVerifyClassShape: NULL pointer received.'
          call mexErrMsgIdAndTxt (trim(eid), trim(mssg))
      end if

      if (mxIsClass(px, class_name) /= 1) then
          eid = 'MEXAPI:WrongInput'
          mssg = 'fmxVerifyClassShape: A variable of invalid class ' // &
     &     'received when an argument of class "' // class_name //      &
     &     '" is expected.'
          call mexErrMsgIdAndTxt (trim(eid), trim(mssg))
      end if

      ! Check fmxGetDble(px) if px is associated with a double
      if (class_name == 'double') then
          if (fmxGetDble(px) == 0) then
              eid = 'MEXAPI:NULLPointer' 
              mssg = 'fmxVerifyClassShape: NULL pointer returned '//    &
     &         'by fmxGetDble.'
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      end if

      m = mxGetM(px)
      n = mxGetN(px)

      if (shape_type == 'rank0' .or. shape_type == 'scalar') then
          if (m /= 1 .or. n /= 1) then
              eid = 'MEXAPI:WrongInput'
              mssg = 'fmxVerifyClassShape: A variable of invalid ' //   &
     &         'shape received when an array of rank 0 (scalar) is ' // &
     &         'expected.' 
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      else if (shape_type == 'rank1' .or. shape_type == 'vector') then
          if ((m /= 1 .or. n < 1) .and. (m < 1 .or. n /= 1)) then
              eid = 'MEXAPI:WrongInput'
              mssg = 'fmxVerifyClassShape: A variable of invalid ' //   &
     &         'shape received when an array of rank 1 (vector) is ' // &
     &         'expected.' 
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      else if (shape_type == 'rank2' .or. shape_type == 'matrix') then
          if (m < 1 .or. n < 1) then
              eid = 'MEXAPI:WrongInput'
              mssg = 'fmxVerifyClassShape: A variable of invalid ' //   &
     &         'shape received when an array of rank 2 (matrix) is ' // &
     &         'expected.' 
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      else if (shape_type == 'column') then
          if (m < 1 .or. n /= 1) then
              eid = 'MEXAPI:WrongInput'
              mssg = 'fmxVerifyClassShape: A variable of invalid ' //   &
     &         'shape received when a column vector is expected.' 
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      else if (shape_type == 'row') then
          if (m /= 1 .or. n < 1) then
              eid = 'MEXAPI:WrongInput'
              mssg = 'fmxVerifyClassShape: A variable of invalid ' //   &
     &         'shape received when a row vector is expected.' 
              call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
          end if
      else
          eid = 'MEXAPI:WrongShapeType'
          mssg = 'fmxVerifyClassShape: An invalid shape type "' //      &
     &     shape_type // '" received.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if
      end subroutine fmxVerifyClassShape


      subroutine alloc_rvector_sp(x, n)
      ! ALLOC_RVECTOR_SP allocates the space for an allocatable
      ! single-precision vector X, whose size is N after allocation. 
      use consts_mod, only : SP, IK, MSSGLEN
      implicit none
      
      ! Input
      integer(IK), intent(in) :: n

      ! Output
      real(SP), allocatable, intent(out) :: x(:)

      ! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'MEXAPI:AllocateFailed'
          mssg = 'ALLOC_RVECTOR_SP: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      x = 0.0_SP  ! Set X = 0; otherwise, compilers may complain.
      
      end subroutine alloc_rvector_sp

      
      subroutine alloc_rmatrix_sp(x, m, n)
      ! ALLOC_RMATRIX_SP allocates the space for a single-precision
      ! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : SP, IK, MSSGLEN
      implicit none
      
      ! Input
      integer(IK), intent(in) :: m, n

      ! Output
      real(SP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'MEXAPI:AllocateFailed'
          mssg = 'ALLOC_RMATRIX_SP: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if
      
      x = 0.0_SP  ! Set X = 0; otherwise, compilers may complain.
      end subroutine alloc_rmatrix_sp


      subroutine alloc_rvector_dp(x, n)
      ! ALLOC_RVECTOR_DP allocates the space for an allocatable
      ! double-precision vector X, whose size is N after allocation. 
      use consts_mod, only : DP, IK, MSSGLEN
      implicit none
      
      ! Input
      integer(IK), intent(in) :: n

      ! Output
      real(DP), allocatable, intent(out) :: x(:)

      ! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'MEXAPI:AllocateFailed'
          mssg = 'ALLOC_RVECTOR_DP: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      x = 0.0_DP  ! Set = 0; otherieise, compilers may complain.

      end subroutine alloc_rvector_dp


      subroutine alloc_rmatrix_dp(x, m, n)
      ! ALLOC_RMATRIX_DP allocates the space for a double-precision
      ! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : DP, IK, MSSGLEN
      implicit none
      
      ! Input
      integer(IK), intent(in) :: m, n

      ! Output
      real(DP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'MEXAPI:AllocateFailed'
          mssg = 'ALLOC_RMATRIX_DP: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      x = 0.0_DP  ! Set = 0; otherieise, compilers may complain.

      end subroutine alloc_rmatrix_dp


#if __USE_QP__ == 1
      subroutine alloc_rvector_qp(x, n)
      ! ALLOC_RVECTOR_QP allocates the space for an allocatable
      ! double-precision vector X, whose size is N after allocation. 
      use consts_mod, only : QP, IK, MSSGLEN
      implicit none
      
      ! Input
      integer(IK), intent(in) :: n

      ! Output
      real(QP), allocatable, intent(out) :: x(:)

      ! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'MEXAPI:AllocateFailed'
          mssg = 'ALLOC_RVECTOR_QP: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      x = 0.0_QP  ! Set = 0; otherieise, compilers may complain.

      end subroutine alloc_rvector_qp


      subroutine alloc_rmatrix_qp(x, m, n)
      ! ALLOC_RMATRIX_DP allocates the space for a double-precision
      ! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : QP, IK, MSSGLEN
      implicit none
      
      ! Input
      integer(IK), intent(in) :: m, n

      ! Output
      real(QP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variable
      integer :: alloc_status
      character(len = MSSGLEN) :: eid, mssg
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          eid = 'MEXAPI:AllocateFailed'
          mssg = 'ALLOC_RMATRIX_QP: Memory allocation fails.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      x = 0.0_QP  ! Set = 0; otherieise, compilers may complain.

      end subroutine alloc_rmatrix_qp
#endif

      
      subroutine read_rscalar(px, x) 
      ! READ_RSCALAR reads the double scalar associated with an
      ! mwPointer PX and saves the data in X, which is a REAL(RP)
      ! scalar.
      use consts_mod, only : RP, ONE, MSSGLEN
      implicit none

      ! Input
      mwPointer, intent(in) :: px

      ! Output
      real(RP), intent(out) :: x

      ! Intermediate variable
      real(DP) :: x_dp(1)
      character(len = MSSGLEN) :: eid, mssg

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'scalar')

      ! Read the input 
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

      ! Convert the input to the type expected by the Fortran code
      x = real(x_dp(1), kind(x))

      ! Check whether the type conversion is proper
      if (abs(x-x_dp(1)) > convTol*max(abs(x), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'READ_RSCALAR: Large error occurs when converting ' // &
     &     'REAL(DP) to REAL(RP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if
      end subroutine read_rscalar 


      subroutine read_rvector(px, x) 
      ! READ_RVECTOR reads the double vector associated with an 
      ! mwPointer PX and saves the data in X, which is a REAL(RP) 
      ! vector.
      use consts_mod, only : RP, IK, ONE, MSSGLEN
      implicit none

      ! Input
      mwPointer, intent(in) :: px
      
      ! Output
      real(RP), allocatable, intent(out) :: x(:)

      ! Intermediate variables
      real(DP), allocatable  :: x_dp(:)
      integer(IK) :: n
      mwSize :: n_mw 
      character(len = MSSGLEN) :: eid, mssg

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'vector')

      ! Get size
      n_mw = int(mxGetM(px)*mxGetN(px), kind(n_mw))
      n = int(n_mw, kind(n))
      
      ! Allocate memory for X and X_DP
      call fmxAllocate(x, n)
      call fmxAllocate(x_dp, n)

      ! Copy input to X_DP
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, n_mw)

      ! Convert X_DP to the type expected by the Fortran code
      x = real(x_dp, kind(x))

      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'READ_RVECTOR: Large error occurs when converting ' // &
     &     'REAL(DP) to REAL(RP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      ! Deallocate X_DP
      deallocate(x_dp)

      end subroutine read_rvector


      subroutine read_rmatrix(px, x) 
      ! READ_RMATRIX reads the double matrix associated with an
      ! mwPointer PX and saves the data in X, which is a REAL(RP)
      ! matrix.
      use consts_mod, only : RP, IK, ONE, MSSGLEN
      implicit none

      ! Input
      mwPointer, intent(in) :: px
      
      ! Output
      real(RP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variables
      real(DP), allocatable  :: x_dp(:, :)
      integer(IK) :: m, n
      mwSize :: xsize
      character(len = MSSGLEN) :: eid, mssg

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'matrix')

      ! Get size
      m = int(mxGetM(px), kind(m))
      n = int(mxGetN(px), kind(n))
      xsize = int(m*n, kind(xsize))
      

      ! Allocate memory for X and X_DP
      call fmxAllocate(x, m, n)
      call fmxAllocate(x_dp, m, n)

      ! Copy input to X_DP
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, xsize)

      ! Convert X_DP to the type expected by the Fortran code
      x = real(x_dp, kind(x))

      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'READ_RMATRIX: Large error occurs when converting ' // &
     &     'REAL(DP) to REAL(RP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      ! Deallocate X_DP
      deallocate(x_dp)

      end subroutine read_rmatrix


      subroutine read_iscalar(px, x) 
      ! READ_ISCALAR reads a MEX input X that is a double scalar with
      ! an integer value. Such a value will be passed to the Fortran
      ! code as an integer. In MEX, data is passed by pointers, but
      ! there are only very limited functions that can read an integer
      ! value from a pointer or write an integer value to a pointer
      ! (mxCopyPtrToInteger1, mxCopyInteger1ToPtr, mxCopyPtrToInteger2,
      ! mxCopyInteger2ToPtr, mxCopyPtrToInteger4, mxCopyInteger4ToPtr;
      ! no function for INTEGER*8). This makes it impossible to pass
      ! integer data properly unless we know the kind of the integer.
      ! Therefore, in general, it is recommended to pass integers as
      ! double variables and then cast them back to integers befor using
      ! them in the Fortran code. Indeed, in MATLAB, even if we define
      ! X = 1000, the class of X is double! To get an integer X, we
      ! would have to define X = INT32(1000) or X = INT64(1000)!
      use consts_mod, only : IK, MSSGLEN
      implicit none
      
      ! Input
      mwPointer, intent(in) :: px

      ! Output
      integer(IK), intent(out) :: x

      ! Intermediate variable
      real(DP) :: x_dp(1)
      character(len = MSSGLEN) :: eid, mssg

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'scalar')

      ! Read the input 
      call mxCopyPtrToReal8(fmxGetDble(px), x_dp, mwOne)

      ! Convert the input to the type expected by the Fortran code
      x = int(x_dp(1), kind(x))

      ! Check whether the type conversion is proper
      if (abs(x-x_dp(1)) > 0.5_DP) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'READ_ISCALAR: Large error occurs when converting ' // &
     &     'REAL(DP) to INTEGER(IK).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if
      end subroutine read_iscalar 


      subroutine write_rscalar(x, px)
      ! WRITE_RSCALAR associates a REAL(RP) scalar X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : RP, ONE, MSSGLEN
      implicit none

      ! Input
      real(RP), intent(in) :: x

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp
      character(len = MSSGLEN) :: eid, mssg

      x_dp = real(x, kind(x_dp))
      
      if (abs(x-x_dp) > convTol*max(abs(x), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'WRITE_RSCALAR: Large error occurs when converting ' //&
     &     'REAL(RP) to REAL(DP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      px = mxCreateDoubleScalar(x_dp)

      end subroutine write_rscalar


      subroutine write_rvector(x, px)
      ! WRITE_RVECTOR associates a REAL(RP) vector X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : RP, IK, ONE, MSSGLEN
      implicit none

      ! Input
      real(RP), intent(in) :: x(:)

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp(size(x))
      integer(IK) :: n
      mwSize :: n_mw
      character(len = MSSGLEN) :: eid, mssg

      ! Get size of X
      n_mw = int(size(x), kind(n_mw))
      n = int(n_mw, kind(n))

      ! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
      
      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'WRITE_RVECTOR: Large error occurs when converting ' //&
     &     'REAL(RP) to REAL(DP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      ! Create a MATLAB matrix using the data in X_DP
      px = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
      call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), n_mw)

      end subroutine write_rvector


      subroutine write_rvector_rowcol(x, px, rowcol)
      ! WRITE_RVECTOR_ROWCOL associates a REAL(RP) vector X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      ! If ROWCOL = 'row', then the vector is passed as a row vector,
      ! otherwise, it will be a column vector.
      use consts_mod, only : RP, IK, ONE, MSSGLEN
      implicit none

      ! Input
      real(RP), intent(in) :: x(:)
      character(len = *), intent(in) :: rowcol

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp(size(x))
      integer(IK) :: n
      mwSize :: n_mw
      character(len = MSSGLEN) :: eid, mssg

      ! Get size of X
      n_mw = int(size(x), kind(n_mw))
      n = int(n_mw, kind(n))

      ! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
      
      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'WRITE_RVECTOR: Large error occurs when converting ' //&
     &     'REAL(RP) to REAL(DP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      ! Create a MATLAB matrix using the data in X_DP
      if (rowcol == 'row' .or. rowcol == 'ROW' .or. rowcol =='Row') then
          px = mxCreateDoubleMatrix(mwOne, n_mw, notComplex)
      else
          px = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
      end if
      call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), n_mw)

      end subroutine write_rvector_rowcol


      subroutine write_rmatrix(x, px)
      ! WRITE_RMATRIX associates a REAL(RP) matrix X with an
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : RP, IK, ONE, MSSGLEN
      implicit none

      ! Input
      real(RP), intent(in) :: x(:, :)

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp(size(x, 1), size(x, 2))
      integer(IK) :: m, n
      mwSize :: m_mw, n_mw
      character(len = MSSGLEN) :: eid, mssg

      ! Get size of X
      m = int(size(x, 1), kind(m))
      n = int(size(x, 2), kind(n))
      m_mw = int(m, kind(m_mw))
      n_mw = int(n, kind(n_mw))

      ! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
      
      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'WRITE_RMATRIX: Large error occurs when converting ' //&
     &     'REAL(RP) to REAL(DP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      ! Create a MATLAB matrix using the data in X_DP
      px = mxCreateDoubleMatrix(m_mw, n_mw, notComplex)
      call mxCopyReal8ToPtr(x_dp, fmxGetDble(px), m_mw*n_mw)

      end subroutine write_rmatrix


      subroutine write_iscalar(x, px)
      ! WRITE_RSCALAR associates an INTEGER(IK) scalar X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : DP, IK, MSSGLEN
      implicit none

      ! Input
      integer(IK), intent(in) :: x

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp
      character(len = MSSGLEN) :: eid, mssg

      x_dp = real(x, kind(x_dp))
      
      if (abs(x - x_dp) > 0.5_DP) then
          eid = 'MEXAPI:ConversionError'
          mssg = 'WRITE_ISCALAR: Large error occurs when converting ' //&
     &     'INTEGER(IK) to REAL(DP).'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      px = mxCreateDoubleScalar(x_dp)

      end subroutine write_iscalar


      subroutine fmxCallMATLAB(fun_ptr, pin, pout)
      ! fmxCallMATLAB executes matlab command
      ! output = feval(fun, input),
      ! where fun_ptr is an mwPointer pointing to the function handle of
      ! fun, while pin and pout are mwPointer arrays associated 
      ! with the inputs and outputs
      use consts_mod, only : INT32, MSSGLEN
      implicit none
      
      mwPointer, intent(in) :: fun_ptr      
      mwPointer, intent(in) :: pin(:)
      mwPointer, intent(out) :: pout(:)
     
      mwPointer :: aug_pin(size(pin) + 1)
      integer(INT32) :: nin
      integer(INT32) :: nout
      character(5), parameter :: FEVAL = 'feval'

      integer(INT32) :: r
      character(len = MSSGLEN) :: eid, mssg
      
      ! Augment the input to include FUN_PTR
      aug_pin = (/ fun_ptr, pin /)  

      ! Get number of inputs and number of outputs
      nin = int(size(aug_pin), kind(nin))
      nout = int(size(pout), kind(nout))
      
      ! integer*4 mexCallMATLAB(integer*4 nout, mwPointer pout, integer*4 nin, mwPointer pin, character*(*) functionName)
      ! If mexCallMATLAB returns 0, the execution is successful.
      r = mexCallMATLAB(nout, pout, nin, aug_pin, FEVAL)
      if (r /= 0) then
          eid = 'MEXAPI:UnsuccessfulCall'
          mssg = 'fmxCallMATLAB: MEX fails to call a MATLAB function.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      if (any(pout == 0)) then
          eid = 'MEXAPI:NULLPointer'
          mssg = 'fmxCallMATLAB: NULL pointer returned when MEX ' //    &
     &     'calls a MATLAB function.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

      end subroutine fmxCallMATLAB

      end module mexapi_mod

      ! Remark: What is the distinction between mx and mex prefixes?
      ! 1. According to 
      ! matlab.izmiran.ru/help/techdoc/matlab_external/ch03cre5.html
      ! "Routines in the API that are prefixed with mx allow you to 
      ! create, access, manipulate, and destroy mxArrays. Routines 
      ! prefixed with mex perform operations back in the MATLAB 
      ! environment."
      ! 2. We use "fmx" prefix for all subroutines defined by us, e.g.,
      ! fmxReadMPtr. Here "fmx" indicates "Fortran" and "MEX".
