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
!    substitution), some lines may become too long and hence get
!    truncated. 

! Coded by Zaikun ZHANG, July 2020


#include "fintrf.h"

! Use mxGetDoubles instead of mxGetPr if possible, the former being
! available since MATLAB R2018b. The following macros should be put
! after the header "fintrf.h" is included, because mxGetDoubles is
! defined in it. MX_GET_DBLE gets the pointer pointing to a real array.
#ifdef MX_GET_DBLE
#undef MX_GET_DBLE
#endif
#ifdef mxGetDoubles
#define MX_GET_DBLE mxGetDoubles
#else
#define MX_GET_DBLE mxGetPr
#endif


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
          module procedure write_rscalar, write_rvector, write_rmatrix
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
          function mexCallMATLAB(nargout, poutput, nargin,              &
     &     pinput, functionName)
              use consts_mod, only : INT32
              implicit none
              integer(INT32) :: mexCallMATLAB
              integer(INT32), intent(in) :: nargout, nargin
              ! N.B.: 
              ! Segmentation Fault will occur if we write poutput(:)
              ! or pinput(:) 
              mwPointer, intent(in) :: pinput(*)
              mwPointer, intent(out) ::  poutput(*) 
              character*(*), intent(in) :: functionName
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

      
      subroutine fmxVerifyNArgin(nargin, expected_nargin)
      ! fmxVerifyNArgin verifies that nargin = expected_nargin.
      implicit none
      integer, intent(in) :: nargin  ! NARGIN is of type INTEGER
      integer, intent(in) :: expected_nargin

      if (nargin /= expected_nargin) then
         call mexErrMsgIdAndTxt ('MEXAPI_MOD:nInput',                   &
     &    'Incorrect number of input arguments received.')
      end if
      end subroutine fmxVerifyNArgin

      subroutine fmxVerifyNArgout(nargout, expected_nargout)
      ! fmxVerifyNArgout verifies that nargout <= expected_nargout.
      implicit none
      integer, intent(in) :: nargout  ! NARGOUT is of type INTEGER
      integer, intent(in) :: expected_nargout

      if (nargout > expected_nargout) then
         call mexErrMsgIdAndTxt ('MEXAPI_MOD:nOutput',                  &
     &    'Too many output arguments.')
      end if
      end subroutine fmxVerifyNArgout

      
      subroutine fmxVerifyClassShape(px, class_name, shape_type)
      ! fmxVerifyClassShape verifies the class and shape of the data
      ! associated with mwPointer px. 
      mwPointer, intent(in) :: px
      character(len = *), intent(in) :: class_name
      character(len = *), intent(in) :: shape_type

      if (px == 0) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',               &
     &     'NULL pointer received when verifying the class and shape' //&
     &     ' of an argument.')
      end if

      if (mxIsClass(px, class_name) /= 1) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',               &
     &     'A variable of invalid class received when an argument ' //  &
     &     'of class "' // class_name // '" is expected.')
      end if

      ! Check MX_GET_DBLE(px) if px is associated with a double
      if (class_name == 'double') then
          if (MX_GET_DBLE(px) == 0) then
              call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',           &
     &         'MX_GET_DBLE returns a NULL pointer when verifying' //   &
     &         ' the class and shape of an argument.')
          end if
      end if

      if (shape_type == 'rank0' .or. shape_type == 'scalar') then
          if (mxGetM(px) /= 1 .or. mxGetN(px) /= 1) then
              call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',           &
     &         'A variable of invalid shape received when an array '  //&
     &         'of rank 0 (scalar) is expected.') 
          end if
      else if (shape_type == 'rank1' .or. shape_type == 'vector') then
          if ((mxGetM(px) /= 1 .or. mxGetN(px) < 1) .and.               &
     &     (mxGetM(px) < 1 .or. mxGetN(px) /= 1)) then
              call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',           &
     &         'A variable of invalid shape received when an array ' // &
     &         'of rank 1 (vector) is expected.') 
          end if
      else if (shape_type == 'rank2' .or. shape_type == 'matrix') then
          if (mxGetM(px) < 1 .or. mxGetN(px) < 1) then
              call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',           &
     &         'A varible of invalid shape received when an array ' //  &
     &         'of rank 2 (matrix) is expected.') 
          end if
      else if (shape_type == 'column') then
          if (mxGetM(px) < 1 .or. mxGetN(px) /= 1) then
              call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',           &
     &         'A variable of invalid shape received when a column ' // &
     &         'vector is expected.')
          end if
      else if (shape_type == 'row') then
          if (mxGetM(px) /= 1 .or. mxGetN(px) < 1) then
              call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongInput',           &
     &         'A variable of invalid shape received when a row ' //    &
     &         'vector is expected.')
          end if
      else
          call mexErrMsgIdAndTxt('MEXAPI_MOD:WrongShapetype',           &
     &         'An invalid shape type "' // shape_type //               &
     &         '" received when verifying the shape of an argument.')
      end if
      end subroutine fmxVerifyClassShape


      subroutine alloc_rvector_sp(x, n)
      ! ALLOC_RVECTOR_SP allocates the space for an allocatable
      ! single-precision vector X, whose size is N after allocation. 
      use consts_mod, only : SP, IK     
      implicit none
      
      ! Input
      integer(IK), intent(in) :: n

      ! Output
      real(SP), allocatable, intent(out) :: x(:)

      ! Intermediate variable
      integer(IK) :: alloc_status
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:AllocateFailed',           &
     &    'Memory allocation fails.')
      end if

      x = 0.0_SP  ! Set X = 0; otherwise, compilers may complain.
      
      end subroutine alloc_rvector_sp

      
      subroutine alloc_rmatrix_sp(x, m, n)
      ! ALLOC_RMATRIX_SP allocates the space for a single-precision
      ! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : SP, IK     
      implicit none
      
      ! Input
      integer(IK), intent(in) :: m, n

      ! Output
      real(SP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variable
      integer(IK) :: alloc_status
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:AllocateFailed',           &
     &    'Memory allocation fails.')
      end if
      
      x = 0.0_SP  ! Set X = 0; otherwise, compilers may complain.
      end subroutine alloc_rmatrix_sp
      

      subroutine alloc_rvector_dp(x, n)
      ! ALLOC_RVECTOR_DP allocates the space for an allocatable
      ! double-precision vector X, whose size is N after allocation. 
      use consts_mod, only : DP, IK     
      implicit none
      
      ! Input
      integer(IK), intent(in) :: n

      ! Output
      real(DP), allocatable, intent(out) :: x(:)

      ! Intermediate variable
      integer(IK) :: alloc_status
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(n), stat = alloc_status)
      if (alloc_status /= 0) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:AllocateFailed',           &
     &    'Memory allocation fails.')
      end if

      x = 0.0_DP  ! Set = 0; otherieise, compilers may complain.

      end subroutine alloc_rvector_dp


      subroutine alloc_rmatrix_dp(x, m, n)
      ! ALLOC_RMATRIX_DP allocates the space for a double-precision
      ! matrix X, whose size is (M, N) after allocation.
      use consts_mod, only : DP, IK     
      implicit none
      
      ! Input
      integer(IK), intent(in) :: m, n

      ! Output
      real(DP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variable
      integer(IK) :: alloc_status
      
      ! Deallocate X if X is allocated. This is necessary. 
      if (allocated(x)) then
          deallocate (x)
      end if

      ! Allocate memory for X
      allocate (x(m, n), stat = alloc_status)
      if (alloc_status /= 0) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:AllocateFailed',           &
     &    'Memory allocation fails.')
      end if

      x = 0.0_DP  ! Set = 0; otherieise, compilers may complain.

      end subroutine alloc_rmatrix_dp

      
      subroutine read_rscalar(px, x) 
      ! READ_RSCALAR reads the double scalar associated with an
      ! mwPointer PX and saves the data in X, which is a REAL(RP)
      ! scalar.
      use consts_mod, only : RP, ONE
      implicit none

      ! Input
      mwPointer, intent(in) :: px

      ! Output
      real(RP), intent(out) :: x

      ! Intermediate variable
      real(DP) :: x_dp(1)

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'scalar')

      ! Read the input 
      call mxCopyPtrToReal8(MX_GET_DBLE(px), x_dp, mwOne)

      ! Convert the input to the type expected by the Fortran code
      x = real(x_dp(1), kind(x))

      ! Check whether the type conversion is proper
      if (abs(x-x_dp(1)) > convTol*max(abs(x), ONE)) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',          &
     &     'Large error occurs when converting REAL(DP) to REAL(RP).')
      end if
      end subroutine read_rscalar 


      subroutine read_rvector(px, x) 
      ! READ_RVECTOR reads the double vector associated with an 
      ! mwPointer PX and saves the data in X, which is a REAL(RP) 
      ! vector.
      use consts_mod, only : RP, IK, ONE
      implicit none

      ! Input
      mwPointer, intent(in) :: px
      
      ! Output
      real(RP), allocatable, intent(out) :: x(:)

      ! Intermediate variables
      real(DP), allocatable  :: x_dp(:)
      integer(IK) :: n
      mwSize :: n_mw 

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'vector')

      ! Get size
      n_mw = int(mxGetM(px)*mxGetN(px), kind(n_mw))
      n = int(n_mw, kind(n))
      
      ! Allocate memory for X and X_DP
      call fmxAllocate(x, n)
      call fmxAllocate(x_dp, n)

      ! Copy input to X_DP
      call mxCopyPtrToReal8(MX_GET_DBLE(px), x_dp, n_mw)

      ! Convert X_DP to the type expected by the Fortran code
      x = real(x_dp, kind(x))

      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',          &
     &     'Large error occurs when converting REAL(DP) to REAL(RP).')
      end if

      ! Deallocate X_DP
      deallocate(x_dp)

      end subroutine read_rvector


      subroutine read_rmatrix(px, x) 
      ! READ_RMATRIX reads the double matrix associated with an
      ! mwPointer PX and saves the data in X, which is a REAL(RP)
      ! matrix.
      use consts_mod, only : RP, IK, ONE
      implicit none

      ! Input
      mwPointer, intent(in) :: px
      
      ! Output
      real(RP), allocatable, intent(out) :: x(:, :)

      ! Intermediate variables
      real(DP), allocatable  :: x_dp(:, :)
      integer(IK) :: m, n
      mwSize :: xsize

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
      call mxCopyPtrToReal8(MX_GET_DBLE(px), x_dp, xsize)

      ! Convert X_DP to the type expected by the Fortran code
      x = real(x_dp, kind(x))

      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',          &
     &     'Large error occurs when converting REAL(DP) to REAL(RP).')
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
      use consts_mod, only : IK
      implicit none
      
      ! Input
      mwPointer, intent(in) :: px

      ! Output
      integer(IK), intent(out) :: x

      ! Intermediate variable
      real(DP) :: x_dp(1)

      ! Check input type and size
      call fmxVerifyClassShape(px, 'double', 'scalar')

      ! Read the input 
      call mxCopyPtrToReal8(MX_GET_DBLE(px), x_dp, mwOne)

      ! Convert the input to the type expected by the Fortran code
      x = int(x_dp(1), kind(x))

      ! Check whether the type conversion is proper
      if (abs(x-x_dp(1)) > 0.5_DP) then
         call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',           &
     &    'Large error occurs when converting REAL(DP) to INTEGER(IK).')
      end if
      end subroutine read_iscalar 


      subroutine write_rscalar(x, px)
      ! WRITE_RSCALAR associates a REAL(RP) scalar X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : RP, ONE
      implicit none

      ! Input
      real(RP), intent(in) :: x

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp

      x_dp = real(x, kind(x_dp))
      
      if (abs(x-x_dp) > convTol*max(abs(x), ONE)) then
         call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',           &
     &    'Large error occurs when converting REAL(RP) to REAL(DP).')
      end if

      px = mxCreateDoubleScalar(x_dp)

      end subroutine write_rscalar


      subroutine write_rvector(x, px)
      ! WRITE_RSCALAR associates a REAL(RP) vector X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : RP, IK, ONE
      implicit none

      ! Input
      real(RP), intent(in) :: x(:)

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp(size(x))
      integer(IK) :: n
      mwSize :: n_mw

      ! Get size of X
      n_mw = int(size(x), kind(n_mw))
      n = int(n_mw, kind(n))

      ! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
      
      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
         call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',           &
     &    'Large error occurs when converting REAL(RP) to REAL(DP).')
      end if

      ! Create a MATLAB matrix using the data in X_DP
      px = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
      call mxCopyReal8ToPtr(x_dp, MX_GET_DBLE(px), n_mw)

      end subroutine write_rvector


      subroutine write_rmatrix(x, px)
      ! WRITE_RSCALAR associates a REAL(RP) matrix X with an
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : RP, IK, ONE
      implicit none

      ! Input
      real(RP), intent(in) :: x(:, :)

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp(size(x, 1), size(x, 2))
      integer(IK) :: m, n
      mwSize :: m_mw, n_mw

      ! Get size of X
      m = int(size(x, 1), kind(m))
      n = int(size(x, 2), kind(n))
      m_mw = int(m, kind(m_mw))
      n_mw = int(n, kind(n_mw))

      ! Convert X to REAL(DP), which is expected by mxCopyReal8ToPtr
      x_dp = real(x, kind(x_dp))
      
      ! Check whether the type conversion is proper
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
         call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',           &
     &    'Large error occurs when converting REAL(RP) to REAL(DP).')
      end if

      ! Create a MATLAB matrix using the data in X_DP
      px = mxCreateDoubleMatrix(m_mw, n_mw, notComplex)
      call mxCopyReal8ToPtr(x_dp, MX_GET_DBLE(px), m_mw*n_mw)

      end subroutine write_rmatrix


      subroutine write_iscalar(x, px)
      ! WRITE_RSCALAR associates an INTEGER(IK) scalar X with an 
      ! mwPointer PX, after which X can be passed to MATLAB either as an
      ! output of mexFunction or an input of mexCallMATLAB.
      use consts_mod, only : DP, IK
      implicit none

      ! Input
      integer(IK), intent(in) :: x

      ! Output
      mwPointer, intent(out) :: px
      
      ! Intermediate variable
      real(DP) :: x_dp

      x_dp = real(x, kind(x_dp))
      
      if (abs(x - x_dp) > 0.5_DP) then
         call mexErrMsgIdAndTxt('MEXAPI_MOD:ConversionError',           &
     &    'Large error occurs when converting INTEGER(IK) to REAL(RP).')
      end if

      px = mxCreateDoubleScalar(x_dp)

      end subroutine write_iscalar


      subroutine fmxCallMATLAB(fun_ptr, pinput, poutput)
      ! fmxCallMATLAB executes matlab command
      ! output = feval(fun, input),
      ! where fun_ptr is an mwPointer pointing to the function handle of
      ! fun, while pinput and poutput are mwPointer arrays associated 
      ! with the inputs and outputs
      use consts_mod, only : INT32
      implicit none
      
      mwPointer, intent(in) :: fun_ptr      
      mwPointer, intent(in) :: pinput(:)
      mwPointer, intent(out) :: poutput(:)
     
      mwPointer :: aug_pinput(size(pinput) + 1)
      integer(INT32) :: nargin
      integer(INT32) :: nargout
      character(5), parameter :: funFeval = 'feval'
      
      ! Augment the input to include FUN_PTR
      aug_pinput = (/ fun_ptr, pinput /)  

      ! Get number of inputs and number of outputs
      nargin = int(size(aug_pinput), kind(nargin))
      nargout = int(size(poutput), kind(nargout))
      
      ! integer*4 mexCallMATLAB(integer*4 nargout, mwPointer poutput, integer*4 nargin, mwPointer pinput, character*(*) functionName)
      ! If mexCallMATLAB returns 0, the execution is successful.
      if (0 /= mexCallMATLAB(nargout, poutput, nargin,                  &
     & aug_pinput, funFeval)) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:UnsuccessfulCall',         &
     &     'MEX fails to call a MATLAB function.')
      end if

      if (any(poutput == 0)) then
          call mexErrMsgIdAndTxt('MEXAPI_MOD:UnsuccessfulCall',         &
     &    'NULL pointer returned when MEX calls a MATLAB function.')
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
