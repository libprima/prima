C GETHUGE subroutine
C
C **********************************************************************
C   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
C               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
C               Department of Applied Mathematics,
C               The Hong Kong Polytechnic University
C **********************************************************************

C Remarks: 
C
C 1. Be careful with the "kind" and storage size for integer-type
C (integer, mwSize, mwIndex) variables/functions. Some of them may be
C 32bit, while they others may be 64bit, depending on the machine, the
C version of matlab, and the compilation option of mex. Do NOT assume
C any two of them to be the same. 
C If ever a Segmentation Fault occurs, check these variables first.
C
C 2. Be careful with the line width limit. After preprocessing (macro
C substitution), some lines may become too long and hence get truncated

#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C     data_huge = gethuge(data_type)
C     This function returns the largest value of data_type on the
C     current platform, the possible values of data_type being
C     'integer', 'float', 'real', 'single', 'double', 'mwSize', 
C     'mwIndex', 'fun', 'function', 'con', 'constraint'.

      use pdfoconst ! See pdfoconst.F. It defines sp, dp, and int4
      implicit none

C     mexFunction arguments:
      integer, intent(in) :: nlhs, nrhs
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(out) :: plhs(nlhs)

C     Function declarations:
      mwPointer, external :: mxCreateDoubleScalar 
C mwPointer mxCreateDoubleScalar(real*8 value)
      mwSize, external :: mxGetM, mxGetN 
C mwPointer mxGetM(mwPointer pm)
      integer(int4), external :: mxGetString  
C integer*4 mxGetString(mwPointer pm, character*(*) str, mwSize strlen)
      integer(int4), external :: mxIsChar 
C integer*4 mxIsChar(mwPointer pm)

C     Variables
      mwSize :: rows, cols ! Size of the input
      integer, parameter :: maxlen = 20 
C The largest length of the input string. It cannot be larger than 20
C given the possible inputs. 
      character(len=maxlen) :: data_type 
C The input string, which specifies the data type
      integer, parameter :: intZero = 0 
C Integer zero of the default kind on the current platform  
      real, parameter :: floatZero = 0.0 
C Floating-point zero of the default real kind on the current platform
      real(sp), parameter :: singleZero = 0.0 
C Single-precision floating-point zero
      real(dp), parameter :: doubleZero = 0.0D0 
C Double-precision floating-point zero
      mwSize, parameter :: msZero = 0 
C Integer zero of type mwSize on the current platform
      mwIndex, parameter :: miZero = 0 
C Integer zero of type mwIndex on the current platform
      integer, parameter :: intHuge = huge(intZero)
      real, parameter :: floatHuge = huge(floatZero)
      real(sp), parameter :: singleHuge = huge(singleZero)
      real(dp), parameter :: doubleHuge = huge(doubleZero)
      mwSize, parameter :: msHuge = huge(msZero)
      mwIndex, parameter :: miHuge = huge(miZero)
      real(dp) :: hugeValue

C     Validate the number of arguments.
      if (nrhs .ne. 1) then
          call mexErrMsgIdAndTxt ('gethuge:nInput',
     &    'gethuge: There should be 1 input argument.')
      end if

      if (nlhs .gt. 1) then
          call mexErrMsgIdAndTxt ('gethuge:nOutput',
     &    'gethuge: At most 1 output argument allowed.')
      end if

C     The input must be a string.
      rows = mxGetM(prhs(1))
      cols = mxGetN(prhs(1))
      if (mxIsChar(prhs(1)) .ne. 1
     &   .or. rows .gt. 1 .or. cols .gt. maxlen) then
          call mexErrMsgIdAndTxt ('gethuge:WrongInput',
     &   'gethuge: The input can only be a string.')
      end if

C     Get the input string.
      if (mxGetString(prhs(1), data_type, cols) .ne. 0) then
          call mexErrMsgIdAndTxt ('gethuge:GetInputFail',
     &    'gethuge: Fail to get the input.')
      end if

C     Define hugeValue.
      if (data_type .eq. 'integer' .or. data_type .eq. 'Integer'
     &    .or. data_type .eq. 'INTEGER') then
          hugeValue = real(intHuge, dp)
      else if (data_type .eq. 'float' .or. data_type .eq. 'real' 
     &    .or. data_type .eq. 'Float' .or. data_type .eq. 'Real' 
     &    .or. data_type .eq. 'FLOAT' .or. data_type .eq. 'REAL') then
          hugeValue = real(floatHuge, dp)
      else if (data_type .eq. 'single' .or. data_type .eq. 'Single'
     &    .or. data_type .eq. 'SINGLE') then
          hugeValue = singleHuge
      else if (data_type .eq. 'double' .or. data_type .eq. 'Double'
     &    .or. data_type .eq. 'DOUBLE') then
          hugeValue = doubleHuge
      else if (data_type .eq. 'mwSize' .or. data_type .eq. 'mwsize' 
     &     .or. data_type .eq. 'MWSIZE') then
          hugeValue = real(msHuge, dp)
      else if (data_type .eq. 'mwIndex' .or. data_type .eq. 'mwindex' 
     &     .or. data_type .eq. 'MWINDEX') then
          hugeValue = real(miHuge, dp)
      else if (data_type .eq. 'fun' .or. data_type .eq. 'Fun' .or. 
     &    data_type .eq. 'FUN' .or. data_type .eq. 'function' .or. 
     &    data_type .eq. 'Function' .or. data_type .eq. 'FUNCTION') 
     &then
          hugeValue = real(HUGEFUN, dp) ! HUGEFUN defined in pdfoconst
      else if (data_type .eq. 'con' .or. data_type .eq. 'Con' .or. 
     &    data_type .eq. 'CON' .or. data_type .eq. 'constraint' .or. 
     &    data_type .eq. 'Constraint' .or. data_type .eq. 'CONSTRAINT') 
     &then
          hugeValue = real(HUGECON, dp) ! HUGECON defined in pdfoconst
      else
          call mexErrMsgIdAndTxt ('gethuge:WrongInput',
     &    'gethuge: Unrecognized input.')
      end if

C     Write output.
      plhs(1) = mxCreateDoubleScalar(hugeValue)

      return
      end subroutine mexFunction
