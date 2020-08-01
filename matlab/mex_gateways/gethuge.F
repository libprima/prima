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

      subroutine mexFunction(nargout, poutput, nargin, pinput)
C     data_huge = gethuge(data_type)
C     This function returns the largest value of data_type on the
C     current platform, the possible values of data_type being
C     'integer', 'float', 'real', 'single', 'double', 'mwSize',
C     'mwIndex', 'fun', 'function', 'con', 'constraint'.

      use consts_mod, only : IK, SP, DP, RP, HUGEFUN, HUGECON
      use mexapi_mod, only : mxGetN, mxGetString
      use mexapi_mod, only : mexErrMsgIdAndTxt
      use mexapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use mexapi_mod, only : fmxVerifyClassShape, fmxWriteMPtr

      implicit none

C     mexFunction arguments:
      integer, intent(in) :: nargout, nargin
      mwPointer, intent(in) :: pinput(nargin)
      mwPointer, intent(out) :: poutput(nargout)

C     Variables
      mwSize :: cols ! Size of the input
      integer, parameter :: maxlen = 100 
C The largest length of the input string. It cannot be larger than 20
C given the possible inputs.
      character(len=maxlen) :: data_type
C The input string, which specifies the data type
      integer(IK), parameter :: intZero = 0
C Integer zero used by the Fortran code
      mwSize, parameter :: msZero = 0
C Integer zero of type mwSize on the current platform
      mwIndex, parameter :: miZero = 0
C Integer zero of type mwIndex on the current platform
      real(DP) :: hugeValue

C     Validate the number of arguments.
      call fmxVerifyNArgin(nargin, 1)
      call fmxVerifyNArgout(nargout, 1)

C     Validate the class and shape of input.
      call fmxVerifyClassShape(pinput(1), 'char', 'row')


C     Get the input string.
      cols = int(mxGetN(pinput(1)), kind(cols))
      if (mxGetString(pinput(1), data_type, cols) /= 0) then
          call mexErrMsgIdAndTxt ('gethuge:GetInputFail',
     &    'gethuge: Fail to get the input.')
      end if

C     Define hugeValue.
C     Note that the REAL values passed to MATLAB via MEX can only be
C     doubles. Therefore, we may need to cap the huge values by taking
C     min to ensure that they will not overflow when cast to doubles.
      if (data_type == 'float' .or. data_type == 'Float'
     &    .or. data_type == 'FLOAT') then
          hugeValue = real(huge(0.0), DP)
      else if (data_type == 'single' .or. data_type == 'Single'
     &    .or. data_type == 'SINGLE') then
          hugeValue = real(huge(0.0_SP), DP)
      else if (data_type == 'double' .or. data_type == 'Double'
     &    .or. data_type == 'DOUBLE') then
          hugeValue = huge(0.0_DP)
      else if (data_type == 'real' .or. data_type == 'Real'
     &    .or. data_type == 'REAL') then
          hugeValue = huge(0.0_DP)
          if (huge(0.0_RP) < hugeValue) then
              hugeValue = real(huge(0.0_RP), DP)
          end if
      else if (data_type == 'integer' .or. data_type == 'Integer'
     &    .or. data_type == 'INTEGER') then
          hugeValue = real(huge(intZero), DP)
      else if (data_type == 'mwSize' .or. data_type == 'mwsize'
     &     .or. data_type == 'MWSIZE') then
          hugeValue = real(huge(msZero), DP)
      else if (data_type == 'mwIndex' .or. data_type == 'mwindex'
     &     .or. data_type == 'MWINDEX') then
          hugeValue = real(huge(miZero), DP)
      else if (data_type == 'fun' .or. data_type == 'Fun' .or.
     &    data_type == 'FUN' .or. data_type == 'function' .or.
     &    data_type == 'Function' .or. data_type == 'FUNCTION')
     &then
          hugeValue = huge(0.0_DP)
          if (HUGEFUN < hugeValue) then
              hugeValue = real(HUGEFUN, DP)
          end if
      else if (data_type == 'con' .or. data_type == 'Con' .or.
     &    data_type == 'CON' .or. data_type == 'constraint' .or.
     &    data_type == 'Constraint' .or. data_type == 'CONSTRAINT')
     &then
          hugeValue = huge(0.0_DP)
          if (HUGECON < hugeValue) then
              hugeValue = real(HUGECON, DP)
          end if
      else
          call mexErrMsgIdAndTxt ('gethuge:WrongInput',
     &    'gethuge: Unrecognized input.')
      end if

C     Write output.
      call fmxWriteMPtr(hugeValue, poutput(1))

      return
      end subroutine mexFunction
