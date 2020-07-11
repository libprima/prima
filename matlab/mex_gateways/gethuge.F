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

      use consts_mod, only : SP, DP, RP, QP, HUGEFUN, HUGECON
      use mex_mod
      implicit none

C     mexFunction arguments:
      integer, intent(in) :: nlhs, nrhs
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(out) :: plhs(nlhs)

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
      real(SP), parameter :: singleZero = 0.0_SP
C Single-precision floating-point zero
      real(DP), parameter :: doubleZero = 0.0_DP 
C Double-precision floating-point zero
      real(RP), parameter :: realZero = 0.0_RP 
C Floating-point zero of the real precision in use (RP)
      mwSize, parameter :: msZero = 0 
C Integer zero of type mwSize on the current platform
      mwIndex, parameter :: miZero = 0 
C Integer zero of type mwIndex on the current platform
      real(DP) :: hugeValue
      real(DP) :: intHuge, msHuge, miHuge
      real(DP) :: floatHuge, singleHuge, doubleHuge, realHuge
      real(QP) :: doubleHuge_qp

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

C     Get the huge values.
      doubleHuge = huge(doubleZero)
      doubleHuge_qp = real(doubleHuge, QP)
      floatHuge = real(min(real(huge(floatZero),QP), doubleHuge_qp), DP)
      singleHuge = real(min(real(huge(singleZero),QP),doubleHuge_qp),DP)
      realHuge = real(min(real(huge(realZero), QP), doubleHuge_qp), DP)
      intHuge = real(min(real(huge(intZero), QP), doubleHuge_qp), DP)
      msHuge = real(min(real(huge(msZero), QP), doubleHuge_qp), DP)
      miHuge = real(min(real(huge(miZero), QP), doubleHuge_qp), DP)

C     Define hugeValue.
      if (data_type .eq. 'float' .or. data_type .eq. 'Float' 
     &    .or. data_type .eq. 'FLOAT') then
          hugeValue= floatHuge
      else if (data_type .eq. 'single' .or. data_type .eq. 'Single'
     &    .or. data_type .eq. 'SINGLE') then
          hugeValue = singleHuge
      else if (data_type .eq. 'double' .or. data_type .eq. 'Double'
     &    .or. data_type .eq. 'DOUBLE') then
          hugeValue = doubleHuge
      else if (data_type .eq. 'real' .or. data_type .eq. 'Real' 
     &    .or. data_type .eq. 'REAL') then
          hugeValue= realHuge
      else if (data_type .eq. 'integer' .or. data_type .eq. 'Integer'
     &    .or. data_type .eq. 'INTEGER') then
          hugeValue = intHuge
      else if (data_type .eq. 'mwSize' .or. data_type .eq. 'mwsize' 
     &     .or. data_type .eq. 'MWSIZE') then
          hugeValue = msHuge
      else if (data_type .eq. 'mwIndex' .or. data_type .eq. 'mwindex' 
     &     .or. data_type .eq. 'MWINDEX') then
          hugeValue = miHuge
      else if (data_type .eq. 'fun' .or. data_type .eq. 'Fun' .or. 
     &    data_type .eq. 'FUN' .or. data_type .eq. 'function' .or. 
     &    data_type .eq. 'Function' .or. data_type .eq. 'FUNCTION') 
     &then
          hugeValue = real(min(real(HUGEFUN,QP),real(doubleHuge,QP)),DP)
      else if (data_type .eq. 'con' .or. data_type .eq. 'Con' .or. 
     &    data_type .eq. 'CON' .or. data_type .eq. 'constraint' .or. 
     &    data_type .eq. 'Constraint' .or. data_type .eq. 'CONSTRAINT') 
     &then
          hugeValue = real(min(real(HUGECON,QP),real(doubleHuge,QP)),DP)
      else
          call mexErrMsgIdAndTxt ('gethuge:WrongInput',
     &    'gethuge: Unrecognized input.')
      end if

C     Write output.
      plhs(1) = mxCreateDoubleScalar(hugeValue)

      return
      end subroutine mexFunction
