!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of gethuge.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 09-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! GETHUGE subroutine
!
! **********************************************************************
!   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!               Department of Applied Mathematics,
!               The Hong Kong Polytechnic University
! **********************************************************************

! Remarks:
!
! 1. Be careful with the "kind" and storage size for integer-type
! (integer, mwSize, mwIndex) variables/functions. Some of them may be
! 32bit, while they others may be 64bit, depending on the machine, the
! version of matlab, and the compilation option of mex. Do NOT assume
! any two of them to be the same.
! If ever a Segmentation Fault occurs, check these variables first.
!
! 2. Be careful with the line width limit. After preprocessing (macro
! substitution), some lines may become too long and hence get truncated.

#include "fintrf.h"

      subroutine mexFunction(nargout, poutput, nargin, pinput)
! Usage: data_huge = gethuge(data_type)
! This function returns the largest value of data_type on the
! current platform, the possible values of data_type being
! 'integer', 'float', 'real', 'single', 'double', 'mwSI',
! 'fun', 'function', 'con', 'constraint'.
!
! In previous versions, GETHUGE accepted 'mwSize' or 'mwIndex' as
! inputs, but it is not the case anymore. This is because mwSize
! and mwIndex are macros defined in fintrf.h, and they will be
! replaced by other strings after preprocessing.
! Instead of 'mwSize' and 'mwIndex', we now use 'mwSI' to get the
! smaller value between huge(msZero) and huge(miZero).

      use consts_mod, only : IK, SP, DP, RP, HUGEFUN, HUGECON, MSSGLEN
      use mexapi_mod, only : mxGetN, mxGetString
      use mexapi_mod, only : mexErrMsgIdAndTxt
      use mexapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use mexapi_mod, only : fmxVerifyClassShape, fmxWriteMPtr

      implicit none

! mexFunction arguments:
      integer, intent(in) :: nargout, nargin
      mwPointer, intent(in) :: pinput(nargin)
      mwPointer, intent(out) :: poutput(nargout)

! Variables
       mwSize :: cols ! Size of the input
! The largest length of the input string.
      integer, parameter :: maxlen = 50
! The input string, which specifies the data type
      character(len=maxlen) :: data_type
! Integer zero used by the Fortran code
      integer(IK), parameter :: intZero = 0
! Integer zero of type mwSize on the current platform
       mwSize, parameter :: msZero = 0
! Integer zero of type mwIndex on the current platform
       mwIndex, parameter :: miZero = 0
! The huge value that will be returned
      real(RP) :: hugeValue
      character(len = MSSGLEN) :: eid, mssg

! Check inputs
      call fmxVerifyNArgin(nargin, 1)
      call fmxVerifyNArgout(nargout, 1)
      call fmxVerifyClassShape(pinput(1), 'char', 'row')

! Get the input string.
      cols = int(mxGetN(pinput(1)), kind(cols))
      if (cols > maxlen) then
          eid = 'gethuge:InvalidInput'
          mssg = 'gethuge: The input is too long.'
          call mexErrMsgIdAndTxt (eid, mssg)
      end if
      if (mxGetString(pinput(1), data_type, cols) /= 0) then
          eid = 'gethuge:GetInputFail'
          mssg = 'Fail to get the input.'
          call mexErrMsgIdAndTxt (eid, mssg)
      end if

! Define hugeValue.
! Note that the REAL values passed to MATLAB via MEX can only be
! doubles. Therefore, we may need to cap the huge values by taking
! min to ensure that they will not overflow when cast to doubles.
      if (data_type == 'float' .or. data_type == 'Float'.or. data_type =&
     &= 'FLOAT') then
          hugeValue = real(huge(0.0), RP)
      else if (data_type == 'single' .or. data_type == 'Single'.or. data&
     &_type == 'SINGLE') then
          hugeValue = real(huge(0.0_SP), RP)
      else if (data_type == 'double' .or. data_type == 'Double'.or. data&
     &_type == 'DOUBLE') then
          hugeValue = real(huge(0.0_DP), RP)
      else if (data_type == 'real' .or. data_type == 'Real' .or. data_ty&
     &pe == 'REAL') then
          hugeValue = real(huge(0.0_DP), RP)
          if (huge(0.0_RP) < hugeValue) then
              hugeValue = huge(0.0_RP)
          end if
      else if (data_type == 'integer' .or. data_type == 'Integer' .or. d&
     &ata_type == 'INTEGER') then
          hugeValue = real(huge(intZero), RP)
      else if (data_type == 'mwSI' .or. data_type == 'mwsi' .or. data_ty&
     &pe == 'MWSI') then
          hugeValue = min(real(huge(msZero),RP), real(huge(miZero),RP))
      else if (data_type == 'fun' .or. data_type == 'Fun' .or. data_type&
     & == 'FUN' .or. data_type == 'function' .or. data_type == 'Function&
     &' .or. data_type == 'FUNCTION') then
          hugeValue = real(huge(0.0_DP), RP)
          if (HUGEFUN < hugeValue) then
              hugeValue = real(HUGEFUN, RP)
          end if
      else if (data_type == 'con' .or. data_type == 'Con' .or. data_type&
     & == 'CON' .or. data_type == 'constraint' .or. data_type == 'Constr&
     &aint' .or. data_type == 'CONSTRAINT') then
          hugeValue = real(huge(0.0_DP), RP)
          if (HUGECON < hugeValue) then
              hugeValue = real(HUGECON, RP)
          end if
      else
          eid = 'GETHUGE:WrongInput'
          mssg = 'GETHUGE: Unrecognized input.'
          call mexErrMsgIdAndTxt (eid, mssg)
      end if

! Write output.
      call fmxWriteMPtr(hugeValue, poutput(1))

      return
      end subroutine mexFunction