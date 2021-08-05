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
! on 06-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! GETHUGE subroutine
!
! **********************************************************************
!   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!               Department of Applied Mathematics,
!               The Hong Kong Polytechnic University
! **********************************************************************
!
! Last Modified: Sunday, May 23, 2021 PM09:27:33

! Remarks:
!
! 1. Be careful with the "kind" and storage size for integer-type
! (integer, mwSize, mwIndex) variables/functions. Some of them may be
! 32bit, while the others may be 64bit, depending on the machine, the
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
      use fmxapi_mod, only : mxGetN, mxGetString
      use fmxapi_mod, only : mexErrMsgIdAndTxt
      use fmxapi_mod, only : mxCreateDoubleScalar
      use fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use fmxapi_mod, only : fmxVerifyClassShape

      implicit none

! mexFunction arguments:
      integer, intent(in) :: nargout, nargin
      mwPointer, intent(in) :: pinput(nargin)
      mwPointer, intent(out) :: poutput(nargout)

! Variables
      mwSize :: cols
      ! Size of the input
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
      character(len=MSSGLEN) :: eid, mssg

! Check inputs
      call fmxVerifyNArgin(nargin, 1)
      call fmxVerifyNArgout(nargout, 1)
      call fmxVerifyClassShape(pinput(1), 'char', 'row')

! Get the input string.
      cols = int(mxGetN(pinput(1)), kind(cols))
      if (cols > maxlen) then
          eid = 'gethuge:InvalidInput'
          mssg = 'gethuge: The input is too long.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if
      if (mxGetString(pinput(1), data_type, cols) /= 0) then
          eid = 'gethuge:GetInputFail'
          mssg = 'Fail to get the input.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Define hugeValue.
! Note that the REAL values passed to MATLAB via MEX can only be
! doubles. Therefore, we may need to cap the huge values by taking
! min to ensure that they will not overflow when cast to doubles.
      if (data_type == 'float' .or. data_type == 'Float' .or. data_type &
     &== 'FLOAT') then
! No overflow in this case
          hugeValue = real(huge(0.0), DP)
      else if (data_type == 'single' .or. data_type == 'Single' .or. dat&
     &a_type == 'SINGLE') then
! No overflow in this case
          hugeValue = real(huge(0.0_SP), DP)
      else if (data_type == 'double' .or. data_type == 'Double' .or. dat&
     &a_type == 'DOUBLE') then
! No overflow in this case
          hugeValue = huge(0.0_DP)
      else if (data_type == 'real' .or. data_type == 'Real' .or. data_ty&
     &pe == 'REAL') then
! HUGE(0.0_RP) may be bigger thant HUGE(0.0_DP) in this case
          if (RP == DP) then
              hugeValue = huge(0.0_DP)
          else
              hugeValue = 10.0_DP**(min(real(log10(huge(0.0_RP)), DP), l&
     &og10(huge(0.0_DP))) - 1.0_DP)
          end if
      else if (data_type == 'integer' .or. data_type == 'Integer' .or. d&
     &ata_type == 'INTEGER') then
! No overflow in this case
          hugeValue = real(huge(intZero), DP)
      else if (data_type == 'mwSI' .or. data_type == 'mwsi' .or. data_ty&
     &pe == 'MWSI') then
! No overflow in this case
          hugeValue = min(real(huge(msZero), DP), real(huge(miZero), DP)&
     &)
      else if (data_type == 'fun' .or. data_type == 'Fun' .or. data_type&
     & == 'FUN' .or. data_type == 'function' .or. data_type == 'Function&
     &' .or. data_type == 'FUNCTION') then
! HUGEFUN < huge(0.0_DP) according to the definition in CONSTS_MOD.
          hugeValue = real(HUGEFUN, DP)
      else if (data_type == 'con' .or. data_type == 'Con' .or. data_type&
     & == 'CON' .or. data_type == 'constraint' .or. data_type == 'Constr&
     &aint' .or. data_type == 'CONSTRAINT') then
! HUGECON < huge(0.0_DP) according to the definition in CONSTS_MOD.
          hugeValue = real(HUGECON, DP)
      else
          eid = 'GETHUGE:WrongInput'
          mssg = 'GETHUGE: Unrecognized input.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Write output.
! Do NOT use fmxWriteMPtr; when DP /= RP, there is no proper
! implementation of fmxWriteMPtr available in FMXAPI_MOD.
      poutput(1) = mxCreateDoubleScalar(hugeValue)

      return
      end subroutine mexFunction