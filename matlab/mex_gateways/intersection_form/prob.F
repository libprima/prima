!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of prob.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 14-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! PROB_MOD is a module defining the optimization problem. In particular,
! it implements CALFUN and CALCFC.  CALFUN evaluates the objective function
! for unconstrained, bound constrained, and linearly constrained problems;
! CALCFC evaluates the objective function and constraint for nonlinearly
! constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.

#include "fintrf.h"

      module prob_mod

      implicit none
      private
      public :: fun_ptr, calfun

      mwPointer :: fun_ptr ! Pointer to objective function


      contains


! The Fortran subroutine that evaluates the objective function.
      subroutine calfun(x, funval)

! Generic modules
      use consts_mod, only : RP, HUGEFUN, INT32, MSSGLEN
      use infnan_mod, only : is_nan
      use fmxapi_mod, only : mxGetM, mxGetN, mxIsDouble
      use fmxapi_mod, only : mxDestroyArray
      use fmxapi_mod, only : mexErrMsgIdAndTxt
      use fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr, fmxCallMATLAB

      implicit none

! Inputs
      real(RP), intent(in) :: x(:)

! Output
      real(RP), intent(out) :: funval

! Intermediate variables
      mwPointer :: pinput(1), poutput(1)
      mwSize :: row, col
      integer(INT32) :: isdble
      character(len = MSSGLEN) :: eid, mssg

! Associate X with INPUT(1)
      call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function that evaluates the objective function
      call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Verify the class and shape of outputs. Indeed, fmxReadMPtr
! does also the verification. We do it here in order to print a
! more informative error message when the verification fails.
      row = mxGetM(poutput(1))
      col = mxGetN(poutput(1))
      isdble = mxIsDouble(poutput(1))
      if (row*col /= 1 .or. isdble /= 1) then
          eid = 'PROBLEM:ObjectiveNotScalar'
          mssg = 'PROBLEM: Objective function does not return a scalar.'
          call mexErrMsgIdAndTxt(trim(eid), trim(mssg))
      end if

! Read the data in OUTPUT
      call fmxReadMPtr(poutput(1), funval)

! Use extreme barrier to cope with 'hidden constraints'
      if (funval > HUGEFUN .or. is_nan(funval)) then
          funval = HUGEFUN
      end if

! Destroy the matrix created by fmxWriteMPtr for X.
! This must be done. Otherwise, the matrix will be destroyed only
! when the MEX function terminates. However, this subroutine will
! be called maybe thousands of times before that.
      call mxDestroyArray(pinput(1))

      end subroutine calfun


      end module prob_mod