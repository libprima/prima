!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of cbfun.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun ZHANG (www.zhangzk.net)
! on 12-Feb-2022.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "fintrf.h"

      module cbfun_mod
!--------------------------------------------------------------------------------------------------!
! This module evaluates callback functions received from MATLAB.
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Friday, February 11, 2022 AM11:24:45
!--------------------------------------------------------------------------------------------------!

      implicit none
      private
      public :: evalcb

      interface evalcb
          module procedure evalcb_f, evalcb_fc
      end interface evalcb


      contains


      subroutine evalcb_f(fun_ptr, x, f)
!--------------------------------------------------------------------------------------------------!
! This subroutine evaluates a MATLAB function F = FUN(X). Here, FUN is represented by a mwPointer
! FUN_PTR pointing to FUN, with mwPointer being a type defined in fintrf.h.
!--------------------------------------------------------------------------------------------------!

! Generic modules
      use, non_intrinsic :: consts_mod, only : RP
      use, non_intrinsic :: debug_mod, only : validate

! Fortran MEX API modules
      use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
      use, non_intrinsic :: fmxapi_mod, only : fmxIsDoubleScalar
      use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr&
     &, fmxCallMATLAB

      implicit none

! Inputs
      mwPointer, intent(in) :: fun_ptr
      real(RP), intent(in) :: x(:)

! Outputs
      real(RP), intent(out) :: f

! Local variables
      character(len=*), parameter :: srname = 'EVALCB_F'
      integer :: i
      mwPointer :: pinput(1), poutput(1)

! Associate the input with INPUT.
      call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function.
      call fmxCallMATLAB(fun_ptr, pinput, poutput)

! Destroy the arrays in PINPUT(:).
! This must be done. Otherwise, the array created for X by fmxWriteMPtr will be destroyed only when
! the MEX function terminates, but this subroutine will be called maybe thousands of times before that.
      do i = 1, size(pinput)
          call mxDestroyArray(pinput(i))
      end do

! Read the data in POUTPUT.
! First, verify the class & shape of outputs (even not debugging). Indeed, fmxReadMPtr does also the
! verification. We do it here in order to print a more informative error message in case of failure.
      call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function r&
     &eturns a scalar', srname)
! Second, copy the data.
      call fmxReadMPtr(poutput(1), f)
! Third, destroy the arrays in POUTPUT.
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
      do i = 1, size(poutput)
          call mxDestroyArray(poutput(i))
      end do

      end subroutine evalcb_f


      subroutine evalcb_fc(funcon_ptr, x, f, constr)
!--------------------------------------------------------------------------------------------------!
! This subroutine evaluates a MATLAB function [F, CONSTR] = FUNCON(X). Here, FUN is represented by
! a mwPointer FUNCON_PTR pointing to FUN, with mwPointer being a type defined in fintrf.h.
!--------------------------------------------------------------------------------------------------!

! Generic modules
      use, non_intrinsic :: consts_mod, only : RP
      use, non_intrinsic :: debug_mod, only : validate

! Fortran MEX API modules
      use, non_intrinsic :: fmxapi_mod, only : mxDestroyArray
      use, non_intrinsic :: fmxapi_mod, only : fmxIsDoubleScalar, fmxIsD&
     &oubleVector
      use, non_intrinsic :: fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr&
     &, fmxCallMATLAB

      implicit none

! Inputs
      mwPointer, intent(in) :: funcon_ptr
      real(RP), intent(in) :: x(:)

! Outputs
      real(RP), intent(out) :: f
      real(RP), intent(out) :: constr(:)

! Local variables
      character(len=*), parameter :: srname = 'EVALCB_FC'
      integer :: i
      mwPointer :: pinput(1), poutput(2)
      real(RP), allocatable :: constr_loc(:)

! Associate the input with PINPUT.
      call fmxWriteMPtr(x, pinput(1))

! Call the MATLAB function.
      call fmxCallMATLAB(funcon_ptr, pinput, poutput)

! Destroy the arrays in PINPUT.
! This must be done. Otherwise, the array created for X by fmxWriteMPtr will be destroyed only when
! the MEX function terminates, but this subroutine will be called maybe thousands of times before that.
      do i = 1, size(pinput)
          call mxDestroyArray(pinput(i))
      end do

! Read the data in POUTPUT.
! First, verify the class & shape of outputs (even not debugging). Indeed, fmxReadMPtr does also the
! verification. We do it here in order to print a more informative error message in case of failure.
      call validate(fmxIsDoubleScalar(poutput(1)), 'Objective function r&
     &eturns a real scalar', srname)
      call validate(fmxIsDoubleVector(poutput(2)), 'Constriant function &
     &returns a real vector', srname)
! Second, copy the data.
      call fmxReadMPtr(poutput(1), f)
      call fmxReadMPtr(poutput(2), constr_loc)
! Third, destroy the arrays in POUTPUT.
! MATLAB allocates dynamic memory to store the arrays in plhs (i.e., poutput) for mexCallMATLAB.
! MATLAB automatically deallocates the dynamic memory when you exit the MEX file. However, this
! subroutine will be called maybe thousands of times before that.
! See https://www.mathworks.com/help/matlab/apiref/mexcallmatlab_fortran.html
      do i = 1, size(poutput)
          call mxDestroyArray(poutput(i))
      end do

! Copy CONSTR_LOC to CONSTR.
! Before copying, check that the size of CONSTR_LOC is correct (even if not debugging).
!--------------------------------------------------------------------------------------------------!
! N.B.: We allow SIZE(CONSTR_LOC) == 1 < SIZE(CONSTR). In this case, we set CONSTR = CONSTR_LOC(1).
! The motivation is to allow the MATLAB function to return a scalar when the evaluation fails, where
! the scalar indicates the failure (e.g., NaN, Inf, or a value with an extremely large magnitude).
!--------------------------------------------------------------------------------------------------!
      call validate(size(constr_loc) == size(constr) .or. (size(constr_l&
     &oc) == 1 .and. size(constr) > 1), 'SIZE(CONSTR_LOC) == SIZE(CONSTR&
     &), or SIZE(CONSTR_LOC) == 1 and SIZE(CONSTR) > 0', srname)
      if (size(constr_loc) == size(constr)) then
          constr = constr_loc
      else
          constr = constr_loc(1)
      end if
! Deallocate CONSTR_LOC, allocated by fmxReadMPtr. Indeed, it would be deallocated automatically.
      deallocate (constr_loc)

      end subroutine evalcb_fc


      end module cbfun_mod