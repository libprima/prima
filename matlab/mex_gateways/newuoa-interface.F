! The mex gateway for NEWUOA
!
! Coded by Zaikun ZHANG, July 2020


#include "fintrf.h"

! Fortran code starts

      module fnewuoa_mod

      use consts_mod, only : DP, IK
      implicit none
      private
      public :: fun_ptr, nf, fhist

      ! Some global variables
      mwPointer, save :: fun_ptr ! Objective function
      integer(IK), save :: nf ! Number of function evaluations
      real(DP), allocatable, save :: fhist(:) ! Func. evaluation history

      end module fnewuoa_mod


      subroutine mexFunction(nargout, poutput, nargin, pinput)
      ! Entry point to Fortran MEX function
      ! If the binary MEX file is named as FUNCTION_NAME.mex*** 
      ! (file-name extension depends on the platform), then the 
      ! following function is callable in matlab:
      ! [xopt, fopt, info, nf, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, npt, iprint)

      ! General modules
      use consts_mod, only : DP, RP, IK
      use newuoa_mod, only : newuoa
      use mexapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use mexapi_mod, only : fmxVerifyClassShape
      use mexapi_mod, only : fmxAllocate, fmxReadMPtr, fmxWriteMPtr

      ! Solver-specific module
      use fnewuoa_mod, only : fun_ptr, nf, fhist

      implicit none

      ! mexFunction arguments (dummy variables)
      ! nargout and nargin are of type INTEGER in MATLAB 2019a documents
      integer, intent(in) :: nargout, nargin
      mwPointer, intent(in) :: pinput(nargin)
      mwPointer, intent(out) :: poutput(nargout)

      ! Intermediate variables
      real(RP), allocatable :: x(:)
      real(RP) :: f
      real(RP) :: rhobeg, rhoend
      real(RP) :: eta1, eta2
      real(RP) :: gamma1, gamma2
      real(RP) :: ftarget
      integer(IK) :: maxfun, npt, iprint, info

      ! Validate number of arguments
      call fmxVerifyNArgin(nargin, 12)
      call fmxVerifyNArgout(nargout, 5)

      ! Verify that input 1 is a function handle; the other inputs will
      ! be verified when read.
      call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')
      
      ! Read inputs (there are 12)
      fun_ptr = pinput(1)  ! FUN_PTR is a pointer to the function handle
      call fmxReadMPtr(pinput(2), x)
      call fmxReadMPtr(pinput(3), rhobeg)
      call fmxReadMPtr(pinput(4), rhoend)
      call fmxReadMPtr(pinput(5), eta1)
      call fmxReadMPtr(pinput(6), eta2)
      call fmxReadMPtr(pinput(7), gamma1)
      call fmxReadMPtr(pinput(8), gamma2)
      call fmxReadMPtr(pinput(9), ftarget)
      call fmxReadMPtr(pinput(10), maxfun)
      call fmxReadMPtr(pinput(11), npt)
      call fmxReadMPtr(pinput(12), iprint)

      ! Initialize global variables
      nf = 0
      call fmxAllocate(fhist, maxfun)
      fhist = huge(0.0_DP)

      ! Call the Fortran code.
      call newuoa(x, f, rhobeg, rhoend, eta1, eta2, gamma1, gamma2,     &
     & ftarget, npt, maxfun, iprint, info)

      ! Write outputs
      call fmxWriteMPtr(x, poutput(1))
      call fmxWriteMPtr(f, poutput(2))
      call fmxWriteMPtr(info, poutput(3))
      call fmxWriteMPtr(nf, poutput(4))
      call fmxWriteMPtr(fhist(1:nf), poutput(5))

      ! Free memory
      deallocate (x)
      deallocate (fhist)

      return
      end subroutine mexFunction


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! The Fortran subroutine that evaluates the objective function
      subroutine calfun(x, funval)

      ! General modules
      use consts_mod, only : RP, HUGEFUN
      use infnan_mod, only : is_nan
      use mexapi_mod, only : mxGetM, mxGetN, mxIsDouble
      use mexapi_mod, only : mxDestroyArray
      use mexapi_mod, only : mexErrMsgIdAndTxt
      use mexapi_mod, only : fmxReadMPtr, fmxWriteMPtr, fmxCallMATLAB

      ! Solver-specific module
      use fnewuoa_mod, only : fun_ptr, nf, fhist

      implicit none

      ! Inputs
      real(RP), intent(in) :: x(:)

      ! Output
      real(RP), intent(out) :: funval

      ! Intermediate variables
      mwPointer :: pinput(1), poutput(1) 

      ! Associate X with INPUT(1)
      call fmxWriteMPtr(x, pinput(1))

      ! Call the MATLAB function that evaluates the objective function
      call fmxCallMATLAB(fun_ptr, pinput, poutput)

      ! Verify the class and shape of outputs. Indeed, fmxReadMPtr
      ! does also the verification. We do it here in order to print a 
      ! more informative error message when the verification fails.
      if (mxGetM(poutput(1))*mxGetN(poutput(1)) /= 1                    &
     &    .or. mxIsDouble(poutput(1)) /= 1) then
           call mexErrMsgIdAndTxt('fnewuoa:ObjectiveNotScalar',         &
     & 'fnewuoa: The objective function should return a scalar value.')
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

      ! Update global variables
      nf = int(nf + 1, kind(nf))
      fhist(nf) = real(funval, kind(fhist))  ! FHIST is REAL(DP)
      return
      end subroutine calfun
