!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of newuoa-interface.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 12-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! The mex gateway for NEWUOA
!
! Coded by Zaikun Zhang in July 2020.


#include "fintrf.h"

! Entry point to Fortran MEX function
      subroutine mexFunction(nargout, poutput, nargin, pinput)
! If the binary MEX file is named as FUNCTION_NAME.mex***
! (file-name extension depends on the platform), then the
! following function is callable in matlab:
! [xopt, fopt, info, nf, xhist, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, maxhist, npt, iprint)

! Generic modules
      use consts_mod, only : RP, IK
      use newuoa_mod, only : newuoa
      use fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use fmxapi_mod, only : fmxVerifyClassShape
      use fmxapi_mod, only : fmxAllocate
      use fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
      use prob_mod, only : fun_ptr, calfun

      implicit none

! mexFunction arguments (dummy variables)
! nargout and nargin are of type INTEGER in MATLAB 2019a documents
      integer, intent(in) :: nargout, nargin
      mwPointer, intent(in) :: pinput(nargin)
      mwPointer, intent(out) :: poutput(nargout)

! Intermediate variables
      integer(IK) :: info
      integer(IK) :: iprint
      integer(IK) :: maxfun
      integer(IK) :: maxfhist
      integer(IK) :: maxxhist
      integer(IK) :: nf
      integer(IK) :: npt
      real(RP) :: eta1
      real(RP) :: eta2
      real(RP) :: f
      real(RP) :: ftarget
      real(RP) :: gamma1
      real(RP) :: gamma2
      real(RP) :: rhobeg
      real(RP) :: rhoend
      real(RP), allocatable :: fhist(:)
      real(RP), allocatable :: x(:)
      real(RP), allocatable :: xhist(:, :)

! Validate number of arguments
      call fmxVerifyNArgin(nargin, 14)
      call fmxVerifyNArgout(nargout, 6)

! Verify that input 1 is a function handle; the other inputs will
! be verified when read.
      call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read inputs (there are 14)
      fun_ptr = pinput(1) ! FUN_PTR is a pointer to the function handle
      call fmxReadMPtr(pinput(2), x)
      call fmxReadMPtr(pinput(3), rhobeg)
      call fmxReadMPtr(pinput(4), rhoend)
      call fmxReadMPtr(pinput(5), eta1)
      call fmxReadMPtr(pinput(6), eta2)
      call fmxReadMPtr(pinput(7), gamma1)
      call fmxReadMPtr(pinput(8), gamma2)
      call fmxReadMPtr(pinput(9), ftarget)
      call fmxReadMPtr(pinput(10), maxfun)
      call fmxReadMPtr(pinput(11), maxxhist)
      call fmxReadMPtr(pinput(12), maxfhist)
      call fmxReadMPtr(pinput(13), npt)
      call fmxReadMPtr(pinput(14), iprint)

! Call the Fortran code.
! There are various cases because XHIST/FHIST may or may not be passed to
! the Fortran code.
      if (maxxhist > 0 .and. maxfhist > 0) then
          call newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun,&
     & npt, iprint, eta1, eta2, gamma1, gamma2, xhist = xhist, fhist = f&
     &hist, maxhist = maxxhist, info = info)
      else if (maxxhist > 0 .and. maxfhist <= 0) then
          call newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun,&
     & npt, iprint, eta1, eta2, gamma1, gamma2, xhist = xhist, maxhist =&
     & maxxhist, info = info)
      else if (maxxhist <= 0 .and. maxfhist > 0) then
          call newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun,&
     & npt, iprint, eta1, eta2, gamma1, gamma2, fhist = fhist, maxhist =&
     & maxfhist, info = info)
      else
          call newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun,&
     & npt, iprint, eta1, eta2, gamma1, gamma2, info = info)
      end if

! After the Fortran code, XHIST or FHIST may not be allocated, because
! they may not have been passed to the Fortran code. We allocate them
! here. Otherwise, fmxWriteMPtr will fail.
      if (.not. allocated(xhist)) then
          call fmxAllocate(xhist, int(size(x), IK), 0_IK)
      end if
      if (.not. allocated(fhist)) then
          call fmxAllocate(fhist, 0_IK)
      end if

! Write outputs
      call fmxWriteMPtr(x, poutput(1))
      call fmxWriteMPtr(f, poutput(2))
      call fmxWriteMPtr(info, poutput(3))
      call fmxWriteMPtr(nf, poutput(4))
      call fmxWriteMPtr(xhist(:, 1 : min(nf, maxxhist)), poutput(5))
      call fmxWriteMPtr(fhist(1 : min(nf, maxfhist)), poutput(6), 'row')

! Free memory
! X was allocated by fmxReadMPtr.
      deallocate(x)
! XHIST was allocated by the Fortran code or fmxAllocate.
      deallocate(xhist)
! FHIST was allocated by the Fortran code or fmxAllocate.
      deallocate(fhist)

      return
      end subroutine mexFunction