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
! on 04-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! The mex gateway for NEWUOA
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Friday, June 25, 2021 PM07:39:47


#include "fintrf.h"

! Entry point to Fortran MEX function
      subroutine mexFunction(nargout, poutput, nargin, pinput)
! If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name extension depends on the
! platform), then the following function is callable in matlab:
! [xopt, fopt, info, nf, xhist, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, npt, iprint, maxhist, output_xhist)

! Generic modules
      use consts_mod, only : RP, IK
      use fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use fmxapi_mod, only : fmxVerifyClassShape
      use fmxapi_mod, only : fmxAllocate
      use fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
      use newuoa_mod, only : newuoa
      use prob_mod, only : fun_ptr, calfun

      implicit none

! mexFunction arguments nargout and nargin are of type INTEGER in MATLAB 2019a documents.
      integer, intent(in) :: nargout, nargin
      mwPointer, intent(in) :: pinput(nargin)
      mwPointer, intent(out) :: poutput(nargout)

! Intermediate variables
      integer(IK) :: info
      integer(IK) :: iprint
      integer(IK) :: maxfun
      integer(IK) :: maxhist
      integer(IK) :: nf
      integer(IK) :: npt
      integer(IK) :: output_xhist
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

! Verify that input 1 is a function handle; the other inputs will be verified when read.
      call fmxVerifyClassShape(pinput(1), 'function_handle', 'rank0')

! Read inputs (there are 14)
      fun_ptr = pinput(1)
      ! FUN_PTR is a pointer to the function handle
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
      call fmxReadMPtr(pinput(13), maxhist)
      call fmxReadMPtr(pinput(14), output_xhist)

! Call the Fortran code.
! There are different cases because XHIST may or may not be passed to the Fortran code.
      if (output_xhist > 0) then
          call newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun,&
     & npt, iprint, eta1, eta2, gamma1, gamma2, xhist=xhist, fhist=fhist&
     &, maxhist=maxhist, info=info)
      else
          call newuoa(calfun, x, f, nf, rhobeg, rhoend, ftarget, maxfun,&
     & npt, iprint, eta1, eta2, gamma1, gamma2, fhist=fhist, maxhist=max&
     &hist, info=info)
      end if

! After the Fortran code, XHIST may not be allocated, because it may not have been passed to the
! Fortran code. We allocate it here. Otherwise, fmxWriteMPtr will fail.
      if (.not. allocated(xhist)) then
          call fmxAllocate(xhist, int(size(x), IK), 0_IK)
      end if

! Write outputs
      call fmxWriteMPtr(x, poutput(1))
      call fmxWriteMPtr(f, poutput(2))
      call fmxWriteMPtr(info, poutput(3))
      call fmxWriteMPtr(nf, poutput(4))
      call fmxWriteMPtr(xhist(:, 1:min(int(nf), size(xhist, 2))), poutpu&
     &t(5))
      call fmxWriteMPtr(fhist(1:min(int(nf), size(fhist))), poutput(6), &
     &'row')
! N.B.:
! 1. INT(NF) converts NF to the default integer type; if not, MIN may complain.
! 2. It can happen that 0 < SIZE(XHIST, 2) < MAXHIST or 0 < SIZE(FHIST) < MAXHIST due to the memory
! limit in the Fortran code.

! Free memory. Indeed, automatic deallocation would take place.
      deallocate (x)
      ! Allocated by fmxReadMPtr.
      deallocate (xhist)
      deallocate (fhist)

      return
      end subroutine mexFunction