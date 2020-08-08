!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of classical-newuoa-interface.F90.
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


! The mex gateway for NEWUOA (classical version)
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


! Fortran code starts

      module classical_newuoa_mod

      use consts_mod, only : DP
      implicit none
      private
      public :: fun_ptr, nf, fhist
      public :: NEWUOA

! Some global variables
      mwPointer :: fun_ptr ! Objective function
      mwSize :: nf ! Number of function evaluations
! nf should be of type mwSize, because it will be later used as an
! argument of mxCreateDoubleMatrix
      real(DP), allocatable :: fhist(:)
! History of function evaluation

      interface
          subroutine NEWUOA(n, npt, x, rhobeg, rhoend, iprint, maxfun, w&
     &, f, info, ftarget)
          use consts_mod, only : DP
          implicit none
          integer, intent(in) :: n, npt, iprint, maxfun
          integer, intent(out) :: info
          real(DP), intent(in) :: rhobeg, rhoend, ftarget
          real(DP), intent(out) :: f
          real(DP), intent(inout) :: x(n), w(*)
! x(:) or w(:) does not  work !!!
          end subroutine NEWUOA
      end interface

      end module classical_newuoa_mod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
! subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
! Entry point to Fortran MEX function
! If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
! extension depends on the platform), then the following function is
! callable in matlab:
! [xopt, fopt, info, nf, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, ftarget, maxfun, npt, iprint)

      use consts_mod, only : DP, MSSGLEN
      use classical_newuoa_mod, only : fun_ptr, nf, fhist, newuoa
      use mexapi_mod, only : mwOne, notComplex
      use mexapi_mod, only : mxGetM, mxGetN
      use mexapi_mod, only : mexErrMsgIdAndTxt
      use mexapi_mod, only : mxCopyReal8ToPtr, mxCopyPtrToReal8
      use mexapi_mod, only : mxCreateDoubleScalar
      use mexapi_mod, only : mxCreateDoubleMatrix
      use mexapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
      use mexapi_mod, only : fmxVerifyClassShape
      use mexapi_mod, only : fmxGetDble

      implicit none

! mexFunction arguments (dummy variables):
      integer, intent(in) :: nlhs, nrhs
! nlhs and nrhs are of type INTEGER according to the documents of MATLAB R2019a
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(inout) :: plhs(nlhs)

! workspace and intermediate variables
      mwSize :: n
! n is the problem dimension; it is of type mwSize because it is
! obtained by mxGetM
      integer :: n_int
! n_int is the INTEGER cast of n, i.e., n_int = int(n); used when a
! variable of type INTEGER is needed (e.g., when calling NEWUOA).
      real(DP), allocatable :: x(:), w(:)
      real(DP) :: f, rhobeg(1), rhoend(1), ftarget(1)
      real(DP) :: maxfun_r(1), npt_r(1), iprint_r(1)
      integer :: info, iprint, allocate_status
! info and iprint are inputs to NEWUOA, which requires them to be
! INTEGERs (not necessarily the same as mwSize)
      integer :: npt, maxfun, nw
! maxfun and npt  will also be inputs to NEWUOA, which requires them
! to be INTEGERs (not necessarily the same as mwSize).
! nw is the size of workspace. nw = (npt+13)*(npt+n)+3*n*(n+3)/2 + 1 (see below).
! nw is the largest size of memory allocated by this code. If the
! value assigned  to nw is so large that overflow occurs, then there
! will be a Segmentation Fault!!!
! The largest value of nw depends on the storage size of nw (~2E9 for
! kind 4, and ~9E18 for kind 8). This imposes an upper limit on the size
! of problem solvable by this code. If nw is INTEGER*4, assuming that
! npt=2n+1, the largest value of n is ~16000. NEWUOA is not designed for
! so large problems.
      character(len = MSSGLEN) :: eid, mssg

! Check inputs
      call fmxVerifyNArgin(nrhs, 8)
      call fmxVerifyNArgout(nlhs, 5)
      call fmxVerifyClassShape(prhs(1), 'function_handle', 'rank0')
      call fmxVerifyClassShape(prhs(2), 'double', 'vector')
      call fmxVerifyClassShape(prhs(3), 'double', 'scalar')
      call fmxVerifyClassShape(prhs(4), 'double', 'scalar')
      call fmxVerifyClassShape(prhs(5), 'double', 'scalar')
      call fmxVerifyClassShape(prhs(6), 'double', 'scalar')
      call fmxVerifyClassShape(prhs(7), 'double', 'scalar')
      call fmxVerifyClassShape(prhs(8), 'double', 'scalar')

! Read inputs (there are 8)
      fun_ptr = prhs(1)
      n = mxGetM((prhs(2)))*mxGetN(prhs(2)) ! n is of type mwSize
      n_int = int(n, kind(n_int))
! n_int is uUsed when a variable of type INTEGER is needed
      if (n .ne. n_int) then
          eid = 'FNEWUOA:IntError'
          mssg = 'FNEWUOA: n does not equal n_int.'
          call mexErrMsgIdAndTxt(eid, mssg)
      end if

      if (allocated(x)) deallocate (x)
      allocate (x(n_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          eid = 'FNEWUOA:InsufficientMemory'
          mssg = 'FNEWUOA: allocate(x) failed.'
          call mexErrMsgIdAndTxt(eid, mssg)
      end if
      call mxCopyPtrToReal8(fmxGetDble(prhs(2)), x(1:n), n)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
! Note the mwOne is of type mwSize; should not use literal constant 1
! NEVER use literal constants in Fortran mex.
      call mxCopyPtrToReal8(fmxGetDble(prhs(3)), rhobeg, mwOne)
      call mxCopyPtrToReal8(fmxGetDble(prhs(4)), rhoend, mwOne)
      call mxCopyPtrToReal8(fmxGetDble(prhs(5)), ftarget, mwOne)
      call mxCopyPtrToReal8(fmxGetDble(prhs(6)), maxfun_r, mwOne)
      maxfun = int(maxfun_r(1), kind(maxfun))
      call mxCopyPtrToReal8(fmxGetDble(prhs(7)), npt_r, mwOne)
      npt = int(npt_r(1), kind(npt))
      call mxCopyPtrToReal8(fmxGetDble(prhs(8)), iprint_r, mwOne)
      iprint = int(iprint_r(1), kind(iprint))
! maxfun, npt, and iprint will be inputs to NEWUOA, which requires
! them to be INTEGERs (not necessary the same as mwSize)

! Allocate workspace
      if (allocated(w)) deallocate (w)
! nw = (npt+13)*(npt+n)+3*n*(n+3)/2 + 1
      nw = (npt+13)*(npt+n_int)+3*n_int*(n_int+3)/2 + 1
      allocate (w(nw), stat=allocate_status)
      if (allocate_status .ne. 0) then
          eid = 'FNEWUOA:InsufficientMemory'
          mssg = 'FNEWUOA: allocate(w) failed.'
          call mexErrMsgIdAndTxt(eid, mssg)
      end if

! Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          eid = 'FNEWUOA:InsufficientMemory'
          mssg = 'FNEWUOA: allocate(fhist) failed.'
          call mexErrMsgIdAndTxt(eid, mssg)
      end if
      fhist = huge(0.0_DP)

! Call NEWUOA
      call newuoa(n_int, npt, x, rhobeg(1), rhoend(1), iprint, maxfun, w&
     &, f, info, ftarget(1))
! Note that n is of type mwSize, yet NEWUOA expects input 1 to be
! of type INTEGER. Therefore, we should use n_int instead of n.

! Write outputs
      plhs(1) = mxCreateDoubleMatrix(n, mwOne, notComplex)
      call mxCopyReal8ToPtr(x(1:n), fmxGetDble(plhs(1)), n)
      plhs(2) = mxCreateDoubleScalar(f)
! Although info and nf are integers logically, they are passed as double
      plhs(3) = mxCreateDoubleScalar(real(info, DP))
      plhs(4) = mxCreateDoubleScalar(real(nf, DP))
      plhs(5) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(fhist(1:nf), fmxGetDble(plhs(5)), nf)

! Free memory
      deallocate (x)
      deallocate (w)
      deallocate (fhist)

      return
      end subroutine mexFunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The Fortran subroutine that evaluates the objective function
      subroutine calfun(n, x, funval)
      use classical_newuoa_mod, only : fun_ptr, nf, fhist
      use consts_mod, only : DP, INT32, MSSGLEN
      use mexapi_mod, only : mwOne, intOne, intTwo, notComplex
      use mexapi_mod, only : mexErrMsgIdAndTxt, mexCallMATLAB
      use mexapi_mod, only : mxCopyReal8ToPtr, mxCopyPtrToReal8
      use mexapi_mod, only : mxCreateDoubleScalar
      use mexapi_mod, only : mxCreateDoubleMatrix
      use mexapi_mod, only : mxDestroyArray
      use mexapi_mod, only : fmxGetDble

      implicit none

! Dummy variables
      integer, intent(in) :: n
! The type of n is derived from the code of NEWUOA. Thus n should
! be of type INTEGER instead of mwSize
      real(DP), intent(in) :: x(n)
      real(DP), intent(out) :: funval


! intermediate variables
      mwSize :: n_mw
! n_mw is the mwSize cast of n: n_mw = int(n, kind(n_mw));
! used when a variable of type mwSize is needed
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      integer(INT32) :: r
      character(5), parameter :: FEVAL = 'feval'
      real(DP) :: funval_array(1)
      mwPointer :: pf
      character(len = MSSGLEN) :: eid, mssg

! Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      if (n .ne. n_mw) then
          eid = 'FNEWUOA:IntError'
          mssg = 'FNEWUOA: n does not equal n_mw.'
          call mexErrMsgIdAndTxt(eid, mssg)
      end if

      plhs(1) = mxCreateDoubleScalar(huge(0.0_DP))
! Output of f_value = feval(fun, x); see below
      prhs(1) = fun_ptr
! First input of f_value = feval(fun, x); see below; fun_ptr is a global variable
      prhs(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
! second input of f_value = feval(fun, x); see below
      call mxCopyReal8ToPtr(x(1:n), fmxGetDble(prhs(2)), n_mw)
! subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)

      r = mexCallMATLAB(intOne, plhs, intTwo, prhs, FEVAL)
! Execute matlab command: f_value = feval(fun, x)
! integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
      if (r /= 0) then
          eid = 'FNEWUOA:UnsuccessfulCall'
          mssg = 'FNEWUOA: mex fails to call fun.'
          call mexErrMsgIdAndTxt(eid, mssg)
      end if

      if (plhs(1) .eq. 0 .or. fmxGetDble(plhs(1)) .eq. 0) then
          eid = 'FNEWUOA:UnsuccessfulCall'
          mssg = 'FNEWUOA: NULL pointer returned when mex calls fun.'
          call mexErrMsgIdAndTxt(eid, mssg)
      else
          pf = fmxGetDble(plhs(1))
          call mxCopyPtrToReal8(pf, funval_array, mwOne)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
          funval = funval_array(1)
      end if

! Free memory; note that plhs and prhs are just temporary variables in
! this subroutine. We are NOT in mexFunction!
      call mxDestroyArray(plhs(1))
      call mxDestroyArray(prhs(2))
! subroutine mxDestroyArray(mwPointer pm)

! Update global variables
      nf = nf + int(1, kind(nf))
! Some compiler (e.g., g95) may complain about implicit conversion if
! written as nf = nf+1
      fhist(nf) = funval
      return
      end subroutine calfun