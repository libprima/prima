C The mex gateway for NEWUOA
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
C 2. Note that MEX generally use double precesion for real values. It
C is not necessarily the case in the Fortran code. Therefore, explicit
C type conversion is necessary whenever real values are exchanged
C between Fortran and MATLAB. Type mismatch will lead to errors like
C Segmentation Fault.
C
C 3. Be careful with the line width limit. After preprocessing (macro
C substitution), some lines may become too long and hence get truncated.
C This is also why we can NOT put a continuation character & at the 73rd
C column --- it may not be the 73rd column after preprocessing!

#include "fintrf.h"

C Use mxGetDoubles instead of mxGetPr if possible, the former being
C available since MATLAB R2018b. The following macros should be put
C after the header "fintrf.h" is included, because mxGetDoubles is
C defined in it. MX_GET_DBLE gets the pointer pointing to a real array.
#ifdef MX_GET_DBLE
#undef MX_GET_DBLE
#endif
#ifdef mxGetDoubles
#define MX_GET_DBLE mxGetDoubles
#else
#define MX_GET_DBLE mxGetPr
#endif

C Fortran code starts

      module fnewuoa_mod

      use consts_mod, only : DP
      implicit none
      private
      public :: fun_ptr, nf, fhist
      public :: NEWUOA

C Some global variables
      mwPointer, save :: fun_ptr ! Objective function
      mwSize, save :: nf ! Number of function evaluations
C nf should be of type mwSize, because it will be later used as an
C argument of mxCreateDoubleMatrix
      real(DP), allocatable, save :: fhist(:)
C History of function evaluation

      interface
          subroutine NEWUOA(npt, x, rhobeg, rhoend, iprint, maxfun,
     &        f, info, ftarget)
          use consts_mod, only : RP, IK
          implicit none
          integer(IK), intent(in) :: npt, iprint, maxfun
          integer(IK), intent(out) :: info
          real(RP), intent(in) :: rhobeg, rhoend, ftarget
          real(RP), intent(out) :: f
          real(RP), intent(inout) :: x(:)  ! x(*) does not work !!!
          end subroutine NEWUOA
      end interface

      end module fnewuoa_mod

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!      interface
!          subroutine mexFunction(nlhs, plhs, nrhs, prhs)
!              integer, intent(in) :: nlhs, nrhs
!              ! N.B.: Do NOT write plhs(:) or prhs(:) 
!              mwPointer, intent(in) :: prhs(*)
!              mwPointer, intent(out) :: plhs(*)
!          end subroutine mexFunction
!      end interface 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
C Entry point to Fortran MEX function
C If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
C extension depends on the platform), then the following function is
C callable in matlab:
C [xopt, fopt, info, nf, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, maxfun, npt, ftarget)

      use consts_mod, only : DP, RP, IK, ONE
      use fnewuoa_mod, only : fun_ptr, nf, fhist, NEWUOA
      use mex_mod
      implicit none

C mexFunction arguments (dummy variables):
      integer, intent(in) :: nlhs, nrhs
C nlhs and nrhs are of type INTEGER according to the documents of MATLAB R2019a
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(out) :: plhs(nlhs)

C workspace and intermediate variables
      mwSize :: n_mw
C n_mw is the problem dimension; it is of type mwSize because it is
C obtained by mxGetM
      integer(IK) :: n
C n is the INTEGER(IK) cast of n_mw, i.e., n = int(n_mw, IK); used when
C calling the Fortran code
      real(RP), allocatable :: x(:)
C x is the variable vector in the Fortran code
      real(DP), allocatable :: x_dp(:)
C x_dp is the REAL(DP) cast of x; used in the MEX subroutines
      real(RP) :: f, rhobeg, rhoend, ftarget
      real(DP) :: f_dp, rhobeg_dp(1), rhoend_dp(1), ftarget_dp(1)
C f, rhobeg, rhoend, and ftarget are inputs to the Fortran code, which
C expects them to be REAL(RP), not necessarily the same as REAL(DP).
C f_dp, rhobeg_dp, rhoend_dp, and ftarget_dp are REAL(DP) the versions.
C rhoend_dp, rhoend_dp, and ftarget_dp need to be arrays as required by
C mxCopyPtrToReal8.
      integer(IK) :: info, iprint, npt, maxfun
      real(DP) :: maxfun_dp(1), npt_dp(1)
C info, iprint, npt, and maxfun are inputs to the Fortran code, which
C expects them to be INTEGER(IK), not necessarily the same as INTEGER
C or mwSize. maxfun_dp and npt_dp are the values of maxfun and npt
C received from MATLAB, which are double arrays.
      integer(IK) :: allocate_status

C Validate number of arguments
      if (nrhs /= 7) then
         call mexErrMsgIdAndTxt ('fnewuoa:nInput',
     &   'fnewuoa: 7 input argument required.')
C subroutine mexErrMsgIdAndTxt(character*(*) errorid, character*(*) errormsg)
      end if
      if (nlhs > 5) then
         call mexErrMsgIdAndTxt ('fnewuoa:nOutput',
     &   'fnewuoa: At most 5 output arguments.')
      end if

C Validate inputs
C Input 1: fun (function handle)
      if (mxIsClass(prhs(1), 'function_handle') /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 1 should be a function handle.')
      end if
C Input 2: x0 (double column)
      if (mxIsDouble(prhs(2)) /= 1 .or. mxGetM(prhs(2)) < 1
     &    .or. mxGetN(prhs(2)) /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 2 should be a column vector of doubles.')
      end if
C Input 3: rhobeg (double scalar)
      if (mxIsDouble(prhs(3)) /= 1 .or. mxGetM(prhs(3)) /= 1
     &    .or. mxGetN(prhs(3)) /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 3 should be a double.')
      end if
C Input 4: rhoend (double scalar)
      if (mxIsDouble(prhs(4)) /= 1 .or. mxGetM(prhs(4)) /= 1
     &    .or. mxGetN(prhs(4)) /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 4 should be a double.')
      end if
C Input 5: maxfun (double scalar)
      if (mxIsDouble(prhs(5)) /= 1 .or. mxGetM(prhs(5)) /= 1
     &    .or. mxGetN(prhs(5)) /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     & 'fnewuoa: Input 5 should be a double (with an integer value).')
      end if
C Input 6: npt (double scalar)
      if (mxIsDouble(prhs(6)) /= 1 .or. mxGetM(prhs(6)) /= 1
     &    .or. mxGetN(prhs(6)) /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     & 'fnewuoa: Input 6 should be a double (with an integer value).')
      end if
C Although inputs 5 and 6 (maxfun and npt) are integers logically,
C they have to be passed to the mexified code as double variables. In
C mex, data is passed by pointers, but there are only very limited
C functions that can read an integer value from a pointer or write
C an interger value to a pointer (mxCopyPtrToInteger1,
C mxCopyInteger1ToPtr, mxCopyPtrToInteger2, mxCopyInteger2ToPtr,
C mxCopyPtrToInteger4, mxCopyInteger4ToPtr; no function for
C INTEGER*8). This makes it impossible to pass integer data properly
C unless we know the kind of the integer. Therefore, in general, it
C is recommended to pass integers as double variables and then cast
C them back to integers when needed. Indeed, in MATLAB, even if we
C define maxfun = 1000, the class of maxfun is double! To get an
C integer maxfun, we would have to define maxfun = INT32(1000) or
C maxfun = INT64(1000)!

C Input 7: ftarget (double scalar)
      if (mxIsDouble(prhs(7)) /= 1 .or. mxGetM(prhs(7)) /= 1
     &    .or. mxGetN(prhs(7)) /= 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 7 should be a double.')
      end if

C Read inputs (there are 7)
      fun_ptr = prhs(1)
      n_mw = mxGetM(prhs(2))
      n = int(n_mw, kind(n))
      if (n /= n_mw) then
          call mexErrMsgIdAndTxt('fnewuoa:IntError',
     &    'fnewuoa: n does not equal n_mw.')
      end if
      if (allocated(x_dp)) deallocate (x_dp)
      allocate (x_dp(n), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fnewuoa:InsufficientMemory',
     &    'fnewuoa: allocate(x_dp) failed.')
      end if
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(2)), x_dp(1:n_mw), n_mw)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(3)), rhobeg_dp, mwOne)
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(4)), rhoend_dp, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
C Note the mwOne is of type mwSize; should not use literal constant 1
C NEVER use literal constants in Fortran mex.
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(5)),
     &     maxfun_dp, mwOne)
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(6)),
     &     npt_dp, mwOne)
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(7)),
     &     ftarget_dp, mwOne)

C Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(int(maxfun_dp(1))), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fnewuoa:InsufficientMemory',
     &    'fnewuoa: allocate(fhist) failed.')
      end if
      fhist = huge(0.0_DP)

C Cast the inputs to the types expected by the Fortran code.
      maxfun = int(maxfun_dp(1), kind(maxfun))
      npt = int(npt_dp(1), kind(npt))
      if (abs(real(maxfun, DP) - maxfun_dp(1)) >= 0.5_DP .or.
     &    abs(real(npt, DP) - npt_dp(1)) >= 0.5_DP) then
          call mexErrMsgIdAndTxt('fnewuoa:RealError',
     &    'fnewuoa: Error during type conversion.')
      end if

      if (allocated(x)) deallocate (x)
      allocate (x(n), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fnewuoa:InsufficientMemory',
     &    'fnewuoa: allocate(x) failed.')
      end if
      x = real(x_dp, RP)
      rhobeg = real(rhobeg_dp(1), kind(rhobeg))
      rhoend = real(rhoend_dp(1), kind(rhoend))
      ftarget = real(ftarget_dp(1), kind(ftarget))
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE) .or.
     &    abs(rhobeg-rhobeg_dp(1)) > convTol*max(abs(rhobeg), ONE) .or.
     &    abs(rhoend-rhoend_dp(1)) > convTol*max(abs(rhoend), ONE) .or.
     &    abs(ftarget-ftarget_dp(1)) > convTol*max(abs(ftarget), ONE))
     &    then
          call mexErrMsgIdAndTxt('fnewuoa:RealError',
     &    'fnewuoa: Error during type conversion.')
      end if

C Call the Fortran code.
      iprint = 0
      call NEWUOA(npt, x, rhobeg, rhoend, iprint, maxfun,f,info,ftarget)

C Cast the outputs to the types expected by the MEX subroutines.
      x_dp = real(x, DP)
      f_dp = real(f, DP)
      if (maxval(abs(x - x_dp)) > convTol*max(maxval(abs(x)), ONE)
     &    .or. abs(f - f_dp) > convTol*max(abs(f), ONE)) then
          call mexErrMsgIdAndTxt('fnewuoa:RealError',
     &    'fnewuoa: Error during type conversion.')
      end if
C Write outputs
      plhs(1) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
      call mxCopyReal8ToPtr(x_dp(1:n_mw), MX_GET_DBLE(plhs(1)), n_mw)
      plhs(2) = mxCreateDoubleScalar(f_dp)
C Although info and nf are integers logically, they are passed as double
      plhs(3) = mxCreateDoubleScalar(real(info, DP))
      plhs(4) = mxCreateDoubleScalar(real(nf, DP))
      plhs(5) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(fhist(1:nf),
     &     MX_GET_DBLE(plhs(5)), nf)

C Free memory
      deallocate (x)
      deallocate (x_dp)
      deallocate (fhist)

      return
      end subroutine mexFunction

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!      interface
!          subroutine calfun(n, x, funval)
!          use consts_mod, only : RP, IK
!          implicit none
!          integer(IK), intent(in) :: n
!          real(RP), intent(in) :: x(n)
!          real(RP), intent(out) :: funval
!          end subroutine calfun
!      end interface
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The Fortran subroutine that evaluates the objective function
      subroutine calfun(n, x, funval)
      use fnewuoa_mod, only : fun_ptr, nf, fhist
      use consts_mod, only : DP, RP, IK, HUGEFUN, ONE
      use infnan_mod, only : is_nan
      use mex_mod
      implicit none

C Dummy variables
C The types of n, x, and funval are derived from the code of NEWUOA.
C Thus n should be of type INTEGER(IK) instead of mwSize, while x and
C funval should be of type REAL(RP) instead of REAL(DP).
      integer(IK), intent(in) :: n
      real(RP), intent(in) :: x(n)
      real(RP), intent(out) :: funval

C intermediate variables
      mwSize :: n_mw
C n_mw is the mwSize cast of n: n_mw = int(n, kind(n_mw));
C used when a variable of type mwSize is needed
      real(DP) :: x_dp(n), funval_dp(1)
C x_dp is the REAL(DP) cast of x: x_dp = real(x, DP);
C used when a variable of type REAL(DP) is needed
C funval_dp is the REAL(DP) cast of x: funval_dp = real(funval, DP);
C used when a variable of type REAL(DP) is needed
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      character(5), parameter :: funFeval = 'feval'

C Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      if (n_mw /= n) then
          call mexErrMsgIdAndTxt('fnewuoa:IntError',
     &    'fnewuoa: n_mw does not equal n.')
      end if
      x_dp = real(x(1:n), DP)  ! cast x to REAL(DP)
      if (maxval(abs(x-x_dp)) > convTol*max(maxval(abs(x)), ONE)) then
          call mexErrMsgIdAndTxt('fnewuoa:RealError',
     &    'fnewuoa: x_dp does not equal x.')
      end if


      plhs(1) = mxCreateDoubleScalar(huge(0.0_DP))
C Output of f_value = feval(fun, x); see below
      prhs(1) = fun_ptr
C First input of f_value = feval(fun, x); see below
      prhs(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
C second input of f_value = feval(fun, x); see below
      call mxCopyReal8ToPtr(x_dp(1:n), MX_GET_DBLE(prhs(2)), n_mw)
C subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)

      if (0 /= mexCallMATLAB(intOne, plhs, intTwo,
     &    prhs, funFeval)) then
C Execute matlab command: f_value = feval(fun, x)
C integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
          call mexErrMsgIdAndTxt('fnewuoa:UnsuccessfulCall',
     &    'fnewuoa: mex fails to call fun.')
      end if

      if (plhs(1) .eq. 0 .or. MX_GET_DBLE(plhs(1)) .eq. 0) then
          call mexErrMsgIdAndTxt('fnewuoa:UnsuccessfulCall',
     &    'fnewuoa: NULL pointer returned when mex calls fun.')
      end if

      if (mxGetM(plhs(1))*mxGetN(plhs(1)) /= 1
     &    .or. mxIsDouble(plhs(1)) /= 1) then
           call mexErrMsgIdAndTxt('fnewuoa:ObjectiveNotScalar',
     & 'fnewuoa: The objective function should return a scalar value.')
      end if

      call mxCopyPtrToReal8(MX_GET_DBLE(plhs(1)), funval_dp, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

C Use extreme barrier to cope with 'hidden constraints'
      if (funval_dp(1) > HUGEFUN .or. is_nan(funval_dp(1))) then
          funval_dp(1) = HUGEFUN 
      end if

      funval = real(funval_dp(1), RP)  ! cast funval_dp to REAL(RP)
      if (abs(funval-funval_dp(1)) > convTol*max(abs(funval), ONE)) then
          call mexErrMsgIdAndTxt('fnewuoa:RealError',
     &    'fnewuoa: funval does not equal funval_dp.')
      end if

C Free memory; note that plhs and prhs are just temporary variables in
C this subroutine. We are NOT in mexFunction!
      call mxDestroyArray(plhs(1))
      call mxDestroyArray(prhs(2))
C subroutine mxDestroyArray(mwPointer pm)

C Update global variables
      nf = nf + int(1, kind(nf))
C Some compiler (e.g., g95) may complain about implicit conversion if
C written as nf = nf+1
      fhist(nf) = funval_dp(1)  ! Note that fhist is REAL(DP)
      return
      end subroutine calfun
