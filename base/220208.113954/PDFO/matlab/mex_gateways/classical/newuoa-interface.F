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
C 2. Be careful with the line width limit. After preprocessing (macro
C substitution), some lines may become too long and hence get truncated

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

      module fnewuoa
      use pdfoconst ! See pdfoconst.F, which defines dp, int4, HUGENUM
      implicit none

      integer(int4), parameter :: notComplex = 0
C notComplex is used in mxCreateDoubleMatrix, whose signature is
C mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwSize, parameter :: mwOne = 1 ! Integer 1 with type mwSize
C mwOne is used in mxCreateDoubleMatrix and mxCopyPtrToReal8,
C whose signatures are
C mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

C Some global variables
      mwPointer, save :: fun_ptr ! Objective function
      mwSize, save :: nf ! Number of function evaluations
C nf should be of type mwSize, because it will be later used as an
C argument of mxCreateDoubleMatrix
      real(dp), allocatable, save :: fhist(:)
C History of function evaluation
      end module fnewuoa

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
C Entry point to Fortran MEX function
C If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
C extension depends on the platform), then the following function is
C callable in matlab:
C [xopt, fopt, info, nf, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, maxfun, npt, ftarget)

      use fnewuoa ! module 'pdfoconst' is used by 'fnewuoa'
      implicit none

C mexFunction arguments (dummy variables):
      integer, intent(in) :: nlhs, nrhs
C nlhs and nrhs are of type INTEGER according to the documents of MATLAB R2019a
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(inout) :: plhs(nlhs)

C function declarations:
      mwPointer, external :: mxCreateDoubleMatrix
C mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwPointer, external :: mxCreateDoubleScalar
C mwPointer mxCreateDoubleScalar(real*8 value)
      mwPointer, external :: MX_GET_DBLE
C mwPointer MX_GET_DBLE(mwPointer pm);
C can be replaced by mxGetDoubles since MATLAB R2018b
      mwSize, external :: mxGetM, mxGetN
C mwPointer mxGetM(mwPointer pm), mxGetN(mwPointer pm)
      integer(int4), external :: mxIsDouble, mxIsClass
C integer*4 mxIsDouble(mwPointer pm);
C integer*4 mxIsClass(mwPointer pm, character*(*) classname)

C workspace and intermediate variables
      mwSize :: n
C n is the problem dimension; it is of type mwSize because it is
C obtained by mxGetM
      integer :: n_int
C n_int is the INTEGER cast of n, i.e., n_int = int(n); used when a
C variable of type INTEGER is needed (e.g., when calling NEWUOA).
      real(dp), allocatable :: x(:), w(:)
      real(dp) :: f, rhobeg, rhoend, maxfun_r, npt_r, ftarget
      integer :: info, iprint, allocate_status
C info and iprint are inputs to NEWUOA, which requires them to be
C INTEGERs (not necessarily the same as mwSize)
      integer :: npt, maxfun, nw
C maxfun and npt  will also be inputs to NEWUOA, which requires them
C to be INTEGERs (not necessarily the same as mwSize).
C nw is the size of workspace. nw = (npt+13)*(npt+n)+3*n*(n+3)/2 + 1 (see below).
C nw is the largest size of memory allocated by this code. If the
C value assigned  to nw is so large that overflow occurs, then there
C will be a Segmentation Fault!!!
C The largest value of nw depends on the storage size of nw (~2E9 for
C kind 4, and ~9E18 for kind 8). This imposes an upper limit on the size
C of problem solvable by this code. If nw is INTEGER*4, assuming that
C npt=2n+1, the largest value of n is ~16000. NEWUOA is not designed for
C so large problems.

C Validate number of arguments
      if (nrhs .ne. 7) then
         call mexErrMsgIdAndTxt ('fnewuoa:nInput',
     &   'fnewuoa: 7 input argument required.')
C subroutine mexErrMsgIdAndTxt(character*(*) errorid, character*(*) errormsg)
      end if
      if (nlhs .gt. 5) then
         call mexErrMsgIdAndTxt ('fnewuoa:nOutput',
     &   'fnewuoa: At most 5 output arguments.')
      end if

C Validate inputs
C Input 1: fun (function handle)
      if (mxIsClass(prhs(1), 'function_handle') .ne. 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 1 should be a function handle.')
      end if
C Input 2: x0 (double column)
      if (mxIsDouble(prhs(2)) .ne. 1 .or. mxGetM(prhs(2)) .lt. 1
     &    .or. mxGetN(prhs(2)) .ne. 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 2 should be a column vector of doubles.')
      end if
C Input 3: rhobeg (double scalar)
      if (mxIsDouble(prhs(3)) .ne. 1 .or. mxGetM(prhs(3)) .ne. 1
     &    .or. mxGetN(prhs(3)) .ne. 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 3 should be a double.')
      end if
C Input 4: rhoend (double scalar)
      if (mxIsDouble(prhs(4)) .ne. 1 .or. mxGetM(prhs(4)) .ne. 1
     &    .or. mxGetN(prhs(4)) .ne. 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 4 should be a double.')
      end if
C Input 5: maxfun (double scalar)
      if (mxIsDouble(prhs(5)) .ne. 1 .or. mxGetM(prhs(5)) .ne. 1
     &    .or. mxGetN(prhs(5)) .ne. 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     & 'fnewuoa: Input 5 should be a double (with an integer value).')
      end if
C Input 6: npt (double scalar)
      if (mxIsDouble(prhs(6)) .ne. 1 .or. mxGetM(prhs(6)) .ne. 1
     &    .or. mxGetN(prhs(6)) .ne. 1) then
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
C integer maxfun, we would have to define maxfun = int32(1000) or
C maxfun = int64(1000)!

C Input 7: ftarget (double scalar)
      if (mxIsDouble(prhs(7)) .ne. 1 .or. mxGetM(prhs(7)) .ne. 1
     &    .or. mxGetN(prhs(7)) .ne. 1) then
          call mexErrMsgIdAndTxt('fnewuoa:WrongInput',
     &    'fnewuoa: Input 7 should be a double.')
      end if

C Read inputs (there are 7)
      fun_ptr = prhs(1)
      n = mxGetM(prhs(2)) ! This is why n should be of type mwSize
      n_int = int(n, kind(n_int))
C n_int is uUsed when a variable of type INTEGER is needed
      if (n .ne. n_int) then
          call mexErrMsgIdAndTxt('fnewuoa:IntError',
     &    'fnewuoa: n does not equal n_int.')
      end if

      if (allocated(x)) deallocate (x)
      allocate (x(n_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fnewuoa:InsufficientMemory',
     &    'fnewuoa: allocate(x) failed.')
      end if
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(2)), x(1:n), n)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(3)), rhobeg, mwOne)
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(4)), rhoend, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
C Note the mwOne is of type mwSize; should not use literal constant 1
C NEVER use literal constants in Fortran mex.

C Check the values of rhobeg and rhoend. We do not check the values of
C other inputs (e.g., n, maxfun, npt) because the Fortran code does it
      if (rhobeg .le. ZERO .or. rhobeg .lt. rhoend .or.
     &rhoend .lt. ZERO) then
          call mexErrMsgIdAndTxt('fnewuoa:InvalidRhobegRhoend',
     &'fnewuoa: rhobeg and rhoend do not satisfy rhobeg >= rhobeg > 0.')
      end if

      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(5)),
     &     maxfun_r, mwOne)
      maxfun = int(maxfun_r, kind(maxfun))
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(6)), npt_r, mwOne)
      npt = int(npt_r, kind(npt))
C maxfun and npt will be inputs to NEWUOA, which requires them to be
C INTEGERs (not necessary the same as mwSize)
      call mxCopyPtrToReal8(MX_GET_DBLE(prhs(7)), ftarget, mwOne)

C Allocate workspace
      if (allocated(w)) deallocate (w)
C      nw = (npt+13)*(npt+n)+3*n*(n+3)/2 + 1
      nw = (npt+13)*(npt+n_int)+3*n_int*(n_int+3)/2 + 1
      allocate (w(nw), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fnewuoa:InsufficientMemory',
     &    'fnewuoa: allocate(w) failed.')
      end if

C Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fnewuoa:InsufficientMemory',
     &    'fnewuoa: allocate(fhist) failed.')
      end if
      fhist = HUGENUM

C Call NEWUOA
      iprint = 0
      call NEWUOA(n_int, npt, x, rhobeg, rhoend, iprint, maxfun, w, f,
     &     info, ftarget)
C Note that n is of type mwSize, yet NEWUOA expects input 1 to be
C of type INTEGER. Therefore, we should use n_int instead of n.

C Write outputs
      plhs(1) = mxCreateDoubleMatrix(n, mwOne, notComplex)
      call mxCopyReal8ToPtr(x(1:n), MX_GET_DBLE(plhs(1)), n)
      plhs(2) = mxCreateDoubleScalar(f)
C Although info and nf are integers logically, they are passed as double
      plhs(3) = mxCreateDoubleScalar(real(info, dp))
      plhs(4) = mxCreateDoubleScalar(real(nf, dp))
      plhs(5) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(fhist(1:nf),
     &     MX_GET_DBLE(plhs(5)), nf)

C Free memory
      deallocate (x)
      deallocate (w)
      deallocate (fhist)

      return
      end subroutine mexFunction

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C The Fortran subroutine that evaluates the objective function
      subroutine calfun(n, x, funval)
      use fnewuoa ! module 'pdfoconst' is used by 'fnewuoa'
      implicit none

C Dummy variables
      integer, intent(in) :: n
C The type of n is derived from the code of NEWUOA. Thus n should
C be of type INTEGER instead of mwSize
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: funval

C Function declarations
      integer(int4), external :: mexCallMATLAB
C integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
      mwPointer, external :: mxCreateDoubleMatrix
C mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwPointer, external :: mxCreateDoubleScalar
C mwPointer mxCreateDoubleScalar(real*8 value)
      mwPointer, external :: MX_GET_DBLE
C mwPointer MX_GET_DBLE(mwPointer pm)

C intermediate variables
      mwSize :: n_mw
C n_mw is the mwSize cast of n: n_mw = int(n, kind(n_mw));
C used when a variable of type mwSize is needed
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      integer(int4), parameter :: intOne=1, intTwo=2
      character(5), parameter :: funFeval = 'feval'
C intOne, intTwo, funFeval are  used when calling mexCallMATLAB

C Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      if (n .ne. n_mw) then
          call mexErrMsgIdAndTxt('fnewuoa:IntError',
     &    'fnewuoa: n does not equal n_mw.')
      end if

      plhs(1) = mxCreateDoubleScalar(HUGENUM)
C Output of f_value = feval(fun, x); see below
      prhs(1) = fun_ptr
C First input of f_value = feval(fun, x); see below; fun_ptr is a global variable
      prhs(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
C second input of f_value = feval(fun, x); see below
      call mxCopyReal8ToPtr(x(1:n), MX_GET_DBLE(prhs(2)), n_mw)
C subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)

      if (0 .ne. mexCallMATLAB(intOne, plhs, intTwo,
     &    prhs, funFeval)) then
C Execute matlab command: f_value = feval(fun, x)
C integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
          call mexErrMsgIdAndTxt('fnewuoa:UnsuccessfulCall',
     &    'fnewuoa: mex fails to call fun.')
      end if

      if (plhs(1) .eq. 0 .or. MX_GET_DBLE(plhs(1)) .eq. 0) then
          call mexErrMsgIdAndTxt('fnewuoa:UnsuccessfulCall',
     &    'fnewuoa: NULL pointer returned when mex calls fun.')
      else
          call mxCopyPtrToReal8(MX_GET_DBLE(plhs(1)),
     &         funval, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
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
      fhist(nf) = funval
      return
      end subroutine calfun
