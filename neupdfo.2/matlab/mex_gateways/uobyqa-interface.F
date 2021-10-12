C The mex gateway for UOBYQA
C
C **********************************************************************
C   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
C               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
C               Department of Applied Mathematics,
C               The Hong Kong Polytechnic University
C
C   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
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
C defined in it. _MGETDB gets the pointer pointing to a real array.
#if defined _MGETDB
#undef _MGETDB
#endif
#if defined mxGetDoubles
#define _MGETDB mxGetDoubles
#else
#define _MGETDB mxGetPr
#endif

C Fortran code starts

      module fuobyqa

      use consts_mod, only : DP, rp, hugefun, zero
      use infnan_mod
      implicit none

      integer*4, parameter :: int32_0 =  0
      integer, parameter :: int32=kind(int32_0)
      integer(int32), parameter :: notComplex = 0
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
      real(rp), allocatable, save :: fhist(:)
C History of function evaluation
      end module fuobyqa

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
C Entry point to Fortran MEX function
C If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
C extension depends on the platform), then the following function is
C callable in matlab:
C [xopt, fopt, info, nf, fhist] = FUNCTION_NAME(fun, x0, rhobeg, rhoend, maxfun, ftarget)

      use fuobyqa ! module 'consts' is used by 'fuobyqa'
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
      mwPointer, external :: _MGETDB
C mwPointer _MGETDB(mwPointer pm);
C can be replaced by mxGetDoubles since MATLAB R2018b
      mwSize, external :: mxGetM, mxGetN
C mwPointer mxGetM(mwPointer pm), mxGetN(mwPointer pm)
      integer(int32), external :: mxIsDouble, mxIsClass
C integer*4 mxIsDouble(mwPointer pm);
C integer*4 mxIsClass(mwPointer pm, character*(*) classname)

C workspace and intermediate variables
      mwSize :: n
C n is the problem dimension; it is of type mwSize because it is
C obtained by mxGetM
      integer :: n_int
C n_int is the INTEGER cast of n, i.e., n_int = int(n); used when a
C variable of type INTEGER is needed (e.g., when calling UOBYQA).
      real(rp), allocatable :: x(:), w(:)
      real(rp) :: f, rhobeg, rhoend, maxfun_r, ftarget
      integer :: info, iprint, allocate_status
C info and iprint are inputs to UOBYQA, which requires them to be
C INTEGERs (not necessarily the same as mwSize)
      integer :: maxfun, nw
C maxfun will also be an input to UOBYQA, which requires it
C to be an INTEGER (not necessarily the same as mwSize).
C nw is the size of workspace.
C nw = (n*(42+n*(23+n*(8+n)))+max(2*n*n+4,18*n))/4 + 1  (see below).
C nw is the largest size of memory allocated by this code. If the
C value assigned to nw is so large that overflow occurs, then there will
C be a Segmentation Fault!!!
C The largest value of nw depends on the storage size of nw (~2E9 for
C INTEGER*4, and ~9E18 for INTEGER*8). This imposes an upper limit on
C the size of problem solvable by this code. If nw is INTEGER*4, the
C largest value of n is ~300. UOBYQA is not designed for so large
C problems. Indeed, when n > 10, NEWUOA/BOBYQA/LINCOA can solve
C unconstrained problems much more efficiently.

C Validate number of arguments
      if (nrhs .ne. 6) then
         call mexErrMsgIdAndTxt ('fuobyqa:nInput',
     &   'fuobyqa: 6 input argument required.')
C subroutine mexErrMsgIdAndTxt(character*(*) errorid, character*(*) errormsg)
      end if
      if (nlhs .gt. 5) then
         call mexErrMsgIdAndTxt ('fuobyqa:nOutput',
     &   'fuobyqa: At most 5 output arguments.')
      end if

C Validate inputs
C Input 1: fun (function handle)
      if (mxIsClass(prhs(1), 'function_handle') .ne. 1) then
          call mexErrMsgIdAndTxt('fuobyqa:WrongInput',
     &    'fuobyqa: Input 1 should be a function handle.')
      end if
C Input 2: x0 (double column)
      if (mxIsDouble(prhs(2)) .ne. 1 .or. mxGetM(prhs(2)) .lt. 1
     &    .or. mxGetN(prhs(2)) .ne. 1) then
          call mexErrMsgIdAndTxt('fuobyqa:WrongInput',
     &    'fuobyqa: Input 2 should be a column vector of doubles.')
      end if
C Input 3: rhobeg (double scalar)
      if (mxIsDouble(prhs(3)) .ne. 1 .or. mxGetM(prhs(3)) .ne. 1
     &    .or. mxGetN(prhs(3)) .ne. 1) then
          call mexErrMsgIdAndTxt('fuobyqa:WrongInput',
     &    'fuobyqa: Input 3 should be a double.')
      end if
C Input 4: rhoend (double scalar)
      if (mxIsDouble(prhs(4)) .ne. 1 .or. mxGetM(prhs(4)) .ne. 1
     &    .or. mxGetN(prhs(4)) .ne. 1) then
          call mexErrMsgIdAndTxt('fuobyqa:WrongInput',
     &    'fuobyqa: Input 4 should be a double.')
      end if
C Input 5: maxfun (double scalar)
      if (mxIsDouble(prhs(5)) .ne. 1 .or. mxGetM(prhs(5)) .ne. 1
     &    .or. mxGetN(prhs(5)) .ne. 1) then
          call mexErrMsgIdAndTxt('fuobyqa:WrongInput',
     & 'fuobyqa: Input 5 should be a double (with an integer value).')
      end if
C Although input 5 (maxfun) is an integer logically,
C it has to be passed to the mexified code as a double variable. In
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

C Input 6: ftarget (double scalar)
      if (mxIsDouble(prhs(6)) .ne. 1 .or. mxGetM(prhs(6)) .ne. 1
     &    .or. mxGetN(prhs(6)) .ne. 1) then
          call mexErrMsgIdAndTxt('fuobyqa:WrongInput',
     &    'fuobyqa: Input 6 should be a double.')
      end if

C Read inputs (there are 6)
      fun_ptr = prhs(1)
      n = mxGetM(prhs(2)) ! This is why n should be of type mwSize
      n_int = int(n, kind(n_int))
C i_int is used when a variable of type INTEGER is needed
      if (n .ne. n_int) then
          call mexErrMsgIdAndTxt('fuobyqa:IntError',
     &    'fuobyqa: n does not equal n_int.')
      end if

      if (allocated(x)) deallocate (x)
      allocate (x(n_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fuobyqa:InsufficientMemory',
     &    'fuobyqa: allocate(x) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(2)), x(1:n), n)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

      call mxCopyPtrToReal8(_MGETDB(prhs(3)), rhobeg, mwOne)
      call mxCopyPtrToReal8(_MGETDB(prhs(4)), rhoend, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
C Note the mwOne is of type mwSize; should not use literal constant 1
C NEVER use literal constants in Fortran mex.

C Check the values of rhobeg and rhoend. We do not check the values of
C other inputs (e.g., n, maxfun, npt) because the Fortran code does it
      if (rhobeg .le. zero .or. rhobeg .lt. rhoend .or.
     &rhoend .lt. zero) then
          call mexErrMsgIdAndTxt('fuobyqa:InvalidRhobegRhoend',
     &'fuobyqa: rhobeg and rhoend do not satisfy rhobeg >= rhobeg > 0.')
      end if

      call mxCopyPtrToReal8(_MGETDB(prhs(5)),
     &     maxfun_r, mwOne)
      maxfun = int(maxfun_r, kind(maxfun))
C maxfun will be an input to UOBYQA, which requires maxfun to be
C an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(6)), ftarget, mwOne)

C Allocate workspace
      if (allocated(w)) deallocate (w)
C      nw = (n*(42+n*(23+n*(8+n)))+max(2*n*n+4,18*n))/4 + 1
      nw = (n_int*(42+n_int*(23+n_int*(8+n_int))) +
     &     max(2*n_int*n_int+4,18*n_int))/4 + 1
      allocate (w(nw), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fuobyqa:InsufficientMemory',
     &    'fuobyqa: allocate(w) failed.')
      end if

C Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fuobyqa:InsufficientMemory',
     &    'fuobyqa: allocate(fhist) failed.')
      end if
      fhist = huge(0.0_DP)

C Call UOBYQA
      iprint = 0
      call UOBYQA(n_int, x, rhobeg, rhoend, iprint, maxfun, w, f, info,
     &     ftarget)
C Note that n is of type mwSize, yet UOBYQA expects input 1 to be
C of type INTEGER. Therefore, we should use n_int instead of n.

C Write outputs
      plhs(1) = mxCreateDoubleMatrix(n, mwOne, notComplex)
      call mxCopyReal8ToPtr(x(1:n), _MGETDB(plhs(1)), n)
      plhs(2) = mxCreateDoubleScalar(f)
C Although info and nf are integers logically, they are passed as double
      plhs(3) = mxCreateDoubleScalar(real(info, rp))
      plhs(4) = mxCreateDoubleScalar(real(nf, rp))
      plhs(5) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(fhist(1:nf),
     &     _MGETDB(plhs(5)), nf)

C Free memory
      deallocate (x)
      deallocate (w)
      deallocate (fhist)

      return
      end subroutine mexFunction

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C The Fortran subroutine that evaluates the objective function
      subroutine calfun(n, x, funval)
      use fuobyqa ! module 'consts' is used by 'fuobyqa'
      implicit none

C Dummy variables
      integer, intent(in) :: n
C The type of n is derived from the code of UOBYQA. Thus n should
C be of type INTEGER instead of mwSize
      real(rp), intent(in) :: x(n)
      real(rp), intent(out) :: funval

C Function declarations
      integer(int32), external :: mexCallMATLAB, mxIsDouble
C integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
C integer*4 mxIsDouble(mwPointer pm);
      mwSize, external :: mxGetM, mxGetN
C mwPointer mxGetM(mwPointer pm), mxGetN(mwPointer pm)
      mwPointer, external :: mxCreateDoubleMatrix
C mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwPointer, external :: mxCreateDoubleScalar
C mwPointer mxCreateDoubleScalar(real*8 value)
      mwPointer, external :: _MGETDB
C mwPointer _MGETDB(mwPointer pm)

C intermediate variables
      mwSize :: n_mw
C n_mw is the mwSize cast of n: n_mw = int(n, kind(n_mw));
C used when a variable of type mwSize is needed
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      integer(int32), parameter :: intOne=1, intTwo=2
      character(5), parameter :: funFeval = 'feval'
C intOne, intTwo, and funFeval are used in mexCallMATLAB

C Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      if (n .ne. n_mw) then
          call mexErrMsgIdAndTxt('fuobyqa:IntError',
     &    'fuobyqa: n does not equal n_mw.')
      end if

      plhs(1) = mxCreateDoubleScalar(huge(0.0_DP))
C Output of f_value = feval(fun, x); see below
      prhs(1) = fun_ptr
C First input of f_value = feval(fun, x); see below; fun_ptr is a global variable
      prhs(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
C Second input of f_value = feval(fun, x); see below
      call mxCopyReal8ToPtr(x(1:n), _MGETDB(prhs(2)), n_mw)
C subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)

      if (0 .ne. mexCallMATLAB(intOne, plhs, intTwo,
     &    prhs, funFeval)) then
C Execute matlab command: f_value = feval(fun, x)
C integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
          call mexErrMsgIdAndTxt('fuobyqa:UnsuccessfulCall',
     &    'fuobyqa: mex fails to call fun.')
      end if

      if (plhs(1) .eq. 0 .or. _MGETDB(plhs(1)) .eq. 0) then
          call mexErrMsgIdAndTxt('fuobyqa:UnsuccessfulCall',
     &    'fuobyqa: NULL pointer returned when mex calls fun.')
      end if

      if (mxGetM(plhs(1))*mxGetN(plhs(1)) .ne. 1
     &    .or. mxIsDouble(plhs(1)) .ne. 1) then
           call mexErrMsgIdAndTxt('fuobyqa:ObjectiveNotScalar',
     & 'fuobyqa: The objective function should return a scalar value.')
      end if

      call mxCopyPtrToReal8(_MGETDB(plhs(1)), funval, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

C Use extreme barrier to cope with 'hidden constraints'
      if (funval .gt. hugefun .or. is_nan(funval)) then
          funval = hugefun ! hugefun is defined in consts
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
