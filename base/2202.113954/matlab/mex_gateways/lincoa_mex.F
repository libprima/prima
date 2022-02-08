C The mex gateway for LINCOA
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

      module flincoa ! Some global variables

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
      real(rp), allocatable, save :: chist(:)
C History of constraint violation
      integer, save :: mresmax, iaresmax
C Size of linear constraint matrix (there are mresmax linear
C constraints, the gradient of each has size iaresmax, which should be
C the same as problem dimension n);
      real(rp), allocatable, save :: Aresmax(:,:), bresmax(:)
C Linear constraint matrix and right-hand side; will be used when calculating constraint violation
      end module flincoa

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
C Entry point to Fortran MEX function
C If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
C extension depends on the platform), then the following function is
C callable in matlab:
C [xopt, fopt, info, nf, fhist, constrviolation, chist] = FUNCTION_NAME(fun, x0, A, b, rhobeg, rhoend, maxfun, npt, ftarget)

      use flincoa ! module 'consts' is used by 'flincoa'
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
      mwSize :: n, m, ia
C n is the problem dimension; m is the number of linear constraints, and
C ia is the length of the gradient of each constraint, which should be n
C they are of type mwSize because they are obtained by mxGetM/mxGetN
      integer :: n_int, m_int, ia_int
C n_int is the INTEGER cast of n, i.e., n_int = int(n); used when a
C variable of type INTEGER is needed (e.g., when calling LINCOA).
      mwIndex :: i, j
C i and j are the index for the entries of A; since m and ia are of type
C mwSize, i and j are of the type mwIndex; see the do loop below
      real(rp), allocatable :: x(:), A(:, :), b(:), w(:)
      real(rp) :: f, rhobeg, rhoend,maxfun_r,npt_r,ftarget,resmax,cval
      integer :: info, iprint, allocate_status
C info and iprint are inputs to LINCOA, which requires them to be
C INTEGERs (not necessarily the same as mwSize)
      integer :: npt, maxfun, nw
C npt and maxfun will also be inputs to LINCOA, which requires
C npt and maxfun to be INTEGERs (not necessarily the same as mwSize).
C nw is the size of workspace.
C nw = m*(2+n)+npt*(4+n+npt)+n*(9+3*n)+max(m+3*n,max(2*m+n,2*npt)) (see
C below) nw is the largest size of memory allocated by this code. If the
C value assigned to nw is so large that overflow occurs, then there will
C be a Segmentation Fault!!!
C The largest value of nw depends on the storage size of nw (~2E9 for
C INTEGER*4, and ~9E18 for INTEGER*8). This imposes an upper limit on
C the size of problem solvable by this code. If nw is INTEGER*4,
C assuming that m=10n and npt=2n+1, the largest value of n is ~10000.
C LINCOA is not designed for so large problems.

C Validate number of arguments
      if (nrhs .ne. 9) then
         call mexErrMsgIdAndTxt ('flincoa:nInput',
     &   'flincoa: 9 input argument required.')
C subroutine mexErrMsgIdAndTxt(character*(*) errorid, character*(*) errormsg)
      end if
      if (nlhs .gt. 7) then
         call mexErrMsgIdAndTxt ('flincoa:nOutput',
     &   'flincoa: At most 7 output arguments.')
      end if

C Validate inputs
C Input 1: fun (function handle)
      if (mxIsClass(prhs(1), 'function_handle') .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &    'flincoa: Input 1 should be a function handle.')
      end if
C Input 2: x0 (double column)
      if (mxIsDouble(prhs(2)) .ne. 1 .or. mxGetM(prhs(2)) .lt. 1
     &    .or. mxGetN(prhs(2)) .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &    'flincoa: Input 2 should be a column vector of doubles.')
      end if
C Input 3: A (double matrix; empty or size(A, 1)=length(x0); different
C from MATLAB by a transpose)
      if (mxIsDouble(prhs(3)) .ne. 1 .or.
     & (mxGetM(prhs(3))*mxGetN(prhs(3)) .gt. 0 .and.
     & mxGetM(prhs(3)) .ne. mxGetM(prhs(2)))) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &    'flincoa: Input 3 should be a matrix of doubles; it is either
     &empty or its number of rows equals the length of input 2.')
      end if
C Input 4: b (double column)
      if (mxIsDouble(prhs(4)) .ne. 1 .or. mxGetM(prhs(4)) .ne.
     &    mxGetN(prhs(3)) .or. mxGetN(prhs(4)) .gt. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &      'flincoa: Input 4 should be a column of doubles; its length
     &equals the number of columns of input 3.')
      end if
C Input 5: rhobeg (double scalar)
      if (mxIsDouble(prhs(5)) .ne. 1 .or. mxGetM(prhs(5)) .ne. 1
     &    .or. mxGetN(prhs(5)) .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &    'flincoa: Input 5 should be a double.')
      end if
C Input 6: rhoend (double scalar)
      if (mxIsDouble(prhs(6)) .ne. 1 .or. mxGetM(prhs(6)) .ne. 1
     &    .or. mxGetN(prhs(6)) .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &    'flincoa: Input 6 should be a double.')
      end if
C Input 7: maxfun (double scalar)
      if (mxIsDouble(prhs(7)) .ne. 1 .or. mxGetM(prhs(7)) .ne. 1
     &    .or. mxGetN(prhs(7)) .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     & 'flincoa: Input 7 should be a double (with an integer value).')
      end if
C Input 8: npt (double scalar)
      if (mxIsDouble(prhs(8)) .ne. 1 .or. mxGetM(prhs(8)) .ne. 1
     &    .or. mxGetN(prhs(8)) .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     & 'flincoa: Input 8 should be a double (with an integer value).')
      end if
C Although inputs 7 and 8 (maxfun and npt) are integers logically,
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

C Input 9: ftarget (double scalar)
      if (mxIsDouble(prhs(9)) .ne. 1 .or. mxGetM(prhs(9)) .ne. 1
     &    .or. mxGetN(prhs(9)) .ne. 1) then
          call mexErrMsgIdAndTxt('flincoa:WrongInput',
     &    'flincoa: Input 9 should be a double.')
      end if

C Read inputs (there are 9)
      fun_ptr = prhs(1)
      n = mxGetM(prhs(2)) ! This is why n should be mwSize
      n_int = int(n, kind(n_int))
C n_int is used when a variable of type INTEGER is needed
      m = mxGetN(prhs(3)) ! This is why m should be mwSize
      m_int = int(m, kind(m_int))
C m_int is used when a variable of type INTEGER is needed
      ia = mxGetM(prhs(3)) ! This is why ia should be mwSize
      ia_int = int(ia, kind(ia_int))
C ia_int is used when a variable of type INTEGER is needed
      if (n .ne. n_int) then
          call mexErrMsgIdAndTxt('flincoa:IntError',
     &    'flincoa: n does not equal n_int.')
      end if
      if (m .ne. m_int) then
          call mexErrMsgIdAndTxt('flincoa:IntError',
     &    'flincoa: m does not equal m_int.')
      end if
      if (ia .ne. ia_int) then
          call mexErrMsgIdAndTxt('flincoa:IntError',
     &    'flincoa: ia does not equal ia_int.')
      end if

      if (allocated(x)) deallocate (x)
      allocate (x(n_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(x) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(2)), x(1:n), n)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

      if (allocated(A)) deallocate (A)
      allocate (A(ia_int, m_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(A) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(3)),
     &     A(1:ia,1:m), ia*m)

      if (allocated(b)) deallocate (b)
      allocate (b(m_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(b) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(4)), b(1:m), m)

      call mxCopyPtrToReal8(_MGETDB(prhs(5)), rhobeg, mwOne)
      call mxCopyPtrToReal8(_MGETDB(prhs(6)), rhoend, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
C Note the mwOne is of type mwSize; should not use literal constant 1
C NEVER use literal constants in Fortran mex.

C Check the values of rhobeg and rhoend. We do not check the values of
C other inputs (e.g., n, maxfun, npt) because the Fortran code does it
      if (rhobeg .le. zero .or. rhobeg .lt. rhoend .or.
     &rhoend .lt. zero) then
          call mexErrMsgIdAndTxt('flincoa:InvalidRhobegRhoend',
     &'flincoa: rhobeg and rhoend do not satisfy rhobeg >= rhobeg > 0.')
      end if

      call mxCopyPtrToReal8(_MGETDB(prhs(7)),
     &     maxfun_r, mwOne)
      maxfun = int(maxfun_r, kind(maxfun))
C maxfun will be an input to LINCOA, which requires maxfun to be
C an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(8)), npt_r, mwOne)
      npt = int(npt_r, kind(npt))
C npt will be an input to LINCOA, which requires maxfun to be
C an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(9)), ftarget, mwOne)

C Allocate workspace
      if (allocated(w)) deallocate (w)
C      nw = m*(2+n)+npt*(4+n+npt)+n*(9+3*n)+max(m+3*n,max(2*m+n,2*npt))
C nw is of type INTEGER
      nw = m_int*(2+n_int)+npt*(4+n_int+npt)+n_int*(9+3*n_int)+
     &      max(m_int+3*n_int,max(2*m_int+n_int,2*npt))
      allocate (w(nw), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(w) failed.')
      end if

C Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(fhist) failed.')
      end if
      fhist = huge(0.0_DP)

      if (allocated(chist)) deallocate (chist)
      allocate (chist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(chist) failed.')
      end if
      chist = huge(0.0_DP)

      mresmax = m_int
      iaresmax = ia_int
      if (allocated(Aresmax)) deallocate (Aresmax)
      allocate (Aresmax(ia_int, m_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(Aresmax) failed.')
      end if
      Aresmax(1:ia, 1:m) = A(1:ia, 1:m)

      if (allocated(bresmax)) deallocate (bresmax)
      allocate (bresmax(m_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('flincoa:InsufficientMemory',
     &    'flincoa: allocate(bresmax) failed.')
      end if
      bresmax(1:m) = b(1:m)

C Call LINCOA
      iprint = 0
      call LINCOA(n_int, npt, m_int, A, ia_int, b, x, rhobeg, rhoend,
     &     iprint, maxfun, w, f, info, ftarget)
C Note that n, m, and ia are of type mwSize, yet LINCOA expects input 1,
C 3, and 5 to be of type INTEGER. Therefore, we should use n_int, m_int,
C and ia_int instead.

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

      resmax = zero
C Constraint violation (named 'RESMAX' in Powell's COBYLA code)
C zero is defined in module 'consts'.
C Linear constraints: A'*x <= b (NOTE the transpose)
      do j = 1, m
          cval = -b(j)
          do i = 1, ia
              cval = cval + A(i,j)*x(i)
          end do
          if (is_nan(cval)) then
              resmax = cval ! Set resmax=NaN if constraint contains NaN
              exit
          else
              resmax = max(resmax, cval) ! Absolute containt violation
          end if
      end do
      plhs(6) = mxCreateDoubleScalar(resmax)
      plhs(7) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(chist(1:nf),
     &     _MGETDB(plhs(7)), nf)

C Free memory
      deallocate (x)
      deallocate (A)
      deallocate (b)
      deallocate (w)
      deallocate (fhist)
      deallocate (chist)
      deallocate (Aresmax)
      deallocate (bresmax)

      end subroutine mexFunction

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C The Fortran subroutine that evaluates the objective function
      subroutine calfun(n, x, funval)
      use flincoa ! module 'consts' is used by 'flincoa'
      implicit none

C dummy variables
      integer, intent(in) :: n
C The type of n is derived from the code of LINCOA. Thus n should
C be of type INTEGER instead of mwSize
      real(rp), intent(in) :: x(n)
      real(rp), intent(out) :: funval

C function declarations
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
      integer :: i, j
C i and j are the index for the entries of Aresmax; since mresmax and
C iaresmax are of type INTEGER, i and j are of type INTEGER. See the do
C loop below
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      integer(int32), parameter :: intOne=1, intTwo=2
      character(5), parameter :: funFeval = 'feval'
C intOne, intTwo, and funFeval are used when calling mexCallMATLAB
      real(rp) :: resmax, cval ! constraint violation/value

C Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      if (n .ne. n_mw) then
          call mexErrMsgIdAndTxt('flincoa:IntError',
     &    'flincoa: n does not equal n_mw.')
      end if

      plhs(1) = mxCreateDoubleScalar(huge(0.0_DP))
C First output of f_value = feval(fun, x); see below
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
          call mexErrMsgIdAndTxt('flincoa:UnsuccessfulCall',
     &    'flincoa: mex fails to call fun.')
      end if

      if (plhs(1) .eq. 0 .or. _MGETDB(plhs(1)) .eq. 0) then
          call mexErrMsgIdAndTxt('flincoa:UnsuccessfulCall',
     &    'flincoa: NULL pointer returned when mex calls fun.')
      end if

      if (mxGetM(plhs(1))*mxGetN(plhs(1)) .ne. 1
     &    .or. mxIsDouble(plhs(1)) .ne. 1) then
           call mexErrMsgIdAndTxt('flincoa:ObjectiveNotScalar',
     & 'flincoa: The objective function should return a scalar value.')
      end if

      call mxCopyPtrToReal8(_MGETDB(plhs(1)),funval, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

C Use extreme barrier to cope with 'hidden constraints'
      if (funval .gt. hugefun .or. is_nan(funval)) then
          funval = hugefun ! hugefun is defined in consts
      end if

      resmax = zero
C Constraint violation (named 'RESMAX' in Powell's COBYLA code)
C zero is defined in module 'consts'.
C Linear constraints: Aresmax'*x <= bresmax (NOTE the transpose)
      do j = 1, mresmax
          cval = -bresmax(j)
          do i = 1, iaresmax
              cval = cval + aresmax(i,j)*x(i)
          end do
          if (is_nan(cval)) then
              resmax = cval ! Set resmax=NaN if constraint contains NaN
              exit
          else
              resmax = max(resmax, cval) !Absolute constraint violation
          end if
      end do

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
      chist(nf) = resmax

      end subroutine calfun
