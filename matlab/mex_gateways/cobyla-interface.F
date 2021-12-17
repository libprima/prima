C The mex gateway for COBYLA
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

      module fcobyla ! Some global variables

      use consts_mod, only : DP, rp, hugecon, hugefun, zero
      use infnan_mod
      use pintrf_mod, only : FUNCON
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
      mwPointer, save :: fun_ptr, constr_ptr
C Objective function and constraints (pointers)
      mwSize, save :: nf ! Number of function evaluations
C nf should be of type mwSize, because it will be later used as an
C argument of mxCreateDoubleMatrix
      real(rp), allocatable, save :: fhist(:)
C History of function evaluation
      real(rp), allocatable, save :: chist(:)
C History of constraint violation
      real(rp), allocatable, save :: conval_x0(:)
C Constraint value at x0
      end module fcobyla

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
C Entry point to Fortran MEX function
C If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
C extension depends on the platform), then the following function is
C callable in matlab:
C [xopt, fopt, info, nf, fhist, convalopt, constrviolation, chist] = FUNCTION_NAME(fun, constr, x0, rhobeg, rhoend, maxfun, m, ftarget, conval_x0)

      use fcobyla
      use cobyla_mod, only : cobyla
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
      mwSize :: n, m
      integer :: n_int, m_int
C n is the problem dimension; it is of type mwSize because it is
C obtained by mxGetM. m is the number of constraints. It is of type
C mwSize because it will be used in mxCopyPtrToReal8.
C n_int and m_int are the INTEGER cast of n and m, i.e., n_int = int(n),
C m_int = int(m), which will be used when variables of type INTEGER are
C needed (e.g., when calling COBYLA).
      real(rp) :: f_x0, ctol
      real(rp), allocatable :: x(:), conval(:)
      real(rp) :: f, rhobeg, rhoend, maxfun_r, m_r, ftarget, resmax
      integer :: info, iprint, allocate_status
C info and iprint are inputs to COBYLA, which requires them to be
C INTEGERs (not necessarily the same as mwSize)
      integer :: maxfun, nw
      integer :: nf_loc
      integer :: maxhist
      real(RP), allocatable :: fhist_loc(:), chist_loc(:)
C maxfun will also be an input to COBYLA, which requires
C maxfun to be INTEGERs (not necessarily the same as mwSize).
C nw is the size of workspace. nw = n*(3*n+2*m+11)+4*m+6 (see below).
C nw is the largest size of memory allocated by this code. If the value
C assigned to nw is so large that overflow occurs, then there will be
C a Segmentation Fault!!!
C The largest value of nw depends on the storage size of nw (~2E9 for
C kind 4, and ~9E18 for kind 8). This imposes an upper limit on the size
C of problem solvable by this code. If nw is INTEGER*4, assuming that
C m=10n, the largest value of n is ~9600. COBYLA is not designed for so
C large problems.
       procedure(FUNCON) :: calcfc

C Validate number of arguments.
      if (nrhs .ne. 11) then
         call mexErrMsgIdAndTxt ('fcobyla:nInput',
     &   'fcobyla: 11 input argument required.')
C subroutine mexErrMsgIdAndTxt(character*(*) errorid, character*(*) errormsg)
      end if
      if (nlhs .gt. 8) then
         call mexErrMsgIdAndTxt ('fcobyla:nOutput',
     &   'fcobyla: At most 8 outputs.')
      end if

C Validate inputs
C Input 1: fun (function handle)
      if (mxIsClass(prhs(1), 'function_handle') .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 1 should be a function handle.')
      end if
C Input 2: constr (function handle)
      if (mxIsClass(prhs(2), 'function_handle') .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 2 should be a function handle.')
      end if
C Input 3: x0 (double column)
      if (mxIsDouble(prhs(3)) .ne. 1 .or. mxGetM(prhs(3)) .lt. 1
     &    .or. mxGetN(prhs(3)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 3 should be a column vector of doubles.')
      end if
C Input 4: rhobeg (double scalar)
      if (mxIsDouble(prhs(4)) .ne. 1 .or. mxGetM(prhs(4)) .ne. 1
     &    .or. mxGetN(prhs(4)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 4 should be a double.')
      end if
C Input 5: rhobend (double scalar)
      if (mxIsDouble(prhs(5)) .ne. 1 .or. mxGetM(prhs(5)) .ne. 1
     &    .or. mxGetN(prhs(5)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 5 should be a double.')
      end if
C Input 6: maxfun (double scalar)
      if (mxIsDouble(prhs(6)) .ne. 1 .or. mxGetM(prhs(6)) .ne. 1
     &    .or. mxGetN(prhs(6)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     & 'fcobyla: Input 6 should be a double (with an integer value).')
      end if
C Input 7: m (double scalar)
      if (mxIsDouble(prhs(7)) .ne. 1 .or. mxGetM(prhs(7)) .ne. 1
     &    .or. mxGetN(prhs(7)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     & 'fcobyla: Input 7 should be a double (with an integer value).')
      end if
C Although inputs 6 and 7 (maxfun and m) are integers logically,
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

C Input 8: ftarget (double scalar)
      if (mxIsDouble(prhs(8)) .ne. 1 .or. mxGetM(prhs(8)) .ne. 1
     &    .or. mxGetN(prhs(8)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 8 should be a double.')
      end if
C Input 9: conval_x0 (double column, can be empty)
      if (mxIsDouble(prhs(9)) .ne. 1 .or.
     &(mxGetM(prhs(9)) .gt. 0 .and. mxGetN(prhs(9)) .gt. 1)) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 9 should be a column vector of doubles.')
      end if
C Input 10: f_x0 (double scalar)
      if (mxIsDouble(prhs(10)) .ne. 1 .or. mxGetM(prhs(10)) .ne. 1
     &    .or. mxGetN(prhs(10)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 10 should be a double.')
      end if
C Input 11: ctol (double scalar)
      if (mxIsDouble(prhs(11)) .ne. 1 .or. mxGetM(prhs(11)) .ne. 1
     &    .or. mxGetN(prhs(11)) .ne. 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Input 11 should be a double.')
      end if

C Read the inputs (there are 11)
      fun_ptr = prhs(1)
      constr_ptr = prhs(2)
      n = mxGetM(prhs(3)) ! This is why n should be of type mwSize
      n_int = int(n, kind(n_int))
C n_int is used when a variable of type INTEGER is needed
      if (n .ne. n_int) then
          call mexErrMsgIdAndTxt('fcobyla:IntError',
     &    'fcobyla: n does not equal n_int.')
      end if
      if (allocated(x)) deallocate (x)
      allocate (x(n_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory',
     &    'fcobyla: allocate(x) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(3)), x(1:n), n)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

      call mxCopyPtrToReal8(_MGETDB(prhs(4)), rhobeg, mwOne)
      call mxCopyPtrToReal8(_MGETDB(prhs(5)), rhoend, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
C Note the mwOne is of type mwSize; should not use literal constant 1
C NEVER use literal constants in Fortran mex.

C Check the values of rhobeg and rhoend. We do not check the values of
C other inputs (e.g., n, maxfun, npt) because the Fortran code does it
      if (rhobeg .le. zero .or. rhobeg .lt. rhoend .or.
     &rhoend .lt. zero) then
          call mexErrMsgIdAndTxt('fcobyla:InvalidRhobegRhoend',
     &'fcobyla: rhobeg and rhoend do not satisfy rhobeg >= rhobeg > 0.')
      end if

      call mxCopyPtrToReal8(_MGETDB(prhs(6)),
     &     maxfun_r, mwOne)
      maxfun = int(maxfun_r, kind(maxfun))
C maxfun will be an input to COBYLA, which requires maxfun to be
C an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(7)), m_r, mwOne)
      m = int(m_r, kind(m))
      m_int = int(m_r, kind(m_int))
      if (m .ne. m_int) then
          call mexErrMsgIdAndTxt('fcobyla:IntError',
     &    'fcobyla: m does not equal m_int.')
      end if
C m will be used in mxCopyPtrToReal8, requiring it to be of type mwSize
C m_int will be an input to COBYLA, which requires it to be
C an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(8)), ftarget, mwOne)
      call mxCopyPtrToReal8(_MGETDB(prhs(10)), f_x0, mwOne)
      call mxCopyPtrToReal8(_MGETDB(prhs(11)), ctol, mwOne)

      if (m .gt. 0 .and. m .ne. mxGetM(prhs(9))) then
C m is number of constraints
          call mexErrMsgIdAndTxt('fcobyla:WrongInput',
     &    'fcobyla: Length of input 9 should be m (input 7).')
      end if
      if (allocated(conval_x0)) deallocate (conval_x0)
      allocate (conval_x0(m_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory',
     &    'fcobyla: allocate(conval_x0) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(9)),
     &     conval_x0(1:m), m)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

C     Allocate workspace
      if (allocated(conval)) deallocate (conval)
      allocate (conval(m_int), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory',
     &    'fcobyla: allocate(conval) failed.')
      end if


C     Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory',
     &    'fcobyla: allocate(fhist) failed.')
      end if
      fhist = huge(0.0_DP)

      if (allocated(chist)) deallocate (chist)
      allocate (chist(maxfun), stat=allocate_status)
      if (allocate_status .ne. 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory',
     &    'fcobyla: allocate(chist) failed.')
      end if
      chist = huge(0.0_DP)

C     Call COBYLA
      iprint = 0
      ctol = epsilon(0.0_DP)
      maxhist = maxfun
!      call COBYLA(calcfc, n_int, m_int, x, rhobeg, rhoend, iprint,
!     & maxfun, f, info, ftarget, resmax, conval)
      call COBYLA(calcfc, x, f,
     & resmax, conval,
     & m_int, f_x0, conval_x0,
     & nf_loc, rhobeg, rhoend, ftarget, ctol, maxfun, iprint,
     & fhist = fhist_loc, chist = chist_loc, maxhist=maxhist,info=info)
      nf = nf_loc
      fhist(1:min(nf_loc, size(fhist), size(fhist_loc))) =
     & fhist_loc(1:min(nf_loc, size(fhist), size(fhist_loc)))
      chist(1:min(nf_loc, size(chist), size(chist_loc))) =
     & chist_loc(1:min(nf_loc, size(chist), size(chist_loc)))
C Note that n is of type mwSize, yet COBYLA expects input 1 to be
C of type INTEGER. Therefore, we should use n_int instead of n. Similar
C fo m/m_int.

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
C Set conval to the value of constr(x)
      plhs(6) = mxCreateDoubleMatrix(m, mwOne, notComplex)
      call mxCopyReal8ToPtr(conval(1:m), _MGETDB(plhs(6)), m)
      plhs(7) = mxCreateDoubleScalar(resmax)
      plhs(8) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(chist(1:nf),
     &     _MGETDB(plhs(8)), nf)

C     Free memory
      deallocate (x)
      deallocate (conval_x0)
      deallocate (fhist)
      deallocate (chist)
      deallocate (conval)

      return
      end subroutine mexFunction

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C The Fortran subroutine that evaluates the objective function
      subroutine calcfc(x, funval, conval)
      use fcobyla
      implicit none

C dummy variables
      integer :: n, m
C The types of n and m are derived from the code of COBYLA. Thus n and m
C should be of type INTEGER instead of mwSize
      real(rp), intent(in) :: x(:)
      real(rp), intent(out) :: funval, conval(:)

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
      mwSize :: n_mw, m_mw
C n_mw is the mwSize cast of n: n_mw = int(n, kind(n_mw));
C used when a variable of type mwSize is needed
      integer :: k
C k is the index for the constraints; since m is of type INTEGER,
C k is of the same type; see the (unique) do lopp below
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      integer(int32), parameter :: intOne=1, intTwo=2
      character(5), parameter :: funFeval = 'feval'
C intOne, intTwo, and funFeval are used when calling mexCallMATLAB
      real(rp) :: resmax ! constraint violation

      n = size(x)
      m = size(conval)
C Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      m_mw = int(m, kind(m_mw)) ! cast m to type mwSize
      if (n .ne. n_mw) then
          call mexErrMsgIdAndTxt('fcobyla:IntError',
     &    'fcobyla: n does not equal n_mw.')
      end if
      if (m .ne. m_mw) then
          call mexErrMsgIdAndTxt('fcobyla:IntError',
     &    'fcobyla: m does not equal m_mw.')
      end if

C Evaluate the objective function (fun_ptr) at x
      plhs(1) = mxCreateDoubleScalar(huge(0.0_DP))
C Output of f_value = feval(fun, x); see below
      prhs(1) = fun_ptr
C First input of f_value = feval(fun, x); see below; fun_ptr is a global variable
      prhs(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
C Second input of f_value = feval(fun, x); see below
      call mxCopyReal8ToPtr(x(1:n), _MGETDB(prhs(2)), n_mw)
C subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)
      if (0 .ne. mexCallMATLAB(intOne, plhs,
     &   intTwo, prhs, funFeval)) then
C Execute matlab command: f_value = feval(fun, x)
C integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
          call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall',
     &    'fcobyla: mex fails to call fun.')
      end if

      if (plhs(1) .eq. 0 .or. _MGETDB(plhs(1)) .eq. 0) then
          call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall',
     &    'fcobyla: NULL pointer returned when mex calls fun.')
      end if

      if (mxGetM(plhs(1))*mxGetN(plhs(1)) .ne. 1
     &    .or. mxIsDouble(plhs(1)) .ne. 1) then
           call mexErrMsgIdAndTxt('fcobyla:ObjectiveNotScalar',
     & 'fcobyla: The objective function should return a scalar value.')
      end if

      call mxCopyPtrToReal8(_MGETDB(plhs(1)), funval, mwOne)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

C Use extreme barrier to cope with 'hidden constraints'
      if (funval .gt. hugefun .or. is_nan(funval)) then
          funval = hugefun ! hugefun is defined in consts
      end if
      !funval = max(-hugefun, funval)

C Free memory; note that plhs and prhs are just temporary variables in
C this subroutine. We are NOT in mexFunction!
      call mxDestroyArray(plhs(1))
C Not yet to free prhs(2), which will be used when evaluating the constraint

C Evaluate the constraint (constr_ptr) at x
      if (nf .eq. 0) then
C The very first iteration needs constr(x0), which was already evaluated in
C the matlab code (to get the value of m) and saved in fcobyla.mod as
C conval_x0. Copy the value directly without calling constr.
          conval(1:m) = conval_x0(1:m)
      else
          plhs(1) = mxCreateDoubleMatrix(m_mw, mwOne, notComplex)
C Output of c_value = feval(constr, x); see below
          prhs(1) = constr_ptr
C First input of c_value = feval(constr, x); see below; constr_ptr is a global variable
C prhs(2) was already set to x when evaluating fun
          if (0 .ne. mexCallMATLAB(intOne, plhs, intTwo,
     &        prhs, funFeval)) then
C Execute matlab command: c_value = feval(constr, x)
              call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall',
     &        'fcobyla: mex fails to call constr.')
          end if
          if (plhs(1) .eq. 0 .or. (m .gt. 0 .and.
     &        _MGETDB(plhs(1)) .eq. 0)) then
              call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall',
     &        'fcobyla: NULL pointer returned when mex calls constr.')
          end if
          if (m .gt. 0 .and. (mxGetM(plhs(1)) .ne. m .or.
     &        mxGetN(plhs(1)) .ne. 1 .or.
     &        mxIsDouble(plhs(1)) .ne. 1)) then
           call mexErrMsgIdAndTxt('fcobyla:ConstrNotScalarVector',
     &'fcobyla: The constraint function should return a scalar vector of
     &size mx1.')
          end if
          call mxCopyPtrToReal8(_MGETDB(plhs(1)),
     &         conval(1:m), m_mw)
C subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
      end if
C Calculate the constraint violation (named 'RESMAX' in Powell's COBYLA code)
      resmax = zero ! zero is defined in module 'consts'.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Code without extreme barrier:
C Constraint: constr(x) >= 0
C      do k = 1, m
C          if (conval(k) .ne. conval(k)) then
C              resmax = conval(k) ! Set resmax=NaN if conval contains NaN
C              exit
C          else
C              resmax = max(resmax, -conval(k))
C          end if
C      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Use extreme barrier to cope with 'hidden constraints'
C Constraint: constr(x) >= 0
      do k = 1, m
          if (conval(k) .lt. -hugecon .or. is_nan(conval(k))) then
              conval(k) = -hugecon ! hugecon is defined in consts
          end if

C This part is NOT extreme barrier. We replace extremely negative values
C of cineq (which leads to no constraint violation) by -hugecon. Otherwise,
C NaN or Inf may occur in the interpolation models.
          if (conval(k) .gt. hugecon) then
              conval(k) = hugecon ! hugecon is defined in consts
          end if

          resmax = max(resmax, -conval(k))
      end do
C
C Free memory; note that plhs and prhs are just temporary variables in
C this subroutine. We are NOT in mexFunction!
      if (nf .gt. 0) call mxDestroyArray(plhs(1)) ! Only if nf >= 1
      call mxDestroyArray(prhs(2))

C Update global variables
      nf = nf + int(1, kind(nf))
C Some compiler (e.g., g95) may complain about implicit conversion if
C written as nf = nf+1
      fhist(nf) = funval
      chist(nf) = resmax

      return
      end subroutine calcfc
