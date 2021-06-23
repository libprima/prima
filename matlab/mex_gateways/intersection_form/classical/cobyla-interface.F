!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of cobyla-interface.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 23-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! The mex gateway for COBYLA
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
! substitution), some lines may become too long and hence get truncated

#include "fintrf.h"

! Use mxGetDoubles instead of mxGetPr if possible, the former being
! available since MATLAB R2018b. The following macros should be put
! after the header "fintrf.h" is included, because mxGetDoubles is
! defined in it. _MGETDB gets the pointer pointing to a real array.
#if defined _MGETDB
#undef _MGETDB
#endif
#if defined mxGetDoubles
#define _MGETDB mxGetDoubles
#else
#define _MGETDB mxGetPr
#endif

! Fortran code starts

      module fcobyla ! Some global variables

      use consts_mod, only : DP, rp, zero
      use infnan_mod
      implicit none

      integer*4, parameter :: int32_0 = 0
      integer, parameter :: INT32 = kind(int32_0)
      integer(INT32), parameter :: notComplex = 0
! notComplex is used in mxCreateDoubleMatrix, whose signature is
! mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwSize, parameter :: mwOne = 1 ! Integer 1 with type mwSize
! mwOne is used in mxCreateDoubleMatrix and mxCopyPtrToReal8,
! whose signatures are
! mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

! Some global variables
      mwPointer, save :: fun_ptr, con_ptr
! Objective function and constraints (pointers)
      mwSize, save :: nf ! Number of function evaluations
! nf should be of type mwSize, because it will be later used as an
! argument of mxCreateDoubleMatrix
      real(rp), allocatable, save :: fhist(:)
! History of function evaluation
      real(rp), allocatable, save :: chist(:)
! History of constraint violation
      real(rp), allocatable, save :: conval_x0(:)
! Constraint value at x0
      end module fcobyla

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
! subroutine mexFunction(integer nlhs, mwPointer plhs, integer nrhs, mwPointer prhs)
! Entry point to Fortran MEX function
! If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name
! extension depends on the platform), then the following function is
! callable in matlab:
! [xopt, fopt, info, nf, fhist, convalopt, constrviolation, chist] = FUNCTION_NAME(fun, con, x0, rhobeg, rhoend, maxfun, m, ftarget, conval_x0)

      use fcobyla
      implicit none

! mexFunction arguments (dummy variables):
      integer, intent(in) :: nlhs, nrhs
! nlhs and nrhs are of type INTEGER according to the documents of MATLAB R2019a
      mwPointer, intent(in) :: prhs(nrhs)
      mwPointer, intent(inout) :: plhs(nlhs)

! function declarations:
      mwPointer, external :: mxCreateDoubleMatrix
! mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwPointer, external :: mxCreateDoubleScalar
! mwPointer mxCreateDoubleScalar(real*8 value)
      mwPointer, external :: _MGETDB
! mwPointer _MGETDB(mwPointer pm);
! can be replaced by mxGetDoubles since MATLAB R2018b
      mwSize, external :: mxGetM, mxGetN
! mwPointer mxGetM(mwPointer pm), mxGetN(mwPointer pm)
      integer(INT32), external :: mxIsDouble, mxIsClass
! integer*4 mxIsDouble(mwPointer pm);
! integer*4 mxIsClass(mwPointer pm, character*(*) classname)

! workspace and intermediate variables
      mwSize :: n, m
      integer :: n_int, m_int
! n is the problem dimension; it is of type mwSize because it is
! obtained by mxGetM. m is the number of constraints. It is of type
! mwSize because it will be used in mxCopyPtrToReal8.
! n_int and m_int are the INTEGER cast of n and m, i.e., n_int = int(n),
! m_int = int(m), which will be used when variables of type INTEGER are
! needed (e.g., when calling COBYLA).
      real(rp), allocatable :: x(:), conval(:), w(:)
      integer, allocatable :: iact(:)
      real(rp) :: f, rhobeg, rhoend, maxfun_r, m_r, ftarget, resmax
      integer :: info, iprint, allocate_status
! info and iprint are inputs to COBYLA, which requires them to be
! INTEGERs (not necessarily the same as mwSize)
      integer :: maxfun, nw
! maxfun will also be an input to COBYLA, which requires
! maxfun to be INTEGERs (not necessarily the same as mwSize).
! nw is the size of workspace. nw = n*(3*n+2*m+11)+4*m+6 (see below).
! nw is the largest size of memory allocated by this code. If the value
! assigned to nw is so large that overflow occurs, then there will be
! a Segmentation Fault!!!
! The largest value of nw depends on the storage size of nw (~2E9 for
! kind 4, and ~9E18 for kind 8). This imposes an upper limit on the size
! of problem solvable by this code. If nw is INTEGER*4, assuming that
! m=10n, the largest value of n is ~9600. COBYLA is not designed for so
! large problems.

! Validate number of arguments.
      if (nrhs /= 9) then
          call mexErrMsgIdAndTxt('fcobyla:nInput', 'fcobyla: 9 input arg&
     &ument required.')
! subroutine mexErrMsgIdAndTxt(character*(*) errorid, character*(*) errormsg)
      end if
      if (nlhs > 8) then
          call mexErrMsgIdAndTxt('fcobyla:nOutput', 'fcobyla: At most 8 &
     &outputs.')
      end if

! Validate inputs
! Input 1: fun (function handle)
      if (mxIsClass(prhs(1), 'function_handle') /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 1&
     & should be a function handle.')
      end if
! Input 2: con (function handle)
      if (mxIsClass(prhs(2), 'function_handle') /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 2&
     & should be a function handle.')
      end if
! Input 3: x0 (double column)
      if (mxIsDouble(prhs(3)) /= 1 .or. mxGetM(prhs(3)) < 1 .or. mxGetN(&
     &prhs(3)) /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 3&
     & should be a column vector of doubles.')
      end if
! Input 4: rhobeg (double scalar)
      if (mxIsDouble(prhs(4)) /= 1 .or. mxGetM(prhs(4)) /= 1 .or. mxGetN&
     &(prhs(4)) /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 4&
     & should be a double.')
      end if
! Input 5: rhobend (double scalar)
      if (mxIsDouble(prhs(5)) /= 1 .or. mxGetM(prhs(5)) /= 1 .or. mxGetN&
     &(prhs(5)) /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 5&
     & should be a double.')
      end if
! Input 6: maxfun (double scalar)
      if (mxIsDouble(prhs(6)) /= 1 .or. mxGetM(prhs(6)) /= 1 .or. mxGetN&
     &(prhs(6)) /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 6&
     & should be a double (with an integer value).')
      end if
! Input 7: m (double scalar)
      if (mxIsDouble(prhs(7)) /= 1 .or. mxGetM(prhs(7)) /= 1 .or. mxGetN&
     &(prhs(7)) /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 7&
     & should be a double (with an integer value).')
      end if
! Although inputs 6 and 7 (maxfun and m) are integers logically,
! they have to be passed to the mexified code as double variables. In
! mex, data is passed by pointers, but there are only very limited
! functions that can read an integer value from a pointer or write
! an interger value to a pointer (mxCopyPtrToInteger1,
! mxCopyInteger1ToPtr, mxCopyPtrToInteger2, mxCopyInteger2ToPtr,
! mxCopyPtrToInteger4, mxCopyInteger4ToPtr; no function for
! INTEGER*8). This makes it impossible to pass integer data properly
! unless we know the kind of the integer. Therefore, in general, it
! is recommended to pass integers as double variables and then cast
! them back to integers when needed. Indeed, in MATLAB, even if we
! define maxfun = 1000, the class of maxfun is double! To get an
! integer maxfun, we would have to define maxfun = int32(1000) or
! maxfun = int64(1000)!

! Input 8: ftarget (double scalar)
      if (mxIsDouble(prhs(8)) /= 1 .or. mxGetM(prhs(8)) /= 1 .or. mxGetN&
     &(prhs(8)) /= 1) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 8&
     & should be a double.')
      end if
! Input 9: conval_x0 (double column, can be empty)
      if (mxIsDouble(prhs(9)) /= 1 .or. (mxGetM(prhs(9)) > 0 .and. mxGet&
     &N(prhs(9)) > 1)) then
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 9&
     & should be a column vector of doubles.')
      end if

! Read the inputs (there are 9)
      fun_ptr = prhs(1)
      con_ptr = prhs(2)
      n = mxGetM(prhs(3)) ! This is why n should be of type mwSize
      n_int = int(n, kind(n_int))
! n_int is used when a variable of type INTEGER is needed
      if (n /= n_int) then
          call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: n does no&
     &t equal n_int.')
      end if
      if (allocated(x)) deallocate (x)
      allocate (x(n_int), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(x) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(3)), x(1:n), n)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

      call mxCopyPtrToReal8(_MGETDB(prhs(4)), rhobeg, mwOne)
      call mxCopyPtrToReal8(_MGETDB(prhs(5)), rhoend, mwOne)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
! Note the mwOne is of type mwSize; should not use literal constant 1
! NEVER use literal constants in Fortran mex.

! Check the values of rhobeg and rhoend. We do not check the values of
! other inputs (e.g., n, maxfun, npt) because the Fortran code does it
      if (rhobeg <= zero .or. rhobeg < rhoend .or. rhoend < zero) then
          call mexErrMsgIdAndTxt('fcobyla:InvalidRhobegRhoend', 'fcobyla&
     &: rhobeg and rhoend do not satisfy rhobeg >= rhobeg > 0.')
      end if

      call mxCopyPtrToReal8(_MGETDB(prhs(6)), maxfun_r, mwOne)
      maxfun = int(maxfun_r, kind(maxfun))
! maxfun will be an input to COBYLA, which requires maxfun to be
! an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(7)), m_r, mwOne)
      m = int(m_r, kind(m))
      m_int = int(m_r, kind(m_int))
      if (m /= m_int) then
          call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: m does no&
     &t equal m_int.')
      end if
! m will be used in mxCopyPtrToReal8, requiring it to be of type mwSize
! m_int will be an input to COBYLA, which requires it to be
! an INTEGER (not necessary the same as mwSize)
      call mxCopyPtrToReal8(_MGETDB(prhs(8)), ftarget, mwOne)

      if (m /= mxGetM(prhs(9))) then ! m is number of constraints
          call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Length &
     &of input 9 should be m (input 7).')
      end if
      if (allocated(conval_x0)) deallocate (conval_x0)
      allocate (conval_x0(m_int), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(conval_x0) failed.')
      end if
      call mxCopyPtrToReal8(_MGETDB(prhs(9)), conval_x0(1:m), m)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

!     Allocate workspace
      if (allocated(conval)) deallocate (conval)
      allocate (conval(m_int), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(conval) failed.')
      end if

      if (allocated(w)) deallocate (w)
      nw = n_int * (3 * n_int + 2 * m_int + 11) + 4 * m_int + 6
      allocate (w(nw), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(w) failed.')
      end if

      if (allocated(iact)) deallocate (iact)
      nw = m_int + 1
      allocate (iact(nw), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(iact) failed.')
      end if

!     Initialize global variables
      nf = 0
      if (allocated(fhist)) deallocate (fhist)
      allocate (fhist(maxfun), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(fhist) failed.')
      end if
      fhist = huge(0.0_DP)

      if (allocated(chist)) deallocate (chist)
      allocate (chist(maxfun), stat=allocate_status)
      if (allocate_status /= 0) then
          call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla:&
     & allocate(chist) failed.')
      end if
      chist = huge(0.0_DP)

!     Call COBYLA
      iprint = 0
      call COBYLA(n_int, m_int, x, rhobeg, rhoend, iprint, maxfun, w, ia&
     &ct, f, info, ftarget, resmax, conval)
! Note that n is of type mwSize, yet COBYLA expects input 1 to be
! of type INTEGER. Therefore, we should use n_int instead of n. Similar
! fo m/m_int.

! Write outputs
      plhs(1) = mxCreateDoubleMatrix(n, mwOne, notComplex)
      call mxCopyReal8ToPtr(x(1:n), _MGETDB(plhs(1)), n)
      plhs(2) = mxCreateDoubleScalar(f)
! Although info and nf are integers logically, they are passed as double
      plhs(3) = mxCreateDoubleScalar(real(info, rp))
      plhs(4) = mxCreateDoubleScalar(real(nf, rp))
      plhs(5) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(fhist(1:nf), _MGETDB(plhs(5)), nf)
! Set conval to the value of con(x)
      plhs(6) = mxCreateDoubleMatrix(m, mwOne, notComplex)
      call mxCopyReal8ToPtr(conval(1:m), _MGETDB(plhs(6)), m)
      plhs(7) = mxCreateDoubleScalar(resmax)
      plhs(8) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
      call mxCopyReal8ToPtr(chist(1:nf), _MGETDB(plhs(8)), nf)

!     Free memory
      deallocate (x)
      deallocate (conval_x0)
      deallocate (w)
      deallocate (iact)
      deallocate (fhist)
      deallocate (chist)
      deallocate (conval)

      return
      end subroutine mexFunction

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! The Fortran subroutine that evaluates the objective function
      subroutine calcfc(n, m, x, funval, conval)
      use fcobyla
      implicit none

! dummy variables
      integer, intent(in) :: n, m
! The types of n and m are derived from the code of COBYLA. Thus n and m
! should be of type INTEGER instead of mwSize
      real(rp), intent(in) :: x(n)
      real(rp), intent(out) :: funval, conval(m)

! function declarations
      integer(INT32), external :: mexCallMATLAB
! integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
      mwPointer, external :: mxCreateDoubleMatrix
! mwPointer mxCreateDoubleMatrix(mwSize m, mwSize n, integer*4 ComplexFlag)
      mwPointer, external :: mxCreateDoubleScalar
! mwPointer mxCreateDoubleScalar(real*8 value)
      mwPointer, external :: _MGETDB
! mwPointer _MGETDB(mwPointer pm)

! intermediate variables
      mwSize :: n_mw, m_mw
! n_mw is the mwSize cast of n: n_mw = int(n, kind(n_mw));
! used when a variable of type mwSize is needed
      integer :: k
! k is the index for the constraints; since m is of type INTEGER,
! k is of the same type; see the (unique) do lopp below
      mwPointer :: plhs(1), prhs(2) ! used in mexCallMATLAB
      integer(INT32), parameter :: intOne = 1, intTwo = 2
      character(5), parameter :: funFeval = 'feval'
! intOne, intTwo, and funFeval are used when calling mexCallMATLAB
      real(rp) :: resmax ! constraint violation

! Start the real business
      n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
      m_mw = int(m, kind(m_mw)) ! cast m to type mwSize
      if (n /= n_mw) then
          call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: n does no&
     &t equal n_mw.')
      end if
      if (m /= m_mw) then
          call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: m does no&
     &t equal m_mw.')
      end if

! Evaluate the objective function (fun_ptr) at x
      plhs(1) = mxCreateDoubleScalar(huge(0.0_DP))
! Output of f_value = feval(fun, x); see below
      prhs(1) = fun_ptr
! First input of f_value = feval(fun, x); see below; fun_ptr is a global variable
      prhs(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
! Second input of f_value = feval(fun, x); see below
      call mxCopyReal8ToPtr(x(1:n), _MGETDB(prhs(2)), n_mw)
! subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)
      if (0 /= mexCallMATLAB(intOne, plhs, intTwo, prhs, funFeval)) then
! Execute matlab command: f_value = feval(fun, x)
! integer*4 mexCallMATLAB(integer*4 nlhs, mwPointer plhs, integer*4 nrhs, mwPointer prhs, character*(*) functionName)
          call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyla: m&
     &ex fails to call fun.')
      end if
      if (plhs(1) == 0 .or. _MGETDB(plhs(1)) == 0) then
          call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyla: N&
     &ULL pointer returned when mex calls fun.')
      else
          call mxCopyPtrToReal8(_MGETDB(plhs(1)), funval, mwOne)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
      end if

! Free memory; note that plhs and prhs are just temporary variables in
! this subroutine. We are NOT in mexFunction!
      call mxDestroyArray(plhs(1))
! Not yet to free prhs(2), which will be used when evaluating the constraint

! Evaluate the constraint (con_ptr) at x
      if (nf == 0) then
! The very first iteration needs con(x0), which was already evaluated in
! the matlab code (to get the value of m) and saved in fcobyla.mod as
! conval_x0. Copy the value directly without calling con.
          conval(1:m) = conval_x0(1:m)
      else
          plhs(1) = mxCreateDoubleMatrix(m_mw, mwOne, notComplex)
! Output of c_value = feval(con, x); see below
          prhs(1) = con_ptr
! First input of c_value = feval(con, x); see below; con_ptr is a global variable
! prhs(2) was already set to x when evaluating fun
          if (0 /= mexCallMATLAB(intOne, plhs, intTwo, prhs, funFeval)) &
     &then
! Execute matlab command: c_value = feval(con, x)
              call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyl&
     &a: mex fails to call con.')
          end if
          if (plhs(1) == 0 .or. (m > 0 .and. _MGETDB(plhs(1)) == 0)) the&
     &n
              call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyl&
     &a: NULL pointer returned when mex calls con.')
          else
              call mxCopyPtrToReal8(_MGETDB(plhs(1)), conval(1:m), m_mw)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
          end if
      end if
! Calculate the constraint violation (named 'RESMAX' in Powell's COBYLA code)
      resmax = zero ! zero is defined in module 'consts'.
! Constraint: con(x) >= 0
      do k = 1, m
          if (is_nan(conval(k))) then
              resmax = conval(k) ! Set resmax=NaN if conval contains NaN
              exit
          else
              resmax = max(resmax, -conval(k))
          end if
      end do
!
! Free memory; note that plhs and prhs are just temporary variables in
! this subroutine. We are NOT in mexFunction!
      if (nf > 0) call mxDestroyArray(plhs(1)) ! Only if nf >= 1
      call mxDestroyArray(prhs(2))

! Update global variables
      nf = nf + int(1, kind(nf))
! Some compiler (e.g., g95) may complain about implicit conversion if
! written as nf = nf+1
      fhist(nf) = funval
      chist(nf) = resmax

      return
      end subroutine calcfc