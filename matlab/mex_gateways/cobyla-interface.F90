! The mex gateway for COBYLA
!
! Coded by Zaikun Zhang in June 2021.
!
! Last Modified: Friday, June 25, 2021 PM07:39:02


#include "fintrf.h"

! Entry point to Fortran MEX function
subroutine mexFunction(nargout, poutput, nargin, pinput)
! If the binary MEX file is named as FUNCTION_NAME.mex*** (file-name extension depends on the
! platform), then the following function is callable in matlab:
! [xopt, fopt, info, nf, fhist, convalopt, constrviolation, chist] = FUNCTION_NAME(fun, con, x0, rhobeg, rhoend, maxfun, m, ftarget, conval_x0)


! Generic modules
use consts_mod, only : RP, IK
use fmxapi_mod, only : fmxVerifyNArgin, fmxVerifyNArgout
use fmxapi_mod, only : fmxVerifyClassShape
use fmxapi_mod, only : fmxAllocate
use fmxapi_mod, only : fmxReadMPtr, fmxWriteMPtr

! Solver-specific module
use cobyla_mod, only : cobyla 

implicit none

! mexFunction arguments nargout and nargin are of type INTEGER in MATLAB 2019a documents.
integer, intent(in) :: nargout, nargin
mwPointer, intent(in) :: pinput(nargin)
mwPointer, intent(out) :: poutput(nargout)

! Intermediate variables
integer(IK) :: info
integer(IK) :: iprint
integer(IK) :: maxfun
integer(IK) :: nf
real(RP) :: f
real(RP) :: ftarget
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

! Validate number of arguments
call fmxVerifyNArgin(nargin, 9)
call fmxVerifyNArgout(nargout, 8)

! Read inputs (there are 14)
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
call fmxReadMPtr(pinput(13), maxhist)
call fmxReadMPtr(pinput(14), output_xhist)

! Validate inputs
! Input 1: fun (function handle)
if (mxIsClass(pinput(1), 'function_handle') /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 1 should be a function handle.')
end if
! Input 2: con (function handle)
if (mxIsClass(pinput(2), 'function_handle') /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 2 should be a function handle.')
end if
! Input 3: x0 (double column)
if (mxIsDouble(pinput(3)) /= 1 .or. mxGetM(pinput(3)) < 1 .or. mxGetN(pinput(3)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 3 should be a column vector of doubles.')
end if
! Input 4: rhobeg (double scalar)
if (mxIsDouble(pinput(4)) /= 1 .or. mxGetM(pinput(4)) /= 1 .or. mxGetN(pinput(4)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 4 should be a double.')
end if
! Input 5: rhobend (double scalar)
if (mxIsDouble(pinput(5)) /= 1 .or. mxGetM(pinput(5)) /= 1 .or. mxGetN(pinput(5)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 5 should be a double.')
end if
! Input 6: maxfun (double scalar)
if (mxIsDouble(pinput(6)) /= 1 .or. mxGetM(pinput(6)) /= 1 .or. mxGetN(pinput(6)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 6 should be a double (with an integer value).')
end if
! Input 7: m (double scalar)
if (mxIsDouble(pinput(7)) /= 1 .or. mxGetM(pinput(7)) /= 1 .or. mxGetN(pinput(7)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 7 should be a double (with an integer value).')
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
if (mxIsDouble(pinput(8)) /= 1 .or. mxGetM(pinput(8)) /= 1 .or. mxGetN(pinput(8)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 8 should be a double.')
end if
! Input 9: conval_x0 (double column, can be empty)
if (mxIsDouble(pinput(9)) /= 1 .or. (mxGetM(pinput(9)) > 0 .and. mxGetN(pinput(9)) > 1)) then
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Input 9 should be a column vector of doubles.')
end if

! Read the inputs (there are 9)
fun_ptr = pinput(1)
con_ptr = pinput(2)
n = mxGetM(pinput(3)) ! This is why n should be of type mwSize
n_int = int(n, kind(n_int))
! n_int is used when a variable of type INTEGER is needed
if (n /= n_int) then
    call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: n does not equal n_int.')
end if
if (allocated(x)) deallocate (x)
allocate (x(n_int), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(x) failed.')
end if
call mxCopyPtrToReal8(_MGETDB(pinput(3)), x(1:n), n)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

call mxCopyPtrToReal8(_MGETDB(pinput(4)), rhobeg, mwOne)
call mxCopyPtrToReal8(_MGETDB(pinput(5)), rhoend, mwOne)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
! Note the mwOne is of type mwSize; should not use literal constant 1
! NEVER use literal constants in Fortran mex.

! Check the values of rhobeg and rhoend. We do not check the values of
! other inputs (e.g., n, maxfun, npt) because the Fortran code does it
if (rhobeg <= zero .or. rhobeg < rhoend .or. rhoend < zero) then
    call mexErrMsgIdAndTxt('fcobyla:InvalidRhobegRhoend', 'fcobyla: rhobeg and rhoend do not satisfy rhobeg >= rhobeg > 0.')
end if

call mxCopyPtrToReal8(_MGETDB(pinput(6)), maxfun_r, mwOne)
maxfun = int(maxfun_r, kind(maxfun))
! maxfun will be an input to COBYLA, which requires maxfun to be
! an INTEGER (not necessary the same as mwSize)
call mxCopyPtrToReal8(_MGETDB(pinput(7)), m_r, mwOne)
m = int(m_r, kind(m))
m_int = int(m_r, kind(m_int))
if (m /= m_int) then
    call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: m does not equal m_int.')
end if
! m will be used in mxCopyPtrToReal8, requiring it to be of type mwSize
! m_int will be an input to COBYLA, which requires it to be
! an INTEGER (not necessary the same as mwSize)
call mxCopyPtrToReal8(_MGETDB(pinput(8)), ftarget, mwOne)

if (m > 0 .and. m /= mxGetM(pinput(9))) then
! m is number of constraints
    call mexErrMsgIdAndTxt('fcobyla:WrongInput', 'fcobyla: Length of input 9 should be m (input 7).')
end if
if (allocated(conval_x0)) deallocate (conval_x0)
allocate (conval_x0(m_int), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(conval_x0) failed.')
end if
call mxCopyPtrToReal8(_MGETDB(pinput(9)), conval_x0(1:m), m)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

!     Allocate workspace
if (allocated(conval)) deallocate (conval)
allocate (conval(m_int), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(conval) failed.')
end if

if (allocated(w)) deallocate (w)
nw = n_int * (3 * n_int + 2 * m_int + 11) + 4 * m_int + 6
allocate (w(nw), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(w) failed.')
end if

if (allocated(iact)) deallocate (iact)
nw = m_int + 1
allocate (iact(nw), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(iact) failed.')
end if

!     Initialize global variables
nf = 0
if (allocated(fhist)) deallocate (fhist)
allocate (fhist(maxfun), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(fhist) failed.')
end if
fhist = huge(0.0_DP)

if (allocated(chist)) deallocate (chist)
allocate (chist(maxfun), stat=allocate_status)
if (allocate_status /= 0) then
    call mexErrMsgIdAndTxt('fcobyla:InsufficientMemory', 'fcobyla: allocate(chist) failed.')
end if
chist = huge(0.0_DP)

!     Call COBYLA
iprint = 0
call COBYLA(n_int, m_int, x, rhobeg, rhoend, iprint, maxfun, iact, f, info, ftarget, resmax, conval)
! Note that n is of type mwSize, yet COBYLA expects input 1 to be
! of type INTEGER. Therefore, we should use n_int instead of n. Similar
! fo m/m_int.

! Write outputs
poutput(1) = mxCreateDoubleMatrix(n, mwOne, notComplex)
call mxCopyReal8ToPtr(x(1:n), _MGETDB(poutput(1)), n)
poutput(2) = mxCreateDoubleScalar(f)
! Although info and nf are integers logically, they are passed as double
poutput(3) = mxCreateDoubleScalar(real(info, rp))
poutput(4) = mxCreateDoubleScalar(real(nf, rp))
poutput(5) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
call mxCopyReal8ToPtr(fhist(1:nf), _MGETDB(poutput(5)), nf)
! Set conval to the value of con(x)
poutput(6) = mxCreateDoubleMatrix(m, mwOne, notComplex)
call mxCopyReal8ToPtr(conval(1:m), _MGETDB(poutput(6)), m)
poutput(7) = mxCreateDoubleScalar(resmax)
poutput(8) = mxCreateDoubleMatrix(mwOne, nf, notComplex)
call mxCopyReal8ToPtr(chist(1:nf), _MGETDB(poutput(8)), nf)

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
integer(INT32), external :: mexCallMATLAB, mxIsDouble
! integer*4 mexCallMATLAB(integer*4 nargout, mwPointer poutput, integer*4 nargin, mwPointer pinput, character*(*) functionName)
! integer*4 mxIsDouble(mwPointer pm);
mwSize, external :: mxGetM, mxGetN
! mwPointer mxGetM(mwPointer pm), mxGetN(mwPointer pm)
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
mwPointer :: poutput(1), pinput(2) ! used in mexCallMATLAB
integer(INT32), parameter :: intOne = 1, intTwo = 2
character(5), parameter :: funFeval = 'feval'
! intOne, intTwo, and funFeval are used when calling mexCallMATLAB
real(rp) :: resmax ! constraint violation

! Start the real business
n_mw = int(n, kind(n_mw)) ! cast n to type mwSize
m_mw = int(m, kind(m_mw)) ! cast m to type mwSize
if (n /= n_mw) then
    call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: n does not equal n_mw.')
end if
if (m /= m_mw) then
    call mexErrMsgIdAndTxt('fcobyla:IntError', 'fcobyla: m does not equal m_mw.')
end if

! Evaluate the objective function (fun_ptr) at x
poutput(1) = mxCreateDoubleScalar(huge(0.0_DP))
! Output of f_value = feval(fun, x); see below
pinput(1) = fun_ptr
! First input of f_value = feval(fun, x); see below; fun_ptr is a global variable
pinput(2) = mxCreateDoubleMatrix(n_mw, mwOne, notComplex)
! Second input of f_value = feval(fun, x); see below
call mxCopyReal8ToPtr(x(1:n), _MGETDB(pinput(2)), n_mw)
! subroutine mxCopyReal8ToPtr(real*8 y(n), mwPointer px, mwSize n)
if (0 /= mexCallMATLAB(intOne, poutput, intTwo, pinput, funFeval)) then
! Execute matlab command: f_value = feval(fun, x)
! integer*4 mexCallMATLAB(integer*4 nargout, mwPointer poutput, integer*4 nargin, mwPointer pinput, character*(*) functionName)
    call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyla: mex fails to call fun.')
end if

if (poutput(1) == 0 .or. _MGETDB(poutput(1)) == 0) then
    call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyla: NULL pointer returned when mex calls fun.')
end if

if (mxGetM(poutput(1)) * mxGetN(poutput(1)) /= 1 .or. mxIsDouble(poutput(1)) /= 1) then
    call mexErrMsgIdAndTxt('fcobyla:ObjectiveNotScalar', 'fcobyla: The objective function should return a scalar value.')
end if

call mxCopyPtrToReal8(_MGETDB(poutput(1)), funval, mwOne)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)

! Use extreme barrier to cope with 'hidden constraints'
if (funval > hugefun .or. is_nan(funval)) then
    funval = hugefun ! hugefun is defined in consts
end if

! Free memory; note that poutput and pinput are just temporary variables in
! this subroutine. We are NOT in mexFunction!
call mxDestroyArray(poutput(1))
! Not yet to free pinput(2), which will be used when evaluating the constraint

! Evaluate the constraint (con_ptr) at x
if (nf == 0) then
! The very first iteration needs con(x0), which was already evaluated in
! the matlab code (to get the value of m) and saved in fcobyla.mod as
! conval_x0. Copy the value directly without calling con.
    conval(1:m) = conval_x0(1:m)
else
    poutput(1) = mxCreateDoubleMatrix(m_mw, mwOne, notComplex)
! Output of c_value = feval(con, x); see below
    pinput(1) = con_ptr
! First input of c_value = feval(con, x); see below; con_ptr is a global variable
! pinput(2) was already set to x when evaluating fun
    if (0 /= mexCallMATLAB(intOne, poutput, intTwo, pinput, funFeval)) then
! Execute matlab command: c_value = feval(con, x)
        call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyla: mex fails to call con.')
    end if
    if (poutput(1) == 0 .or. (m > 0 .and. _MGETDB(poutput(1)) == 0)) then
        call mexErrMsgIdAndTxt('fcobyla:UnsuccessfulCall', 'fcobyla: NULL pointer returned when mex calls con.')
    end if
    if (m > 0 .and. (mxGetM(poutput(1)) /= m .or. mxGetN(poutput(1)) /= 1 .or. mxIsDouble(poutput(1)) /= 1)) then
        call mexErrMsgIdAndTxt('fcobyla:ConstrNotScalarVector', 'fcobyla: The constraint function should return a scalar vector of &
            &size mx1.')
    end if
    call mxCopyPtrToReal8(_MGETDB(poutput(1)), conval(1:m), m_mw)
! subroutine mxCopyPtrToReal8(mwPointer px, real*8 y(n), mwSize n)
end if
! Calculate the constraint violation (named 'RESMAX' in Powell's COBYLA code)
resmax = zero ! zero is defined in module 'consts'.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Code without extreme barrier:
! Constraint: con(x) >= 0
!do k = 1, m
!    if (conval(k) .ne. conval(k)) then
!        resmax = conval(k) ! Set resmax=NaN if conval contains NaN
!        exit
!    else
!        resmax = max(resmax, -conval(k))
!    end if
!end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use extreme barrier to cope with 'hidden constraints'
! Constraint: con(x) >= 0
do k = 1, m
    if (conval(k) < -hugecon .or. is_nan(conval(k))) then
        conval(k) = -hugecon ! hugecon is defined in consts
    end if

! This part is NOT extreme barrier. We replace extremely negative values
! of cineq (which leads to no constraint violation) by -hugecon. Otherwise,
! NaN or Inf may occur in the interpolation models.
    if (conval(k) > hugecon) then
        conval(k) = hugecon ! hugecon is defined in consts
    end if

    resmax = max(resmax, -conval(k))
end do
!
! Free memory; note that poutput and pinput are just temporary variables in
! this subroutine. We are NOT in mexFunction!
if (nf > 0) call mxDestroyArray(poutput(1)) ! Only if nf >= 1
call mxDestroyArray(pinput(2))

! Update global variables
nf = nf + int(1, kind(nf))
! Some compiler (e.g., g95) may complain about implicit conversion if
! written as nf = nf+1
fhist(nf) = funval
chist(nf) = resmax

return
end subroutine calcfc
