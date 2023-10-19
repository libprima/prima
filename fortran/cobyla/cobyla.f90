module cobyla_mod
!--------------------------------------------------------------------------------------------------!
! COBYLA_MOD is a module providing the reference implementation of of Powell's COBYLA algorithm in
!
! M. J. D. Powell, A direct search optimization method that models the objective and constraint
! functions by linear interpolation, In Advances in Optimization and Numerical Analysis, eds. S.
! Gomez and J. P. Hennart, pages 51--67, Springer Verlag, Dordrecht, Netherlands, 1994
!
! COBYLA approximately solves
!
!   min F(X) subject to NLCONSTR(X) <= 0, Aineq*X <= Bineq, Aeq*x = Beq, XL <= X <= XU,
!
! where X is a vector of variables that has N components, F is a real-valued objective function, and
! NLCONSTR(X) is an M_NLCON-dimensional vector-valued mapping representing the nonlinear constraints,
! Aineq is an Mineq-by-N matrix, Bineq is an Mineq-dimensional real vector, Aeq is an Meq-by-N
! matrix, Beq is an Meq-dimensional real vector, XL is an N-dimensional real vector, and XU is an
! N-dimensional real vector.
!
! The algorithm employs linear approximations to the objective and nonlinear constraint functions,
! the approximations being formed by linear interpolation at N + 1 points in the space of the
! variables. We regard these interpolation points as vertices of a simplex. The parameter RHO
! controls the size of the simplex and it is reduced automatically from RHOBEG to RHOEND. For each
! RHO the subroutine tries to achieve a good vector of variables for the current size, and then RHO
! is reduced until the value RHOEND is reached. Therefore RHOBEG and RHOEND should be set to
! reasonable initial changes to and the required accuracy in the variables respectively, but this
! accuracy should be viewed as a subject for experimentation because it is not guaranteed. The
! subroutine has an advantage over many of its competitors, however, which is that it treats each
! constraint individually when calculating a change to the variables, instead of lumping the
! constraints together into a single penalty function. The name of the subroutine is derived from
! the phrase Constrained Optimization BY Linear Approximations.
!
! N.B.:
! 1. In Powell's implementation, the constraints are in the form of CONSTR(X) >= 0, whereas we
! consider CONSTR(X) <= 0, where CONSTR is a vector-valued function that wraps all constraints.
! 2. Powell's implementation does not deal with bound and linear constraints explicitly, but we do.
! 3. Our formulation of constraints is consistent with FMINCON of MATLAB.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on the COBYLA paper and Powell's code, with
! modernization, bug fixes, and improvements.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Thursday, October 19, 2023 AM10:39:21
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobyla


contains


subroutine cobyla(calcfc, m_nlcon, x, &
    & f, cstrv, nlconstr, &
    & Aineq, bineq, &
    & Aeq, beq, &
    & xl, xu, &
    & f0, nlconstr0, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, chist, nlchist, maxhist, maxfilt, info)
!--------------------------------------------------------------------------------------------------!
! Among all the arguments, only CALCFC, M_NLCON, and X are obligatory. The others are OPTIONAL and
! you can neglect them unless you are familiar with the algorithm. Any unspecified optional input
! will take the default value detailed below. For instance, we may invoke the solver as follows.
!
! ! First define CALCFC, M_NLCON, and X, and then do the following.
! call cobyla(calcfc, m_nlcon, x, f, cstrv)
!
! or
!
! ! First define CALCFC, M_NLCON, and X, and then do the following.
! call cobyla(calcfc, m_nlcon, x, f, cstrv, rhobeg = 0.5D0, rhoend = 1.0D-3, maxfun = 100)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! IMPORTANT NOTICE: The user must set M_NLCON correctly to the number of nonlinear constraints,
! ! namely the size of NLCONSTR introduced below. Set it to 0 if there is no nonlinear constraint.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! See examples/cobyla_exmp.f90 for a concrete example.
!
! A detailed introduction to the arguments is as follows.
! N.B.: RP and IK are defined in the module CONSTS_MOD. See consts.F90 under the directory named
! "common". By default, RP = kind(0.0D0) and IK = kind(0), with REAL(RP) being the double-precision
! real, and INTEGER(IK) being the default integer. For ADVANCED USERS, RP and IK can be defined by
! setting PRIMA_REAL_PRECISION and PRIMA_INTEGER_KIND in common/ppf.h. Use the default if unsure.
!
! CALCFC
!   Input, subroutine.
!   CALCFC(X, F, NLCONSTR) should evaluate the objective function and nonlinear constraints at the
!   given REAL(RP) vector X; it should set the objective function value to the REAL(RP) scalar F
!   and the nonlinear constraint value to the REAL(RP) vector NLCONSTR. It must be provided by the
!   user, and its definition must conform to the following interface:
!   !-------------------------------------------------------------------------!
!    subroutine calcfc(x, f, nlconstr)
!    real(RP), intent(in) :: x(:)
!    real(RP), intent(out) :: f
!    real(RP), intent(out) :: nlconstr(:)
!    end subroutine calcfc
!   !-------------------------------------------------------------------------!
!   Besides, the subroutine should NOT access NLCONSTR beyond NLCONSTR(1:M_NLCON), where M_NLCON
!   is the second compulsory argument (see below), signifying the number of nonlinear constraints.
!
! M_NLCON
!   Input, INTEGER(IK) scalar.
!   M_NLCON must be set to the number of nonlinear constraints, namely the size of NLCONSTR(X).
!   N.B.:
!   1. M_NLCON must be specified correctly, or the program will crash!!!
!   2. Why don't we define M_NLCON as optional and default it to 0 when it is absent? This is
!   because we need to allocate memory for CONSTR_LOC using M_NLCON. To ensure that the size of
!   CONSTR_LOC is correct, we require the user to specify M_NLCON explicitly.
!
! X
!   Input and output, REAL(RP) vector.
!   As an input, X should be an N-dimensional vector that contains the starting point, N being the
!   dimension of the problem. As an output, X will be set to an approximate minimizer.
!
! F
!   Output, REAL(RP) scalar.
!   F will be set to the objective function value of X at exit.
!
! CSTRV
!   Output, REAL(RP) scalar.
!   CSTRV will be set to the constraint violation of X at exit, i.e.,
!   MAXVAL([0, XL - X, X - XU, Aineq*X - Bineq, ABS(Aeq*X -Beq), NLCONSTR(X)]).
!
! NLCONSTR
!   Output, REAL(RP) vector.
!   NLCONSTR should be an M_NLCON-dimensional vector and will be set to the nonlinear constraint
!   value of X at exit.
!
! Aineq, Bineq
!   Input, REAL(RP) matrix of size [Mineq, N] and REAL vector of size Mineq unless they are both
!   empty, default: [] and [].
!   Aineq and Bineq represent the linear inequality constraints: Aineq*X <= Bineq.
!
! Aeq, Beq
!   Input, REAL(RP) matrix of size [Meq, N] and REAL vector of size Meq unless they are both
!   empty, default: [] and [].
!   Aeq and Beq represent the linear equality constraints: Aeq*X = Beq.
!
! XL, XU
!   Input, REAL(RP) vectors of size N unless they are both empty, default: [] and [].
!   XL and XU represent the lower and upper bounds of the variables: XL <= X <= XU.
!
! F0
!   Input, REAL(RP) scalar.
!   F0, if present, should be set to the objective function value of the starting X.
!
! NLCONSTR0
!   Input, REAL(RP) vector.
!   NLCONSTR0, if present, should be set to the nonlinear constraint value at the starting X; in
!   addition, SIZE(NLCONSTR0) must be M_NLCON, or the solver will abort.
!
! NF
!   Output, INTEGER(IK) scalar.
!   NF will be set to the number of calls of CALCFC at exit.
!
! RHOBEG, RHOEND
!   Inputs, REAL(RP) scalars, default: RHOBEG = 1, RHOEND = 10^-6. RHOBEG and RHOEND must be set to
!   the initial and final values of a trust-region radius, both being positive and RHOEND <= RHOBEG.
!   Typically RHOBEG should be about one tenth of the greatest expected change to a variable, and
!   RHOEND should indicate the accuracy that is required in the final values of the variables.
!
! FTARGET
!   Input, REAL(RP) scalar, default: -Inf.
!   FTARGET is the target function value. The algorithm will terminate when a feasible point with a
!   function value <= FTARGET is found.
!
! CTOL
!   Input, REAL(RP) scalar, default: machine epsilon.
!   CTOL is the tolerance of constraint violation. X is considered feasible if CSTRV(X) <= CTOL.
!   N.B.: 1. CTOL is absolute, not relative. 2. CTOL is used only when selecting the returned X.
!   It does not affect the iterations of the algorithm.
!
! CWEIGHT
!   Input, REAL(RP) scalar, default: CWEIGHT_DFT defined in the module CONSTS_MOD in common/consts.F90.
!   CWEIGHT is the weight that the constraint violation takes in the selection of the returned X.
!
! MAXFUN
!   Input, INTEGER(IK) scalar, default: MAXFUN_DIM_DFT*N with MAXFUN_DIM_DFT defined in the module
!   CONSTS_MOD (see common/consts.F90). MAXFUN is the maximal number of calls of CALCFC.
!
! IPRINT
!   Input, INTEGER(IK) scalar, default: 0.
!   The value of IPRINT should be set to 0, 1, -1, 2, -2, 3, or -3, which controls how much
!   information will be printed during the computation:
!   0: there will be no printing;
!   1: a message will be printed to the screen at the return, showing the best vector of variables
!      found and its objective function value;
!   2: in addition to 1, each new value of RHO is printed to the screen, with the best vector of
!      variables so far and its objective function value; each new value of CPEN is also printed;
!   3: in addition to 2, each function evaluation with its variables will be printed to the screen;
!   -1, -2, -3: the same information as 1, 2, 3 will be printed, not to the screen but to a file
!      named COBYLA_output.txt; the file will be created if it does not exist; the new output will
!      be appended to the end of this file if it already exists.
!   Note that IPRINT = +/-3 can be costly in terms of time and/or space.
!
! ETA1, ETA2, GAMMA1, GAMMA2
!   Input, REAL(RP) scalars, default: ETA1 = 0.1, ETA2 = 0.7, GAMMA1 = 0.5, and GAMMA2 = 2.
!   ETA1, ETA2, GAMMA1, and GAMMA2 are parameters in the updating scheme of the trust-region radius
!   detailed in the subroutine TRRAD in trustregion.f90. Roughly speaking, the trust-region radius
!   is contracted by a factor of GAMMA1 when the reduction ratio is below ETA1, and enlarged by a
!   factor of GAMMA2 when the reduction ratio is above ETA2. It is required that 0 < ETA1 <= ETA2
!   < 1 and 0 < GAMMA1 < 1 < GAMMA2. Normally, ETA1 <= 0.25. It is NOT advised to set ETA1 >= 0.5.
!
! XHIST, FHIST, CHIST, NLCHIST, MAXHIST
!   XHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   FHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   CHIST: Output, ALLOCATABLE rank 1 REAL(RP) array;
!   NLCHIST: Output, ALLOCATABLE rank 2 REAL(RP) array;
!   MAXHIST: Input, INTEGER(IK) scalar, default: MAXFUN
!   XHIST, if present, will output the history of iterates; FHIST, if present, will output the
!   history function values; CHIST, if present, will output the history of constraint violations;
!   NLCHIST, if present, will output the history of nonlinear constraint values; MAXHIST should be
!   a nonnegative integer, and XHIST/FHIST/CHIST/NLCHIST will output only the history of the last
!   MAXHIST iterations. Therefore, MAXHIST= 0 means XHIST/FHIST/NLCHIST/CHIST will output nothing,
!   while setting MAXHIST = MAXFUN requests XHIST/FHIST/CHIST/NLCHIST to output all the history.
!   If XHIST is present, its size at exit will be (N, min(NF, MAXHIST)); if FHIST is present, its
!   size at exit will be min(NF, MAXHIST); if CHIST is present, its size at exit will be
!   min(NF, MAXHIST); if NLCHIST is present, its size at exit will be (M_NLCON, min(NF, MAXHIST)).
!
!   IMPORTANT NOTICE:
!   Setting MAXHIST to a large value can be costly in terms of memory for large problems.
!   MAXHIST will be reset to a smaller value if the memory needed exceeds MAXHISTMEM defined in
!   CONSTS_MOD (see consts.F90 under the directory named "common").
!   Use *HIST with caution!!! (N.B.: the algorithm is NOT designed for large problems).
!
! MAXFILT
!   Input, INTEGER(IK) scalar.
!   MAXFILT is a nonnegative integer indicating the maximal length of the filter used for selecting
!   the returned solution; default: MAXFILT_DFT (a value lower than MIN_MAXFILT is not recommended);
!   see common/consts.F90 for the definitions of MAXFILT_DFT and MIN_MAXFILT.
!
! INFO
!   Output, INTEGER(IK) scalar.
!   INFO is the exit flag. It will be set to one of the following values defined in the module
!   INFOS_MOD (see common/infos.f90):
!   SMALL_TR_RADIUS: the lower bound for the trust region radius is reached;
!   FTARGET_ACHIEVED: the target function value is reached;
!   MAXFUN_REACHED: the objective function has been evaluated MAXFUN times;
!   MAXTR_REACHED: the trust region iteration has been performed MAXTR times (MAXTR = 2*MAXFUN);
!   NAN_INF_X: NaN or Inf occurs in X;
!   DAMAGING_ROUNDING: rounding errors are becoming damaging.
!   !--------------------------------------------------------------------------!
!   The following case(s) should NEVER occur unless there is a bug.
!   NAN_INF_F: the objective function returns NaN or +Inf;
!   NAN_INF_MODEL: NaN or Inf occurs in the model;
!   TRSUBP_FAILED: a trust region step failed to reduce the model
!   !--------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, TEN, TENTH, EPS, REALMAX
use, non_intrinsic :: debug_mod, only : assert, errstop, warning
use, non_intrinsic :: evaluate_mod, only : evaluate, moderatex, moderatec
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
use, non_intrinsic :: infos_mod, only : INVALID_INPUT
use, non_intrinsic :: linalg_mod, only : trueloc, matprod
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJCON
use, non_intrinsic :: selectx_mod, only : isbetter
use, non_intrinsic :: preproc_mod, only : preproc
use, non_intrinsic :: string_mod, only : num2str

! Solver-specific modules
use, non_intrinsic :: cobylb_mod, only : cobylb

implicit none

! Compulsory arguments
procedure(OBJCON) :: calcfc ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
real(RP), intent(inout) :: x(:)     ! X(N)
integer(IK), intent(in) :: m_nlcon  ! Number of constraints defined in CALCFC

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfilt
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
real(RP), intent(in), optional :: Aeq(:, :)  ! Aeq(Meq, N)
real(RP), intent(in), optional :: Aineq(:, :)  ! Aineq(Mineq, N)
real(RP), intent(in), optional :: beq(:)  ! Beq(Meq)
real(RP), intent(in), optional :: bineq(:)  ! Bineq(Mineq)
real(RP), intent(in), optional :: nlconstr0(:)    ! NLCONSTR0(M_NLCON)
real(RP), intent(in), optional :: ctol
real(RP), intent(in), optional :: cweight
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: f0
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend
real(RP), intent(in), optional :: xl(:)  ! XL(N)
real(RP), intent(in), optional :: xu(:)  ! XU(N)

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), optional :: cstrv
real(RP), intent(out), optional :: f
real(RP), intent(out), optional :: nlconstr(:) ! NLCONSTR(M_NLCON)
real(RP), intent(out), optional, allocatable :: chist(:)    ! CHIST(MAXCHIST)
real(RP), intent(out), optional, allocatable :: fhist(:)    ! FHIST(MAXFHIST)
real(RP), intent(out), optional, allocatable :: nlchist(:, :)     ! NLCHIST(M_NLCON, MAXCONHIST)
real(RP), intent(out), optional, allocatable :: xhist(:, :) ! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'COBYLA'
character(len=*), parameter :: srname = 'COBYLA'
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: m
integer(IK) :: maxfilt_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: meq
integer(IK) :: mineq
integer(IK) :: mxl
integer(IK) :: mxu
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: nhist
integer(IK), allocatable :: ixl(:)
integer(IK), allocatable :: ixu(:)
real(RP) :: cstrv_loc
real(RP) :: ctol_loc
real(RP) :: cweight_loc
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: f_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP), allocatable :: Aeq_loc(:, :)  ! Aeq_LOC(Meq, N)
real(RP), allocatable :: Aineq_loc(:, :)  ! Aineq_LOC(Mineq, N)
real(RP), allocatable :: amat(:, :)  ! AMAT(N, M_LCON); each column corresponds to a linear constraint
real(RP), allocatable :: beq_loc(:)  ! Beq_LOC(Meq)
real(RP), allocatable :: bineq_loc(:)  ! Bineq_LOC(Mineq)
real(RP), allocatable :: bvec(:)  ! BVEC(M_LCON)
real(RP), allocatable :: chist_loc(:)  ! CHIST_LOC(MAXCHIST)
real(RP), allocatable :: conhist_loc(:, :)  ! CONHIST_LOC(M_NLCON, MAXCONHIST)
real(RP), allocatable :: constr_loc(:)  ! CONSTR_LOC(M_NLCON)
real(RP), allocatable :: fhist_loc(:)   ! FHIST_LOC(MAXFHIST)
real(RP), allocatable :: xhist_loc(:, :)  ! XHIST_LOC(N, MAXXHIST)
real(RP), allocatable :: xl_loc(:)  ! XL_LOC(N)
real(RP), allocatable :: xu_loc(:)  ! XU_LOC(N)

! Sizes
if (present(bineq)) then
    mineq = int(size(bineq), kind(mineq))
else
    mineq = 0
end if
if (present(beq)) then
    meq = int(size(beq), kind(meq))
else
    meq = 0
end if
if (present(xl)) then
    mxl = int(count(xl > -REALMAX), kind(mxl))
else
    mxl = 0
end if
if (present(xu)) then
    mxu = int(count(xu < REALMAX), kind(mxu))
else
    mxu = 0
end if
m = mxu + mxl + 2_IK * meq + mineq + m_nlcon
n = int(size(x), kind(n))


! Preconditions
if (DEBUGGING) then
    call assert(m_nlcon >= 0, 'M_NLCON >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    if (present(nlconstr)) then
        call assert(size(nlconstr) == m_nlcon, 'SIZE(NLCONSTR) == M_NLCON', srname)
    end if
    call assert(present(Aineq) .eqv. present(bineq), 'Aineq and Bineq are both present or both absent', srname)
    if (present(Aineq)) then
        call assert((size(Aineq, 1) == mineq .and. size(Aineq, 2) == n) &
            & .or. (size(Aineq, 1) == 0 .and. size(Aineq, 2) == 0 .and. mineq == 0), &
            & 'SIZE(Aineq) == [Mineq, N] unless Aineq and Bineq are both empty', srname)
    end if
    call assert(present(Aeq) .eqv. present(beq), 'Aeq and Beq are both present or both absent', srname)
    if (present(Aeq)) then
        call assert((size(Aeq, 1) == meq .and. size(Aeq, 2) == n) &
            & .or. (size(Aeq, 1) == 0 .and. size(Aeq, 2) == 0 .and. meq == 0), &
            & 'SIZE(Aeq) == [Meq, N] unless Aeq and Beq are both empty', srname)
    end if
    if (present(xl)) then
        call assert(size(xl) == n .or. size(xl) == 0, 'SIZE(XL) == N unless XL is empty', srname)
    end if
    if (present(xu)) then
        call assert(size(xu) == n .or. size(xu) == 0, 'SIZE(XU) == N unless XU is empty', srname)
    end if
    call assert(present(f0) .eqv. present(nlconstr0), 'F0 and NLCONSTR0 are both present or both absent', srname)
end if

! Exit if the size of NLCONSTR0 is inconsistent with M_NLCON.
if (present(nlconstr0)) then
    if (size(nlconstr0) /= m_nlcon) then
        if (DEBUGGING) then
            call errstop(srname, 'SIZE(NLCONSTR0) /= M_NLCON. Exiting', INVALID_INPUT)
        else
            call warning(srname, 'SIZE(NLCONSTR0) /= M_NLCON. Exiting')
            return  ! This may be problematic, as outputs like F are undefined.
        end if
    end if
end if


! Read the inputs.

call safealloc(Aineq_loc, mineq, n)  ! NOT removable even in F2003, as Aineq may be absent or of size 0-by-0.
if (present(Aineq) .and. mineq > 0) then
    ! We must check Mineq > 0. Otherwise, the size of Aineq_LOC may be changed to 0-by-0 due to
    ! automatic (re)allocation if that is the size of Aineq; we allow Aineq to be 0-by-0, but
    ! Aineq_LOC should be n-by-0.
    Aineq_loc = Aineq
end if

call safealloc(bineq_loc, mineq)  ! NOT removable even in F2003, as Bineq may be absent.
if (present(bineq)) then
    bineq_loc = bineq
end if

call safealloc(Aeq_loc, meq, n)  ! NOT removable even in F2003, as Aeq may be absent or of size 0-by-0.
if (present(Aeq) .and. meq > 0) then
    ! We must check Meq > 0. Otherwise, the size of Aeq_LOC may be changed to 0-by-0 due to
    ! automatic (re)allocation if that is the size of Aeq; we allow Aeq to be 0-by-0, but
    ! Aeq_LOC should be n-by-0.
    Aeq_loc = Aeq
end if

call safealloc(beq_loc, meq)  ! NOT removable even in F2003, as Beq may be absent.
if (present(beq)) then
    beq_loc = beq
end if

call safealloc(xl_loc, n)  ! NOT removable even in F2003, as XL may be absent.
if (present(xl)) then
    xl_loc = xl
else
    xl_loc = -REALMAX
end if
call safealloc(ixl, mxl)
ixl = trueloc(xl_loc > -REALMAX)

call safealloc(xu_loc, n)  ! NOT removable even in F2003, as XU may be absent.
if (present(xu)) then
    xu_loc = xu
else
    xu_loc = REALMAX
end if
call safealloc(ixu, mxu)
ixu = trueloc(xu_loc < REALMAX)

! Wrap the linear and bound constraints into a single constraint: AMAT^T*X <= BVEC.
call get_lincon(Aeq_loc, Aineq_loc, beq_loc, bineq_loc, xl_loc, xu_loc, amat, bvec)

! Allocate memory for CONSTR_LOC.
call safealloc(constr_loc, m)  ! NOT removable even in F2003!

! Set [F_LOC, CONSTR_LOC] to [F(X0), CONSTR(X0)] after evaluating the latter if needed. In this way,
! COBYLB only needs one interface.
constr_loc(1:m - m_nlcon) = moderatec(matprod(x, amat) - bvec)
if (present(f0) .and. present(nlconstr0) .and. all(is_finite(x))) then
    f_loc = f0
    constr_loc(m - m_nlcon + 1:m) = nlconstr0
else
    x = moderatex(x)
    call evaluate(calcfc, x, f_loc, constr_loc(m - m_nlcon + 1:m))
    ! N.B.: Do NOT call FMSG, SAVEHIST, or SAVEFILT for the function/constraint evaluation at X0.
    ! They will be called during the initialization, which will read the function/constraint at X0.
end if
cstrv_loc = maxval([ZERO, constr_loc])

! If RHOBEG is present, then RHOBEG_LOC is a copy of RHOBEG; otherwise, RHOBEG_LOC takes the default
! value for RHOBEG, taking the value of RHOEND into account. Note that RHOEND is considered only if
! it is present and it is VALID (i.e., finite and positive). The other inputs are read similarly.
if (present(rhobeg)) then
    rhobeg_loc = rhobeg
elseif (present(rhoend)) then
    ! Fortran does not take short-circuit evaluation of logic expressions. Thus it is WRONG to
    ! combine the evaluation of PRESENT(RHOEND) and the evaluation of IS_FINITE(RHOEND) as
    ! "IF (PRESENT(RHOEND) .AND. IS_FINITE(RHOEND))". The compiler may choose to evaluate the
    ! IS_FINITE(RHOEND) even if PRESENT(RHOEND) is false!
    if (is_finite(rhoend) .and. rhoend > 0) then
        rhobeg_loc = max(TEN * rhoend, RHOBEG_DFT)
    else
        rhobeg_loc = RHOBEG_DFT
    end if
else
    rhobeg_loc = RHOBEG_DFT
end if

if (present(rhoend)) then
    rhoend_loc = rhoend
elseif (rhobeg_loc > 0) then
    rhoend_loc = max(EPS, min(TENTH * rhobeg_loc, RHOEND_DFT))
else
    rhoend_loc = RHOEND_DFT
end if

if (present(ctol)) then
    ctol_loc = ctol
else
    ctol_loc = CTOL_DFT
end if

if (present(cweight)) then
    cweight_loc = cweight
else
    cweight_loc = CWEIGHT_DFT
end if

if (present(ftarget)) then
    ftarget_loc = ftarget
else
    ftarget_loc = FTARGET_DFT
end if

if (present(maxfun)) then
    maxfun_loc = maxfun
else
    maxfun_loc = MAXFUN_DIM_DFT * n
end if

if (present(iprint)) then
    iprint_loc = iprint
else
    iprint_loc = IPRINT_DFT
end if

if (present(eta1)) then
    eta1_loc = eta1
elseif (present(eta2)) then
    if (eta2 > 0 .and. eta2 < 1) then
        eta1_loc = max(EPS, eta2 / 7.0_RP)
    end if
else
    eta1_loc = TENTH
end if

if (present(eta2)) then
    eta2_loc = eta2
elseif (eta1_loc > 0 .and. eta1_loc < 1) then
    eta2_loc = (eta1_loc + TWO) / 3.0_RP
else
    eta2_loc = 0.7_RP
end if

if (present(gamma1)) then
    gamma1_loc = gamma1
else
    gamma1_loc = HALF
end if

if (present(gamma2)) then
    gamma2_loc = gamma2
else
    gamma2_loc = TWO
end if

if (present(maxhist)) then
    maxhist_loc = maxhist
else
    maxhist_loc = maxval([maxfun_loc, n + 2_IK, MAXFUN_DIM_DFT * n])
end if

if (present(maxfilt)) then
    maxfilt_loc = maxfilt
else
    maxfilt_loc = MAXFILT_DFT
end if

! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & m=m, is_constrained=(m > 0), ctol=ctol_loc, cweight=cweight_loc, eta1=eta1_loc, &
    & eta2=eta2_loc, gamma1=gamma1_loc, gamma2=gamma2_loc, maxfilt=maxfilt_loc)

! Further revise MAXHIST_LOC according to MAXHISTMEM, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! CHIST = NaN(1, MAXFUN), NLCHIST = NaN(M_NLCON, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST =
! NaN(N, MAXFUN) if they are requested; replace MAXFUN with 0 for the history not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc, &
    & present(chist), chist_loc, m, present(nlchist), conhist_loc)


!-------------------- Call COBYLB, which performs the real calculations. --------------------------!
!call cobylb(calcfc_internal, iprint_loc, maxfilt_loc, maxfun_loc, ctol_loc, cweight_loc, eta1_loc, eta2_loc, &
call cobylb(calcfc, iprint_loc, maxfilt_loc, maxfun_loc, amat, bvec, ctol_loc, cweight_loc, &
    & eta1_loc, eta2_loc, ftarget_loc, gamma1_loc, gamma2_loc, rhobeg_loc, rhoend_loc, constr_loc, &
    & f_loc, x, nf_loc, chist_loc, conhist_loc, cstrv_loc, fhist_loc, xhist_loc, info_loc)
!--------------------------------------------------------------------------------------------------!

! Deallocate variables not needed any more. Indeed, automatic allocation will take place at exit.
deallocate (Aineq_loc, Aeq_loc, amat, bineq_loc, beq_loc, bvec, xl_loc, xu_loc)


! Write the outputs.

if (present(f)) then
    f = f_loc
end if

if (present(cstrv)) then
    cstrv = cstrv_loc
end if

if (present(nlconstr)) then
    nlconstr = constr_loc(m - m_nlcon + 1:m)
end if
deallocate (constr_loc)

if (present(nf)) then
    nf = nf_loc
end if

if (present(info)) then
    info = info_loc
end if

! Copy XHIST_LOC to XHIST if needed.
if (present(xhist)) then
    nhist = min(nf_loc, int(size(xhist_loc, 2), IK))
    !----------------------------------------------------!
    call safealloc(xhist, n, nhist)  ! Removable in F2003.
    !----------------------------------------------------!
    xhist = xhist_loc(:, 1:nhist)
    ! N.B.:
    ! 0. Allocate XHIST as long as it is present, even if the size is 0; otherwise, it will be
    ! illegal to enquire XHIST after exit.
    ! 1. Even though Fortran 2003 supports automatic (re)allocation of allocatable arrays upon
    ! intrinsic assignment, we keep the line of SAFEALLOC, because some very new compilers (Absoft
    ! Fortran 21.0) are still not standard-compliant in this respect.
    ! 2. NF may not be present. Hence we should NOT use NF but NF_LOC.
    ! 3. When SIZE(XHIST_LOC, 2) > NF_LOC, which is the normal case in practice, XHIST_LOC contains
    ! GARBAGE in XHIST_LOC(:, NF_LOC + 1 : END). Therefore, we MUST cap XHIST at NF_LOC so that
    ! XHIST contains only valid history. For this reason, there is no way to avoid allocating
    ! two copies of memory for XHIST unless we declare it to be a POINTER instead of ALLOCATABLE.
end if
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet we prefer to deallocate
! them immediately when they finish their jobs.
deallocate (xhist_loc)

! Copy FHIST_LOC to FHIST if needed.
if (present(fhist)) then
    nhist = min(nf_loc, int(size(fhist_loc), IK))
    !--------------------------------------------------!
    call safealloc(fhist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    fhist = fhist_loc(1:nhist)  ! The same as XHIST, we must cap FHIST at NF_LOC.
end if
deallocate (fhist_loc)

! Copy CHIST_LOC to CHIST if needed.
if (present(chist)) then
    nhist = min(nf_loc, int(size(chist_loc), IK))
    !--------------------------------------------------!
    call safealloc(chist, nhist)  ! Removable in F2003.
    !--------------------------------------------------!
    chist = chist_loc(1:nhist)  ! The same as XHIST, we must cap CHIST at NF_LOC.
end if
deallocate (chist_loc)

! Copy CONHIST_LOC to NLCHIST if needed.
! N.B.: We need only the nonlinear part of the history. Therefore, one may modify COBYLB so that it
! records only the history of nonlinear constraints rather than that of all constraints, which is
! what CONHIST_LOC records in the current implementation. This will save memory and time, but will
! complicate a bit the code of COBYLB. We prefer to keep the code simple, assuming that such a
! difference in memory and time is not a problem in the context of derivative-free optimization.
! A similar comment can be made on CONSTR_LOC and NLCONSTR, which are related to CONFILT in COBYLB.
if (present(nlchist)) then
    nhist = min(nf_loc, int(size(conhist_loc, 2), IK))
    !---------------------------------------------------------------!
    call safealloc(nlchist, m_nlcon, nhist)  ! Removable in F2003.
    !---------------------------------------------------------------!
    nlchist = conhist_loc(m - m_nlcon + 1:m, 1:nhist)  ! The same as XHIST, we must cap NLCHIST at NF_LOC.
end if
deallocate (conhist_loc)

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist) .or. present(chist) .or. present(nlchist)) .and. maxhist_loc < nf_loc) then
    call warning(solver, 'Only the history of the last '//num2str(maxhist_loc)//' iteration(s) is recorded')
end if

! Postconditions
if (DEBUGGING) then
    call assert(nf_loc <= maxfun_loc, 'NF <= MAXFUN', srname)
    call assert(size(x) == n .and. .not. any(is_nan(x)), 'SIZE(X) == N, X does not contain NaN', srname)
    nhist = min(nf_loc, maxhist_loc)
    if (present(xhist)) then
        call assert(size(xhist, 1) == n .and. size(xhist, 2) == nhist, 'SIZE(XHIST) == [N, NHIST]', srname)
        call assert(.not. any(is_nan(xhist)), 'XHIST does not contain NaN', srname)
    end if
    if (present(fhist)) then
        call assert(size(fhist) == nhist, 'SIZE(FHIST) == NHIST', srname)
        call assert(.not. any(is_nan(fhist) .or. is_posinf(fhist)), 'FHIST does not contain NaN/+Inf', srname)
    end if
    if (present(chist)) then
        call assert(size(chist) == nhist, 'SIZE(CHIST) == NHIST', srname)
        call assert(.not. any(chist < 0 .or. is_nan(chist) .or. is_posinf(chist)), &
            & 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
    end if
    if (present(nlchist)) then
        call assert(size(nlchist, 1) == m_nlcon .and. size(nlchist, 2) == nhist, 'SIZE(NLCHIST) == [M_NLCON, NHIST]', srname)
        call assert(.not. any(is_nan(nlchist) .or. is_posinf(nlchist)), 'NLCHIST does not contain NaN/+Inf', srname)
    end if
    if (present(fhist) .and. present(chist)) then
        call assert(.not. any(isbetter(fhist(1:nhist), chist(1:nhist), f_loc, cstrv_loc, ctol_loc)), &
            & 'No point in the history is better than X', srname)
    end if
end if

end subroutine cobyla


subroutine get_lincon(Aeq, Aineq, beq, bineq, xl, xu, amat, bvec)
!--------------------------------------------------------------------------------------------------!
! This subroutine wraps the linear and bound constraints into a single constraint: AMAT^T*X <= BVEC.
! N.B.:
! 1. The linear inequality constraints received by COBYLA is AMAT^T * X <= BVEC. Note that Each
! column of AMAT corresponds to a constraint. This is different from Aineq and Aeq, whose rows
! correspond to constraints. AMAT is defined in this way because it is accessed in columns during
! the computation, and because Fortran saves arrays in the column-major order. In Python/C
! implementations, AMAT should be transposed.
! 2. LINCOA normalizes the linear constraints so that each constraint has a gradient of norm 1.
! However, COBYLA does not do this.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : eye, trueloc
use, non_intrinsic :: memory_mod, only : safealloc

implicit none

! Inputs
real(RP), intent(in) :: Aeq(:, :)
real(RP), intent(in) :: Aineq(:, :)
real(RP), intent(in) :: beq(:)
real(RP), intent(in) :: bineq(:)
real(RP), intent(in) :: xl(:)
real(RP), intent(in) :: xu(:)

! Outputs
real(RP), intent(out), allocatable :: amat(:, :)
real(RP), intent(out), allocatable :: bvec(:)

! Local variables
character(len=*), parameter :: srname = 'GET_LINCON'
integer(IK) :: m_lcon
integer(IK) :: meq
integer(IK) :: mineq
integer(IK) :: mxl
integer(IK) :: mxu
integer(IK) :: n
integer(IK), allocatable :: ixl(:)
integer(IK), allocatable :: ixu(:)
real(RP) :: idmat(size(xl), size(xl))

! Sizes
n = int(size(xl), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(Aineq, 1) == size(bineq) .and. size(Aineq, 2) == n, 'SIZE(AINEQ) == [SIZE(BINEQ), N]', srname)
    call assert(size(Aeq, 1) == size(beq) .and. size(Aeq, 2) == n, 'SIZE(AEQ) == [SIZE(BEQ), N]', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == SIZE(XU) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

! Decide the number of nontrivial constraints.
mxl = int(count(xl > -REALMAX), kind(mxl))
mxu = int(count(xu < REALMAX), kind(mxu))
meq = int(size(beq), kind(meq))
mineq = int(size(bineq), kind(mineq))
m_lcon = mxl + mxu + 2_IK * meq + mineq  ! The final number of linear inequality constraints.

! Allocate memory. Removable in F2003.
call safealloc(ixl, mxl)
call safealloc(ixu, mxu)
call safealloc(amat, n, m_lcon)
call safealloc(bvec, m_lcon)

! Define the indices of the nontrivial bound constraints.
ixl = trueloc(xl > -REALMAX)
ixu = trueloc(xu < REALMAX)

! Wrap the linear constraints.
! The bound constraint XL <= X <= XU is handled as two constraints -X <= -XL, X <= XU.
! The equality constraint Aeq*X = Beq is handled as two constraints -Aeq*X <= -Beq, Aeq*X <= Beq.
! N.B.:
! 1. The treatment of the equality constraints is naive. One may choose to eliminate them instead.
! 2. The code below is quite inefficient in terms of memory, but we prefer readability.
idmat = eye(n, n)
amat = reshape(shape=shape(amat), source= &
    & [-idmat(:, ixl), idmat(:, ixu), -transpose(Aeq), transpose(Aeq), transpose(Aineq)])
bvec = [-xl(ixl), xu(ixu), -beq, beq, bineq]
!!MATLAB code:
!!amat = [-idmat(:, ixl), idmat(:, ixu), -Aeq', Aeq', Aineq'];
!!bvec = [-xl(ixl); xu(ixu); -beq; beq; bineq];

! Deallocate memory.
deallocate (ixl, ixu)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(amat, 1) == size(xl) .and. size(amat, 2) == size(bvec), &
        & 'SIZE(AMAT) == [SIZE(X), SIZE(BVEC)]', srname)
end if
end subroutine get_lincon


end module cobyla_mod
