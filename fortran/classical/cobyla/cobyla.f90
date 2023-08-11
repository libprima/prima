module cobyla_mod
!--------------------------------------------------------------------------------------------------!
! Classical mode. Not maintained. Strongly discouraged. Please use the modernized version instead.
!
! The usage is the same as the modernized version.
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: cobyla


contains


subroutine cobyla(calcfc, m_nlcon, x, f, &
    & cstrv, nlconstr, &
    & Aineq, bineq, &
    & Aeq, beq, &
    & xl, xu, &
    & f0, nlconstr0, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, iprint, &
    & eta1, eta2, gamma1, gamma2, xhist, fhist, chist, nlchist, maxhist, maxfilt, info)

! Common modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS, REALMAX
use, non_intrinsic :: debug_mod, only : assert, errstop, warning
use, non_intrinsic :: evaluate_mod, only : evaluate, moderatex, moderatec
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : trueloc, matprod
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJCON
use, non_intrinsic :: preproc_mod, only : preproc
use, non_intrinsic :: string_mod, only : num2str

! Solver-specific modules
use, non_intrinsic :: cobylb_mod, only : cobylb

implicit none

! Compulsory arguments
procedure(OBJCON) :: calcfc
real(RP), intent(inout) :: x(:)
real(RP), intent(out) :: f
integer(IK), intent(in) :: m_nlcon

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
real(RP), intent(out), optional :: nlconstr(:)
real(RP), intent(out), allocatable, optional :: chist(:)
real(RP), intent(out), allocatable, optional :: nlchist(:, :)
real(RP), intent(out), allocatable, optional :: fhist(:)
real(RP), intent(out), allocatable, optional :: xhist(:, :)
real(RP), intent(out), optional :: cstrv

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
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
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP), allocatable :: Aeq_loc(:, :)  ! Aeq_LOC(Meq, N)
real(RP), allocatable :: Aineq_loc(:, :)  ! Aineq_LOC(Mineq, N)
real(RP), allocatable :: beq_loc(:)  ! Beq_LOC(Meq)
real(RP), allocatable :: bineq_loc(:)  ! Bineq_LOC(Mineq)
real(RP), allocatable :: constr_loc(:)
real(RP), allocatable :: chist_loc(:)
real(RP), allocatable :: conhist_loc(:, :)
real(RP), allocatable :: fhist_loc(:)
real(RP), allocatable :: xhist_loc(:, :)
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
    call assert(mineq >= 0, 'Mineq >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
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
    call assert(present(f0) .eqv. present(nlconstr0), 'F0 and NLCONSTR0 are both present or both absent', srname)
end if

! Exit if the size of NLCONSTR0 is inconsistent with M_NLCON.
if (present(nlconstr0)) then
    if (size(nlconstr0) /= m_nlcon) then
        if (DEBUGGING) then
            call errstop(srname, 'SIZE(NLCONSTR0) /= M_NLCON. Exiting')
        else
            call warning(srname, 'SIZE(NLCONSTR0) /= M_NLCON. Exiting')
            return
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

! Allocate memory for CONSTR_LOC.
call safealloc(constr_loc, m)  ! NOT removable even in F2003!
!! If the user provides the function & constraint value at X0, then set up F_X0 and NLCONSTR_X0.
if (present(f0) .and. present(nlconstr0) .and. all(is_finite(x))) then
    f = f0
    constr_loc = moderatec(-[x(ixl) - xl_loc(ixl), xu_loc(ixu) - x(ixu), &
    & matprod(Aeq_loc, x) - beq_loc, beq_loc - matprod(Aeq_loc, x), &
    & bineq_loc - matprod(Aineq_loc, x), -nlconstr0])
    constr_loc = -constr_loc
    cstrv_loc = maxval([ZERO, -constr_loc])
else
    ! Replace any NaN in X by ZERO and Inf/-Inf in X by REALMAX/-REALMAX.
    x = moderatex(x)
    call evaluate(calcfc_internal, x, f, constr_loc)
    constr_loc = -constr_loc
    cstrv_loc = maxval([ZERO, -constr_loc])
end if

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
    if (eta2 > ZERO .and. eta2 < ONE) then
        eta1_loc = max(EPS, eta2 / 7.0_RP)
    end if
else
    eta1_loc = TENTH
end if

if (present(eta2)) then
    eta2_loc = eta2
elseif (eta1_loc > ZERO .and. eta1_loc < ONE) then
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
    & m=m, ctol=ctol_loc, cweight=cweight_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, &
    & gamma2=gamma2_loc, maxfilt=maxfilt_loc)

! Further revise MAXHIST_LOC according to MAXMEMORY, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! CHIST = NaN(1, MAXFUN), NLCHIST = NaN(M_NLCON, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc, &
    & present(chist), chist_loc, m, present(nlchist), conhist_loc)

!--------------------------------------------------------------------------------------------------!
!-------------------- Call COBYLB, which performs the real calculations. --------------------------!
!!!! ETA1, ETA2, GAMMA1, GAMMA2, MAXFILT, CTOL, CWEIGHT are not used in the classical mode. !!!!
call cobylb(calcfc_internal, iprint_loc, maxfun_loc, rhobeg_loc, rhoend_loc, constr_loc, x, cstrv_loc, f, info_loc, &
    & nf_loc, xhist_loc, fhist_loc, chist_loc, conhist_loc)
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

! Write the outputs.

! Copy CONSTR_LOC to NLCONSTR if needed.
if (present(nlconstr)) then
    nlconstr = constr_loc(m - m_nlcon + 1:m)
end if
deallocate (constr_loc)

if (present(cstrv)) then
    cstrv = cstrv_loc
end if

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
if (present(nlchist)) then
    nhist = min(nf_loc, int(size(conhist_loc, 2), IK))
    !----------------------------------------------------------!
    call safealloc(nlchist, m_nlcon, nhist)  ! Removable in F2003.
    !----------------------------------------------------------!
    nlchist = conhist_loc(m - m_nlcon + 1:m, 1:nhist)  ! The same as XHIST, we must cap NLCHIST at NF_LOC.
end if
deallocate (conhist_loc)

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist) .or. present(chist) .or. present(nlchist)) .and. maxhist_loc < nf_loc) then
    call warning(solver, 'Only the history of the last '//num2str(maxhist_loc)//' iteration(s) is recorded')
end if


contains


subroutine calcfc_internal(x_internal, f_internal, constr_internal)
implicit none
! Inputs
real(RP), intent(in) :: x_internal(:)
! Outputs
real(RP), intent(out) :: f_internal
real(RP), intent(out) :: constr_internal(:)
! Local variables
real(RP) :: constr_nlc(m_nlcon)

call calcfc(x_internal, f_internal, constr_nlc)
constr_internal = -[x_internal(ixl) - xl_loc(ixl), xu_loc(ixu) - x_internal(ixu), &
    & matprod(Aeq_loc, x_internal) - beq_loc, beq_loc - matprod(Aeq_loc, x_internal), &
    & bineq_loc - matprod(Aineq_loc, x_internal), -constr_nlc]

end subroutine


end subroutine cobyla


end module cobyla_mod
