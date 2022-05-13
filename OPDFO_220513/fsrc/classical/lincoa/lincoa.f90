module lincoa_mod
!--------------------------------------------------------------------------------------------------!
! Classical mode. Not maintained. Not recommended. Please use the modernized version instead.
!
! The usage is the same as the modernized version.
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lincoa


contains


subroutine lincoa(calfun, x, f, &
    & cstrv, &
    & A, b, &
    & nf, rhobeg, rhoend, ftarget, ctol, cweight, maxfun, npt, iprint, &
    & eta1, eta2, gamma1, gamma2, xhist, fhist, chist, maxhist, maxfilt, info)

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, MAXFILT_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, CTOL_DFT, CWEIGHT_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, errstop, warning
use, non_intrinsic :: evaluate_mod, only : evaluate, moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: preproc_mod, only : preproc

implicit none

! Compulsory arguments
procedure(OBJ) :: calfun
real(RP), intent(inout) :: x(:)
real(RP), intent(out) :: f

! Optional inputs
integer(IK), intent(in), optional :: iprint
integer(IK), intent(in), optional :: maxfilt
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
integer(IK), intent(in), optional :: npt
real(RP), intent(in), optional :: A(:, :)
real(RP), intent(in), optional :: b(:)
real(RP), intent(in), optional :: ctol
real(RP), intent(in), optional :: cweight
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), allocatable, optional :: chist(:)
real(RP), intent(out), allocatable, optional :: fhist(:)
real(RP), intent(out), allocatable, optional :: xhist(:, :)
real(RP), intent(out), optional :: cstrv

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'LINCOA'
character(len=MSGLEN) :: wmsg
integer(IK) :: i
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfilt_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: npt_loc
integer(IK) :: nhist
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
real(RP), allocatable :: A_loc(:, :)
real(RP), allocatable :: b_loc(:)
real(RP), allocatable :: chist_loc(:)
real(RP), allocatable :: fhist_loc(:)
real(RP), allocatable :: xhist_loc(:, :)
real(RP) :: smallx
real(RP) :: sum
logical :: constr_modified

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Working variables
real(RP) :: temp
real(RP), allocatable :: w(:)
integer(IK) :: iw, iamat, ib, ndim, ixb, ixp, ifv, ixs, ixo, igo, ihq, ipq, ibmat, izmat, istp, isp,&
    & ixn, iac, irc, iqf, irf, ipqw, j, m
integer(IK), allocatable :: iact(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sizes
if (present(b)) then
    m = int(size(b), kind(m))
else
    m = 0_IK
end if
n = int(size(x), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(present(A) .eqv. present(b), 'A and B are both present or both absent', srname)
    if (present(A)) then
        call assert((size(A, 1) == n .and. size(A, 2) == m) &
            & .or. (size(A, 1) == 0 .and. size(A, 2) == 0 .and. m == 0), &
            & 'SIZE(A) == [N, M] unless A and B are both empty', srname)
    end if
end if

! Read the inputs

x = moderatex(x)

call safealloc(A_loc, n, m)
if (present(A)) then
    A_loc = A
end if

call safealloc(b_loc, m)
if (present(b)) then
    b_loc = b
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
    if (is_finite(rhoend) .and. rhoend > ZERO) then
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

if (present(npt)) then
    npt_loc = npt
elseif (maxfun_loc >= 1) then
    npt_loc = max(n + 2_IK, min(maxfun_loc - 1_IK, 2_IK * n + 1_IK))
else
    npt_loc = 2_IK * n + 1_IK
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
    maxhist_loc = maxval([maxfun_loc, n + 3_IK, MAXFUN_DIM_DFT * n])
end if

if (present(maxfilt)) then
    maxfilt_loc = maxfilt
else
    maxfilt_loc = MAXFILT_DFT
end if

! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & npt=npt_loc, ctol=ctol_loc, cweight=cweight_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, &
    & gamma2=gamma2_loc, maxfilt=maxfilt_loc)

! Further revise MAXHIST_LOC according to MAXMEMORY, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! CHIST = NaN(1, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc, present(chist), chist_loc)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Working space
call safealloc(w, m * (2_IK + n) + npt_loc * (4_IK + n + npt_loc) + n * (9_IK + 3_IK * n) + &
    & max(m + 3_IK * n, 2_IK * m + n, 2_IK * npt_loc))

! Normalize the constraints, and copy the resultant constraint matrix and right hand sides into
! working space, after increasing the right hand sides if necessary so that the starting point
! is feasible.
iamat = max(m + 3 * n, 2 * m + n, 2 * npt_loc) + 1
ib = iamat + m * n
constr_modified = .false.; 
smallx = 1.0E-6_RP * rhoend_loc
if (m > 0) then
    iw = iamat - 1
    do j = 1, m
        sum = ZERO
        temp = ZERO
        do i = 1, n
            sum = sum + A_loc(i, j) * x(i)
            temp = temp + A_loc(i, j)**2
        end do
        if (temp <= 0) then
            if (present(info)) then
                info = 12
            end if
            return
        end if
        temp = sqrt(temp)
        constr_modified = constr_modified .or. (sum - b_loc(j) > smallx * temp)
        w(ib + j - 1) = max(b_loc(j), sum) / temp
        do i = 1, n
            iw = iw + 1
            w(iw) = A_loc(i, j) / temp
        end do
    end do
end if

ndim = npt_loc + n
ixb = ib + m
ixp = ixb + n
ifv = ixp + n * npt_loc
ixs = ifv + npt_loc
ixo = ixs + n
igo = ixo + n
ihq = igo + n
ipq = ihq + (n * (n + 1)) / 2
ibmat = ipq + npt_loc
izmat = ibmat + ndim * n
istp = izmat + npt_loc * (npt_loc - n - 1)
isp = istp + n
ixn = isp + npt_loc + npt_loc
iac = ixn + n
irc = iac + n
iqf = irc + m
irf = iqf + n * n
ipqw = irf + (n * (n + 1)) / 2

call safealloc(iact, m)
call lincob(calfun, n, npt_loc, m, w(iamat), w(ib), x, rhobeg_loc, rhoend_loc, iprint_loc, &
& maxfun_loc, w(ixb), w(ixp), w(ifv), w(ixs), w(ixo), w(igo), w(ihq), &
& w(ipq), w(ibmat), w(izmat), ndim, w(istp), w(isp), w(ixn), iact, &
& w(irc), w(iqf), w(irf), w(ipqw), w, f, info_loc, ftarget_loc, &
& A_loc, b_loc, cstrv_loc, nf_loc, &
& xhist_loc, int(size(xhist_loc, 2), kind=IK), fhist_loc, int(size(fhist_loc), kind=IK), chist_loc, &
& int(size(chist_loc), kind=IK))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Deallocate A_LOC and B_LOC. Indeed, automatic allocation will take place at exit.
deallocate (A_loc)
deallocate (b_loc)


! Write the outputs.

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
    ! XHIST cointains only valid history. For this reason, there is no way to avoid allocating
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

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist) .or. present(chist)) .and. maxhist_loc < nf_loc) then
    write (wmsg, ifmt) maxhist_loc
    call warning(solver, 'Only the history of the last '//trim(wmsg)//' iteration(s) is recoreded')
end if

end subroutine lincoa


end module lincoa_mod
