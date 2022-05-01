module bobyqa_mod
!--------------------------------------------------------------------------------------------------!
! Classical mode. Not maintained. Not recommended. Please use the modernized version instead.
!
! The usage is the same as the modernized version.
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: bobyqa


contains


subroutine bobyqa(calfun, x, f, &
    & xl, xu, &
    & nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, eta1, eta2, gamma1, gamma2, &
    & xhist, fhist, maxhist, honour_x0, info)

! Generic modules
use, non_intrinsic :: consts_mod, only : DEBUGGING
use, non_intrinsic :: consts_mod, only : MAXFUN_DIM_DFT, IPRINT_DFT
use, non_intrinsic :: consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, EPS, HUGENUM, MSGLEN
use, non_intrinsic :: debug_mod, only : assert, warning
use, non_intrinsic :: evaluate_mod, only : moderatex
use, non_intrinsic :: history_mod, only : prehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
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
integer(IK), intent(in), optional :: maxfun
integer(IK), intent(in), optional :: maxhist
integer(IK), intent(in), optional :: npt
logical, intent(in), optional :: honour_x0
real(RP), intent(in), optional :: eta1
real(RP), intent(in), optional :: eta2
real(RP), intent(in), optional :: ftarget
real(RP), intent(in), optional :: gamma1
real(RP), intent(in), optional :: gamma2
real(RP), intent(in), optional :: rhobeg
real(RP), intent(in), optional :: rhoend
real(RP), intent(in), optional :: xl(:)
real(RP), intent(in), optional :: xu(:)

! Optional outputs
integer(IK), intent(out), optional :: info
integer(IK), intent(out), optional :: nf
real(RP), intent(out), allocatable, optional :: fhist(:)
real(RP), intent(out), allocatable, optional :: xhist(:, :)

! Local variables
character(len=*), parameter :: ifmt = '(I0)'  ! I0: use the minimum number of digits needed to print
character(len=*), parameter :: solver = 'LINCOA'
character(len=*), parameter :: srname = 'LINCOA'
character(len=MSGLEN) :: wmsg
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfun_loc
integer(IK) :: maxhist_loc
integer(IK) :: n
integer(IK) :: nf_loc
integer(IK) :: npt_loc
integer(IK) :: nhist
real(RP) :: eta1_loc
real(RP) :: eta2_loc
real(RP) :: ftarget_loc
real(RP) :: gamma1_loc
real(RP) :: gamma2_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP) :: xl_loc(size(x))
real(RP) :: xu_loc(size(x))
real(RP), allocatable :: fhist_loc(:)
real(RP), allocatable :: xhist_loc(:, :)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Working variables
real(RP) :: temp
real(RP), allocatable :: w(:)
integer(IK) :: ndim, ixb, ixp, ifv, ixo, igo, ihq, ipq, ibmat, izmat, isl, isu, ixn, ixa, id, ivl, iw, &
    & j, jsl, jsu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sizes
n = int(size(x), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    if (present(xl)) then
        call assert(size(xl) == n .or. size(xl) == 0, 'SIZE(XL) == N unless XL is empty', srname)
    end if
    if (present(xu)) then
        call assert(size(xu) == n .or. size(xu) == 0, 'SIZE(XU) == N unless XU is empty', srname)
    end if
end if

! Read the inputs

if (present(xl)) then
    xl_loc = xl
else
    xl_loc = -HUGENUM
end if
where (is_nan(xl_loc))
    xl_loc = -HUGENUM
end where

if (present(xu)) then
    xu_loc = xu
else
    xu_loc = HUGENUM
end if
where (is_nan(xu_loc))
    xu_loc = HUGENUM
end where

x = max(xl_loc, min(xu_loc, moderatex(x)))

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

! Preprocess the inputs in case some of them are invalid. It does nothing if all inputs are valid.
call preproc(solver, n, iprint_loc, maxfun_loc, maxhist_loc, ftarget_loc, rhobeg_loc, rhoend_loc, &
    & npt=npt_loc, eta1=eta1_loc, eta2=eta2_loc, gamma1=gamma1_loc, gamma2=gamma2_loc, xl = xl_loc, xu = xu_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!! Revise X (see below) and RHOBEG, RHOEND
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! XU and XL may change a bit (due to rounding ???). It was oberved in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOBEG = MIN(XU-XL)/2 in the interfaces, then it may happen
! that RHOBEG > MIN(XU-XL)/2. That is why we do the following. After
! this, INFO=6 should never occur.
rhobeg_loc = min(0.5D0 * (1.0D0 - 1.0D-5) * minval(xu_loc(1:n) - xl_loc(1:n)), rhobeg_loc)
! For the same reason, we ensure RHOEND <= RHOBEG by the following.
rhoend_loc = min(rhobeg_loc, rhoend_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Further revise MAXHIST_LOC according to MAXMEMORY, and allocate memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
call prehist(maxhist_loc, n, present(xhist), xhist_loc, present(fhist), fhist_loc)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Working space
call safealloc(w, int((npt_loc + 5) * (npt_loc + n) + 3 * n * (n + 5) / 2, IK))

!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
ndim = npt + n
ixb = 1
ixp = ixb + n
ifv = ixp + n * npt
ixo = ifv + npt
igo = ixo + n
ihq = igo + n
ipq = ihq + (n * (n + 1)) / 2
ibmat = ipq + npt
izmat = ibmat + ndim * n
isl = izmat + npt * (npt - n - 1)
isu = isl + n
ixn = isu + n
ixa = ixn + n
id = ixa + n
ivl = id + n
iw = ivl + ndim
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
do j = 1, n
    temp = xu_loc(j) - xl_loc(j)
    if (temp < rhobeg_loc + rhobeg_loc) then
        if (present(info)) then
            info = 6
        end if
        return
    end if
    jsl = isl + j - 1
    jsu = jsl + n
    w(jsl) = xl_loc(j) - x(j)
    w(jsu) = xu_loc(j) - x(j)
    if (w(jsl) >= -rhobeg_loc) then
        if (w(jsl) >= ZERO) then
            x(j) = xl_loc(j)
            w(jsl) = ZERO
            w(jsu) = temp
        else
            x(j) = xl_loc(j) + rhobeg_loc
            w(jsl) = -rhobeg_loc
            w(jsu) = max(xu_loc(j) - x(j), rhobeg_loc)
        end if
    else if (w(jsu) <= rhobeg_loc) then
        if (w(jsu) <= ZERO) then
            x(j) = xu_loc(j)
            w(jsl) = -temp
            w(jsu) = ZERO
        else
            x(j) = xu_loc(j) - rhobeg_loc
            w(jsl) = min(xl_loc(j) - x(j), -rhobeg_loc)
            w(jsu) = rhobeg_loc
        end if
    end if
end do
!
!     Make the call of BOBYQB.
!
call bobyqb(calfun, n, npt_loc, x, xl_loc, xu_loc, rhobeg_loc, rhoend_loc, iprint_loc, maxfun_loc, &
    & w(ixb), w(ixp), w(ifv), w(ixo), w(igo), w(ihq), w(ipq), w(ibmat), w(izmat), &
& ndim, w(isl), w(isu), w(ixn), w(ixa), w(id), w(ivl), w(iw), f, info_loc, ftarget_loc, &
& nf_loc, xhist_loc, int(size(xhist_loc, 2), kind=IK), fhist_loc, int(size(fhist_loc), kind=IK))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate (w)

! Write the outputs.

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

! If NF_LOC > MAXHIST_LOC, warn that not all history is recorded.
if ((present(xhist) .or. present(fhist)) .and. maxhist_loc < nf_loc) then
    write (wmsg, ifmt) maxhist_loc
    call warning(solver, 'Only the history of the last '//trim(wmsg)//' iteration(s) is recorded')
end if

end subroutine bobyqa


end module bobyqa_mod
