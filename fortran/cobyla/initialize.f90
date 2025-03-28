module initialize_cobyla_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines for initialization.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Saturday, March 16, 2024 PM04:48:11
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initxfc, initfilt


contains


subroutine initxfc(calcfc, iprint, maxfun, constr0, amat, bvec, ctol, f0, ftarget, rhobeg, x0, nf, chist, &
    & conhist, conmat, cval, fhist, fval, sim, simi, xhist, evaluated, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine does the initialization concerning X, function values, and constraints.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TENTH, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : eye, inv, isinv, maximum
use, non_intrinsic :: message_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJCON

implicit none

! Inputs
procedure(OBJCON) :: calcfc ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: constr0(:)  ! CONSTR0(M)
real(RP), intent(in) :: amat(:, :)
real(RP), intent(in) :: bvec(:)
real(RP), intent(in) :: ctol
real(RP), intent(in) :: f0
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)   ! X0(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
logical, intent(out) :: evaluated(:)    ! EVALUATED(N+1)
real(RP), intent(out) :: chist(:)   ! CHIST(MAXCHIST)
real(RP), intent(out) :: conhist(:, :)  ! CONHIST(M, MAXCONHIST)
real(RP), intent(out) :: conmat(:, :)   ! CONMAT(M, N+1)
real(RP), intent(out) :: cval(:)    ! CVAL(N+1)
real(RP), intent(out) :: fhist(:)   ! FHIST(MAXFHIST)
real(RP), intent(out) :: fval(:)    ! FVAL(N+1)
real(RP), intent(out) :: sim(:, :)  ! SIM(N, N+1)
real(RP), intent(out) :: simi(:, :) ! SIMI(N, N)
real(RP), intent(out) :: xhist(:, :)! XHIST(N, MAXXHIST)

! Local variables
character(len=*), parameter :: solver = 'COBYLA'
character(len=*), parameter :: srname = 'INITIALIZE'
integer(IK) :: j
integer(IK) :: k
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: subinfo
real(RP) :: constr(size(conmat, 1))
real(RP) :: cstrv
real(RP) :: f
real(RP) :: x(size(x0))
real(RP), parameter :: itol = TENTH

! Sizes
m = int(size(conmat, 1), kind(m))
n = int(size(sim, 1), kind(n))
maxchist = int(size(chist), kind(maxchist))
maxconhist = int(size(conhist, 2), kind(maxconhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxhist = int(max(maxchist, maxconhist, maxfhist, maxxhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(size(cval) == n + 1, 'SIZE(CVAL) == N+1', srname)
    call assert(size(fval) == n + 1, 'SIZE(FVAL) == N+1', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(size(evaluated) == n + 1, 'SIZE(EVALUATED) == N + 1', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    call assert(size(conhist, 1) == m .and. maxconhist * (maxconhist - maxhist) == 0, &
        & 'SIZE(CONHIST, 1) == M, SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(size(x0) == n .and. all(is_finite(x0)), 'SIZE(X0) == N, X0 is finite', srname)
    call assert(rhobeg > 0, 'RHOBEG > 0', srname)
end if

!====================!
! Calculation starts !
!====================!

! Initialize INFO to the default value. At return, an INFO different from this value will indicate
! an abnormal return.
info = INFO_DFT

! Initialize the simplex. It will be revised during the initialization.
sim = eye(n, n + 1_IK) * rhobeg
sim(:, n + 1) = x0

! Initialize the matrix SIMI. This initial value will be discarded at the end of the initialization.
! If we do not do this, compilers may complain if we return due to CHECKEXIT before SIMI is set.
simi = eye(n) / rhobeg

! EVALUATED(J) = TRUE iff the function/constraint of SIM(:, J) has been evaluated.
evaluated = .false.

! Initialize XHIST, FHIST, CHIST, CONHIST, FVAL, CVAL, and CONMAT. Otherwise, compilers may complain
!that they are not (completely) initialized if the initialization aborts due to abnormality (see
!CHECKEXIT).
! N.B.: 1. Initializing them to NaN would be more reasonable (NaN is not available in Fortran).
! 2. Do not initialize the models if the current initialization aborts due to abnormality. Otherwise,
! errors or exceptions may occur, as FVAL and XPT etc are uninitialized.
xhist = -REALMAX
fhist = REALMAX
chist = REALMAX
conhist = REALMAX
fval = REALMAX
cval = REALMAX
conmat = REALMAX

do k = 1, n + 1_IK
    x = sim(:, n + 1)
    ! We will evaluate F corresponding to SIM(:, J).
    if (k == 1) then
        j = n + 1_IK
        f = f0
        constr = constr0
    else
        j = k - 1_IK
        x(j) = x(j) + rhobeg
        call evaluate(calcfc, x, f, constr, amat, bvec)
    end if
    cstrv = maximum([ZERO, constr])

    ! Print a message about the function/constraint evaluation according to IPRINT.
    call fmsg(solver, 'Initialization', iprint, k, rhobeg, f, x, cstrv, constr)
    ! Save X, F, CONSTR, CSTRV into the history.
    call savehist(k, x, xhist, f, fhist, cstrv, chist, constr, conhist)

    ! Save F, CONSTR, and CSTRV to FVAL, CONMAT, and CVAL respectively. This must be done before
    ! checking whether to exit. If exit, FVAL, CONMAT, and CVAL will define FFILT, CONFILT, and
    ! CFILT, which will define the returned X, F, CONSTR, and CSTRV.
    evaluated(j) = .true.
    fval(j) = f
    conmat(:, j) = constr
    cval(j) = cstrv

    ! Check whether to exit.
    subinfo = checkexit(maxfun, k, cstrv, ctol, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if

    ! Exchange the new vertex of the initial simplex with the optimal vertex if necessary.
    ! This is the ONLY part that is essentially non-parallel.
    if (j <= n .and. fval(j) < fval(n + 1)) then
        fval([j, n + 1_IK]) = fval([n + 1_IK, j])
        cval([j, n + 1_IK]) = cval([n + 1_IK, j])
        conmat(:, [j, n + 1_IK]) = conmat(:, [n + 1_IK, j])
        sim(:, n + 1) = x
        sim(j, 1:j) = -rhobeg  ! SIM(:, 1:N) is lower triangular.
    end if
end do

nf = int(count(evaluated), kind(nf))

if (all(evaluated)) then
    ! Initialize SIMI to the inverse of SIM(:, 1:N).
    simi = inv(sim(:, 1:n))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nf <= maxfun, 'NF <= MAXFUN', srname)
    call assert(size(evaluated) == n + 1, 'SIZE(EVALUATED) == N + 1', srname)
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) &
        & .or. is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain negative values or NaN/+Inf', srname)
    call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxconhist, &
        & 'SIZE(CONHIST) == [M, MAXCONHIST]', srname)
    call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
        & is_posinf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/+Inf', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
        & 'FHIST does not contain NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL does not contain NaN/+Inf', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(sum(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol) .or. any(.not. evaluated), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

end subroutine initxfc


subroutine initfilt(conmat, ctol, cweight, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)
!--------------------------------------------------------------------------------------------------!
! This subroutine initializes the filters (XFILT, etc) that will be used when selecting X at the
! end of the solver.
! N.B.:
! 1. Why not initialize the filters using XHIST, etc? Because the history is empty if the user
! chooses not to output it.
! 2. We decouple INITXFC and INITFILT so that it is easier to parallelize the former if needed.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: selectx_mod, only : savefilt
implicit none

! Inputs
real(RP), intent(in) :: conmat(:, :)
real(RP), intent(in) :: ctol
real(RP), intent(in) :: cweight
real(RP), intent(in) :: cval(:)
real(RP), intent(in) :: fval(:)
real(RP), intent(in) :: sim(:, :)
logical, intent(in) :: evaluated(:)

! In-outputs
integer(IK), intent(inout) :: nfilt
real(RP), intent(inout) :: cfilt(:)
real(RP), intent(inout) :: confilt(:, :)
real(RP), intent(inout) :: ffilt(:)
real(RP), intent(inout) :: xfilt(:, :)

! Local variables
character(len=*), parameter :: srname = 'INITFILT'
integer(IK) :: i
integer(IK) :: m
integer(IK) :: maxfilt
integer(IK) :: n
real(RP) :: x(size(sim, 1))

! Sizes
m = int(size(conmat, 1), kind(m))
n = int(size(sim, 1), kind(n))
maxfilt = int(size(ffilt), kind(maxfilt))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(maxfilt >= 1, 'MAXFILT >= 1', srname)
    call assert(size(confilt, 1) == m .and. size(confilt, 2) == maxfilt, 'SIZE(CONFILT) == [M, MAXFILT]', srname)
    call assert(size(cfilt) == maxfilt, 'SIZE(CFILT) == MAXFILT', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == maxfilt, 'SIZE(XFILT) == [N, MAXFILT]', srname)
    call assert(size(ffilt) == maxfilt, 'SIZE(FFILT) == MAXFILT', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL does not contain NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(sum(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(evaluated) == n + 1, 'SIZE(EVALUATED) == N + 1', srname)
end if

!====================!
! Calculation starts !
!====================!

nfilt = 0
do i = 1, n + 1_IK
    if (evaluated(i)) then
        if (i <= n) then
            x = sim(:, i) + sim(:, n + 1)
        else
            x = sim(:, i)  ! I == N+1
        end if
        call savefilt(cval(i), ctol, cweight, fval(i), x, nfilt, cfilt, ffilt, xfilt, conmat(:, i), confilt)
    end if
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nfilt <= maxfilt, 'NFILT <= MAXFILT', srname)
    call assert(size(confilt, 1) == m .and. size(confilt, 2) == maxfilt, 'SIZE(CONFILT) == [M, MAXFILT]', srname)
    call assert(.not. any(is_nan(confilt(:, 1:nfilt)) .or. is_posinf(confilt(:, 1:nfilt))), &
        & 'CONFILT does not contain NaN/+Inf', srname)
    call assert(size(cfilt) == maxfilt, 'SIZE(CFILT) == MAXFILT', srname)
    call assert(.not. any(cfilt(1:nfilt) < 0 .or. is_nan(cfilt(1:nfilt)) .or. is_posinf(cfilt(1:nfilt))), &
        & 'CFILT does not contain negative values or NaN/Inf', srname)
    call assert(size(xfilt, 1) == n .and. size(xfilt, 2) == maxfilt, 'SIZE(XFILT) == [N, MAXFILT]', srname)
    call assert(.not. any(is_nan(xfilt(:, 1:nfilt))), 'XFILT does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(size(ffilt) == maxfilt, 'SIZE(FFILT) == MAXFILT', srname)
    call assert(.not. any(is_nan(ffilt(1:nfilt)) .or. is_posinf(ffilt(1:nfilt))), &
        & 'FFILT does not contain NaN/+Inf', srname)
end if
end subroutine initfilt


end module initialize_cobyla_mod
