! INITIALIZE_MOD is a module containing subroutine(s) for initialization.
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and the COBYLA paper.
!
! Last Modified: Friday, August 27, 2021 PM03:40:44

module initialize_mod

implicit none
private
public :: initxfc, initfilt


contains

subroutine initxfc(iprint, maxfun, ctol, ftarget, rho, x0, nf, chist, conhist, conmat, cval, fhist, &
    & fval, sim, xhist, evaluated, info)

! Generic modules
use consts_mod, only : RP, IK, ZERO, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, NAN_X, NAN_INF_F
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop, verisize
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : eye

! Solver-specific modules
use history_mod, only : savehist

implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ctol
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rho
real(RP), intent(in) :: x0(:)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: chist(:)
real(RP), intent(out) :: conhist(:, :)
real(RP), intent(out) :: conmat(:, :)
real(RP), intent(out) :: cval(:)
real(RP), intent(out) :: fhist(:)
real(RP), intent(out) :: fval(:)
real(RP), intent(out) :: sim(:, :)
real(RP), intent(out) :: xhist(:, :)
logical, intent(out) :: evaluated(:)

! Local variables
integer(IK) :: j
integer(IK) :: k
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
real(RP) :: con(size(conmat, 1))
real(RP) :: cstrv
real(RP) :: f
real(RP) :: x(size(x0))
character(len=*), parameter :: srname = 'INITIALIZE'

! Get and verify the sizes.
m = size(conmat, 1)
n = size(sim, 1)
maxchist = size(chist)
maxconhist = size(conhist, 2)
maxfhist = size(fhist)
maxxhist = size(xhist, 2)
maxhist = max(maxchist, maxconhist, maxfhist, maxxhist)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(SIM) is invalid')
    end if
    call verisize(conmat, m, n + 1)
    call verisize(cval, n + 1)
    call verisize(fval, n + 1)
    call verisize(sim, n, n + 1)
    call verisize(evaluated, n + 1)
    if (maxchist > 0) then
        call verisize(chist, maxhist)
    end if
    if (maxconhist > 0) then
        call verisize(conhist, m, maxhist)
    end if
    if (maxfhist > 0) then
        call verisize(fhist, maxhist)
    end if
    if (maxxhist > 0) then
        call verisize(xhist, n, maxhist)
    end if
end if

sim = rho * eye(n, n + 1)
sim(:, n + 1) = x0
fval = HUGENUM
cval = HUGENUM
info = 0_IK
! EVALUATED(J) = TRUE iff the function/constraint of SIM(:, J) has been evaluated.
evaluated = .false.

do k = 1, n + 1
    x = sim(:, n + 1)
    ! We will evaluate F corresponding to SIM(:, J).
    if (k == 1) then
        j = n + 1
    else
        j = k - 1
        x(j) = x(j) + rho
    end if

    if (any(is_nan(x))) then
        ! Set F and CON to NaN. This is necessary if the initial X contains NaN.
        f = sum(x)
        con = f
    else
        call calcfc(n, m, x, f, con)  ! Evaluate F and CON.
    end if
    evaluated(j) = .true.
    if (any(is_nan(con))) then
        cstrv = sum(con)  ! Set CSTRV to NaN.
    else
        cstrv = maxval([-con, ZERO])  ! Constraint violation for constraints CON(X) >= 0.
    end if
    fval(j) = f
    conmat(:, j) = con
    cval(j) = cstrv
    ! Save X, F, CON, CSTRV into the history.
    call savehist(k, con, cstrv, f, x, chist, conhist, fhist, xhist)
    ! Check whether to exit.
    if (any(is_nan(x))) then
        info = NAN_X
        exit
    end if
    if (is_nan(f) .or. is_posinf(f) .or. is_nan(cstrv) .or. is_posinf(cstrv)) then
        info = NAN_INF_F
        exit
    end if
    if (f <= ftarget .and. cstrv <= ctol) then
        info = FTARGET_ACHIEVED
        exit
    end if
    if (k >= maxfun) then
        info = MAXFUN_REACHED
        exit
    end if

    ! Exchange the new vertex of the initial simplex with the optimal vertex if necessary.
    ! This is the ONLY part that is essentially non-parallel.
    if (j <= n .and. fval(j) < fval(n + 1)) then
        fval([j, n + 1]) = fval([n + 1, j])
        cval([j, n + 1]) = cval([n + 1, j])
        conmat(:, [j, n + 1]) = conmat(:, [n + 1, j])
        sim(:, n + 1) = x
        sim(j, 1:j) = -rho
    end if

end do

nf = int(count(evaluated), kind(nf))

! It is unnecessary to call UPDATEPOLE in the end, as long as we ensure the following.
! 1. UPDATEPOLE is called at the beginning of a trust-region iteration.
! 2. SELECTX is called before the possible exit after initialization (due to errors like NAN_X).

end subroutine initxfc


subroutine initfilt(conmat, ctol, cval, fval, sim, evaluated, nfilt, cfilt, confilt, ffilt, xfilt)
! This subroutine initializes the filters (XFILT, etc) that will be used when selecting X at the
! end of the solver.
! N.B.:
! 1. Why not initialize the filters using XHIST, etc? Because the history is empty if the user
! chooses not to output it.
! 2. We decouple INITFILT with INITXFC so that it is easier to parallelize the latter if needed.
use consts_mod, only : RP, IK, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
use hist_mod, only : savefilt
implicit none

! Inputs
real(RP), intent(in) :: conmat(:, :)
real(RP), intent(in) :: ctol
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
integer(IK) :: i
integer(IK) :: m
integer(IK) :: maxfilt
integer(IK) :: n
character(len=*), parameter :: srname = 'INITFILT'

! Get and verify the sizes.
m = size(conmat, 1)
n = size(sim, 1)
maxfilt = size(ffilt)
if (DEBUGGING) then
    if (n == 0) then
        call errstop(srname, 'SIZE(SIM, 1) == 0')
    end if
    if (maxfilt == 0) then
        call errstop(srname, 'SIZE(FFILT) == 0')
    end if
    call verisize(conmat, m, n + 1)
    call verisize(cval, n + 1)
    call verisize(fval, n + 1)
    call verisize(sim, n, n + 1)
    call verisize(evaluated, n + 1)
    call verisize(confilt, m, maxfilt)
    call verisize(cfilt, maxfilt)
    call verisize(xfilt, n, maxfilt)
end if

nfilt = 0_IK
do i = 1, n
    if (evaluated(i)) then
        call savefilt(conmat(:, i), cval(i), ctol, fval(i), sim(:, i) + sim(:, n + 1), nfilt, cfilt, confilt, ffilt, xfilt)
    end if
end do
if (evaluated(n + 1)) then
    call savefilt(conmat(:, n + 1), cval(n + 1), ctol, fval(n + 1), sim(:, n + 1), nfilt, cfilt, confilt, ffilt, xfilt)
end if

end subroutine initfilt


end module initialize_mod
