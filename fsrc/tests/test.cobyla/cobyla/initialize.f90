! INITIALIZE_MOD is a module containing subroutine(s) for initialization.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Started: July 2021
!
! Last Modified: Wednesday, September 22, 2021 AM11:55:41

module initialize_mod

implicit none
private
public :: initxfc, initfilt


contains

subroutine initxfc(calcfc, iprint, maxfun, ctol, ftarget, rho, x0, nf, chist, conhist, conmat, cval, fhist, &
    & fval, sim, xhist, evaluated, info)

! Generic modules
use, non_intrinsic :: pintrf_mod, only : FUNCON
use, non_intrinsic :: evaluate_mod, only : evalfc
use, non_intrinsic :: consts_mod, only : RP, IK, HUGENUM, DEBUGGING
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: debug_mod, only : errstop, verisize
use, non_intrinsic :: output_mod, only : retmssg, rhomssg, fmssg
use, non_intrinsic :: linalg_mod, only : eye
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: checkexit_mod, only : checkexit

implicit none

! Inputs
procedure(FUNCON) :: calcfc
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
integer(IK) :: subinfo
real(RP) :: constr(size(conmat, 1))
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
info = INFO_DFT
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
    call evalfc(calcfc, x, f, constr, cstrv)
    evaluated(j) = .true.
    ! Save X, F, CONSTR, CSTRV into the history.
    call savehist(k, constr, cstrv, f, x, chist, conhist, fhist, xhist)
    ! Save F, CONSTR, and CSTRV to FVAL, CONMAT, and CVAL respectively. This must be done before
    ! checking whether to exit. If exit, FVAL, CONMAT, and CVAL will define FFILT, CONFILT, and
    ! CFILT, which will define the returned X, F, CONSTR, and CSTRV.
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
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop, verisize
use, non_intrinsic :: selectx_mod, only : savefilt
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
