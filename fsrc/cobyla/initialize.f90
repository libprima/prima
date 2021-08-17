! INITIALIZE_MOD is a module containing subroutines for initializing
! FVAL, XBASE, XPT, GQ, HQ, PQ, IDZ, ZMAT, and BMAT.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.
!
! Last Modified: Wednesday, August 18, 2021 AM02:42:43

module initialize_mod

implicit none
private
public :: initialize


contains

subroutine initialize(iprint, maxfun, ctol, ftarget, rho, x0, nf, conmat, cval, fval, sim, simi, info)

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, TENTH, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, NAN_X, NAN_INF_F, DAMAGING_ROUNDING
use infnan_mod, only : is_nan, is_posinf, is_neginf
use debug_mod, only : errstop, verisize
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : matprod, inprod, eye

! Solver-specific modules
use update_mod, only : updatepole

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
real(RP), intent(out) :: conmat(:, :)
real(RP), intent(out) :: cval(:)
real(RP), intent(out) :: fval(:)
real(RP), intent(out) :: sim(:, :)
real(RP), intent(out) :: simi(:, :)

! Local variables
integer(IK) :: i
integer(IK) :: j
integer(IK) :: jopt
integer(IK) :: k
integer(IK) :: m
integer(IK) :: n
integer(IK) :: subinfo
real(RP) :: con(size(conmat, 1))
real(RP) :: erri(size(x0), size(x0))
real(RP) :: x(size(x0))
logical :: evaluated(size(x0) + 1)
character(len=SRNLEN), parameter :: srname = 'INITIALIZE'

! Get and verify the sizes
m = size(conmat, 1)
n = size(conmat, 2) - 1
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(CONMAT) is invalid')
    end if
    call verisize(cval, n + 1)
    call verisize(fval, n + 1)
    call verisize(sim, n, n + 1)
    call verisize(simi, n, n)
end if

sim = rho * eye(n, n + 1)
sim(:, n + 1) = x0
simi = eye(n, n) / rho

info = 0_IK

! EVALUATED(J) = TRUE iff the function/constraint of SIM(:, J) has been evaluated.
evaluated = .false.

do k = 1, n + 1
    x = sim(:, n + 1)
    if (k == 1) then
        j = n + 1
    else
        j = k - 1
        x(j) = x(j) + rho
    end if

    ! Exit if X contains NaN.
    if (any(is_nan(x))) then
        fval(j) = sum(x)  ! Set F to NaN.
        cval(j) = fval(j)  ! Set constraint violation to NaN.
        con = fval(j)  ! Set constraint values to NaN.
        info = NAN_X
        exit
    end if

    call calcfc(n, m, x, fval(j), con)
    evaluated(j) = .true.
    cval(j) = maxval([-con, ZERO])  ! Constraint violation for constraints con(x) >= 0.
    conmat(:, j) = con

!if (k == IPRINT - 1 .or. IPRINT == 3) then
!    print 70, k, F, RESMAX, (X(I), I=1, IPTEM)
!70  format(/3X, 'k =', I5, 3X, 'F =', 1PE13.6, 4X, 'MAXCV =', 1PE13.6 / 3X, 'X =', 1PE13.6, 1P4E15.6)
!    if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!80  format(1PE19.6, 1P4E15.6)
!end if

    ! Exchange the new vertex of the initial simplex with the optimal vertex if necessary.
    ! This is the ONLY part that is essentially non-parallel.
    if (j <= n .and. fval(j) < fval(n + 1)) then
        fval([j, n + 1]) = fval([n + 1, j])
        cval([j, n + 1]) = cval([n + 1, j])
        conmat(:, [j, n + 1]) = conmat(:, [n + 1, j])
        sim(:, n + 1) = x
        sim(j, 1:j) = -rho
        simi(j, 1:j) = -sum(simi(:, 1:j), dim=1)
    end if

    ! Exit if the objective function value or the constraints contain NaN/Inf.
    if (is_nan(fval(j)) .or. is_posinf(fval(j))) then
        info = NAN_INF_F
        exit
    end if
    if (any(is_nan(con) .or. is_neginf(con))) then
        cval(j) = sum(abs(con))  ! Set constraint violation to NaN or Inf
        info = NAN_INF_F
        exit
    end if
    ! Exit if FTARGET is achieved by a feasible point.
    if (fval(j) <= ftarget .and. cval(j) <= ctol) then
        info = FTARGET_ACHIEVED
        exit
    end if
    ! Exit if MAXFUN is reached.
    if (k >= maxfun) then
        info = MAXFUN_REACHED
        exit
    end if
end do

nf = int(count(evaluated), kind(nf))

! It is unnecessary to call UPDATEPOLE in the end, as long as we ensure the following.
! 1. UPDATEPOLE is called at the beginning of a trust-region iteration.
! 2. SELECTX is called before the possible exit after initialization (due to errors like NAN_X).

return
end subroutine initialize


end module initialize_mod
