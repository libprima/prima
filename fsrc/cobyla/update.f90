module update_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Started: July 2021
!
! Last Modified: Friday, November 19, 2021 PM06:08:34
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updatexfc, updatepole, findpole


contains


subroutine updatexfc(jdrop, constr, cstrv, d, f, conmat, cval, fval, sim, simi)
!--------------------------------------------------------------------------------------------------!
! This subroutine revises the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and CVAL.
! N.B.: UPDATEXFC does NOT manipulate the simplex so that SIM(:, N+1) is the best vertex;
! that is the job of UPDATEPOLE, which is called before each trust-region/geometry step.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, TENTH, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_neginf, is_posinf, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, inprod, outprod, eye
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
integer(IK), intent(in) :: jdrop
real(RP), intent(in) :: constr(:)
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: d(:)
real(RP), intent(in) :: f

! In-outputs
real(RP), intent(inout) :: conmat(:, :)
real(RP), intent(inout) :: cval(:)
real(RP), intent(inout) :: fval(:)
real(RP), intent(inout) :: sim(:, :)
real(RP), intent(inout) :: simi(:, :)

! Local variables
character(len=*), parameter :: srname = 'UPDATEXFC'
integer(IK) :: m
integer(IK) :: n
real(RP), parameter :: itol = TENTH
real(RP) :: simi_jdrop(size(simi, 2))

! Sizes
m = size(constr)
n = size(sim, 1)

! Preconditions
if (DEBUGGING) then
    call assert(jdrop >= 1 .and. jdrop <= n, '1 <= JDROP <= N', srname)
    call assert(.not. any(is_nan(constr) .or. is_neginf(constr)), 'CONSTR does not contain NaN/-Inf', srname)
    call assert(.not. (is_nan(cstrv) .or. is_posinf(cstrv)), 'CSTRV is not NaN/+Inf', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_neginf(conmat)), 'CONMAT does not contain NaN/-Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL is not NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(all(abs(matprod(simi, sim(:, 1:n)) - eye(n)) <= itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

sim(:, jdrop) = d
simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d)
simi = simi - outprod(matprod(simi, d), simi_jdrop)
simi(jdrop, :) = simi_jdrop
fval(jdrop) = f
conmat(:, jdrop) = constr
cval(jdrop) = cstrv

!====================!
!  Calculation ends  !
!====================!

if (DEBUGGING) then
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_neginf(conmat)), 'CONMAT does not contain NaN/-Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL is not NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(all(abs(matprod(simi, sim(:, 1:n)) - eye(n)) <= itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if
end subroutine updatexfc


subroutine updatepole(cpen, evaluated, conmat, cval, fval, sim, simi, info)
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), namely the
! "Pole Position" as per Powell's code. CONMAT, CVAL, FVAL, and SIMI are updated accordingly.
! N.B.: In precise arithmetic, the following two procedures produce the same results:
! 1. apply UPDATEPOLE to SIM twice the first time with CPEN = CPEN1 and the second with CPEN = CPEN2;
! 2. apply UPDATEPOLE to SIM with CPEN = CPEN2.
! In finite-precision arithmetic, however, they may produce different results unless CPEN1 = CPEN2.

use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, TENTH, DEBUGGING
use, non_intrinsic :: info_mod, only : DAMAGING_ROUNDING
use, non_intrinsic :: debug_mod, only : errstop, verisize
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: linalg_mod, only : matprod, eye, inv

implicit none

! Input
real(RP), intent(in) :: cpen
logical, intent(in) :: evaluated(:)

! In-outputs
real(RP), intent(inout) :: conmat(:, :)
real(RP), intent(inout) :: cval(:)
real(RP), intent(inout) :: fval(:)
real(RP), intent(inout) :: sim(:, :)
real(RP), intent(inout) :: simi(:, :)

! Local variables
integer(IK) :: info
integer(IK) :: jopt
integer(IK) :: m
integer(IK) :: n
real(RP) :: conmat_old(size(conmat, 1), size(conmat, 2))
real(RP) :: cval_old(size(cval))
real(RP) :: fval_old(size(fval))
real(RP) :: erri(size(sim, 1), size(sim, 1))
real(RP), parameter :: itol = TENTH
real(RP) :: sim_jopt(size(sim, 1))
real(RP) :: sim_old(size(sim, 1), size(sim, 2))
real(RP) :: simi_old(size(simi, 1), size(simi, 2))
character(len=*), parameter :: srname = 'UPDATEPOLE'


! Get and verify the sizes.
m = size(conmat, 1)
n = size(sim, 1)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(SIM) is invalid')
    end if
    call verisize(evaluated, n + 1)
    call verisize(conmat, m, n + 1)
    call verisize(cval, n + 1)
    call verisize(fval, n + 1)
    call verisize(sim, n, n + 1)
    call verisize(simi, n, n)
end if

info = 0_IK

! Identify the optimal vertex of the current simplex.
! N.B.: Maybe not all vertex of the simplex are initialized! Use EVALUATED as a mask.
jopt = findpole(cpen, evaluated, cval, fval)

! Switch the best vertex into SIM(:, N+1) if it is not there already. Then update CONMAT etc.
! Before the update, save a copy of CONMAT etc. If the update is unsuccessful due to damaging
! rounding errors, we restore them for COBYLA to extract X/F/C from the data before the damage.
fval_old = fval
conmat_old = conmat
cval_old = cval
sim_old = sim
simi_old = simi
if (jopt <= n) then
    fval([jopt, n + 1]) = fval([n + 1, jopt])
    conmat(:, [jopt, n + 1]) = conmat(:, [n + 1, jopt]) ! Exchange CONMAT(:, JOPT) AND CONMAT(:, N+1)
    cval([jopt, n + 1]) = cval([n + 1, jopt])
    sim(:, n + 1) = sim(:, n + 1) + sim(:, jopt)
    sim_jopt = sim(:, jopt)
    sim(:, jopt) = ZERO
    sim(:, 1:n) = sim(:, 1:n) - spread(sim_jopt, dim=2, ncopies=n)
    ! The above update is equivalent to multiply SIM(:, 1:N) from the right side by a matrix whose
    ! JOPT-th row is [-1, -1, ..., -1], while all the other rows are the same as those of the
    ! identity matrix. It is easy to check that the inverse of this matrix is itself. Therefore,
    ! SIMI should be updated by a multiplication with this matrix (i.e., its inverse) from the left
    ! side, as is done in the following line. The JOPT-th row of the updated SIMI is minus the sum
    ! of all rows of the original SIMI, whereas all the other rows remain unchanged.
    simi(jopt, :) = -sum(simi, dim=1)
end if

! Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
! Do this only if EVALUATED contains only TRUE.
if (all(evaluated)) then
    erri = matprod(simi, sim(:, 1:n)) - eye(n)
    ! Recalculate SIMI if the updated one is damaged by rounding errors.
    if (any(is_nan(erri)) .or. any(abs(erri) > itol)) then
        simi = inv(sim(:, 1:n))
        erri = matprod(simi, sim(:, 1:n)) - eye(n)
    end if
    if (any(is_nan(erri)) .or. any(abs(erri) > itol)) then
        info = DAMAGING_ROUNDING
        fval = fval_old
        conmat = conmat_old
        cval = cval_old
        sim = sim_old
        simi = simi_old
    end if
end if

end subroutine updatepole


function findpole(cpen, evaluated, cval, fval) result(jopt)

use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop

implicit none

! Input
real(RP), intent(in) :: cpen
logical, intent(in) :: evaluated(:)
real(RP), intent(inout) :: cval(:)
real(RP), intent(inout) :: fval(:)

! Output
integer(IK) :: jopt

! Local variables
real(RP) :: phi(size(cval))
real(RP) :: phimin
character(len=*), parameter :: srname = 'FINDPOLE'

if (DEBUGGING) then
    if (size(cval) /= size(fval)) then
        call errstop(srname, 'SIZE(CVAL) /= SIZE(FVAL)')
    end if
end if

! Identify the optimal vertex of the current simplex.
! N.B.: Maybe not all vertex of the simplex are initialized! Use EVALUATED as a mask.
jopt = size(cval) ! We use N + 1 as the default value of JOPT.
phi = fval + cpen * cval
phimin = minval(phi, mask=evaluated)
if (phimin < phi(jopt)) then  ! We keep JOPT = N + 1 unless there is a strictly better choice.
    jopt = int(minloc(phi, mask=evaluated, dim=1), kind(jopt))
end if
if (cpen <= ZERO .and. any(cval < cval(jopt) .and. phi <= phimin .and. evaluated)) then
    ! (CPEN <= ZERO) is indeed (CPEN == ZERO), and (PHI <= PHIMIN) is indeed (PHI == PHIMIN). We
    ! write them in this way to avoid equality comparison of real numbers.
    jopt = int(minloc(cval, mask=(phi <= phimin .and. evaluated), dim=1), kind(jopt))
end if
end function findpole

end module update_mod
