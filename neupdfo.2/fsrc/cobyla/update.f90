module update_mod

implicit none
private
public :: updatexfc, updatepole, findpole


contains


subroutine updatexfc(jdrop, constr, cstrv, d, f, conmat, cval, fval, sim, simi)
! Revise the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and CVAL.
! N.B.: UPDATEXFC does NOT manipulate the simplex so that SIM(:, N+1) is the best vertex;
! that is the job of UPDATEPOLE, which is called before each trust-region/geometry step.

! Generic modules
use consts_mod, only : IK, RP, DEBUGGING
use lina_mod, only : matprod, inprod, outprod

implicit none

! Input
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
real(RP) :: simi_jdrop(size(simi, 2))

sim(:, jdrop) = d
simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d)
simi = simi - outprod(matprod(simi, d), simi_jdrop)
simi(jdrop, :) = simi_jdrop
fval(jdrop) = f
conmat(:, jdrop) = constr
cval(jdrop) = cstrv

end subroutine updatexfc


subroutine updatepole(cpen, evaluated, conmat, cval, fval, sim, simi, info)
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), namely the
! "Pole Position" as per Powell's code. CONMAT, CVAL, FVAL, and SIMI are updated accordingly.
! N.B.: In precise arithmetic, the following two procedures produce the same results:
! 1. apply UPDATEPOLE to SIM twice the first time with CPEN = CPEN1 and the second with CPEN = CPEN2;
! 2. apply UPDATEPOLE to SIM with CPEN = CPEN2.
! In finite-precision arithmetic, however, they may produce different results unless CPEN1 = CPEN2.

use consts_mod, only : IK, RP, ZERO, TENTH, DEBUGGING
use info_mod, only : DAMAGING_ROUNDING
use debug_mod, only : errstop, verisize
use infnan_mod, only : is_nan
use lina_mod, only : matprod, eye

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
    erri = matprod(simi, sim(:, 1:n)) - eye(n, n)
    if (any(is_nan(erri)) .or. maxval(abs(erri)) > TENTH) then
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

use consts_mod, only : IK, RP, ZERO, DEBUGGING
use debug_mod, only : errstop

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
