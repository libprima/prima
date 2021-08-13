module update_mod

contains

subroutine updatepole(cpen, evaluated, datmat, sim, simi, info)
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), namely the
! "Pole Position" as per Powell's code. If necessary, DATMAT of SIMI are updated accordingly.
! N.B.: In precise arithmetic, the following two procedures produce the same results:
! 1. apply UPDATEPOLE to SIM twice the first time with CPEN = CPEN1 and the second with CPEN = CPEN2;
! 2. apply UPDATEPOLE to SIM with CPEN = CPEN2.
! In finite-precision arithmetic, however, they may produce different results unless CPEN1 = CPEN2.

use consts_mod, only : IK, RP, ZERO, TENTH, DEBUGGING, SRNLEN
use info_mod, only : DAMAGING_ROUNDING
use debug_mod, only : errstop, verisize
use infnan_mod, only : is_nan
use lina_mod, only : matprod, eye

implicit none

! Input
real(RP), intent(in) :: cpen
logical, intent(in) :: evaluated(:)

! In-outputs
real(RP), intent(inout) :: datmat(:, :)
real(RP), intent(inout) :: sim(:, :)
real(RP), intent(inout) :: simi(:, :)

! Local variables
integer(IK) :: info
integer(IK) :: jopt
integer(IK) :: m
integer(IK) :: n
real(RP) :: datmat_old(size(datmat, 1), size(datmat, 2))
real(RP) :: erri(size(sim, 1), size(sim, 1))
real(RP) :: phi(size(sim, 2))
real(RP) :: phimin
real(RP) :: sim_jopt(size(sim, 1))
real(RP) :: sim_old(size(sim, 1), size(sim, 2))
real(RP) :: simi_old(size(simi, 1), size(simi, 2))
character(len=SRNLEN), parameter :: srname = 'UPDATEPOLE'


! Get and verify the sizes.
m = size(datmat, 1) - 2
n = size(datmat, 2) - 1
if (DEBUGGING) then
    if (m < 0 .or. n < 1) then
        call errstop(srname, 'SIZE(DATMAT) is invalid')
    end if
    call verisize(evaluated, n + 1)
    call verisize(sim, n, n + 1)
    call verisize(simi, n, n)
end if

info = 0

! Identify the optimal vertex of the current simplex.
! N.B.: Maybe not all vertex of the simplex are initialized! Use EVALUATED as a mask.
jopt = findpole(cpen, evaluated, datmat)

! Switch the best vertex into SIM(:, N+1) if it is not there already. Then update SIMI and DATMAT.
! Before the update, save a copy of DATMAT, SIM, and SIMI. If the update is unsuccessful due to
! damaging rounding errors, we restore them for COBYLA to extract X/F from SIM/DATMAT at exit.
datmat_old = datmat
sim_old = sim
simi_old = simi
if (jopt <= n) then
    datmat(:, [jopt, n + 1]) = datmat(:, [n + 1, jopt]) ! Exchange DATMAT(:, JOPT) AND DATMAT(:, N+1)
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
        datmat = datmat_old
        sim = sim_old
        simi = simi_old
    end if
end if

end subroutine updatepole


function findpole(cpen, evaluated, datmat) result(jopt)

use consts_mod, only : IK, RP, ZERO
use debug_mod, only : errstop, verisize

implicit none

! Input
real(RP), intent(in) :: cpen
logical, intent(in) :: evaluated(:)
real(RP), intent(inout) :: datmat(:, :)

! Output
integer(IK) :: jopt

! Local variables
integer(IK) :: m
real(RP) :: phi(size(datmat, 2))
real(RP) :: phimin

m = size(datmat, 1) - 2

! Identify the optimal vertex of the current simplex.
! N.B.: Maybe not all vertex of the simplex are initialized! Use EVALUATED as a mask.
jopt = size(datmat, 2) ! We use N + 1 as the default value of JOPT.
phi = datmat(m + 1, :) + cpen * datmat(m + 2, :)
phimin = minval(phi, mask=evaluated)
if (phimin < phi(jopt)) then  ! We keep JOPT = N + 1 unless there is a strictly better choice.
    jopt = int(minloc(phi, mask=evaluated, dim=1), kind(jopt))
end if
if (cpen <= ZERO .and. any(datmat(m + 2, :) < datmat(m + 2, jopt) .and. phi <= phimin .and. evaluated)) then
    ! (CPEN <= ZERO) is indeed (CPEN == ZERO), and (PHI <= PHIMIN) is indeed (PHI == PHIMIN). We
    ! write them in this way to avoid equality comparison of real numbers.
    jopt = int(minloc(datmat(m + 2, :), mask=(phi <= phimin .and. evaluated), dim=1), kind(jopt))
end if
end function findpole

end module update_mod
