module update_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Started: July 2021
!
! Last Modified: Thursday, December 02, 2021 PM03:57:11
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
use, non_intrinsic :: linalg_mod, only : matprod, inprod, outprod, isinv
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
    call assert(jdrop >= 0 .and. jdrop <= n, '1 <= JDROP <= N', srname)
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
    !!!call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do nothing when JDROP is 0. This can only happen after a trust-region step.
if (jdrop <= 0) then  ! JDROP < 0 is impossible if the input is correct.
    return
end if

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
    !!!call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if
end subroutine updatexfc


subroutine updatepole(cpen, conmat, cval, fval, sim, simi, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), namely the
! "Pole Position" as per Powell's code. CONMAT, CVAL, FVAL, and SIMI are updated accordingly.
! N.B.: In precise arithmetic, the following two procedures produce the same results:
! 1. apply UPDATEPOLE to SIM twice the first time with CPEN = CPEN1 and the second with CPEN = CPEN2;
! 2. apply UPDATEPOLE to SIM with CPEN = CPEN2.
! In finite-precision arithmetic, however, they may produce different results unless CPEN1 = CPEN2.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, TENTH, DEBUGGING
use, non_intrinsic :: info_mod, only : DAMAGING_ROUNDING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_neginf, is_posinf, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, eye, inv, isinv

implicit none

! Inputs
real(RP), intent(in) :: cpen

! In-outputs
real(RP), intent(inout) :: conmat(:, :)
real(RP), intent(inout) :: cval(:)
real(RP), intent(inout) :: fval(:)
real(RP), intent(inout) :: sim(:, :)
real(RP), intent(inout) :: simi(:, :)

! Local variables
character(len=*), parameter :: srname = 'UPDATEPOLE'
integer(IK) :: info
integer(IK) :: jopt
integer(IK) :: m
integer(IK) :: n
real(RP) :: conmat_old(size(conmat, 1), size(conmat, 2))
real(RP) :: cval_old(size(cval))
real(RP) :: erri(size(sim, 1), size(sim, 1))
real(RP) :: fval_old(size(fval))
real(RP) :: sim_jopt(size(sim, 1))
real(RP) :: sim_old(size(sim, 1), size(sim, 2))
real(RP) :: simi_old(size(simi, 1), size(simi, 2))
real(RP), parameter :: itol = TENTH

! Sizes
m = size(conmat, 1)
n = size(sim, 1)

! Preconditions
if (DEBUGGING) then
    call assert(cpen >= 0, 'CPEN >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
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
    !!!call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

info = 0_IK

! Identify the optimal vertex of the current simplex.
jopt = findpole(cpen, cval, fval)

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

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(findpole(cpen, cval, fval) == n + 1 .or. info == DAMAGING_ROUNDING, &
        & 'The best point is SIM(:, N+1) unless the rounding is damaging', srname)
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
    ! Do not check SIMI = SIM(:, 1:N)^{-1}, as it may not be true due to damaging rounding.
    !call assert(isinv(sim(:, 1:n), simi, itol) .or. , 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

end subroutine updatepole


function findpole(cpen, cval, fval) result(jopt)
!--------------------------------------------------------------------------------------------------!
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf

implicit none

! Inputs
real(RP), intent(in) :: cpen
real(RP), intent(inout) :: cval(:)
real(RP), intent(inout) :: fval(:)

! Outputs
integer(IK) :: jopt

! Local variables
character(len=*), parameter :: srname = 'FINDPOLE'
integer(IK) :: n
real(RP) :: phi(size(cval))
real(RP) :: phimin

! Size
n = int(size(fval) - 1, kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(cpen >= 0, 'CPEN >= 0', srname)
    call assert(size(cval) == n + 1 .and. .not. any(is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL is not NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! Identify the optimal vertex of the current simplex.
jopt = size(fval) ! We use N + 1 as the default value of JOPT.
phi = fval + cpen * cval
phimin = minval(phi)
if (phimin < phi(jopt)) then  ! We keep JOPT = N + 1 unless there is a strictly better choice.
    jopt = int(minloc(phi, dim=1), kind(jopt))
end if
if (cpen <= ZERO .and. any(cval < cval(jopt) .and. phi <= phimin)) then
    ! (CPEN <= ZERO) is indeed (CPEN == ZERO), and (PHI <= PHIMIN) is indeed (PHI == PHIMIN). We
    ! write them in this way to avoid equality comparison of real numbers.
    jopt = int(minloc(cval, mask=(phi <= phimin), dim=1), kind(jopt))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(jopt >= 1 .and. jopt <= n + 1, '1 <= JOPT <= N+1', srname)
end if
end function findpole

end module update_mod
