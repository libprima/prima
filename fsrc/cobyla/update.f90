module update_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Monday, February 28, 2022 PM09:43:48
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updatexfc, updatepole, findpole


contains


subroutine updatexfc(jdrop, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine revises the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and CVAL.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, TENTH, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_neginf, is_posinf, is_finite
use, non_intrinsic :: info_mod, only : INFO_DFT, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : matprod, inprod, outprod, isinv
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
integer(IK), intent(in) :: jdrop
real(RP), intent(in) :: constr(:)   ! CONSTR(M)
real(RP), intent(in) :: cpen
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: f

! In-outputs
real(RP), intent(inout) :: conmat(:, :)     ! CONMAT(M, N+1)
real(RP), intent(inout) :: cval(:)  ! CVAL(N+1)
real(RP), intent(inout) :: fval(:)  ! FVAL(N+1)
real(RP), intent(inout) :: sim(:, :)! SIM(N, N+1)
real(RP), intent(inout) :: simi(:, :)   ! SIMI(N, N)

! Outputs
integer(IK), intent(out) :: info

! Local variables
character(len=*), parameter :: srname = 'UPDATEXFC'
integer(IK) :: m
integer(IK) :: n
real(RP), parameter :: itol = TENTH
real(RP) :: simi_jdrop(size(simi, 2))

! Sizes
m = int(size(constr), kind(m))
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
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
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do nothing when JDROP is 0. This can only happen after a trust-region step.
if (jdrop <= 0) then  ! JDROP < 0 is impossible if the input is correct.
    info = INFO_DFT  ! INFO must be set, as it is an output!!!
    return
end if

sim(:, jdrop) = d
simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d)
simi = simi - outprod(matprod(simi, d), simi_jdrop)
simi(jdrop, :) = simi_jdrop
fval(jdrop) = f
conmat(:, jdrop) = constr
cval(jdrop) = cstrv

! Switch the best vertex to the pole position SIM(:, N+1) if it is not there already.
call updatepole(cpen, conmat, cval, fval, sim, simi, info)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
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
    call assert(isinv(sim(:, 1:n), simi, itol) .or. info == DAMAGING_ROUNDING, &
       & 'SIMI = SIM(:, 1:N)^{-1} unless the rounding is damaging', srname)
end if
end subroutine updatexfc


subroutine updatepole(cpen, conmat, cval, fval, sim, simi, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), namely the
! "Pole Position" as per Powell. CONMAT, CVAL, FVAL, and SIMI are updated accordingly.
!
! N.B. 1: In precise arithmetic, the following two procedures produce the same results:
! 1) apply UPDATEPOLE to SIM twice, first with CPEN = CPEN1 and then with CPEN = CPEN2;
! 2) apply UPDATEPOLE to SIM with CPEN = CPEN2.
! In finite-precision arithmetic, however, they may produce different results unless CPEN1 = CPEN2.
!
! N.B. 2: When JOPT == N+1, the best vertex is already at the pole position, so there is nothing to
! switch. However, as in Powell's code, the code below will check whether SIMI is good enough to
! work as the inverse of SIM(:, 1:N) or not. If not, Powell's code would invoke an error return of
! COBYLB; our implementation, however, will try calculating SIMI from scratch; if the recalculated
! SIMI is still of poor quality, then UPDATEPOLE will return with INFO = DAMAGING_ROUNDING,
! informing COBYLB that SIMI is poor due to damaging rounding errors.
!
! N.B. 3: UPDATEPOLE should be called when and only when FINDPOLE can potentially returns a value
! other than N+1. The value of FINDPOLE is determined by CPEN, CVAL, and FVAL, the latter two being
! decided by SIM. Thus UPDATEPOLE should be called after CPEN or SIM changes. COBYLA updates CPEN at
! only two places: the beginning of each trust-region iteration, and when RESENHANCE is called;
! SIM is updated only by UPDATEXFC, which itself calls UPDATEPOLE internally. Therefore, we only
! need to call UPDATEPOLE after updating CPEN at the beginning of each trust-region iteration and
! after each invocation of RESENHANCE.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, TENTH, DEBUGGING
use, non_intrinsic :: info_mod, only : DAMAGING_ROUNDING, INFO_DFT
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_neginf, is_posinf, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, eye, inv, isinv

implicit none

! Inputs
real(RP), intent(in) :: cpen

! In-outputs
real(RP), intent(inout) :: conmat(:, :) ! CONMAT(M, N+1)
real(RP), intent(inout) :: cval(:)  ! CVAL(N+1)
real(RP), intent(inout) :: fval(:)  ! FVAL(N+1)
real(RP), intent(inout) :: sim(:, :)! SIM(N, N+1)
real(RP), intent(inout) :: simi(:, :)    ! SIMI(N, N)

! Outputs
integer(IK), intent(out) :: info

! Local variables
character(len=*), parameter :: srname = 'UPDATEPOLE'
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
m = int(size(conmat, 1), kind(m))
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(cpen >= 0, 'CPEN >= 0', srname)
    call assert(m >= 0, 'M >= 0', srname)
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
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

! INFO must be set, as it is an output.
info = INFO_DFT

! Identify the optimal vertex of the current simplex.
jopt = findpole(cpen, cval, fval)

! Switch the best vertex to the pole position SIM(:, N+1) if it is not there already. Then update
! CONMAT etc. Before the update, save a copy of CONMAT etc. If the update is unsuccessful due to
! damaging rounding errors, we restore them for COBYLA to extract X/F/C from the undamaged data.
fval_old = fval
conmat_old = conmat
cval_old = cval
sim_old = sim
simi_old = simi
if (jopt >= 1 .and. jopt <= n) then
    ! Unless there is a bug in FINDPOLE, it is guaranteed that JOPT >= 1.
    ! When JOPT == N + 1, there is nothing to switch; in addition, SIMI(JOPT, :) will be illegal.
    fval([jopt, n + 1_IK]) = fval([n + 1_IK, jopt])
    conmat(:, [jopt, n + 1_IK]) = conmat(:, [n + 1_IK, jopt]) ! Exchange CONMAT(:, JOPT) AND CONMAT(:, N+1)
    cval([jopt, n + 1_IK]) = cval([n + 1_IK, jopt])
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
    simi(jopt, :) = -sum(simi, dim=1)  ! Must ensure that 1 <= JOPT <= N!
end if

! Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
erri = matprod(simi, sim(:, 1:n)) - eye(n)

! Calculate SIMI from scratch if the current one is damaged by rounding errors.
if (any(is_nan(erri)) .or. any(abs(erri) > itol)) then
    simi = inv(sim(:, 1:n))
    erri = matprod(simi, sim(:, 1:n)) - eye(n)
end if

! If the recalculated SIMI is still damaged, then restore the data to the version before the update.
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
    call assert(isinv(sim(:, 1:n), simi, itol) .or. info == DAMAGING_ROUNDING, &
       & 'SIMI = SIM(:, 1:N)^{-1} unless the rounding is damaging', srname)
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
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf

implicit none

! Inputs
real(RP), intent(in) :: cpen
real(RP), intent(inout) :: cval(:)  ! CVAL(N+1)
real(RP), intent(inout) :: fval(:)  ! FVAL(N+1)

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
jopt = int(size(fval), kind(jopt))  ! We use N + 1 as the default value of JOPT.
phi = fval + cpen * cval
phimin = minval(phi)
if (phimin < phi(jopt)) then  ! We keep JOPT = N + 1 unless there is a strictly better choice.
    jopt = int(minloc(phi, dim=1), kind(jopt))
end if
if (cpen <= ZERO .and. any(cval < cval(jopt) .and. phi <= phimin)) then
    ! (CPEN <= ZERO) is indeed (CPEN == ZERO), and (PHI <= PHIMIN) is indeed (PHI == PHIMIN).
    ! We code in this way to avoid equality comparison of real numbers.
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
