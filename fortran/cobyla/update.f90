module update_cobyla_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2021
!
! Last Modified: Monday, August 07, 2023 AM03:53:59
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: updatexfc, updatepole, findpole


contains


subroutine updatexfc(jdrop, constr, cpen, cstrv, d, f, conmat, cval, fval, sim, simi, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine revises the simplex by updating the elements of SIM, SIMI, FVAL, CONMAT, and CVAL.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, ONE, TENTH, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : INFO_DFT, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : matprod, inprod, outprod, maximum, eye, inv, isinv
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
real(RP) :: erri
real(RP) :: erri_test
real(RP) :: sim_old(size(sim, 1), size(sim, 2))
real(RP) :: simi_jdrop(size(simi, 2))
real(RP) :: simi_old(size(simi, 1), size(simi, 2))
real(RP) :: simi_test(size(simi, 1), size(simi, 2))
real(RP) :: simid(size(simi, 1))
real(RP) :: sum_simi(size(simi, 2))
real(RP), parameter :: itol = ONE

! Sizes
m = int(size(constr), kind(m))
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(jdrop >= 0 .and. jdrop <= n + 1, '1 <= JDROP <= N+1', srname)
    call assert(.not. any(is_nan(constr) .or. is_posinf(constr)), 'CONSTR does not contain NaN/+Inf', srname)
    call assert(.not. (is_nan(cstrv) .or. is_posinf(cstrv)), 'CSTRV is not NaN/+Inf', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol), 'SIMI = SIM(:, 1:N)^{-1}', srname)
end if

!====================!
! Calculation starts !
!====================!

! Do nothing when JDROP is 0. This can only happen after a trust-region step.
if (jdrop <= 0) then  ! JDROP < 0 is impossible if the input is correct.
    info = INFO_DFT  ! INFO must be set, as it is an output!
    return
end if

sim_old = sim
simi_old = simi
if (jdrop <= n) then
    sim(:, jdrop) = d
    simi_jdrop = simi(jdrop, :) / inprod(simi(jdrop, :), d)
    simi = simi - outprod(matprod(simi, d), simi_jdrop)
    simi(jdrop, :) = simi_jdrop
else  ! JDROP = N+1
    sim(:, n + 1) = sim(:, n + 1) + d
    sim(:, 1:n) = sim(:, 1:n) - spread(d, dim=2, ncopies=n)
    simid = matprod(simi, d)
    sum_simi = sum(simi, dim=1)
    simi = simi + outprod(simid, sum_simi / (ONE - sum(simid)))
end if

! Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
! Calculate SIMI from scratch if the current one is damaged by rounding errors.
erri = maximum(abs(matprod(simi, sim(:, 1:n)) - eye(n)))  ! MAXIMUM(X) returns NaN if X contains NaN
if (erri > TENTH * itol .or. is_nan(erri)) then
    simi_test = inv(sim(:, 1:n))
    erri_test = maximum(abs(matprod(simi_test, sim(:, 1:n)) - eye(n)))
    if (erri_test < erri .or. (is_nan(erri) .and. .not. is_nan(erri_test))) then
        simi = simi_test
        erri = erri_test
    end if
end if

! If SIMI is satisfactory, then update FVAL, CONMAT, CVAL, and the pole position. Otherwise, restore
! SIM and SIMI, and return with INFO = DAMAGING_ROUNDING.
if (erri <= itol) then
    fval(jdrop) = f
    conmat(:, jdrop) = constr
    cval(jdrop) = cstrv
    ! Switch the best vertex to the pole position SIM(:, N+1) if it is not there already.
    call updatepole(cpen, conmat, cval, fval, sim, simi, info)
else  ! ERRI > ITOL or ERRI is NaN
    info = DAMAGING_ROUNDING
    sim = sim_old
    simi = simi_old
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
    call assert(size(simi, 1) == n .and. size(simi, 2) == n, 'SIZE(SIMI) == [N, N]', srname)
    call assert(all(is_finite(simi)), 'SIMI is finite', srname)
    call assert(isinv(sim(:, 1:n), simi, itol) .or. info == DAMAGING_ROUNDING, &
       & 'SIMI = SIM(:, 1:N)^{-1} unless the rounding is damaging', srname)
end if
end subroutine updatexfc


subroutine updatepole(cpen, conmat, cval, fval, sim, simi, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine identifies the best vertex of the current simplex with respect to the merit
! function PHI = F + CPEN * CSTRV, and then switch this vertex to SIM(:, N + 1), which Powell called
! the "pole position" in his comments. CONMAT, CVAL, FVAL, and SIMI are updated accordingly.
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
! only two places: the beginning of each trust-region iteration, and when REDRHO is called;
! SIM is updated only by UPDATEXFC, which itself calls UPDATEPOLE internally. Therefore, we only
! need to call UPDATEPOLE after updating CPEN at the beginning of each trust-region iteration and
! after each invocation of REDRHO.
!--------------------------------------------------------------------------------------------------!

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, ONE, TENTH, DEBUGGING
use, non_intrinsic :: infos_mod, only : DAMAGING_ROUNDING, INFO_DFT
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: linalg_mod, only : matprod, eye, inv, isinv, maximum

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
real(RP) :: erri
real(RP) :: erri_test
real(RP) :: sim_jopt(size(sim, 1))
real(RP) :: sim_old(size(sim, 1), size(sim, 2))
real(RP) :: simi_old(size(simi, 1), size(simi, 2))
real(RP) :: simi_test(size(simi, 1), size(simi, 2))
real(RP), parameter :: itol = ONE

! Sizes
m = int(size(conmat, 1), kind(m))
n = int(size(sim, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(cpen > 0, 'CPEN > 0', srname)
    call assert(size(conmat, 1) == m .and. size(conmat, 2) == n + 1, 'SIZE(CONMAT) = [M, N+1]', srname)
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
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

! Switch the best vertex to the pole position SIM(:, N+1) if it is not there already, and update
! SIMI. Before the update, save a copy of SIM and SIMI. If the update is unsuccessful due to
! damaging rounding errors, we restore them and return with INFO = DAMAGING_ROUNDING.
sim_old = sim
simi_old = simi
if (jopt >= 1 .and. jopt <= n) then
    ! Unless there is a bug in FINDPOLE, it is guaranteed that JOPT >= 1.
    ! When JOPT == N + 1, there is nothing to switch; in addition, SIMI(JOPT, :) will be illegal.
    sim(:, n + 1) = sim(:, n + 1) + sim(:, jopt)
    sim_jopt = sim(:, jopt)
    sim(:, jopt) = ZERO
    sim(:, 1:n) = sim(:, 1:n) - spread(sim_jopt, dim=2, ncopies=n)
    !!MATLAB: sim(:, 1:n) = sim(:, 1:n) - sim_jopt; % sim_jopt should be a column! Implicit expansion
    ! The above update is equivalent to multiply SIM(:, 1:N) from the right side by a matrix whose
    ! JOPT-th row is [-1, -1, ..., -1], while all the other rows are the same as those of the
    ! identity matrix. It is easy to check that the inverse of this matrix is itself. Therefore,
    ! SIMI should be updated by a multiplication with this matrix (i.e., its inverse) from the left
    ! side, as is done in the following line. The JOPT-th row of the updated SIMI is minus the sum
    ! of all rows of the original SIMI, whereas all the other rows remain unchanged.
    simi(jopt, :) = -sum(simi, dim=1)  ! Must ensure that 1 <= JOPT <= N!
end if

! Check whether SIMI is a poor approximation to the inverse of SIM(:, 1:N).
! Calculate SIMI from scratch if the current one is damaged by rounding errors.
erri = maximum(abs(matprod(simi, sim(:, 1:n)) - eye(n)))  ! MAXIMUM(X) returns NaN if X contains NaN
if (erri > TENTH * itol .or. is_nan(erri)) then
    simi_test = inv(sim(:, 1:n))
    erri_test = maximum(abs(matprod(simi_test, sim(:, 1:n)) - eye(n)))
    if (erri_test < erri .or. (is_nan(erri) .and. .not. is_nan(erri_test))) then
        simi = simi_test
        erri = erri_test
    end if
end if

! If SIMI is satisfactory, then update FVAL, CONMAT, and CVAL. Otherwise, restore SIM and SIMI, and
! return with INFO = DAMAGING_ROUNDING.
if (erri <= itol) then
    if (jopt >= 1 .and. jopt <= n) then
        fval([jopt, n + 1_IK]) = fval([n + 1_IK, jopt])
        conmat(:, [jopt, n + 1_IK]) = conmat(:, [n + 1_IK, jopt])
        cval([jopt, n + 1_IK]) = cval([n + 1_IK, jopt])
    end if
else  ! ERRI > ITOL or ERRI is NaN
    info = DAMAGING_ROUNDING
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
    call assert(.not. any(is_nan(conmat) .or. is_posinf(conmat)), 'CONMAT does not contain NaN/+Inf', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
    call assert(size(fval) == n + 1 .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == N+1 and FVAL is not NaN/+Inf', srname)
    call assert(size(sim, 1) == n .and. size(sim, 2) == n + 1, 'SIZE(SIM) == [N, N+1]', srname)
    call assert(all(is_finite(sim)), 'SIM is finite', srname)
    call assert(all(maxval(abs(sim(:, 1:n)), dim=1) > 0), 'SIM(:, 1:N) has no zero column', srname)
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

! Common modules
use, non_intrinsic :: consts_mod, only : IK, RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf

implicit none

! Inputs
real(RP), intent(in) :: cpen
real(RP), intent(in) :: cval(:)  ! CVAL(N+1)
real(RP), intent(in) :: fval(:)  ! FVAL(N+1)

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
    call assert(cpen > 0, 'CPEN > 0', srname)
    call assert(size(cval) == n + 1 .and. .not. any(cval < 0 .or. is_nan(cval) .or. is_posinf(cval)), &
        & 'SIZE(CVAL) == N+1 and CVAL does not contain negative values or NaN/+Inf', srname)
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
! Essentially, JOPT = MINLOC(PHI). However, we keep JOPT = N + 1 unless there is a strictly better
! choice. When there are multiple choices, we choose the JOPT with the smallest value of CVAL.
if (phimin < phi(jopt) .or. any(cval < cval(jopt) .and. phi <= phi(jopt))) then
    jopt = int(minloc(cval, mask=(phi <= phimin), dim=1), kind(jopt))
    !!MATLAB: cmin = min(cval(phi <= phimin)); jopt = find(phi <= phimin & cval <= cmin, 1, 'first');
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(jopt >= 1 .and. jopt <= n + 1, '1 <= JOPT <= N+1', srname)
    call assert(jopt == n + 1 .or. phi(jopt) < phi(n + 1) .or. (phi(jopt) <= phi(n + 1) .and. cval(jopt) < cval(n + 1)), &
        & 'JOPT = N+1 unless PHI(JOPT) < PHI(N+1) or PHI(JOPT) <= PHI(N+1) and CVAL(JOPT) < CVAL(N+1)', srname)
end if
end function findpole


end module update_cobyla_mod
