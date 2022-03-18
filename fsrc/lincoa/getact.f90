module getact_mod
!--------------------------------------------------------------------------------------------------!
! This module provides the GETACT subroutine of LINCOA, which picks the current active set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Friday, March 18, 2022 PM09:09:54
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, psd, vlam)
!--------------------------------------------------------------------------------------------------!
! The main purpose of GETACT is to pick the current active set. It is defined by the property that
! the projection of -G into the space orthogonal to the active constraint normals is as large as
! possible, subject to this projected steepest descent direction moving no closer to the boundary of
! every constraint whose current residual is at most 0.2*SNORM. On return, the settings in NACT,
! IACT, QFAC and RFAC are all appropriate to this choice of active set. The projected steepest
! descent direction itself is returned in PSD, which can be ZERO occasionally.
!
! AMAT, NACT, IACT, QFAC and RFAC are the same as the terms with these names in SUBROUTINE LINCOB.
! The current values must be set on entry. NACT, IACT, QFAC and RFAC are kept up to date when GETACT
! changes the current active set.
!
! SNORM, RESNEW, RESACT, G and PSD are the same as the terms with these names in SUBROUTINE TRSTEP.
! The elements of RESNEW and RESACT are also kept up to date.
!
! VLAM and W are used for working space, the vector VLAM being reserved for the Lagrange multipliers
! of the calculation.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert !, validate
use, non_intrinsic :: linalg_mod, only : matprod, inprod, eye, istriu, isorth

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: snorm

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)

! Outputs
real(RP), intent(out) :: psd(:)  ! PSD(N)
real(RP), intent(out) :: vlam(:)  ! VLAM(N)

! Local variables
character(len=*), parameter :: srname = 'GETACT'
real(RP) :: apsd(size(amat, 2))
real(RP) :: dd
real(RP) :: tol
real(RP) :: vmu(size(g))
real(RP) :: ddsav, dnorm, tdel, temp, violmx, vmult
integer(IK) :: i, ic, j, l

logical :: to_loop
logical :: mask(size(amat, 2))

integer(IK) :: m
integer(IK) :: n

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(g), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(psd) == n, 'SIZE(PSD) == N', srname)
    call assert(size(vlam) == n, 'SIZE(VLAM) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

! Quick return when M = 0?

! Set some constants and a temporary VLAM.
tdel = 0.2_RP * snorm
ddsav = TWO * inprod(g, g)
vlam = ZERO

! Set the initial QFAC to the identity matrix in the case NACT = 0.
if (nact == 0) then
    qfac = eye(n)
end if

! Remove any constraints from the initial active set whose residuals exceed TDEL.
do ic = nact, 1, -1
    if (resact(ic) > tdel) then
        ! Delete the constraint with index IACT(IC) from the active set, and set NACT = NACT - 1.
        call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
    end if
end do

! Remove any constraints from the initial active set whose Lagrange multipliers are nonnegative,
! and set the surviving multipliers.
! The following loop will run for at most NACT times, since each call of DEL_ACT reduces NACT by 1.
to_loop = (nact > 0)  ! Better implementation?
do while (to_loop)
    to_loop = .false.
    do ic = nact, 1, -1
        !temp = ZERO
        !do i = 1, n
        !    temp = temp + qfac(i, ic) * g(i)
        !end do
        temp = inprod(g, qfac(:, ic))
        do j = ic + 1, nact
            temp = temp - rfac(ic, j) * vlam(j)
        end do
        !temp = inprod(g, qfac(:, ic)) - inprod(rfac(ic, ic + 1:nact), vlam(ic + 1:nact))

        if (temp >= ZERO) then
            ! Delete the constraint with index IACT(IC) from the active set, and set NACT = NACT - 1.
            call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
            to_loop = (nact > 0)  ! Restart the loop if any active constraint is removed.
            exit
        else
            vlam(ic) = temp / rfac(ic, ic)
        end if
    end do
end do  ! End of DO WHILE (TO_LOOP)

!do while (nact > 0)
! vlam(1:nact) = lsqr(qfac(:, 1:nact), rfac(1:nact, 1:nact), g)
! if (any(vlam(1:nact) >= 0)) then
!     ic = findloc(vlam(1:nact) >= 0, back = .true.)
!     call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
! else
!     exit
! end if
!end do

! Set the new search direction D. Terminate if the 2-norm of D is ZERO or does not decrease, or if
! NACT=N holds. The situation NACT=N occurs for sufficiently large SNORM if the origin is in the
! convex hull of the constraint gradients.
dd = ZERO  ! Must be set, in case NACT = N at this point.
psd = ZERO  ! Must be set, in case NACT = N at this point.
do while (nact < n)    ! Infinite cycling possible?
    psd(nact + 1:n) = matprod(g, qfac(:, nact + 1:n))
    psd = -matprod(qfac(:, nact + 1:n), psd(nact + 1:n)) ! Projection of -G to range(QFAC(:,NACT+1:N))
    dd = inprod(psd, psd)

    if (dd >= ddsav) then
        psd = ZERO
        dd = ZERO  ! Why???
        exit
    end if
    if (dd == ZERO) exit
    ddsav = dd
    dnorm = sqrt(dd)

    ! Pick L, which is the index of the most violated constraint. Terminate if this violation is 0.
    apsd = matprod(psd, amat)
    mask = (resnew > 0 .and. resnew <= tdel .and. apsd > (dnorm / snorm) * resnew)
    if (any(mask)) then
        l = int(maxloc(apsd, mask=mask, dim=1), IK)
        violmx = apsd(l)
        ! MATLAB: apsd(mask) = -Inf; [violmx , l] = max(apsd);
    else
        exit
    end if

    ! Terminate if a positive value of VIOLMX may be due to computer rounding errors.
    ! N.B.: Theoretically (but not numerically), APSD(IACT(1:NACT)) = 0.
    ! Caution: when NACT = 0, MAXVAL returns -HUGE(APSD)!
    !if (violmx < 0.01_RP * dnorm .and. violmx <= TEN * maxval(abs(apsd(iact(1:nact))))) exit
    if (violmx < 0.01_RP * dnorm .and. violmx <= TEN * maxval([ZERO, abs(apsd(iact(1:nact)))])) exit

    call add_act(l, amat(:, l), iact, nact, qfac, resact, resnew, rfac, vlam)  ! NACT = NACT + 1, VLAM(NACT) = 0.

    ! Set the components of the vector VMU if VIOLMX is positive.
    ! N.B.:
    ! 1. In theory, NACT > 0 is not needed in the condition below, because VIOLMX is necessarily 0
    ! when NACT is 0. We keep NACT > 0 for security: when NACT <= 0, RFAC(NACT, NACT) is invalid.
    ! 2. The loop will run for at most NACT <= N times: if VIOLMX > 0, then IC > 0, and hence
    ! VLAM(IC) = 0, which implies that DEL_ACT will be called to reduce NACT by 1.
    do while (violmx > 0 .and. nact > 0)  ! Infinite cycling possible?
        vmu(nact) = ONE / rfac(nact, nact)**2  ! Here, NACT must be positive! Reason for SIGFAULT?
        do i = nact - 1, 1, -1
            vmu(i) = -inprod(rfac(i, i + 1:nact), vmu(i + 1:nact)) / rfac(i, i)
        end do
        ! vmu(1:nact) = lsqr(qfac(:, 1:nact), rfac(1:nact, 1:nact), qfac(:, nact)) / rfac(nact, nact)

        ! Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
        !vmult = minval([violmx, vlam(1:nact)/vmu(1:nact)])
        !ic = int(minloc([violmx, vlam(1:nact)/vmu(1:nact)], dim = 1) - 1, IK)

        vmult = violmx
        ic = 0
        do j = 1, nact
            if (vlam(j) >= vmult * vmu(j)) then  ! VLAM(1:NACT) <= 0 according to its updates.
                ! What about vlam(j) > vmult * vmu(j) ?
                ic = j
                vmult = vlam(j) / vmu(j)  ! What if VMU(J) = 0?
            end if
        end do
        violmx = max(violmx - vmult, ZERO)

        vlam(1:nact) = vlam(1:nact) - vmult * vmu(1:nact)
        if (ic > 0) then
            vlam(ic) = ZERO
        end if

        ! Reduce the active set if necessary, so that all components of the new VLAM are negative,
        ! with resetting of the residuals of the constraints that become inactive.
        do ic = nact, 1, -1
            !if (vlam(ic) >= 0) then
            if (.not. vlam(ic) < 0) then
                ! Powell's version: IF (.NOT. VLAM(IC) < 0) THEN
                ! Delete the constraint with index IACT(IC) from the active set, and set NACT = NACT - 1.
                call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)  ! NACT = NACT - 1
            end if
        end do
    end do  ! End of DO WHILE (VIOLMX > 0 .AND. NACT > 0)

    if (nact == 0) exit  ! It can only come from DEL_ACT when VLAM(1:NACT) >= 0. Possible at all?
end do  ! End of DO WHILE (NACT < N)

if (nact == n) then  ! The projected steepest gradient descent direction must be ZERO in this case.
    psd = ZERO
end if

! if (nact == 0) then
!   psd = -g
! end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)  ! Is NACT = 0 possible?
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(psd) == n, 'SIZE(PSD) == N', srname)
    call assert(size(vlam) == n, 'SIZE(VLAM) == N', srname)
end if

end subroutine getact


subroutine add_act(l, c, iact, nact, qfac, resact, resnew, rfac, vlam)
!--------------------------------------------------------------------------------------------------!
! This subroutine adds the constraint with index L to the active set as the (NACT+ )-th active
! constriant, updates IACT, QFAC, etc accordingly, and increments NACT to NACT+1. Here, C is the
! gradient of the new active constraint.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : qradd, istriu, isorth
implicit none

! Inputs
integer(IK), intent(in) :: l
real(RP), intent(in) :: c(:)  ! C(N)

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)
real(RP), intent(inout) :: vlam(:)  ! VLAM(N)

! Local variables (debugging only)
character(len=*), parameter :: srname = 'ADD_ACT'
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol

! Sizes
m = size(iact)
n = size(vlam)

! Preconditions
if (DEBUGGING) then
    call assert(m >= 1, 'M >= 1', srname)  ! Should not be called when M == 0.
    call assert(n >= 1, 'N >= 1', srname)
    call assert(nact >= 0 .and. nact <= min(m, n) - 1_IK, '0 <= NACT <= MIN(M, N)-1', srname)
    call assert(l >= 1 .and. l <= m, '1 <= L <= M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(.not. any(iact(1:nact) == l), 'L is not in IACT(1:NACT)', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

!====================!
! Calculation starts !
!====================!

! QRADD applies Givens rotations to the last (N-NACT) columns of QFAC so that the first (NACT+1)
! columns of QFAC are the ones required for the addition of the L-th constraint, and add the
! appropriate column to RFAC.
! N.B.: QRADD always augment NACT by 1. This is different from the strategy in COBYLA. Is it ensured
! that C cannot be linearly represented by the gradients of the existing active constraints?
call qradd(c, qfac, rfac, nact)  ! NACT is increased by 1!

! Indeed, it suffices to pass RFAC(:, 1:NACT+1) to QRADD as follows.
!!call qradd(c, qfac, rfac(:, 1:nact + 1), nact)  ! NACT is increased by 1!

! Update IACT, RESACT, RESNEW, and VLAM. N.B.: NACT has been increased by 1 in QRADD.
iact(nact) = l
resact(nact) = resnew(l)  ! RESACT(NACT) = RESNEW(IACT(NACT))
resnew(l) = ZERO  ! RESNEW(IACT(NACT)) = ZERO  ! Why not TINYCV? See DECACT.
vlam(nact) = ZERO

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact >= 1 .and. nact <= min(m, n), '1 <= NACT <= MIN(M, N)', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

end subroutine add_act


subroutine del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)
!--------------------------------------------------------------------------------------------------!
! This subroutine deletes the constraint with index IACT(IC) from the active set, updates IACT,
! QFAC, etc accordingly, and reduces NACT to NACT-1.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, EPS, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : qrexc, isorth, istriu
implicit none

! Inputs
integer(IK), intent(in) :: ic

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M)
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N)
real(RP), intent(inout) :: resact(:)  ! RESACT(M)
real(RP), intent(inout) :: resnew(:)  ! RESNEW(M)
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N)
real(RP), intent(inout) :: vlam(:)  ! VLAM(N)

! Local variables (debugging only)
character(len=*), parameter :: srname = 'DEL_ACT'
integer(IK) :: l
integer(IK) :: m
integer(IK) :: n
real(RP) :: tol

! Sizes
m = size(iact)
n = size(vlam)

! Preconditions
! Preconditions
if (DEBUGGING) then
    call assert(m >= 1, 'M >= 1', srname)  ! Should not be called when M == 0.
    call assert(n >= 1, 'N >= 1', srname)
    call assert(nact >= 1 .and. nact <= min(m, n), '1 <= NACT <= MIN(M, N)', srname)
    call assert(ic >= 1 .and. ic <= nact, '1 <= IC <= NACT', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(m + 1_IK, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
    l = iact(ic)  ! For debugging only
end if

!====================!
! Calculation starts !
!====================!

! The following instructions rearrange the active constraints so that the new value of IACT(NACT) is
! the old value of IACT(IC). QREXC implements the updates of QFAC and RFAC by sequence of Givens
! rotations. Then NACT is reduced by one.

call qrexc(qfac, rfac(:, 1:nact), ic)  ! QREXC does nothing if IC == NACT.
! Indeed, it suffices to pass QFAC(:, 1:NACT) and RFAC(1:NACT, 1:NACT) to QREXC as follows.
!!call qrexc(qfac(:, 1:nact), rfac(1:nact, 1:nact), ic)

iact(ic:nact) = [iact(ic + 1:nact), iact(ic)]
resact(ic:nact) = [resact(ic + 1:nact), resact(ic)]
resnew(iact(nact)) = max(resact(nact), TINYCV)
vlam(ic:nact) = [vlam(ic + 1:nact), vlam(ic)]
nact = nact - 1_IK

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(nact >= 0 .and. nact <= min(m, n) - 1, '1 <= NACT <= MIN(M, N)-1', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(.not. any(iact(1:nact) == l), 'L is not in IACT(1:NACT)', srname)

    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)

    call assert(size(resact) == m, 'SIZE(RESACT) == M', srname)
    call assert(size(resnew) == m, 'SIZE(RESNEW) == M', srname)
end if

end subroutine del_act


end module getact_mod
