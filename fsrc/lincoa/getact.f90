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
! Last Modified: Thursday, March 17, 2022 PM08:53:26
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: getact


contains


subroutine getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, dd, dw, vlam)
!--------------------------------------------------------------------------------------------------!
! The main purpose of GETACT is to pick the current active set. It is defined by the property that
! the projection of -G into the space orthogonal to the active constraint normals is as large as
! possible, subject to this projected steepest descent direction moving no closer to the boundary of
! every constraint whose current residual is at most 0.2*SNORM. On return, the settings in NACT,
! IACT, QFAC and RFAC are all appropriate to this choice of active set.
! Occasionally this projected direction is ZERO, and then the final value of W(1) is set to ZERO.
! Otherwise, the direction itself is returned in DW, and W(1) is set to the square of the length of
! the direction.
!
! AMAT, NACT, IACT, QFAC and RFAC are the same as the terms with these names in SUBROUTINE LINCOB.
! The current values must be set on entry. NACT, IACT, QFAC and RFAC are kept up to date when GETACT
! changes the current active set.
!
! SNORM, RESNEW, RESACT, G and DW are the same as the terms with these names in SUBROUTINE TRSTEP.
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
real(RP), intent(out) :: dd
real(RP), intent(out) :: dw(:)  ! DW(N)  ; better name?
real(RP), intent(out) :: vlam(:)  ! VLAM(N)

! Local variables
character(len=*), parameter :: srname = 'GETACT'
real(RP) :: dwa(size(amat, 2))
real(RP) :: w(size(g)), tol
real(RP) :: cvtol, ddsav, dnorm, tdel, temp, violmx, vmult
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

    call assert(size(dw) == n, 'SIZE(DW) == N', srname)
    call assert(size(vlam) == n, 'SIZE(VLAM) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

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
        temp = ZERO
        do i = 1, n
            temp = temp + qfac(i, ic) * g(i)
        end do
        if (ic < nact) then
            do j = ic + 1, nact
                temp = temp - rfac(ic, j) * vlam(j)
            end do
        end if
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

! Set the new search direction D. Terminate if the 2-norm of D is ZERO or does not decrease, or if
! NACT=N holds. The situation NACT=N occurs for sufficiently large SNORM if the origin is in the
! convex hull of the constraint gradients.
dd = ZERO  ! Must be set, in case NACT = N at this point.
do while (nact < n)    ! Infinite cycling possible?
    w(nact + 1:n) = matprod(g, qfac(:, nact + 1:n)) 
    dw = -matprod(qfac(:, nact + 1:n), w(nact + 1:n))  ! Projection of -G to range(QFAC(:, NACT+1:N))
    dd = inprod(dw, dw)

    if (dd >= ddsav) then
        dd = ZERO  ! Why???
        exit
    end if
    if (dd == ZERO) exit
    ddsav = dd
    dnorm = sqrt(dd)

    ! Pick the next integer L or terminate, a positive value of L being the index of the most
    ! violated constraint. The purpose of CVTOL below is to estimate whether a positive value of
    ! VIOLMX may be due to computer rounding errors.
    l = 0
    violmx = ZERO  ! Without this, VIOLMX is undefined if M = 0, and compilers will complain.
    cvtol = ZERO  ! Without this, CVTOL is undefined if M = 0, and compilers will complain.
    if (m > 0) then
        dwa = matprod(dw, amat)
        mask = (resnew > 0 .and. resnew <= tdel .and. dwa > (dnorm / snorm) * resnew)
        if (any(mask)) then
            l = int(maxloc(dwa, mask=mask, dim=1), IK)
            violmx = dwa(l)
            ! MATLAB: dwa(mask) = -Inf; [violmx , l] = max(dwa);
        end if

        if (violmx > ZERO .and. violmx < 0.01_RP * dnorm) then
            cvtol = maxval([abs(dwa(iact(1:nact))), cvtol])
        end if
    end if

    w(1) = ONE  ! Why?

    if (l == 0 .or. violmx <= TEN * cvtol) then
        exit
    end if

    call add_act(l, amat(:, l), iact, nact, qfac, resact, resnew, rfac, vlam)  ! NACT = NACT + 1

    ! Set the components of the vector VMU if VIOLMX is positive.
    ! N.B.: In theory, NACT > 0 is not needed in the condition below, because VIOLMX is necessarily
    ! 0 when NACT is 0. We keep NACT > 0 for security: when NACT <= 0, RFAC(NACT, NACT) is invalid.
    do while (violmx > 0 .and. nact > 0)  ! Infinite cycling possible?
        w(nact) = ONE / rfac(nact, nact)**2  ! Here, NACT must be positive! Reason for SIGFAULT?
        do i = nact - 1, 1, -1
            w(i) = -inprod(rfac(i, i + 1:nact), w(i + 1:nact)) / rfac(i, i)
        end do

        ! Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
        !vmult = minval([violmx, vlam(1:nact)/w(1:nact)])
        !ic = int(minloc([violmx, vlam(1:nact)/w(1:nact)], dim = 1) - 1, IK)

        vmult = violmx
        ic = 0
        do j = 1, nact
            if (vlam(j) >= vmult * w(j)) then
                ic = j
                vmult = vlam(j) / w(j)  ! What if W(J) <= 0?
            end if
        end do
        !do j = 1, nact
        !    vlam(j) = vlam(j) - vmult * w(j)
        !end do
        vlam(1:nact) = vlam(1:nact) - vmult * w(1:nact)
        if (ic > 0) then
            vlam(ic) = ZERO
        end if
        violmx = max(violmx - vmult, ZERO)
        !if (ic == 0) violmx = ZERO  ! NACT = 0 ==> IC = 0 ==> VIOLMX = 0. It seems artificial.

        ! Reduce the active set if necessary, so that all components of the new VLAM are negative,
        ! with resetting of the residuals of the constraints that become inactive.
        do ic = nact, 1, -1
            if (.not. vlam(ic) < ZERO) then
                ! Delete the constraint with index IACT(IC) from the active set, and set NACT = NACT - 1.
                call del_act(ic, iact, nact, qfac, resact, resnew, rfac, vlam)  ! NACT = NACT - 1
            end if
        end do

    end do  ! End of DO WHILE (VIOLMX > 0)

    if (nact == 0) exit  ! It can only come from DEL_ACT when VLAM(1:NACT) >= 0. Possible at all?
end do  ! End of DO WHILE (NACT < N)

if (nact == n) then
    ! Why DD should be 0? Because DD is the square of length of the projected steepest descent
    ! direction, projection made to the orthogonal complement of the active constraint gradients.
    ! NACT = N means that the orthogonal complement is {0}.
    dd = ZERO
end if

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

    call assert(size(dw) == n, 'SIZE(DW) == N', srname)
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
