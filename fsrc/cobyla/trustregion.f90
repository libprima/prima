module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region iterations.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: June 2021
!
! Last Modified: Wednesday, November 17, 2021 PM09:15:37
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstlp


contains

function trstlp(A, b, rho) result(d)
! This subroutine calculates an N-component vector D by the following two stages. In the first
! stage, D is set to the shortest vector that minimizes the greatest violation of the constraints
!       dot_product(A(1:N, K), D) >= B(K),  K = 1, 2, 3, ..., M,
! subject to the Euclidean length of D being at most RHO. If its length is strictly less than RHO,
! then we use the resultant freedom in D to minimize the objective function
!       dot_product(-A(1:N, M+1), D)
! subject to no increase in any greatest constraint violation. This notation allows the gradient of
! the objective function to be regarded as the gradient of a constraint. Therefore the two stages
! are distinguished by MCON == M and MCON > M respectively.
!
! It is possible but rare that a degeneracy may prevent D from attaining the target length RHO.
!
! CSTRV is the largest constraint violation of the current D: MAXVAL([B(1:M)-A(:,1:M)^T*D), ZERO]).
! ICON is the index of a most violated constraint if CSTRV is positive.
!
! NACT is the number of constraints in the active set and IACT(1), ...,IACT(NACT) are their indices,
! while the remainder of IACT contains a permutation of the remaining constraint indices.
! N.B.: NACT <= min(M, N). Obviously, NACT <= M. In addition, The constraints in IACT(1, ..., NACT)
! have linearly independent gradients (see the comments above the instructions that delete a
! constraint from the active set to make room for the new active constraint with index IACT(ICON));
! it can also be seen from the update of NACT: starting from 0, NACT is incremented only if NACT < N.
!
! Further, Z is an orthogonal matrix whose first NACT columns can be regarded as the result of
! Gram-Schmidt applied to the active constraint gradients. For J = 1, 2, ..., NACT, the number
! ZDOTA(J) is the scalar product of the J-th column of Z with the gradient of the J-th active
! constraint. D is the current vector of variables and here the residuals of the active constraints
! should be zero. Further, the active constraints have nonnegative Lagrange multipliers that are
! held at the beginning of VMULTC. The remainder of this vector holds the residuals of the inactive
! constraints at d, the ordering of the components of VMULTC being in agreement with the permutation
! of the indices of the constraints that is in IACT. All these residuals are nonnegative, which is
! achieved by the shift CSTRV that makes the least residual zero.

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop, verisize

implicit none

! Inputs
real(RP), intent(in) :: A(:, :)  ! (N, M+1)
real(RP), intent(in) :: b(:)  ! M+1
real(RP), intent(in) :: rho

! Output
real(RP) :: d(size(A, 1))  ! N

! Local variables
integer(IK) :: iact(size(b))
integer(IK) :: m
integer(IK) :: nact
real(RP) :: vmultc(size(b))
real(RP) :: z(size(d), size(d))
character(len=*), parameter :: srname = 'TRSTLP'

! Get and verify the sizes.
m = size(A, 2) - 1_IK

if (DEBUGGING) then
    if (m < 0 .or. size(A, 1) <= 0) then
        call errstop(srname, 'SIZE(A) is invalid')
    end if
    call verisize(b, m + 1)
end if

call trstlp_sub(iact(1:m), nact, 1, A(:, 1:m), b(1:m), rho, d, vmultc(1:m), z)
call trstlp_sub(iact, nact, 2, A, b, rho, d, vmultc, z)

end function trstlp


!--------------------------------------------------------------------------------------------------!
!---------------- QUESTION: What are exactly the objective and algorithm of trstlp_sub? -----------!
! The algorithm was NOT documented in the COBYLA paper. A note should be written to introduce it!
! As a major part of the algorithm, the code maintains and updates the QR factorization of
! A(IACT(1:NACT)), i.e., the gradients of all the active (linear) constraints. The matrix Z is
! indeed Q, and the vector ZDOTA is the diagonal of R. The factorization is updated by Givens
! rotations when an index is added in or removed from IACT.
!--------------------------------------------------------------------------------------------------!


subroutine trstlp_sub(iact, nact, stage, A, b, rho, d, vmultc, z)
! This subroutine does the real calculations for TRSTLP, both stage 1 and stage 2.
! Major differences between stage 1 and stage 2:
! 1. Initialization. Stage 2 inherits the values of some variables from stage 1, so they are
! initialized in stage 1 but not in stage 2.
! 2. CSTRV. CSTRV is updated after at iteration in stage 1, whereas it remains a constant in stage 2.
! 3. SDIRN. See the definition of SDIRN in the code for details.
! 4. OPTNEW. The two stages have different objectives, so OPTNEW is updated differently.
! 5. STEP. STEP <= CSTRV in stage 1.

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, EPS, HUGENUM, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: debug_mod, only : assert, errstop, verisize
use, non_intrinsic :: linalg_mod, only : inprod, matprod, eye, isminor, lsqr, qradd, qrexc
implicit none

! Inputs
integer(IK), intent(in) :: stage
real(RP), intent(in) :: A(:, :)  ! (N, MCON)
real(RP), intent(in) :: b(:)  ! MCON
real(RP), intent(in) :: rho

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! MCON
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: d(:)  ! N
real(RP), intent(inout) :: vmultc(:)  ! MCON
real(RP), intent(inout) :: z(:, :)  ! (N, N)

! Local variables
integer(IK) :: icon
integer(IK) :: iter
integer(IK) :: k
integer(IK) :: m
integer(IK) :: maxiter
integer(IK) :: mcon
integer(IK) :: n
integer(IK) :: nactold
integer(IK) :: nactsav
integer(IK) :: nfail
real(RP) :: cstrv
real(RP) :: cvold
real(RP) :: cvsabs(size(b))
real(RP) :: cvshift(size(b))
real(RP) :: dd
real(RP) :: dnew(size(d))
real(RP) :: dold(size(d))
real(RP) :: frac
real(RP) :: fracmult(size(b))
real(RP) :: optnew
real(RP) :: optold
real(RP) :: sd
real(RP) :: sdirn(size(d))
real(RP) :: ss
real(RP) :: step
real(RP) :: vmultd(size(b))
real(RP) :: zdasav(size(z, 2))
real(RP) :: zdota(size(z, 2))
character(len=*), parameter :: srname = 'TRSTLP_SUB'

! Get and verify the sizes.
n = size(A, 1)
mcon = size(A, 2)

if (DEBUGGING) then
    if (n <= 0 .and. mcon > 0) then
        call errstop(srname, 'SIZE(A, 1) <= 0 while SIZE(A, 2) > 0')
    end if
    call verisize(b, mcon)
    call verisize(iact, mcon)
    call verisize(vmultc, mcon)
    call verisize(d, n)
    call verisize(z, n, n)
end if

! Initialization according to STAGE.
if (stage == 1) then
    iact = [(k, k=1, size(iact))]
    nact = 0_IK
    d = ZERO
    cstrv = maxval([b, ZERO])
    vmultc = cstrv - b
    z = eye(n)
    if (mcon == 0 .or. cstrv <= ZERO) then
        ! Check whether a quick return is possible. Make sure the In-outputs have been initialized.
        return
    end if

    m = mcon
    icon = maxloc(b, dim=1)
    sdirn = ZERO
else
    if (inprod(d, d) >= rho**2) then
        ! Check whether a quick return is possible.
        return
    end if

    iact(mcon) = mcon
    vmultc(mcon) = ZERO
    m = mcon - 1_IK
    icon = mcon

    ! In Powell's code, stage 2 uses the ZDOTA and CSTRV calculated by stage 1. Here we re-calculate
    ! them so that they need not be passed from stage 1 to 2, and hence the coupling is reduced.
    cstrv = maxval([b(1:m) - matprod(d, A(:, 1:m)), ZERO])
    zdota(1:nact) = [(inprod(z(:, k), A(:, iact(k))), k=1, nact)]
end if

! More initialization.
optold = HUGENUM
nactold = nact
nfail = 0_IK

! Zaikun 20211012: Is it true that IACT(NACT) is always the constraint with the largest vioaltion (in stage 1)?

!----------------------------------------------------------------------------------------------!
! Zaikun 20211011: VMULTD is computed from scratch at each iteration, but VMULTC is inherited.
!----------------------------------------------------------------------------------------------!

! Powell's code can encounter infinite cycling, which did happen when testing the following CUTEst
! problems: DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf/NaN appear in D due to extremely large values in A (up to 10^219).
! To resolve this, we set the maximal number of iterations to MAXITER, and terminate in case
! Inf/NaN occurs in D.
maxiter = min(50000_IK, 100_IK * max(m, n))
do iter = 1, maxiter
    if (stage == 1) then
        optnew = cstrv
    else
        optnew = -inprod(d, A(:, mcon))
    end if

    ! End the current stage of the calculation if 3 consecutive iterations have either failed to
    ! reduce the best calculated value of the objective function or to increase the number of active
    ! constraints since the best value was calculated. This strategy prevents cycling, but there is
    ! a remote possibility that it will cause premature termination.
    if (optnew < optold .or. nact > nactold) then
        nactold = nact
        nfail = 0_IK
    else
        nfail = nfail + 1_IK
    end if
    optold = min(optold, optnew)
    if (nfail == 3) then
        exit
    end if

    ! If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to the active set.
    ! Apply Givens rotations so that the last N-NACT-1 columns of Z are orthogonal to the gradient
    ! of the new constraint, a scalar product being set to zero if its nonzero value could be due to
    ! computer rounding errors, which is tested by ISMINOR.
    if (icon > nact) then
        zdasav = zdota
        nactsav = nact
        call qradd(A(:, iact(icon)), z, zdota, nact)  ! QRADD may update NACT tp NACT + 1.

        if (nact == nactsav + 1) then
            ! N.B.: It is problematic to index arrays using [NACT, ICON] when NACT == ICON.
            ! Zaikun 20211012: Why should VMULTC(NACT) = 0?
            if (nact /= icon) then
                vmultc([icon, nact]) = [vmultc(nact), ZERO]
                iact([icon, nact]) = iact([nact, icon])
            else
                vmultc(nact) = ZERO
            end if
        else
            !----------------------------! 1st VMULTD CALCULATION STARTS  !------------------------!
            ! Zaikun 20211011:
            ! 1. VMULTD is calculated from scratch for the first time (out of 2) in one iteration.
            ! 2. The sole purpose of this VMULTD is to compute VMULTC(1 : NACT) and check possible
            ! exit immediately after this loop. Only VMULTD(1 : NACT) is needed.
            ! 3. VMULTD will be computed from scratch again later.
            ! 4. NOTE that IACT has not been updated to replace IACT(NACT) with IACT(ICON). Thus
            ! A(:, IACT(1:NACT)) is the UNUPDATED version before QRADD (Z(:, 1:NACT) remains the
            ! same before and after QRADD). Therefore, if we supply ZDOTA to LSQR (as Rdiag) as
            ! Powell did, we should use the UNUPDATED version, namely ZDASAV.
            vmultd(1:nact) = lsqr(A(:, iact(1:nact)), A(:, iact(icon)), z(:, 1:nact), zdasav(1:nact))
            !----------------------------! 1st VMULTD CALCULATION ENDS  !--------------------------!

            frac = minval(vmultc(1:nact) / vmultd(1:nact), mask=(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m))
            if (frac < ZERO .or. .not. any(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m)) then
                exit  ! This can be triggered by NACT == 0, among other possibilities.
            end if

            ! Revise the Lagrange multipliers. The revision is not applicable to VMULTC(NACT + 1:M).
            vmultc(1:nact) = max(ZERO, vmultc(1:nact) - frac * vmultd(1:nact))

            ! Zaikun 20210811: Powell's code includes the following, which is IMPOSSIBLE TO REACH.
            !--------------------------------------------------------------------------------------!
            !if (icon < nact) then
            !    do k = icon, nact-1
            !        hypt = sqrt(zdota(k+1)**2+inprod(z(:, k), A(:, iact(k+1)))**2)
            !        grot = planerot([zdota(k+1), inprod(z(:, k), A(:, iact(k+1)))])
            !        z(:, [k, k+1]) = matprod(z(:, [k+1, k]), transpose(grot))
            !        zdota([k, k+1]) = [hypt, (zdota(k+1) / hypt) * zdota(k)]
            !    end do
            !    iact(icon:nact) = [iact(icon+1:nact), iact(icon)]
            !    vmultc(icon:nact) = [vmultc(icon+1:nact), vmultc(icon)]
            !end if
            !--------------------------------------------------------------------------------------!

            ! Reorder the active constraints so that the one to be replaced is at the end of the list.
            ! Exit if the new value of ZDOTA(NACT) is not acceptable. Note that the opposite of
            ! 'ABS(ZDOTA(NACT)) > 0' is not 'ABS(ZDOTA(NACT) <= 0)', as ZDOTA(NACT) can be NaN.
            if (abs(zdota(nact)) > 0) then
                vmultc([icon, nact]) = [ZERO, frac]  ! VMULTC([ICON, NACT]) is valid as ICON > NACT.
                iact([icon, nact]) = iact([nact, icon])
            else
                exit
            end if
        end if

        ! Ensure that the objective function continues to be treated as the last active constraint
        ! if stage 2 is in progress.
        ! Zaikun 20211011, 20211111: Is it guaranteed for stage 2 that IACT(NACT-1) = MCON when
        ! IACT(NACT) /= MCON???? If not, then how does the following procedure ensure that MCON is
        ! the last of IACT(1:NACT)?
        if (stage == 2 .and. iact(nact) /= mcon) then
            call qrexc(A(:, iact(1:nact)), Z, zdota(1:nact), nact - 1)
            iact([nact - 1, nact]) = iact([nact, nact - 1])
            !!??zdota(nact-1:nact) = [(inprod(z(:, k), A(:, iact(k))), k = nact-1, nact)]
            vmultc([nact - 1, nact]) = vmultc([nact, nact - 1])
        end if
        ! Zaikun 20211117: It turns out that the last few lines do not guarantee IACT(NACT) == N in
        ! stage 2; the following test cannot be passed. IS THIS A BUG??!!
        !call assert(iact(nact) == mcon .or. stage == 1, 'IACT(NACT) == MCON in stage 2', srname)

        ! Set SDIRN to the direction of the next change to the current vector of variables.
        ! Usually during stage 1 the vector SDIRN gives a search direction that reduces all the
        ! active constraint violations by one simultaneously.
        if (stage == 1) then
            sdirn = sdirn - ((inprod(sdirn, A(:, iact(nact))) - ONE) / zdota(nact)) * z(:, nact)
        else
            sdirn = (ONE / zdota(nact)) * z(:, nact)
            ! SDIRN = Z(:, NACT)/(A(:,IACT(NACT))^T*Z(:, NACT))
            ! SDIRN^T*A(:, IACT(NACT)) = 1, SDIRN is orthogonal to A(:, IACT(1:NACT-1)) and is
            ! parallel to Z(:, NACT).
        end if
    else
        ! Delete the constraint with the index IACT(ICON) from the active set, which is done by
        ! reordering IACT(ICONT:NACT) into [IACT(ICON+1:NACT), IACT(ICON)]. In theory, ICON > 0.
        ! To be safe, the condition below requires ICON > 0, which does not exist in Powell's code.
        if (icon < nact .and. icon > 0) then
            call qrexc(A(:, iact(1:nact)), Z, zdota(1:nact), icon)
            iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
            !!??zdota(icon:nact) = [(inprod(z(:, k), A(:, iact(k))), k = icon, nact)]
            vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(icon)]
        end if
        nact = nact - 1

        ! Set SDIRN to the direction of the next change to the current vector of variables.
        if (stage == 1) then
            sdirn = sdirn - inprod(sdirn, z(:, nact + 1)) * z(:, nact + 1)
            ! SDIRN is orthogonal to Z(:, NACT+1)
        else
            sdirn = (ONE / zdota(nact)) * z(:, nact)
            ! SDIRN = Z(:, NACT)/(A(:,IACT(NACT))^T*Z(:, NACT))
            ! SDIRN^T*A(:, IACT(NACT)) = 1, SDIRN is orthogonal to A(:, IACT(1:NACT-1)) and is
            ! parallel to Z(:, NACT).
        end if
    end if

    ! Calculate the step to the boundary of the trust region or take the step that reduces CSTRV to
    ! zero. The two statements below that include the factor EPS prevent some harmless underflows
    ! that occurred in a test calculation (here, EPS is the machine epsilon; Powell's original code
    ! used 1.0E-6, and Powell's code was written in SINGLE PRECISION). Further, we skip the step if
    ! it could be zero within a reasonable tolerance for computer rounding errors.
    dd = rho**2 - sum(d**2, mask=(abs(d) >= EPS * rho))
    if (dd <= ZERO) then
        exit
    end if
    ss = inprod(sdirn, sdirn)
    sd = inprod(sdirn, d)
    if (abs(sd) >= EPS * sqrt(ss * dd)) then
        step = dd / (sqrt(ss * dd + sd**2) + sd)
    else
        step = dd / (sqrt(ss * dd) + sd)
    end if
    if (stage == 1) then
        if (isminor(cstrv, step)) then
            exit
        end if
        step = min(step, cstrv)
    end if

    ! Set DNEW to the new variables if STEP is the steplength, and reduce CSTRV to the corresponding
    ! maximum residual if stage 1 is being done.
    dnew = d + step * sdirn
    if (stage == 1) then
        cvold = cstrv
        cstrv = maxval([b(iact(1:nact)) - matprod(dnew, A(:, iact(1:nact))), ZERO])
        ! N.B.: CSTRV will be used when calculating VMULTD(NACT+1 : MCON).
    end if

    !--------------------------------! 2nd VMULTD CALCULATION STARTS !-----------------------------!
    ! Zaikun 20211011:
    ! 1. VMULTD is computed from scratch for the second (out of 2) time in one iteration.
    ! 2. VMULTD(1:NACT) and VMULTD(NACT+1:MCON) are calculated separately with no coupling.
    ! 3. VMULTD will be calculated from scratch again in the next iteration.
    !
    ! Set VMULTD to the VMULTC vector that would occur if D became DNEW. A device is included to
    ! force VMULTD(K)=ZERO if deviations from this value can be attributed to computer rounding
    ! errors. First calculate the new Lagrange multipliers.
    vmultd(1:nact) = lsqr(A(:, iact(1:nact)), dnew, z(:, 1:nact), zdota(1:nact))
    if (stage == 2) then
        vmultd(nact) = max(ZERO, vmultd(nact))  ! This seems never activated.
    end if

    ! Complete VMULTD by finding the new constraint residuals. (Powell wrote "Complete VMULTC ...")
    cvshift = matprod(dnew, A(:, iact)) - b(iact) + cstrv  ! Only CVSHIFT(nact+1:mcon) is needed.
    cvsabs = matprod(abs(dnew), abs(A(:, iact))) + abs(b(iact)) + cstrv
    where (isminor(cvshift, cvsabs))
        cvshift = ZERO
    end where
    vmultd(nact + 1:mcon) = cvshift(nact + 1:mcon)
    !--------------------------------! 2nd VMULTD CALCULATION ENDS !-------------------------------!

    ! Calculate the fraction of the step from D to DNEW that will be taken.
    fracmult = vmultc / (vmultc - vmultd)  !
    frac = min(ONE, minval(fracmult, mask=(vmultd < ZERO .and. .not. is_nan(fracmult))))

    if (frac < ONE) then
        icon = minloc(fracmult, mask=(vmultd < ZERO .and. .not. is_nan(fracmult)), dim=1)
    else
        icon = 0_IK  ! This will trigger an exit after the update of D, VMULTC, and CSTRV.
    end if

    ! Update D, VMULTC and CSTRV.
    dold = d
    d = (ONE - frac) * d + frac * dnew
    ! Exit in case of Inf/NaN in D.
    if (.not. is_finite(sum(abs(d)))) then
        d = dold
        exit
    end if

    vmultc = max(ZERO, (ONE - frac) * vmultc + frac * vmultd)
    if (stage == 1) then
        cstrv = (ONE - frac) * cvold + frac * cstrv
        ! In theory, CSTRV = MAXVAL([B(1:M) - MATPROD(D, A(:, 1:M)), ZERO]), yet the CSTRV updated
        ! as above can be quite different from this value if A has huge entries (e.g., > 1E20).
    end if

    if (icon == 0) then
        exit
    end if
end do

end subroutine trstlp_sub

end module trustregion_mod
