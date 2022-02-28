module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of COBYLA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the COBYLA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: June 2021
!
! Last Modified: Monday, February 28, 2022 PM02:44:45
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstlp
public :: trrad


contains


function trstlp(A, b, delta) result(d)
!--------------------------------------------------------------------------------------------------!
! This subroutine calculates an N-component vector D by the following two stages. In the first
! stage, D is set to the shortest vector that minimizes the greatest violation of the constraints
!       dot_product(A(1:N, K), D) >= B(K),  K = 1, 2, 3, ..., M,
! subject to the Euclidean length of D being at most DELTA. If its length is strictly less than DELTA,
! then we use the resultant freedom in D to minimize the objective function
!       dot_product(-A(1:N, M+1), D)
! subject to no increase in any greatest constraint violation. This notation allows the gradient of
! the objective function to be regarded as the gradient of a constraint. Therefore the two stages
! are distinguished by MCON == M and MCON > M respectively.
!
! It is possible but rare that a degeneracy may prevent D from attaining the target length DELTA.
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
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : norm
implicit none

! Inputs
real(RP), intent(in) :: A(:, :)  ! (N, M+1)
real(RP), intent(in) :: b(:)  ! M+1
real(RP), intent(in) :: delta

! Outputs
real(RP) :: d(size(A, 1))  ! N

! Local variables
character(len=*), parameter :: srname = 'TRSTLP'
integer(IK) :: iact(size(b))
integer(IK) :: m
integer(IK) :: nact
real(RP) :: vmultc(size(b))
real(RP) :: z(size(d), size(d))

! Sizes
m = int(size(A, 2) - 1, kind(m))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(size(A, 1) >= 1 .and. size(A, 2) >= 1, 'SIZE(A) >= [1, 1]', srname)
    call assert(size(b) == size(A, 2), 'SIZE(B) == size(A, 2)', srname)
end if

!====================!
! Calculation starts !
!====================!

call trstlp_sub(iact(1:m), nact, 1_IK, A(:, 1:m), b(1:m), delta, d, vmultc(1:m), z)
call trstlp_sub(iact, nact, 2_IK, A, b, delta, d, vmultc, z)

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == size(A, 1), 'SIZE(D) == SIZE(A, 1)', srname)
    call assert(all(is_finite(d)), 'D is finite', srname)
    call assert(norm(d) <= TWO * delta, '|D| <= 2*DELTA', srname)
end if
end function trstlp


!--------------------------------------------------------------------------------------------------!
! N.B.:
! 1. The algorithm was NOT documented in the COBYLA paper. A note should be written to introduce it!
! 2. As a major part of the algorithm, the code maintains and updates the QR factorization of
! A(IACT(1:NACT)), i.e., the gradients of all the active (linear) constraints. The matrix Z is
! indeed Q, and the vector ZDOTA is the diagonal of R. The factorization is updated by Givens
! rotations when an index is added in or removed from IACT.
! 3. There are probably better algorithms available for this trust-region linear programming problem.
!--------------------------------------------------------------------------------------------------!


subroutine trstlp_sub(iact, nact, stage, A, b, delta, d, vmultc, z)
!--------------------------------------------------------------------------------------------------!
! This subroutine does the real calculations for TRSTLP, both stage 1 and stage 2.
! Major differences between stage 1 and stage 2:
! 1. Initialization. Stage 2 inherits the values of some variables from stage 1, so they are
! initialized in stage 1 but not in stage 2.
! 2. CSTRV. CSTRV is updated after at iteration in stage 1, whereas it remains a constant in stage 2.
! 3. SDIRN. See the definition of SDIRN in the code for details.
! 4. OPTNEW. The two stages have different objectives, so OPTNEW is updated differently.
! 5. STEP. STEP <= CSTRV in stage 1.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod, eye, isminor, lsqr, qradd, qrexc, norm
implicit none

! Inputs
integer(IK), intent(in) :: stage
real(RP), intent(in) :: A(:, :)  ! (N, MCON)
real(RP), intent(in) :: b(:)  ! MCON
real(RP), intent(in) :: delta

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! MCON
integer(IK), intent(inout) :: nact
real(RP), intent(inout) :: d(:)  ! N
real(RP), intent(inout) :: vmultc(:)  ! MCON
real(RP), intent(inout) :: z(:, :)  ! (N, N)

! Local variables
character(len=*), parameter :: srname = 'TRSTLP_SUB'
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
!real(RP) :: cvold
real(RP) :: cvsabs(size(b))
real(RP) :: cvshift(size(b))
real(RP) :: dd
real(RP) :: dnew(size(d))
real(RP) :: dold(size(d))
real(RP) :: frac
real(RP) :: fracmult(size(vmultc))
real(RP) :: optnew
real(RP) :: optold
real(RP) :: sd
real(RP) :: sdirn(size(d))
real(RP) :: ss
real(RP) :: step
real(RP) :: vmultd(size(vmultc))
real(RP) :: zdasav(size(z, 2))
real(RP) :: zdota(size(z, 2))

! Sizes
n = int(size(A, 1), kind(n))
mcon = int(size(A, 2), kind(mcon))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(stage == 1 .or. stage == 2, 'STAGE == 1 or 2', srname)
    call assert((mcon >= 0 .and. stage == 1) .or. (mcon >= 1 .and. stage == 2), &
        & 'MCON >= 1 in stage 1 and MCON >= 0 in stage 2', srname)
    call assert(size(b) == mcon, 'SIZE(B) == MCON', srname)
    call assert(size(iact) == mcon, 'SIZE(IACT) == MCON', srname)
    call assert(size(vmultc) == mcon, 'SIZE(VMULTC) == MCON', srname)
    call assert(size(d) == n, 'SIZE(D) == N', srname)
    call assert(size(z, 1) == n .and. size(z, 2) == n, 'SIZE(Z) == [N, N]', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    if (stage == 2) then
        call assert(all(is_finite(d)) .and. norm(d) <= TWO * delta, &
            & 'D is finite and |D| <= 2*DELTA at the beginning of stage 2', srname)
        call assert((nact >= 0 .and. nact <= min(mcon, n)), &
            & '0 <= NACT <= MIN(MCON, N) at the beginning of stage 2', srname)
        call assert(all(vmultc(1:mcon - 1) >= 0), 'VMULTC >= 0 at the beginning of stage 2', srname)
        ! N.B.: Stage 1 defines only VMULTC(1:M); VMULTC(M+1) is undefined!
    end if
end if

!====================!
! Calculation starts !
!====================!

! Initialization according to STAGE.
if (stage == 1) then
    iact = [(k, k=1, int(size(iact), IK))]
    nact = 0_IK
    d = ZERO
    cstrv = maxval([b, ZERO])
    vmultc = cstrv - b
    z = eye(n)
    if (mcon == 0 .or. cstrv <= ZERO) then
        ! Check whether a quick return is possible. Make sure the In-outputs have been initialized.
        return
    end if

    if (all(is_nan(b))) then
        return
    else
        icon = int(maxloc(b, mask=(.not. is_nan(b)), dim=1), kind(icon))
    end if
    m = mcon
    sdirn = ZERO
else
    if (inprod(d, d) >= delta**2) then
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
end if
zdota(1:nact) = [(inprod(z(:, k), A(:, iact(k))), k=1, nact)]

! More initialization.
optold = HUGENUM
nactold = nact
nfail = 0_IK

!----------------------------------------------------------------------------------------------!
! Zaikun 20211011: VMULTD is computed from scratch at each iteration, but VMULTC is inherited.
!----------------------------------------------------------------------------------------------!

! Powell's code can encounter infinite cycling, which did happen when testing the following CUTEst
! problems: DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322. Indeed, in all
! these cases, Inf/NaN appear in D due to extremely large values in A (up to 10^219). To resolve
! this, we set the maximal number of iterations to MAXITER, and terminate if Inf/NaN occurs in D.
! In MATLAB or Python: MAXITER = MIN(10000, 100*MAX(M, N))
maxiter = int(min(int(10_IK**min(4, range(0_IK)), IK), 100_IK * max(m, n)), kind(maxiter))
do iter = 1, maxiter
    if (DEBUGGING) then
        call assert(all(vmultc >= 0), 'VMULTC >= 0', srname)
    end if
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
        zdasav(1:nact) = zdota(1:nact)
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
            ! Zaikun 20211011:
            ! 1. VMULTD is calculated from scratch for the first time (out of 2) in one iteration.
            ! 2. NOTE that IACT has not been updated to replace IACT(NACT) with IACT(ICON). Thus
            ! A(:, IACT(1:NACT)) is the UNUPDATED version before QRADD (Z(:, 1:NACT) remains the
            ! same before and after QRADD). Therefore, if we supply ZDOTA to LSQR (as Rdiag) as
            ! Powell did, we should use the UNUPDATED version, namely ZDASAV.
            vmultd(1:nact) = lsqr(A(:, iact(1:nact)), A(:, iact(icon)), z(:, 1:nact), zdasav(1:nact))
            if (.not. any(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m)) then
                ! N.B.: This can be triggered by NACT == 0 (among other possibilities)! This is
                ! important, because NACT will be used as an index in the sequel.
                exit
            end if
            ! VMULTD(NACT+1:MCON) is not used, but we have to initialize it in Fortran, or compilers
            ! complain about the WHERE construct below (another solution: restrict WHERE to 1:NACT).
            vmultd(nact + 1:mcon) = -ONE  ! SIZE(VMULTD) = MCON

            ! Revise the Lagrange multipliers. The revision is not applicable to VMULTC(NACT + 1:M).
            where (vmultd > 0 .and. iact <= m)
                fracmult = vmultc / vmultd
            elsewhere
                fracmult = HUGENUM
            end where
            ! Only the places where VMULTD > 0 and IACT <= M is relevant blow, if any.
            frac = minval(fracmult(1:nact))  ! FRACMULT(NACT+1:MCON) may contain garbage.
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
            ! N.B.: We cannot arrive here with NACT == 0, which should have triggered an exit above.
            if (abs(zdota(nact)) <= 0 .or. is_nan(zdota(nact))) then
                exit
            end if
            vmultc([icon, nact]) = [ZERO, frac]  ! VMULTC([ICON, NACT]) is valid as ICON > NACT.
            iact([icon, nact]) = iact([nact, icon])
        end if

        ! Ensure that the objective function continues to be treated as the last active constraint
        ! if stage 2 is in progress.
        ! Zaikun 20211011, 20211111: Is it guaranteed for stage 2 that IACT(NACT-1) = MCON when
        ! IACT(NACT) /= MCON??? If not, then how does the following procedure ensure that MCON is
        ! the last of IACT(1:NACT)?
        if (stage == 2 .and. iact(nact) /= mcon) then
            call qrexc(A(:, iact(1:nact)), z, zdota(1:nact), nact - 1_IK)
            iact([nact - 1_IK, nact]) = iact([nact, nact - 1_IK])
            !!??zdota(nact-1:nact) = [(inprod(z(:, k), A(:, iact(k))), k = nact-1_IK, nact)]
            vmultc([nact - 1_IK, nact]) = vmultc([nact, nact - 1_IK])
        end if
        ! Zaikun 20211117: It turns out that the last few lines do not guarantee IACT(NACT) == N in
        ! stage 2; the following test cannot be passed. IS THIS A BUG??!!
        !call assert(iact(nact) == mcon .or. stage == 1, 'IACT(NACT) == MCON in stage 2', srname)

        !------- Powell's code does not include the following ------!
        if (abs(zdota(nact)) <= 0 .or. is_nan(zdota(nact))) then
            exit
        end if
        !-----------------------------------------------------------!

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
            call qrexc(A(:, iact(1:nact)), z, zdota(1:nact), icon)
            iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
            !!??zdota(icon:nact) = [(inprod(z(:, k), A(:, iact(k))), k = icon, nact)]
            vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(icon)]
        end if
        nact = nact - 1_IK

        !------- Powell's code does not include the following ------!
        if (abs(zdota(nact)) <= 0 .or. is_nan(zdota(nact))) then
            exit
        end if
        !-----------------------------------------------------------!

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

    ! Calculate the step to the trust region boundary or take the step that reduces CSTRV to 0.
    !----------------------------------------------------------------------------------------------!
    ! The following calculation of STEP is simpler than Powell's approach and seems to improve the
    ! performance of COBYLA. In our test, SQRT(SS * DD + SD**2) + SD can become 0 due to a tiny SS
    ! or DD or due to a negative SD; hence we reformulated STEP as below and exit when SS becomes 0.
    ! We also found that removing the precaution about underflows is beneficial to the overall
    ! performance of COBYLA --- the underflows are harmless anyway.
    !----------------------------------------------------------------!
    dd = delta**2 - inprod(d, d)
    ss = inprod(sdirn, sdirn)
    if (dd <= 0 .or. ss <= 0) then
        exit
    end if
    sd = inprod(sdirn, d)
    step = (sqrt(ss * dd + sd**2) - sd) / ss
    !----------------------------------------------------------------!
    ! Powell's approach and comments are as follows.
    !----------------------------------------------------------------!
    ! The two statements below that include the factor EPS prevent
    ! some harmless underflows that occurred in a test calculation
    ! (Zaikun: here, EPS is the machine epsilon; Powell's original
    ! code used 1.0E-6, and Powell's code was written in SINGLE
    ! PRECISION). Further, we skip the step if it could be zero within
    ! a reasonable tolerance for computer rounding errors.
    !
    !!dd = delta**2 - sum(d**2, mask=(abs(d) >= EPS * delta))
    !!ss = inprod(sdirn, sdirn)
    !!if (dd  <= 0) then
    !!    exit
    !!end if
    !!sd = inprod(sdirn, d)
    !!if (abs(sd) >= EPS * sqrt(ss * dd)) then
    !!    step = dd / (sqrt(ss * dd + sd**2) + sd)
    !!else
    !!    step = dd / (sqrt(ss * dd) + sd)
    !!end if
    !----------------------------------------------------------------!
    !----------------------------------------------------------------------------------------------!

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
        !cvold = cstrv
        cstrv = maxval([b(iact(1:nact)) - matprod(dnew, A(:, iact(1:nact))), ZERO])
        ! N.B.: CSTRV will be used when calculating VMULTD(NACT+1 : MCON).
    end if

    ! Zaikun 20211011:
    ! 1. VMULTD is computed from scratch for the second (out of 2) time in one iteration.
    ! 2. VMULTD(1:NACT) and VMULTD(NACT+1:MCON) are calculated separately with no coupling.
    ! 3. VMULTD will be calculated from scratch again in the next iteration.
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

    ! Calculate the fraction of the step from D to DNEW that will be taken.
    where (vmultd < 0)
        fracmult = vmultc / (vmultc - vmultd)
    elsewhere
        fracmult = HUGENUM
    end where
    ! Only the places where VMULTD < 0 is relevant below, if any.
    frac = minval([ONE, fracmult])
    icon = int(minloc([ONE, fracmult], dim=1), kind(icon)) - 1_IK

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
        !cstrv = (ONE - frac) * cvold + frac * cstrv  ! Powell's version
        ! In theory, CSTRV = MAXVAL([B(1:M) - MATPROD(D, A(:, 1:M)), ZERO]), yet the CSTRV updated
        ! as above can be quite different from this value if A has huge entries (e.g., > 1E20).
        cstrv = maxval([b(1:m) - matprod(d, A(:, 1:m)), ZERO])
    end if

    if (icon < 1 .or. icon > mcon) then
        ! In Powell's code, the condition is ICON == 0. Indeed, ICON < 0 cannot hold unless
        ! FRACMULT contains only NaN, which should not happen; ICON > MCON can never occur.
        exit
    end if
end do

! Due to rounding errors, it can happen that |D| > DELTA. We brutally scale D down if |D| > 2*DELTA.
if (norm(d) > TWO * delta) then
    d = (d / norm(d)) * delta
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(iact) == mcon, 'SIZE(IACT) == MCON', srname)
    call assert(size(vmultc) == mcon, 'SIZE(VMULTC) == MCON', srname)
    call assert(all(vmultc >= 0), 'VMULTC >= 0', srname)
    call assert(size(d) == n, 'SIZE(D) == N', srname)
    call assert(all(is_finite(d)), 'D is finite', srname)
    call assert(norm(d) <= TWO * delta, '|D| <= 2*DELTA', srname)
    call assert(size(z, 1) == n .and. size(z, 2) == n, 'SIZE(Z) == [N, N]', srname)
    call assert(nact >= 0 .and. nact <= min(mcon, n), '0 <= NACT <= MIN(MCON, N)', srname)
end if

end subroutine trstlp_sub


function trrad(delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio) result(delta)
!--------------------------------------------------------------------------------------------------!
! This function updates the trust region radius according to RATIO and DNORM.
!--------------------------------------------------------------------------------------------------!

! Generic module
use, non_intrinsic :: consts_mod, only : RP, HALF, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Input
real(RP), intent(in) :: delta_in   ! Current trust-region radius
real(RP), intent(in) :: dnorm   ! Norm of current trust-region step
real(RP), intent(in) :: eta1    ! Ratio threshold for contraction
real(RP), intent(in) :: eta2    ! Ratio threshold for expansion
real(RP), intent(in) :: gamma1  ! Contraction factor
real(RP), intent(in) :: gamma2  ! Expansion factor
real(RP), intent(in) :: ratio   ! Reduction ratio

! Outputs
real(RP) :: delta

! Local variables
character(len=*), parameter :: srname = 'TRRAD'

! Preconditions
if (DEBUGGING) then
    call assert(delta_in >= dnorm .and. dnorm > 0, 'DELTA_IN >= DNORM > 0', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(eta1 >= 0 .and. eta1 <= eta2 .and. eta2 < 1, '0 <= ETA1 <= ETA2 < 1', srname)
    call assert(gamma1 > 0 .and. gamma1 < 1 .and. gamma2 > 1, '0 < GAMMA1 < 1 < GAMMA2', srname)
    ! By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
    ! NaN, which should NOT happen due to the moderated extreme barrier.
    call assert(.not. is_nan(ratio), 'RATIO is not NaN', srname)
end if

!====================!
! Calculation starts !
!====================!

if (ratio <= eta1) then
    delta = gamma1 * dnorm
elseif (ratio <= eta2) then
    delta = max(HALF * delta_in, dnorm)
else
    delta = max(HALF * delta_in, gamma2 * dnorm)  ! Powell' version (taken from NEWUOA)
    !delta = max(delta_in, gamma2 * dnorm)  ! Modified version
    ! For noise-free CUTEst problems of <= 100 variables, Powell's version works slightly better
    ! than the modified one.
end if

! For noisy problems, the following may work better.
!!if (ratio <= eta1) then
!!    delta = gamma1 * dnorm
!!elseif (ratio <= eta2) then  ! Ensure DELTA >= DELTA_IN
!!    delta = delta_in
!!else  ! Ensure DELTA > DELTA_IN with a constant factor
!!    delta = max(delta_in * (1.0_RP + gamma2) / 2.0_RP, gamma2 * dnorm)
!!end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(delta > 0, 'DELTA > 0', srname)
end if

end function trrad


end module trustregion_mod
