module trustregion_mod

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

!---------------------------------------------------------------------------!
!-- QUESTION: What are exactly the objective and algorithm of trstlp_sub? --!
! The algorithm was NOT documented in the COBYLA paper. A note should be
! written to introduce it!
! As a major part of the algorithm, the code maintains and updates the QR
! factorization of A(IACT(1:NACT)), i.e., the gradients of all the active
! (linear) constraints. The matrix Z is indeed Q, and the vector ZDOTA
! is the diagonal of R. The factorization is updated by Givens rotations when
! an index is added in or removed from IACT.
! Zaikun 20211011: This subroutine SHOULD BE rewritten to make the update of
! the QR factorization more manifested. Include subroutines that update the
! factorization when IACT is changed by one element !!!!!!!!!!!!!!!!!!!!!!!!!
! This is related to the VMD function (see the end of the module). Also, note
! that in stage 2 the last index in IACT is always M+1, corresponding to the
! (linear) objective function.
! The following functions will be useful, and they can be included into linalg:
! qradd(A, Q, a, R optional, Rdiag optional) : add a new column
! qrdel(A, Q, i, R optional, Rdiag optional) : remove a column
! qrexc(A, Q, i, j, R optional, Rdiag optional) : exchange two columns
! lsqr(A, x, Q optional, R optional, Rdiag options) : linear least squares
!---------------------------------------------------------------------------!

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
!!!!!! use, non_intrinsic:: linalg_mod, only : inprod, matprod, eye, planerot, isminor
use, non_intrinsic :: linalg_mod, only : inprod, matprod, eye, isminor !!!!!!!!!!!!!!!!!!
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
integer(IK) :: nfail
real(RP) :: cgrad(size(d))
real(RP) :: cgz(size(d))
real(RP) :: cgzabs(size(d))
real(RP) :: cstrv
real(RP) :: cvold
real(RP) :: cvsabs(size(b))
real(RP) :: cvshift(size(b))
real(RP) :: dd
real(RP) :: dnew(size(d))
real(RP) :: dold(size(d))
real(RP) :: frac
real(RP) :: frtmp(size(b))
real(RP) :: grot(2, 2)
real(RP) :: hypt
real(RP) :: optnew
real(RP) :: optold
real(RP) :: sd
real(RP) :: sdirn(size(d))
real(RP) :: ss
real(RP) :: step
real(RP) :: vmultd(size(b))
real(RP) :: zda
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
    ! them so that they need to be passed from stage 1 to stage 2, and hence the coupling is reduced.
    cstrv = maxval([b(1:m) - matprod(d, A(:, 1:m)), ZERO])
    zdota(1:nact) = [(inprod(z(:, k), A(:, iact(k))), k=1, nact)]
end if

! More initialization.
optold = HUGENUM
nactold = nact
nfail = 0_IK


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
        !!!!!!!!!! This is QRADD-------------------------------------------------------------------
        cgrad = A(:, iact(icon))
        cgz = matprod(cgrad, z)
        cgzabs = matprod(abs(cgrad), abs(z))
        where (isminor(cgz, cgzabs))  ! Code in MATLAB: CGZ(ISMINOR(CGZ, CGZABS)) = ZERO
            cgz = ZERO
        end where
        do k = n - 1, nact + 1, -1
            ! Apply a 2D Givens rotation to Z(:, [K, K+1]) from the right to zero CGRAD'*Z(:, K+1) out.
            if (abs(cgz(k + 1)) > ZERO) then
                ! Powell wrote CGZ(K+1) /= ZERO instead of ABS(CGZ(K+1)) > ZERO. The two
                ! conditions differ if CGZ(K+1) is NaN.
                grot = PLANEROT_TMP(cgz([k, k + 1]))
                z(:, [k, k + 1]) = matprod(z(:, [k, k + 1]), transpose(grot))
                cgz(k) = sqrt(cgz(k)**2 + cgz(k + 1)**2)
            end if
        end do
        !!!!!!!!!! This is QRADD-------------------------------------------------------------------

        if (nact < n .and. abs(cgz(nact + 1)) > ZERO) then
            ! Add the new constraint if this can be done without a deletion from the active set.
            ! Powell wrote "CGZ(NACT+1) /= ZERO" instead of "ABS(CGZ(NACT+1)) > ZERO".
            ! Zaikun 20211012: Shouldn't this be QREXC? No, because A(:, NACT) is not in the range
            ! of the QR factorization (before the factorization).
            ! Zaikun 20211012: Is it true that IACT(NACT) is always the constraint with the largest
            ! vioaltion (in stage 1)?
            nact = nact + 1
            if (nact /= icon) then
                ! N.B.: It is problematic to index arrays using [NACT, ICON] when NACT == ICON.
                ! Indeed, Sec. 9.5.3.3.3 of the Fortran 2018 standard says: "If a vector subscript
                ! has two or more elements with the same value, an array section with that vector
                ! subscript is NOT definable and shall NOT be defined or become undefined."
                ! Zaikun 20211012: Why should VMULTC(NACT) = 0?
                vmultc([icon, nact]) = [vmultc(nact), ZERO]
                iact([icon, nact]) = iact([nact, icon])
            else
                ! Zaikun 20211012: Is this needed ? Isn't it true naturally?
                vmultc(nact) = ZERO
            end if
            !!!!!!!!!!!!!!!!!!! This is also QRADD------------------------------------------------
            !!!!!!!!!!!!!!!!!!!!!!! Calculate ZDOTA from scratch? Test it.
            zdota(nact) = cgz(nact)  ! Indeed, ZDOTA(NACT) = INPROD(Z(:, NACT), A(:, IACT(NACT)))
            !zdota(nact) = inprod(z(:, nact), a(:, iact(nact)))
            !!!!!!!!!!!!!!!!!!!!!!! Test
            !!!!!!!!!!!!!!!!!!! This is also QRADD------------------------------------------------
        else
            ! The next instruction is reached if a deletion has to be made from the active set in
            ! order to make room for the new active constraint, because the new constraint gradient
            ! is a linear combination of the gradients of the old active constraints. Set the
            ! elements of VMULTD to the multipliers (i.e., coefficients) of the linear combination.
            !
            ! Zaikun 20210811: Powell wrote the following comment, but IOUT is NEVER DEFINED. It
            ! seems that Powell's code always deletes the constraint with index IACT(NACT).
            !--------------------------------------------------------------------------------------!
            ! Further, set IOUT to the index of the constraint to be deleted, but branch if no
            ! suitable index can be found.
            !--------------------------------------------------------------------------------------!

            !----------------------------! 1st VMULTD CALCULATION STARTS  !-------------------------!
            ! Zaikun 20211011:
            ! 1. VMULTD is calculated from scratch for the first time (out of 2) in one iteration.
            ! 2. The sole purpose of this VMULTD is to compute VMULTC(1 : NACT) and check possible
            ! exit immediately after this loop. Only VMULTD(1 : NACT) is needed.
            ! 3. VMULTD will be computed from scratch again later.
            vmultd(1:nact) = vmd(cgrad, A(:, iact(1:nact)), z(:, 1:nact), zdota(1:nact))
            !!!vmultd(1:nact) = vmd(cgrad, A(:, iact(1:nact)), z(:, 1:nact))
            !----------------------------! 1st VMULTD CALCULATION ENDS  !--------------------------!

            frac = minval(vmultc(1:nact) / vmultd(1:nact), mask=(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m))
            if (frac < ZERO .or. .not. any(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m)) then
                exit
            end if

            ! Revise the Lagrange multipliers and reorder the active constraints so that the one to
            ! be replaced is at the end of the list. Also calculate the new value of ZDOTA(NACT) and
            ! branch if it is not acceptable.
            vmultc(1:nact) = max(ZERO, vmultc(1:nact) - frac * vmultd(1:nact))

            ! Zaikun 20210811: Powell's code includes the following, but it is IMPOSSIBLE TO REACH.
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

            zda = inprod(z(:, nact), A(:, iact(icon)))
            if (abs(zda) > ZERO) then
                ! Since ICON /= NACT, it is valid to use [ICON, NACT] to index arrays.
                vmultc([icon, nact]) = [ZERO, frac]
                iact([icon, nact]) = iact([nact, icon])
                zdota(nact) = zda  ! Indeed, ZDOTA(NACT) = INPROD(Z(:, NACT), A(:, IACT(NACT)))
            else
                exit
            end if
        end if

        ! Ensure that the objective function continues to be treated as the last active constraint
        ! if stage 2 is in progress.
        ! Zaikun 20211011: Is it for stage 2 that IACT(NACT-1) = MCON when IACT(NACT) /= MCON?
        if (stage == 2 .and. iact(nact) /= mcon) then
            !!!!!! This is QREXC-----------------------------------------------------------------
            ! HYPT is positive because ZDOTA(NACT) is nonzero.
            hypt = sqrt(zdota(nact)**2 + inprod(z(:, nact - 1), A(:, iact(nact)))**2)
            grot = PLANEROT_TMP([zdota(nact), inprod(z(:, nact - 1), A(:, iact(nact)))])
            z(:, [nact - 1, nact]) = matprod(z(:, [nact, nact - 1]), transpose(grot))

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate ZDOTA from scratch? Test it.
            zdota([nact - 1, nact]) = [hypt, (zdota(nact) / hypt) * zdota(nact - 1)]   ! zdota([nact-1, nact]) = [hypt, grot(1, 1)* zdota(nact-1)]
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Test

            iact([nact - 1, nact]) = iact([nact, nact - 1])

            !!!!!!! Test
            !zdota(nact-1:nact) = [(inprod(z(:, k), A(:, iact(k))), k = nact-1, nact)]
            !!!!!!! Test
            !!!!!! This is QREXC-----------------------------------------------------------------

            vmultc([nact - 1, nact]) = vmultc([nact, nact - 1])
        end if

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
        ! Delete the constraint with the index IACT(ICON) from the active set. In theory, ICON > 0.
        ! To be safe, the condition below requires ICON > 0, which does not exist in Powell's code.

        !!!!!!!!!!!! This is QRDEL-----------------------------------------------------------------
        if (icon < nact .and. icon > 0) then
            do k = icon, nact - 1
                ! Zaikun 20210811: What if HYPT = 0?
                hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:, iact(k + 1)))**2)
                grot = PLANEROT_TMP([zdota(k + 1), inprod(z(:, k), A(:, iact(k + 1)))])
                z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), transpose(grot))

                !!!!!!!!! Calculate ZDOTA from scratch? Test it.
                zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) * zdota(k)]
                !!!!!!!!! Test

            end do
            iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]

            !!!!!!!!!!!!!! Test
            !zdota(icon:nact) = [(inprod(z(:, k), A(:, iact(k))), k = icon, nact)]
            !!!!!!!!!!!!!! Test
        !!!!!!!!!!!! This is QRDEL-----------------------------------------------------------------

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

    ! Set VMULTD to the VMULTC vector that would occur if D became DNEW. A device is included to
    ! force VMULTD(K)=ZERO if deviations from this value can be attributed to computer rounding
    ! errors. First calculate the new Lagrange multipliers.
    vmultd(1:nact) = vmd(dnew, A(:, iact(1:nact)), z(:, 1:nact), zdota(1:nact))
    !!!vmultd(1:nact) = vmd(dnew, A(:, iact(1:nact)), z(:, 1:nact))
    if (stage == 2) then
        vmultd(nact) = max(ZERO, vmultd(nact))  ! This seems a safeguard never activated.
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
    frtmp = vmultc / (vmultc - vmultd)  !
    if (any(vmultd < ZERO .and. .not. is_nan(frtmp))) then
        frac = min(ONE, minval(frtmp, mask=(vmultd < ZERO)))
    else
        frac = ONE
    end if

    if (frac < ONE) then
        icon = minloc(frtmp, mask=(vmultd < ZERO), dim=1)
    else
        icon = 0_IK
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

!--------------------------------------------------------------------------------------------------!
! The following function indeed solves the linear least squares problem
! min |V*VMULTD-u|
! 1. Z is the Q in the QR factorization of V, and zdotv is the diagonal of R.
! 2. In COBYLA, V is of full column rank.
! 3. In COBYLA, it seems that u (CGRAD and DNEW) is in the column space of V. Not sure yet.
! 4. It should be included into linalg_mod as
! function lsqr(A, b, Q = Q)
! maybe with a bit generalization.
!--------------------------------------------------------------------------------------------------!


function vmd(u, V, Z, zdotv) result(vmultd)
!function vmd(u, V, Z) result(vmultd)
!--------------------------------------------------------------------------------------------------!
! This function calculates VMULTD(1:NACT) for a vector U. Here,
! V represents A(:, IACT(1:NACT))
! Z represents Z(:, 1:NACT)
! ZDOTV represents ZDOTA(1:NACT), and it equals DIAG(Z^T*V) so that is should be calculated
! internally!!!
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : isminor, inprod, matprod, eye, norm
implicit none

! Inputs
real(RP), intent(in) :: u(:)
real(RP), intent(in) :: V(:, :)
real(RP), intent(in) :: Z(:, :)
real(RP), intent(in) :: zdotv(:)

! Outputs
real(RP) :: vmultd(size(V, 2))

! Local variables
character(len=*), parameter :: srname = 'VMD'
integer(IK) :: k
real(RP), parameter :: tol = 1.0E4_RP * sqrt(EPS)
real(RP) :: w(size(u))
real(RP) :: wzk
real(RP) :: wzkabs
!!!real(RP) :: zdotv(size(V, 2))

! Preconditions
if (DEBUGGING) then
    call assert(size(u) == size(V, 1), 'SIZE(U) == SIZE(V, 1)', srname)
    call assert(size(Z, 1) == size(V, 1) .and. size(Z, 2) == size(V, 2), 'SIZE(Z) == SIZE(V)', srname)
    call assert(size(zdotv) == size(V, 2), 'SIZE(ZDOTV) == SIZE(V, 2)', srname)
    !call assert(maxval(abs(matprod(transpose(Z), Z) - eye(size(Z, 2)))) &
    !    & <= tol*real(size(Z, 2), RP), 'The columns of Z are orthogonal', srname)
    !call assert(maxval(abs(zdotv - [(inprod(Z(:, k), V(:, k)), k = 1, size(zdotv))])) &
    !    & <= tol*max(maxval(abs(V)), ONE), 'ZDOTV = DIAG(Z^T*V)', srname)
end if

w = u  ! Local copy of U; U is INTENT(IN) and should not be modified.
!!!zdotv = [(inprod(Z(:, k), V(:, k)), k=1, size(zdotv))]

do k = int(size(V, 2), kind(k)), 1, -1
    wzk = inprod(w, Z(:, k))
    wzkabs = inprod(abs(w), abs(Z(:, k)))
    ! Powell's original code sets UZK = 0 when ISMINOR(UZK, UZKABS) = TRUE, and then takes
    ! VMULTD(K) = UZK/ZDOTV(K), which is NaN if UZK = 0 = ZDOTV(K) (in precise arithmetic, this
    ! cannot happen). The following code avoids NaN.
    if (isminor(wzk, wzkabs)) then
        vmultd(k) = ZERO
    else
        vmultd(k) = wzk / zdotv(k)
    end if
    w = w - vmultd(k) * V(:, k)
end do

if (DEBUGGING) then
    !call assert(norm(matprod(u-matprod(V, vmultd), V)) <= tol*max(norm(matprod(u, V)), ONE), &
    !    & 'V*VMULTD is the projection of U to the column space of V', srname)
end if

end function vmd

!!!!!! This is temporary!!! It must be merged with the PLANEROT in linalg_mod.
function PLANEROT_TMP(x) result(G)
! As in MATLAB, PLANEROT(X) returns a 2x2 Givens matrix G for X in R^2 so that Y = G*X has Y(2) = 0.
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, HUGENUM
use, non_intrinsic :: linalg_mod, only : eye
use, non_intrinsic :: debug_mod, only : verisize

implicit none
real(RP), intent(in) :: x(:)
real(RP) :: G(2, 2)

real(RP) :: c, s, r, scaling

call verisize(x, 2_IK)

if (abs(x(2)) > ZERO) then
    !>>>> SHOULD BE
    !r = sqrt(sum(x**2))
    !if (r > ZERO .and. r <= HUGENUM) then
    !    c = x(1) / r
    !    s = x(2) / r
    !else
    !    scaling = maxval(abs(x))
    !    r = sqrt(sum((x/scaling)**2))
    !    c = (x(1) / scaling) / r
    !    s = (x(2) / scaling) / r
    !end if
    !G = reshape([c, -s, s, c], [2, 2])
    !>>>> SHOUD BE

    !<<<< TEMP
    r = sqrt(sum(x**2))
    c = x(1) / r
    s = x(2) / r
    G = reshape([c, -s, s, c], [2, 2])
    !<<<< TEMP

elseif (x(1) < ZERO) then
    ! Setting G = -EYE(2, 2) in this case ensures the continuity of G with respect to X except at 0.
    G = -eye(2_IK)
else
    G = eye(2_IK)
end if
end function PLANEROT_TMP


! A stabilized Givens rotation that avoids over/underflow and keeps continuity (see wikipedia)
!function [c, s, r] = givens_rotation(a, b)
!    if b == 0;
!        c = sign(a);  % Continuity
!        if (c == 0);
!            c = 1.0; % Unlike other languages, MatLab's sign function returns 0 on input 0.
!        end;
!        s = 0;
!        r = abs(a);
!    elseif a == 0;
!        c = 0;
!        s = sign(b);
!        r = abs(b);
!    elseif abs(a) > abs(b);
!        t = b/a;
!        u = sign(a) * sqrt(1+t * t);
!        c = 1/u;
!        s = c*t;
!        r = a*u;
!    else
!        t = a/b;
!        u = sign(b) * sqrt(1+t * t);
!        s = 1/u;
!        c = s*t;
!        r = b*u;
!    end
!end

end module trustregion_mod
