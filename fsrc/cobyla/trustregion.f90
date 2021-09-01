module trustregion_mod

implicit none
private
public :: trstlp


contains

function trstlp(A, b, rho) result(d)
! This subroutine calculates an N-component vector D by the following two stages. In the first
! stage, D is set to the shortest vector that minimizes the greatest violation of the constraints
!       dot_product(A(1:N, K), D) >= B(K),  K= 1, 2, 3, ..., M,
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
! NACT is the number of constraints in the active set and IACT(1),...,IACT(NACT) are their indices,
! while the remainder of IACT contains a permutation of the remaining constraint indices.
! N.B.: NACT <= min(M, N). Obviously, NACT <= M. In addition, The constraints in IACT(1, ..., NACT)
! have linearly independent gradients (see the comments above the instructions that delete a
! constraint from the active set to make room for the new active constraint with index IACT(ICON));
! it can also be seen from the update of NACT: starting from 0, NACT is incremented only if NACT < N.
!
! Further, Z is an orthogonal matrix whose first NACT columns can be regarded as the result of
! Gram-Schmidt applied to the active constraint gradients. For J=1, 2, ..., NACT, the number
! ZDOTA(J) is the scalar product of the J-th column of Z with the gradient of the J-th active
! constraint. D is the current vector of variables and here the residuals of the active constraints
! should be zero. Further, the active constraints have nonnegative Lagrange multipliers that are
! held at the beginning of VMULTC. The remainder of this vector holds the residuals of the inactive
! constraints at d, the ordering of the components of VMULTC being in agreement with the permutation
! of the indices of the constraints that is in IACT. All these residuals are nonnegative, which is
! achieved by the shift CSTRV that makes the least residual zero.

! Generic modules
use consts_mod, only : RP, IK, DEBUGGING
use debug_mod, only : errstop, verisize

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
use consts_mod, only : RP, IK, ZERO, ONE, EPS, HUGENUM, DEBUGGING
use infnan_mod, only : is_nan, is_finite
use debug_mod, only : errstop, verisize
use lina_mod, only : inprod, matprod, eye, planerot, isminor

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
real(RP) :: cgzk
real(RP) :: cgzkabs
real(RP) :: cstrv
real(RP) :: cvold
real(RP) :: cvsabs(size(b))
real(RP) :: cvshift(size(b))
real(RP) :: dd
real(RP) :: dnew(size(d))
real(RP) :: dold(size(d))
real(RP) :: dtmp(size(d))
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
real(RP) :: zdd
real(RP) :: zddabs
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
    z = eye(n, n)
    if (mcon == 0 .or. cstrv <= ZERO) then
        ! Check whether a quick return is possible. Make sure the In-outputs have been initialized.
        return
    end if

    m = mcon
    icon = maxloc(b, dim=1)
    sdirn = ZERO
else
    if (inprod(d, d) >= rho * rho) then
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
        cgrad = A(:, iact(icon))
        cgz = matprod(cgrad, z)
        cgzabs = matprod(abs(cgrad), abs(z))
        where (isminor(cgz, cgzabs))  ! Code in MATLAB: CGZ(ISMINOR(CGZ, CGZABS)) = ZERO
            cgz = ZERO
        end where
        do k = n - 1, nact + 1, -1
            ! Apply a 2D Givens rotation to Z(:, [K,K+1]) from the right to zero CGRAD'*Z(:, K+1) out.
            if (abs(cgz(k + 1)) > ZERO) then
                ! Powell wrote CGZ(K + 1) /= ZERO instead of ABS(CGZ(K + 1)) > ZERO. The two
                ! conditions differ if CGZ(K + 1) is NaN.
                grot = planerot(cgz([k, k + 1]))
                z(:, [k, k + 1]) = matprod(z(:, [k, k + 1]), transpose(grot))
                cgz(k) = sqrt(cgz(k)**2 + cgz(k + 1)**2)
            end if
        end do

        if (nact < n .and. abs(cgz(nact + 1)) > ZERO) then
            ! Add the new constraint if this can be done without a deletion from the active set.
            ! Powell wrote "CGZ(NACT + 1) /= ZERO" instead of "ABS(CGZ(NACT + 1)) > ZERO".
            nact = nact + 1
            vmultc([icon, nact]) = [vmultc(nact), ZERO]
            iact([icon, nact]) = iact([nact, icon])
            zdota(nact) = cgz(nact)  ! Indeed, ZDOTA(NACT) = INPROD(Z(:, NACT), A(:, IACT(NACT)))
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
            do k = nact, 1, -1  ! NACT >= N >= 1
                cgzk = inprod(cgrad, z(:, k))
                cgzkabs = inprod(abs(cgrad), abs(z(:, k)))
                if (isminor(cgzk, cgzkabs)) then
                    vmultd(k) = ZERO
                else
                    vmultd(k) = cgzk / zdota(k)
                end if
                cgrad = cgrad - vmultd(k) * A(:, iact(k))
            end do

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
            !    do k = icon, nact - 1
            !        hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:, iact(k + 1)))**2)
            !        grot = planerot([zdota(k + 1), inprod(z(:, k), A(:, iact(k + 1)))])
            !        z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), transpose(grot))
            !        zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) * zdota(k)]
            !    end do
            !    iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
            !    vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(icon)]
            !end if
            !--------------------------------------------------------------------------------------!

            zda = inprod(z(:, nact), A(:, iact(icon)))
            if (abs(zda) > ZERO) then
                vmultc([icon, nact]) = [ZERO, frac]
                iact([icon, nact]) = iact([nact, icon])
                zdota(nact) = zda  ! Indeed, ZDOTA(NACT) = INPROD(Z(:, NACT), A(:, IACT(NACT)))
            else
                exit
            end if
        end if

        ! Ensure that the objective function continues to be treated as the last active constraint
        ! if stage 2 is in progress.
        if (stage == 2 .and. iact(nact) /= mcon) then
            ! HYPT is positive because ZDOTA(NACT) is nonzero.
            hypt = sqrt(zdota(nact)**2 + inprod(z(:, nact - 1), A(:, iact(nact)))**2)
            grot = planerot([zdota(nact), inprod(z(:, nact - 1), A(:, iact(nact)))])
            z(:, [nact - 1, nact]) = matprod(z(:, [nact, nact - 1]), transpose(grot))
            zdota([nact - 1, nact]) = [hypt, (zdota(nact) / hypt) * zdota(nact - 1)]
            iact([nact - 1, nact]) = iact([nact, nact - 1])
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
        if (icon < nact .and. icon > 0) then
            do k = icon, nact - 1
                ! Zaikun 20210811: What if HYPT = 0?
                hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:, iact(k + 1)))**2)
                grot = planerot([zdota(k + 1), inprod(z(:, k), A(:, iact(k + 1)))])
                z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), transpose(grot))
                zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) * zdota(k)]
            end do
            iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
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
        step = dd / (sqrt(ss * dd + sd * sd) + sd)
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

    ! Set VMULTD to the VMULTC vector that would occur if D became DNEW. A device is included to
    ! force VMULTD(K)=ZERO if deviations from this value can be attributed to computer rounding
    ! errors. First calculate the new Lagrange multipliers.
    dtmp = dnew  ! Use DTMP instead of DNEW for the calculation, retaining DNEW for later usage.
    do k = nact, 1, -1
        zdd = inprod(z(:, k), dtmp)
        zddabs = inprod(abs(z(:, k)), abs(dtmp))
        ! Powell's original code sets ZDD = 0 when ISMINOR(ZDD, ZDDABS) = TRUE, and then takes
        ! VMULTD(K) = ZDD/ZDOTA(K), which is NaN if ZDD = 0 = ZDOTA(K). The following code avoids NaN.
        if (isminor(zdd, zddabs)) then
            vmultd(k) = ZERO
        else
            vmultd(k) = zdd / zdota(k)
        end if
        dtmp = dtmp - vmultd(k) * A(:, iact(k))
    end do
    if (stage == 2) then
        vmultd(nact) = max(ZERO, vmultd(nact))
    end if

    ! Complete VMULTD by finding the new constraint residuals. (Powell wrote "Complete VMULTC ...")
    cvshift = matprod(dnew, A(:, iact)) - b(iact) + cstrv  ! Only CVSHIFT(nact+1:mcon) is needed.
    cvsabs = matprod(abs(dnew), abs(A(:, iact))) + abs(b(iact)) + cstrv
    where (isminor(cvshift, cvsabs))
        cvshift = ZERO
    end where
    vmultd(nact + 1:mcon) = cvshift(nact + 1:mcon)

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

end module trustregion_mod
