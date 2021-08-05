module trustregion_mod

contains

subroutine trstlp(n, m, A, b, rho, d, ifull, iact)

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
use infnan_mod, only : is_nan, is_posinf, is_finite
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : inprod, matprod, eye, planerot, isminor

implicit none

integer(IK), intent(in) :: n
integer(IK), intent(in) :: m
real(RP), intent(in) :: A(:, :)  !(n, m+1)
real(RP), intent(in) :: b(:)
real(RP), intent(in) :: rho
real(RP), intent(inout) :: d(:)
integer(IK), intent(out) :: ifull
integer(IK), intent(out) :: iact(:)

real(RP) :: hypt
real(RP) :: z(n, n)
real(RP) :: zdota(n)
real(RP) :: vmultc(m + 1)
real(RP) :: sdirn(n)
real(RP) :: vmultd(m + 1)
real(RP) :: cgrad(n)
real(RP) :: cgz(n)
real(RP) :: cgzabs(n)
real(RP) :: cgzk
real(RP) :: cgzkabs
real(RP) :: dnew(n)
real(RP) :: tot, cnew(n)

real(RP) :: alpha
real(RP) :: beta
real(RP) :: dd
real(RP) :: grot(2, 2)
real(RP) :: optnew
real(RP) :: optold
real(RP) :: ratio
real(RP) :: cstrv
real(RP) :: cvold
real(RP) :: sd
real(RP) :: sp
real(RP) :: spabs
real(RP) :: ss
real(RP) :: step
real(RP) :: stpful
real(RP) :: summ
real(RP) :: summabs
real(RP) :: summd
real(RP) :: temp
real(RP) :: tempa
real(RP) :: tmpv(m + 1)
real(RP) :: tmpvabs(m + 1)
real(RP) :: vsave
real(RP) :: zdotw
real(RP) :: zdvabs
real(RP) :: zdwabs
real(RP) :: dsav(size(d))  ! N
integer(IK) :: i
integer(IK) :: icon
integer(IK) :: icount
integer(IK) :: isave
integer(IK) :: iter
integer(IK) :: j
integer(IK) :: k
integer(IK) :: kk
integer(IK) :: kl
integer(IK) :: kw
integer(IK) :: maxiter
integer(IK) :: mcon
integer(IK) :: nact
integer(IK) :: nactx

real(RP) :: ACCA, ACCB
integer(IK) :: KP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
! It is possible that a degeneracy may prevent D from attaining the target length RHO. Then the
! value IFULL=0 would be set, but usually IFULL=1 on return.
!
! In general NACT is the number of constraints in the active set and IACT(1),...,IACT(NACT) are
! their indices, while the remainder of IACT contains a permutation of the remaining constraint
! indices.
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
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
iter = 0
maxiter = min(10000_IK, 100_IK * n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IFULL = 1
MCON = M
NACT = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      CSTRV=0.0
CSTRV = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! Question: what is the size of B? M, M+1, or M+2?
! It seems to be M+1. Then when is B(M+1) used?
! What is the size of iact, vmultc? m or m+1?
! What is the upper bound of NACT? M or N, or both?

! Initialize Z and some other variables. The value of CSTRV will be appropriate to D=0, while ICON
! will be the index of a most violated constraint if CSTRV is positive. Usually during the first
! stage the vector SDIRN gives a search direction that reduces all the active constraint violations
! by one simultaneously.
!!!! What is size (number of columns) of Z and related variables (CGZ, ZDOTA, ZDOTW)???????
z = eye(n, n)
d = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cstrv = maxval([b(1:m), ZERO])
icon = maxloc(b(1:m), dim=1)
iact(1:m) = [(k, k=1, m)]
vmultc(1:m) = cstrv - b(1:m)

if (cstrv <= zero) goto 480

sdirn = ZERO

dsav = d   !!!! HOW TO AVOID THIS???


! End the current stage of the calculation if 3 consecutive iterations have either failed to reduce
! the best calculated value of the objective function or to increase the number of active
! constraints since the best value was calculated. This strategy prevents cycling, but there is
! a remote possibility that it will cause premature termination.
!
60 optold = ZERO

icount = 0_IK
70 if (mcon == m) then
    optnew = cstrv
else
    optnew = -inprod(d, A(:, mcon))
end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019
! The original code can still encounter infinite cycling,
! which did happen when testing the following CUTEst problems:
! DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in A (up to 10^219).
! To avoid wasting energy, we do the following.

if (is_finite(sum(abs(d)))) then
    dsav = d
else
    d = dsav
    return
end if
iter = iter + 1
if (iter > maxiter) then
    return
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! What is ICOUNT ????
if (icount == 0 .or. optnew < optold) then
    optold = optnew
    nactx = nact
    icount = 3
else if (nact > nactx) then
    nactx = nact
    icount = 3
else
    icount = icount - 1
    if (icount == 0) goto 490
end if


! If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to the active set. Apply
! Givens rotations so that the last N-NACT-1 columns of Z are orthogonal to the gradient of the new
! constraint, a scalar product being set to zero if its nonzero value could be due to computer
! rounding errors.
if (icon <= nact) goto 260
kk = iact(icon)
cgrad = A(:, iact(icon))
cgz = matprod(cgrad, z)
cgzabs = matprod(abs(cgrad), abs(z))
! The following code in MATLAB: CGZ(ISMINOR(CGZ, CGZABS)) = ZERO
where (isminor(cgz, cgzabs))
    cgz = ZERO
end where
do k = n - 1, nact + 1, -1
    ! Apply a 2x2 Givens rotation to Z(:, [K,K+1]) from the right so that CGRAD'*Z(:, K+1) becomes 0.
    if (abs(cgz(k + 1)) > ZERO) then
        ! Powell wrote "CGZ(K + 1) /= ZERO" instead of "ABS(CGZ(K + 1)) > ZERO". The two conditions
        ! differ if CGZ(K + 1) is NaN.
        grot = planerot(cgz([k, k + 1]))
        z(:, [k, k + 1]) = matprod(z(:, [k, k + 1]), transpose(grot))
        cgz(k) = sqrt(cgz(k)**2 + cgz(k + 1)**2)
    end if
end do

!Add the new constraint if this can be done without a deletion from the active set.
if (nact < n .and. abs(cgz(nact + 1)) > ZERO) then
    ! Powell wrote "CGZ(NACT + 1) /= ZERO" instead of "ABS(CGZ(NACT + 1)) > ZERO", the two
    ! conditions differ if CGZ(NACT + 1) is NaN.
    nact = nact + 1
    zdota(nact) = cgz(nact)
    vmultc(icon) = vmultc(nact)
    vmultc(nact) = ZERO
    goto 210
end if

! The next instruction is reached if a deletion has to be made from the active set in order to make
! room for the new active constraint, because the new constraint gradient is a linear combination of
! the gradients of the old active constraints. Set the elements of VMULTD to the multipliers of the
! linear combination. Further, set IOUT to the index of the constraint to be deleted, but branch if
! no suitable index can be found.

!!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! The following part seems to lead to a memory error when VMULTD(K) is accessed (by why not Z(I, K) and others?????)
do k = nact, 1, -1
    ! Is it possible that NACT = 0???? In Powell's code, the DO loop will be carried out for one time if NACT = 0.
    cgzk = inprod(cgrad, z(:, k))
    cgzkabs = inprod(abs(cgrad), abs(z(:, k)))
    if (isminor(cgzk, cgzkabs)) then
        vmultd(k) = ZERO
    else
        vmultd(k) = cgzk / zdota(k)
    end if
    if (k >= 2) then
        cgrad = cgrad - vmultd(k) * A(:, iact(k))
    end if
end do

ratio = minval(vmultc(1:nact) / vmultd(1:nact), mask=(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m))
if (ratio < ZERO .or. .not. any(vmultd(1:nact) > ZERO .and. iact(1:nact) <= m)) goto 490

! Revise the Lagrange multipliers and reorder the active constraints so that the one to be replaced
! is at the end of the list. Also calculate the new value of ZDOTA(NACT) and branch if it is not
! acceptable.
vmultc(1:nact) = max(ZERO, vmultc(1:nact) - ratio * vmultd(1:nact))

do k = icon, nact - 1
    hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:, iact(k + 1)))**2)  ! What if HYPT = 0???
    grot = planerot([zdota(k + 1), inprod(z(:, k), A(:, iact(k + 1)))])
    z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), transpose(grot))
    zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) * zdota(k)]
end do
if (icon < nact) then
    ! This seems quite rare --- it never happens on CUTEst problems with n <= 50 and m <= 5000.
    iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
    vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(icon)]
end if

if (inprod(z(:, nact), A(:, iact(icon))) == ZERO) goto 490
zdota(nact) = inprod(z(:, nact), A(:, iact(icon)))
vmultc(icon) = ZERO
vmultc(nact) = ratio

! Update IACT and ensure that the objective function continues to be treated as the last active
! constraint when MCON>M.
210 iact([icon, nact]) = iact([nact, icon])
if (mcon > m .and. kk /= mcon) then
    hypt = sqrt(zdota(nact)**2 + inprod(z(:, nact - 1), A(:, kk))**2)  ! What if HYPT = 0 ???
    grot = planerot([zdota(nact), inprod(z(:, nact - 1), A(:, kk))])
    z(:, [nact - 1, nact]) = matprod(z(:, [nact, nact - 1]), transpose(grot))
    zdota([nact - 1, nact]) = [hypt, (zdota(nact) / hypt) * zdota(nact - 1)]
    iact([nact - 1, nact]) = [kk, iact(nact - 1)]
    vmultc([nact - 1, nact]) = vmultc([nact, nact - 1])
end if

! If stage one is in progress, then set SDIRN to the direction of the next change to the current
! vector of variables.
if (mcon > m) goto 320
kk = iact(nact)
sdirn = sdirn - ((inprod(sdirn, A(:, iact(nact))) - ONE) / zdota(nact)) * z(:, nact)
goto 340

! Delete the constraint that has the index IACT(ICON) from the active set.
260 do k = icon, nact - 1
    hypt = sqrt(zdota(k + 1)**2 + inprod(z(:, k), A(:, iact(k + 1)))**2)  ! What if HYPT = 0 ???
    grot = planerot([zdota(k + 1), inprod(z(:, k), A(:, iact(k + 1)))])
    z(:, [k, k + 1]) = matprod(z(:, [k + 1, k]), transpose(grot))
    zdota([k, k + 1]) = [hypt, (zdota(k + 1) / hypt) * zdota(k)]
end do
if (icon < nact) then
    iact(icon:nact) = [iact(icon + 1:nact), iact(icon)]
    vmultc(icon:nact) = [vmultc(icon + 1:nact), vmultc(icon)]
end if
nact = nact - 1

! If stage one is in progress, then set SDIRN to the direction of the next change to the current
! vector of variables.
if (mcon > m) goto 320
sdirn = sdirn - inprod(sdirn, z(:, nact + 1)) * z(:, nact + 1)
goto 340
! Pick the next search direction of stage two.
320 sdirn = (ONE / zdota(nact)) * z(:, nact)

! Calculate the step to the boundary of the trust region or take the step  that reduces CSTRV to
! zero. The two statements below that include the factor 1.0E-6 prevent some harmless underflows
! that occurred in a test calculation. Further, we skip the step if it could be zero within
! a reasonable tolerance for computer rounding errors.
340 dd = rho**2 - sum(d**2, mask=(abs(d) >= 1.0E-6_RP * rho))
sd = inprod(sdirn, d)
ss = inprod(sdirn, sdirn)
if (dd <= ZERO) goto 490
temp = sqrt(ss * dd)
if (abs(sd) >= 1.0E-6_RP * temp) then
    temp = sqrt(ss * dd + sd * sd)
end if

stpful = dd / (temp + sd)
step = stpful
if (mcon == m) then
    if (isminor(cstrv, step)) then
        goto 480
    end if
    step = min(step, cstrv)
end if

! Set DNEW to the new variables if STEP is the steplength, and reduce CSTRV to the corresponding
! maximum residual if stage one is being done. Because DNEW will be changed during the calculation
! of some Lagrange multipliers, it will be restored to the following value later.
dnew = d + step * sdirn
if (mcon == m) then
    cvold = cstrv
    cstrv = maxval([b(iact(1:nact)) - matprod(dnew, A(:, iact(1:nact))), ZERO])
end if

! Set VMULTD to the VMULTC vector that would occur if D became DNEW. A device is included to force
! VMULTD(K)=0.0 if deviations from this value can be attributed to computer rounding errors. First
! calculate the new Lagrange multipliers.
do k = nact, 1, -1
    ! What if NACT = 0? Is it possible? Powell's code will carry out the loop for one time.
    zdotw = inprod(z(:, k), dnew)
    zdwabs = inprod(abs(z(:, k)), abs(dnew))
    if (isminor(zdotw, zdwabs)) then
        vmultd(k) = ZERO
    else
        vmultd(k) = zdotw / zdota(k)
    end if
    if (k >= 2) then
        dnew = dnew - vmultd(k) * A(:, iact(k))
    end if
end do
if (mcon > m) then
    vmultd(nact) = max(ZERO, vmultd(nact))
end if

! Complete VMULTD by finding the new constraint residuals.
dnew = d + step * sdirn
tmpv(nact + 1:mcon) = matprod(dnew, A(:, iact(nact + 1:mcon))) - b(iact(nact + 1:mcon)) + cstrv
tmpvabs(nact + 1:mcon) = matprod(abs(dnew), abs(A(:, iact(nact + 1:mcon)))) + abs(b(iact(nact + 1:mcon))) + cstrv
where (isminor(tmpv, tmpvabs))
    tmpv = ZERO
end where
vmultd(nact + 1:mcon) = tmpv(nact + 1:mcon)

! Calculate the fraction of the step from D to DNEW that will be taken.
tmpv = vmultc / (vmultc - vmultd)
ratio = min(ONE, minval(tmpv(1:mcon), mask=(vmultd(1:mcon) < ZERO)))
if (ratio < ONE) then
    icon = minloc(tmpv(1:mcon), mask=(vmultd(1:mcon) < ZERO), dim=1)
else
    icon = 0
end if

! Update d, VMULTC and CSTRV.
d = (ONE - ratio) * d + ratio * dnew
vmultc(1:mcon) = max(ZERO, (ONE - ratio) * vmultc(1:mcon) + ratio * vmultd(1:mcon))
if (MCON == M) then
    !cstrv = cvold + ratio * (cstrv - cvold)
    cstrv = (ONE - ratio) * cvold + ratio * cstrv
end if

! If the full step is not acceptable then begin another iteration. Otherwise switch to stage two or
! end the calculation.
if (icon > 0) then
    goto 70
end if
if (abs(step - stpful) <= ZERO) then
    goto 500
end if
480 mcon = m + 1
icon = mcon
iact(mcon) = mcon
vmultc(mcon) = ZERO
goto 60

! We employ any freedom that may be available to reduce the objective function before returning a D
! whose length is less than RHO.
490 if (mcon == m) goto 480
ifull = 0
500 return

end subroutine trstlp

end module trustregion_mod
