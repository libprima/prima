module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of LINCOA.
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
! Last Modified: Saturday, April 16, 2022 PM11:53:59
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstep


contains


subroutine trstep(amat, delta, gq, hq, pq, rescon, xpt, iact, nact, qfac, rfac, ngetact, snorm, step)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, HALF, EPS, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, inprod, solve, isorth, istriu, issymmetric
use, non_intrinsic :: powalg_mod, only : hess_mul

! Solver-specific modules
use, non_intrinsic :: getact_mod, only : getact

implicit none

! Inputs
real(RP), intent(in) :: amat(:, :)  ! AMAT(N, M)
real(RP), intent(in) :: delta
real(RP), intent(in) :: gq(:)  ! GQ(N)
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)  ! PQ(NPT)
real(RP), intent(in) :: rescon(:)  ! RESCON(M)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
integer(IK), intent(inout) :: iact(:)  ! IACT(M); Will be updated in GETACT
integer(IK), intent(inout) :: nact  ! Will be updated in GETACT
real(RP), intent(inout) :: qfac(:, :)  ! QFAC(N, N); Will be updated in GETACT
real(RP), intent(inout) :: rfac(:, :)  ! RFAC(N, N); Will be updated in GETACT

! Outputs
integer(IK), intent(out) :: ngetact
real(RP), intent(out) :: snorm
real(RP), intent(out) :: step(:)  ! STEP(N)

! Local variables
character(len=*), parameter :: srname = 'TRSTEP'
real(RP) :: d(size(gq))
real(RP) :: dw(size(gq))
real(RP) :: w(max(size(amat, 2), 2 * size(gq)))
real(RP) :: resact(size(amat, 2))
real(RP) :: resnew(size(amat, 2))
real(RP) :: tol
real(RP) :: g(size(gq))
real(RP) :: vlam(size(gq))
real(RP) :: ad, adw, alpbd, alpha, alphm, alpht, beta, ctest, &
&        dd, dg, dgd, ds, bstep, reduct, resmax, rhs, scaling, snsq, ss, summ, temp, wgd
integer(IK) :: i, icount, j, jsav, k
integer(IK) :: m
integer(IK) :: n
integer(IK) :: npt

! Sizes.
m = int(size(amat, 2), kind(m))
n = int(size(gq), kind(n))
npt = int(size(pq), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(size(amat, 1) == n .and. size(amat, 2) == m, 'SIZE(AMAT) == [N, M]', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)

    call assert(size(rescon) == m, 'SIZE(RESCON) == M', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(iact) == m, 'SIZE(IACT) == M', srname)
    call assert(all(iact(1:nact) >= 1 .and. iact(1:nact) <= m), '1 <= IACT <= M', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    tol = max(1.0E-10_RP, min(1.0E-1_RP, 1.0E8_RP * EPS * real(n, RP)))
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
end if

g = gq

!
!     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
!       are the same as the terms with these names in LINCOB. If RESCON(J)
!       is negative, then |RESCON(J)| must be no less than the trust region
!       radius, so that the J-th constraint can be ignored.
!     SNORM is set to the trust region radius DELTA initially. On the
!       return, however, it is the length of the calculated STEP, which is
!       set to ZERO if the constraints do not allow a long enough step.
!     STEP is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current STEP, except that
!       on return the value of G(1) is set to ONE instead of to ZERO if
!       and only if GETACT is called more than once.
!     RESNEW, RESACT, D, DW and W are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       ZERO value of RESNEW(J) indicates that the J-th constraint is active,
!       and otherwise RESNEW(J) is set to the greater of TINY and the actual
!       residual of the J-th constraint for the current STEP. RESACT holds
!       the residuals of the active constraints, which may be positive.
!       D is the search direction of each line search. DW is either another
!       search direction or the change in gradient along D. The length of W
!       must be at least MAX[M,2*N].
!
!     Set some numbers for the conjugate gradient iterations.
!
ctest = 0.01_RP
!--------------------------------------------------------------------------------------------------!
! Zaikun 20220302: The following line is added so that SNORM is INTENT(OUT) rather than
! INTNT(INOUT). Is SNORM really the norm of step, or just DELTA??? Check also GETACT.
snorm = delta
!--------------------------------------------------------------------------------------------------!
snsq = snorm * snorm

! Set the initial elements of RESNEW, RESACT and STEP.
! 1. RESNEW(J) < 0 indicates that the J-th constraint does not restrict the CG steps of the current
! trust region calculation. In other words, RESCON >= DELTA.
! 2. RESNEW(J) = 0 indicates that the J-th constraint is active, namely B(J) = AMAT(:, J)^T*(XOPT+D).
! 3. RESNEW(J) > 0 means that RESNEW(J) = max(B(J) - AMAT(:, J)^T*(XOPT+D), TINYCV).
where (rescon >= snorm)
    resnew = -ONE
elsewhere(rescon >= 0)
    resnew = max(rescon, TINYCV)
elsewhere
    resnew = rescon
end where
!!MATLAB: resnew = rescon; resnew(rescon >= 0) = max(rescon, TINYCV); resnew(rescon >= snorm) = -1;
resnew(iact(1:nact)) = ZERO
resact(1:nact) = rescon(iact(1:nact))

step = ZERO
ss = ZERO
reduct = ZERO
ngetact = 0

! GETACT picks the active set for the current STEP. It also sets DW to the vector closest to -G that
! is orthogonal to the normals of the active constraints. DW is scaled to have length 0.2*SNORM, as
! then a move of DW from STEP is allowed by the linear constraints.
40 ngetact = ngetact + 1
call getact(amat, g, snorm, iact, nact, qfac, resact, resnew, rfac, dw)
dd = inprod(dw, dw)
if (dd == ZERO) then
    goto 320
end if
scaling = 0.2_RP * snorm / sqrt(dd)
dw = scaling * dw
!do i = 1, n
!    dw(i) = scaling * dw(i)
!end do


! If the modulus of the residual of an active constraint is substantial, then set D to the shortest
! move from STEP to the boundaries of the active constraints.
resmax = ZERO
resmax = maxval([ZERO, resact(1:nact)])
!if (nact > 0) then
!    do k = 1, nact
!        resmax = max(resmax, resact(k))
!    end do
!end if

bstep = ZERO
if (resmax > 1.0D-4 * snorm) then
    !do k = 1, nact
    !    temp = resact(k)
    !    do i = 1, k - 1
    !        temp = temp - rfac(i, k) * vlam(i)
    !    end do
    !    vlam(k) = temp / rfac(k, k)
    !end do
    ! Calculation changed
    vlam(1:nact) = solve(transpose(rfac(1:nact, 1:nact)), resact(1:nact))

    !do i = 1, n
    !    d(i) = ZERO
    !    do k = 1, nact
    !        d(i) = d(i) + vlam(k) * qfac(i, k)
    !    end do
    !end do
    d = matprod(qfac(:, 1:nact), vlam(1:nact))
!
!     The vector D that has just been calculated is also the shortest move
!       from STEP+DW to the boundaries of the active constraints. Set BSTEP
!       to the greatest steplength of this move that satisfies the trust
!       region bound.
!
    !rhs = snsq
    !ds = ZERO
    !dd = ZERO
    !do i = 1, n
    !    summ = step(i) + dw(i)
    !    rhs = rhs - summ * summ
    !    ds = ds + d(i) * summ
    !    dd = dd + d(i)**2
    !end do

    ds = inprod(d, step + dw)
    dd = sum(d**2)
    rhs = snsq - sum((step + dw)**2)  ! Calculation changed

    if (rhs > ZERO) then
        temp = sqrt(ds * ds + dd * rhs)
        if (ds <= ZERO) then
            ! Zaikun 20210925
            ! What if we are at the first iteration? BLEN = DELTA/|D|? See TRSAPP.F90 of NEWUOA.
            bstep = (temp - ds) / dd
        else
            bstep = rhs / (temp + ds)
        end if
    end if
!
!     Reduce the steplength BSTEP if necessary so that the move along D
!       also satisfies the linear constraints.
!
    !j = 0
!110 if (bstep > ZERO) then
    !    j = j + 1
    !    if (resnew(j) > ZERO) then
    !        ad = ZERO
    !        adw = ZERO
    !        do i = 1, n
    !            ad = ad + amat(i, j) * d(i)
    !            adw = adw + amat(i, j) * dw(i)
    !        end do
    !        if (ad > ZERO) then
    !            temp = max((resnew(j) - adw) / ad, ZERO)
    !            bstep = min(bstep, temp)
    !        end if
    !    end if
    !    if (j < m) goto 110
    !end if

    do j = 1, m
        if (.not. bstep > 0) exit
        if (resnew(j) > ZERO) then
            !ad = ZERO
            !adw = ZERO
            !do i = 1, n
            !    ad = ad + amat(i, j) * d(i)
            !    adw = adw + amat(i, j) * dw(i)
            !end do
            ad = inprod(d, amat(:, j))
            adw = inprod(dw, amat(:, j))

            if (ad > ZERO) then
                temp = max((resnew(j) - adw) / ad, ZERO)
                bstep = min(bstep, temp)
            end if
        end if
    end do

    bstep = min(bstep, ONE)
end if
!
!     Set the next direction for seeking a reduction in the model function
!       subject to the trust region bound and the linear constraints.
!
if (bstep <= ZERO) then
    !do i = 1, n
    !    d(i) = dw(i)
    !end do
    d = dw
    icount = nact
else
    !do i = 1, n
    !    d(i) = dw(i) + bstep * d(i)
    !end do
    d = dw + bstep * d
    icount = nact - 1
end if
alpbd = ONE
!
!     Set ALPHA to the steplength from STEP along D to the trust region
!       boundary. Return if the first derivative term of this step is
!       sufficiently small or if no further progress is possible.
!
150 icount = icount + 1
rhs = snsq - ss
if (rhs <= ZERO) goto 320
!dg = ZERO
!ds = ZERO
!dd = ZERO
!do i = 1, n
!    dg = dg + d(i) * g(i)
!    ds = ds + d(i) * step(i)
!    dd = dd + d(i)**2
!end do
dg = inprod(d, g)
ds = inprod(d, step)
dd = inprod(d, d)
if (dg >= ZERO) goto 320
temp = sqrt(rhs * dd + ds * ds)
if (ds <= ZERO) then
    alpha = (temp - ds) / dd
else
    alpha = rhs / (temp + ds)
end if
if (-alpha * dg <= ctest * reduct) goto 320
!
!     Set DW to the change in gradient along D.
!

!----------------------------------------------------!
!do j = 1, n
!    dw(j) = ZERO
!    do i = 1, j
!        if (i < j) dw(j) = dw(j) + hq(i, j) * d(i)
!        dw(i) = dw(i) + hq(i, j) * d(j)
!    end do
!end do
!do k = 1, npt
!    temp = ZERO
!    do j = 1, n
!        temp = temp + xpt(j, k) * d(j)
!    end do
!    temp = pq(k) * temp
!    do i = 1, n
!        dw(i) = dw(i) + temp * xpt(i, k)
!    end do
!end do

dw = hess_mul(d, xpt, pq, hq)
!----------------------------------------------------!

!
!     Set DGD to the curvature of the model along D. Then reduce ALPHA if
!       necessary to the value that minimizes the model.
!
dgd = inprod(d, dw)
!dgd = ZERO
!do i = 1, n
!    dgd = dgd + d(i) * dw(i)
!end do
alpht = alpha
if (dg + alpha * dgd > ZERO) then
    alpha = -dg / dgd
end if
!
!     Make a further reduction in ALPHA if necessary to preserve feasibility,
!       and put some scalar products of D with constraint gradients in W.
!
alphm = alpha
jsav = 0

do j = 1, m
    ad = ZERO
    if (resnew(j) > ZERO) then
        !do i = 1, n
        !    ad = ad + amat(i, j) * d(i)
        !end do
        ad = inprod(d, amat(:, j))
        if (alpha * ad > resnew(j)) then
            alpha = resnew(j) / ad
            jsav = j
        end if
    end if
    w(j) = ad
end do



alpha = max(alpha, alpbd)
alpha = min(alpha, alphm)
if (icount == nact) alpha = min(alpha, ONE)
!
!     Update STEP, G, RESNEW, RESACT and REDUCT.
!
step = step + alpha * d
ss = sum(step**2)
g = g + alpha * dw
!do i = 1, n
!    step(i) = step(i) + alpha * d(i)
!    ss = ss + step(i)**2
!    g(i) = g(i) + alpha * dw(i)
!end do

where (resnew > 0)
    resnew = max(resnew - alpha * w(1:m), TINYCV)
end where
!!MATLAB: mask = (resnew > 0); resnew(mask) = max(resnew(mask) - alpha * w(mask), TINYCV);

!if (m > 0) then
!do j = 1, m
!    if (resnew(j) > ZERO) then
!        resnew(j) = max(resnew(j) - alpha * w(j), TINYCV)
!    end if
!end do
!end if
if (icount == nact .and. nact > 0) then
    !do k = 1, nact
    !    resact(k) = (ONE - bstep) * resact(k)
    !end do
    resact(1:nact) = (ONE - bstep) * resact(1:nact)
end if
reduct = reduct - alpha * (dg + HALF * alpha * dgd)
!
!     Test for termination. Branch to label 40 if there is a new active
!       constraint and if the distance from STEP to the trust region
!       boundary is at least 0.2*SNORM.
!
! Zaikun 2019-08-29: the code can encounter infinite cycling due to NaN
! values. Exit when NGETACT is large or NaN detected.
! Caution: 1. MIN accepts only data with the same KIND; 2. Integer overflow.
if (ngetact > min(10000, 100 * int(m + 1) * int(n)) .or.  &
& alpha /= alpha .or. alpht /= alpht .or. &
& alphm /= alphm .or. dgd /= dgd .or. dg /= dg .or. &
& ss /= ss .or. snsq /= snsq .or. reduct /= reduct) then
    goto 320
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (alpha == alpht) goto 320
temp = -alphm * (dg + HALF * alphm * dgd)
if (temp <= ctest * reduct) goto 320
if (jsav > 0) then
    if (ss <= 0.64_RP * snsq) goto 40
    goto 320
end if
if (icount == n) goto 320
!
!     Calculate the next search direction, which is conjugate to the
!       previous ONE except in the case ICOUNT=NACT.
!
if (nact > 0) then
    !do j = nact + 1, n
    !    w(j) = ZERO
    !    do i = 1, n
    !        w(j) = w(j) + g(i) * qfac(i, j)
    !    end do
    !end do
    !do i = 1, n
    !    temp = ZERO
    !    do j = nact + 1, n
    !        temp = temp + qfac(i, j) * w(j)
    !    end do
    !    w(n + i) = temp
    !end do
    w(n + 1:2 * n) = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
else
    !do i = 1, n
    !    w(n + i) = g(i)
    !end do
    w(n + 1:2 * n) = g
end if
if (icount == nact) then
    beta = ZERO
else
    wgd = inprod(w(n + 1:2 * n), dw)
    !wgd = ZERO
    !do i = 1, n
    !   wgd = wgd + w(n + i) * dw(i)
    !end do
    beta = wgd / dgd
end if
d = -w(n + 1:2 * n) + beta * d
!do i = 1, n
!    d(i) = -w(n + i) + beta * d(i)
!end do
alpbd = ZERO
goto 150
!
!     Return from the subroutine.
!
320 snorm = ZERO
if (reduct > ZERO) snorm = sqrt(ss)

if (DEBUGGING) then
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
end if

end subroutine trstep


end module trustregion_mod
