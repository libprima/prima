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
! Last Modified: Sunday, April 17, 2022 AM12:21:36
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
&        dd, dg, dgd, ds, bstep, reduct, resmax, rhs, scaling, snsq, ss, temp, wgd
integer(IK) :: icount, j, jsav
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
!       set to zero if the constraints prohibit a long enough step.
!     STEP is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current STEP.
!     RESNEW, RESACT, D, DW and W are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       zero value of RESNEW(J) indicates that the J-th constraint is active,
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
if (.not. dd > 0) then
    goto 320
end if
scaling = 0.2_RP * snorm / sqrt(dd)
dw = scaling * dw

! If the modulus of the residual of an active constraint is substantial, then set D to the shortest
! move from STEP to the boundaries of the active constraints.
resmax = ZERO
resmax = maxval([ZERO, resact(1:nact)])

bstep = ZERO
if (resmax > 1.0D-4 * snorm) then
    vlam(1:nact) = solve(transpose(rfac(1:nact, 1:nact)), resact(1:nact))
    d = matprod(qfac(:, 1:nact), vlam(1:nact))

    ! The vector D that has just been calculated is also the shortest move from STEP+DW to the
    ! boundaries of the active constraints. Set BSTEP to the greatest steplength of this move that
    ! satisfies the trust region bound.
    ds = inprod(d, step + dw)
    dd = sum(d**2)
    rhs = snsq - sum((step + dw)**2)  ! Calculation changed

    if (rhs > 0) then
        temp = sqrt(ds * ds + dd * rhs)
        if (ds <= 0) then
            ! Zaikun 20210925
            ! What if we are at the first iteration? BLEN = DELTA/|D|? See TRSAPP.F90 of NEWUOA.
            bstep = (temp - ds) / dd
        else
            bstep = rhs / (temp + ds)
        end if
    end if

    ! Reduce BSTEP if necessary so that the move along D also satisfies the linear constraints.
    do j = 1, m
        if (.not. bstep > 0) exit
        if (resnew(j) > 0) then
            ad = inprod(d, amat(:, j))
            adw = inprod(dw, amat(:, j))
            if (ad > 0) then
                temp = max((resnew(j) - adw) / ad, ZERO)
                bstep = min(bstep, temp)
            end if
        end if
    end do

    bstep = min(bstep, ONE)
end if

! Set the next direction for seeking a reduction in the model function subject to the trust region
! bound and the linear constraints.
if (bstep <= 0) then
    d = dw
    icount = nact
else
    d = dw + bstep * d
    icount = nact - 1
end if
alpbd = ONE

! Set ALPHA to the steplength from STEP along D to the trust region boundary. Return if the first
! derivative term of this step is sufficiently small or if no further progress is possible.
150 icount = icount + 1
rhs = snsq - ss
if (rhs <= 0) goto 320
dg = inprod(d, g)
ds = inprod(d, step)
dd = inprod(d, d)
if (dg >= 0) goto 320
temp = sqrt(rhs * dd + ds * ds)
if (ds <= 0) then
    alpha = (temp - ds) / dd
else
    alpha = rhs / (temp + ds)
end if
if (-alpha * dg <= ctest * reduct) goto 320

dw = hess_mul(d, xpt, pq, hq)

! Set DGD to the curvature of the model along D. Then reduce ALPHA if necessary to the value that
! minimizes the model.
dgd = inprod(d, dw)
alpht = alpha
if (dg + alpha * dgd > 0) then
    alpha = -dg / dgd
end if

! Make a further reduction in ALPHA if necessary to preserve feasibility, and put some scalar
! products of D with constraint gradients in W.
alphm = alpha
jsav = 0

do j = 1, m
    ad = ZERO
    if (resnew(j) > 0) then
        ad = inprod(d, amat(:, j))
        if (ad > 0) then
            if (alpha > resnew(j) / ad) then
                alpha = resnew(j) / ad
                jsav = j
            end if
        end if
    end if
    w(j) = ad
end do


alpha = max(alpha, alpbd)
alpha = min(alpha, alphm)
if (icount == nact) alpha = min(alpha, ONE)

! Update STEP, G, RESNEW, RESACT and REDUCT.
step = step + alpha * d
ss = sum(step**2)
g = g + alpha * dw

where (resnew > 0)
    resnew = max(resnew - alpha * w(1:m), TINYCV)
end where
!!MATLAB: mask = (resnew > 0); resnew(mask) = max(resnew(mask) - alpha * w(mask), TINYCV);

if (icount == nact .and. nact > 0) then
    resact(1:nact) = (ONE - bstep) * resact(1:nact)
end if
reduct = reduct - alpha * (dg + HALF * alpha * dgd)

! Test for termination. Branch to label 40 if there is a new active constraint and if the distance
! from STEP to the trust region boundary is at least 0.2*SNORM.

! Zaikun 2019-08-29: the code can encounter infinite cycling due to NaN values. Exit when NGETACT is
! large or NaN detected.
! Caution: 1. MIN accepts only data with the same KIND; 2. Integer overflow.
if (ngetact > min(10000, 100 * int(m + 1) * int(n)) .or.  &
& alpha /= alpha .or. alpht /= alpht .or. &
& alphm /= alphm .or. dgd /= dgd .or. dg /= dg .or. &
& ss /= ss .or. snsq /= snsq .or. reduct /= reduct) then
    goto 320
end if
if (alpha == alpht) goto 320
temp = -alphm * (dg + HALF * alphm * dgd)
if (temp <= ctest * reduct) goto 320
if (jsav > 0) then
    if (ss <= 0.64_RP * snsq) goto 40
    goto 320
end if
if (icount == n) goto 320

! Calculate the next search direction, which is conjugate to the previous one unless ICOUNT == NACT.
if (nact > 0) then
    w(n + 1:2 * n) = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
else
    w(n + 1:2 * n) = g
end if
if (icount == nact) then
    beta = ZERO
else
    wgd = inprod(w(n + 1:2 * n), dw)
    beta = wgd / dgd
end if
d = -w(n + 1:2 * n) + beta * d
alpbd = ZERO
goto 150

! Return from the subroutine.
320 snorm = ZERO
if (reduct > 0) snorm = sqrt(ss)

if (DEBUGGING) then
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
end if

end subroutine trstep


end module trustregion_mod
