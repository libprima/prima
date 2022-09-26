module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of LINCOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Sunday, August 21, 2022 AM05:33:20
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstep


contains


subroutine trstep(amat, delta, gq, hq, pq, rescon, xpt, iact, nact, qfac, rfac, ngetact, s)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, TWO, HALF, EPS, TINYCV, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : matprod, inprod, norm, solve, isorth, istriu, &
    & issymmetric, trueloc
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
real(RP), intent(out) :: s(:)  ! S(N)

! Local variables
logical :: get_act
character(len=*), parameter :: srname = 'TRSTEP'
real(RP) :: d(size(gq))
real(RP) :: dw(size(gq))
real(RP) :: gw(size(gq))
real(RP) :: resact(size(amat, 2))
real(RP) :: resnew(size(amat, 2))
real(RP) :: tol
real(RP) :: g(size(gq))
real(RP) :: frac(size(amat, 2)), ad(size(amat, 2)), restmp(size(amat, 2)), alpbd, alpha, alphm, alpht, &
& beta, ctest, &
&        dd, dg, dgd, ds, bstep, reduct, delres, scaling, delsq, ss, temp, sold(size(s))
integer(IK) :: maxiter, iter, icount, jsav  ! What is ICOUNT counting? Better name for ICOUNT?
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
    call wassert(all(is_finite(xpt)), 'XPT is finite', srname)
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

!====================!
! Calculation starts !
!====================!

g = gq

! Return if G is not finite. Otherwise, GETACT will fail in the debugging mode.
if (.not. is_finite(sum(abs(g)))) then
    ngetact = 0
    s = ZERO
    return
end if

!
!     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
!       are the same as the terms with these names in LINCOB. If RESCON(J)
!       is negative, then |RESCON(J)| must be no less than the trust region
!       radius, so that the J-th constraint can be ignored.
!     S is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current S.
!     RESNEW, RESACT, D, DW and GW are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       zero value of RESNEW(J) indicates that the J-th constraint is active,
!       and otherwise RESNEW(J) is set to the greater of TINYCV and the actual
!       residual of the J-th constraint for the current S. RESACT holds
!       the residuals of the active constraints, which may be positive.
!       D is the search direction of each line search. DW is either another
!       search direction or the change in gradient along D.
!
!     Set some numbers for the conjugate gradient iterations.
!
ctest = 0.01_RP
delsq = delta * delta

! Set the initial elements of RESNEW, RESACT and S.
! 1. RESNEW(J) < 0 indicates that the J-th constraint does not restrict the CG steps of the current
! trust region calculation. In other words, RESCON >= DELTA.
! 2. RESNEW(J) = 0 indicates that J is an entry of IACT(1:NACT).
! 3. RESNEW(J) > 0 means that RESNEW(J) = max(B(J) - AMAT(:, J)^T*(XOPT+D), TINYCV).
! N.B.: The order of the following lines is important, as the later ones override the earlier.
resnew = rescon
resnew(trueloc(rescon >= 0)) = max(TINYCV, rescon(trueloc(rescon >= 0)))
resnew(trueloc(rescon >= delta)) = -ONE
!!MATLAB:
!!resnew = rescon; resnew(rescon >= 0) = max(TINYCV, rescon(rescon >= 0)); resnew(rescon >= delta) = -1;
resnew(iact(1:nact)) = ZERO
resact(1:nact) = rescon(iact(1:nact))

s = ZERO
ss = ZERO
alpbd = ONE  ! Artificial value. Not used.
reduct = ZERO
ngetact = 0
get_act = .true.  ! Better name?

icount = -1_IK  ! Artificial value. Not used. To entertain compilers. To be removed.

maxiter = min(1000_IK, 10_IK * (m + n))  ! What is the theoretical upper bound of ITER?
do iter = 1, maxiter  ! Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
    if (get_act) then
        ! GETACT picks the active set for the current S. It also sets DW to the vector closest to -G
        ! that is orthogonal to the normals of the active constraints. DW is scaled to have length
        ! 0.2*DELTA. Then a move of DW from S is allowed by the linear constraints: DW reduces
        ! the values of the nearly active constraints; it changes the inactive constraints by at
        ! most 0.2*DELTA, but the residuals of the inactive constraints at no less than 0.2*DELTA.
        ngetact = ngetact + 1
        call getact(amat, delta, g, iact, nact, qfac, resact, resnew, rfac, dw)
        dd = inprod(dw, dw)
        if (dd <= 0 .or. is_nan(dd)) then
            exit
        end if
        scaling = 0.2_RP * delta / sqrt(dd)
        dw = scaling * dw

        bstep = ZERO
        if (any(resact(1:nact) > 1.0E-4_RP * delta)) then
            ! If the modulus of the residual of an active constraint is substantial (namely, is more
            ! than 1.0E-4*DELTA), then set D to the shortest move from S to the boundaries of the
            ! active constraints.
            ! N.B.: We prefer `ANY(X > Y)` to `MAXVAL(X) > Y`, as Fortran standards do not specify
            ! MAXVAL(X) when X contains NaN, and MATLAB/Python/R/Julia behave differently in this
            ! respect. Moreover, MATLAB defines max(X) = [] if X == [], differing from mathematics
            ! and other languages.
            d = matprod(qfac(:, 1:nact), solve(transpose(rfac(1:nact, 1:nact)), resact(1:nact)))
            !!MATLAB: d = qfac(:, 1:nact) * (rfac(1:nact, 1:nact)' \ resact(1:nact))

            ! The vector D that has just been calculated is also the shortest move from S+DW to the
            ! boundaries of the active constraints. Set BSTEP to the greatest steplength of this
            ! move that satisfies both the trust region bound and the linear constraints.
            ds = inprod(d, s + dw)
            dd = sum(d**2)
            delres = delsq - sum((s + dw)**2)  ! Calculation changed

            if (delres > 0) then
                ! Set BSTEP to the greatest value so that S + DW + BSTEP*D satisfies the trust
                ! region bound.
                temp = sqrt(ds * ds + dd * delres)
                if (ds <= 0) then
                    bstep = (temp - ds) / dd
                else
                    bstep = delres / (temp + ds)
                end if
                ! Reduce BSTEP so that the move along D also satisfies the linear constraints.
                ad = -ONE
                ad(trueloc(resnew > 0)) = matprod(d, amat(:, trueloc(resnew > 0)))
                frac = ONE
                restmp(trueloc(ad > 0)) = resnew(trueloc(ad > 0)) - matprod(dw, amat(:, trueloc(ad > 0)))
                frac(trueloc(ad > 0)) = restmp(trueloc(ad > 0)) / ad(trueloc(ad > 0))
                bstep = minval([ONE, bstep, frac])  ! BSTEP = MINVAL([ONES, BSTEP, FRAC(TRUELOC(AD>0))])
            end if
        end if

        ! Set the next direction for seeking a reduction in the model function subject to the trust
        ! region bound and the linear constraints.
        if (bstep > 0) then
            d = dw + bstep * d
            icount = nact - 1
        else
            d = dw
            icount = nact
        end if
        alpbd = ONE
    end if

    ! Set ALPHA to the steplength from S along D to the trust region boundary. Return if the first
    ! derivative term of this step is sufficiently small or if no further progress is possible.
    icount = icount + 1
    delres = delsq - ss
    if (delres <= 0) then
        exit
    end if
    dg = inprod(d, g)
    ds = inprod(d, s)
    dd = inprod(d, d)
    if (dg >= 0 .or. is_nan(dg)) then
        exit
    end if
    temp = sqrt(ds * ds + dd * delres)
    if (ds <= 0) then
        alpha = (temp - ds) / dd
    else
        alpha = delres / (temp + ds)
    end if
    if (-alpha * dg <= ctest * reduct) then
        exit
    end if

    dw = hess_mul(d, xpt, pq, hq)

    ! Set DGD to the curvature of the model along D. Then reduce ALPHA if necessary to the value
    ! that minimizes the model.
    dgd = inprod(d, dw)
    alpht = alpha
    if (dg + alpha * dgd > 0) then
        alpha = -dg / dgd
    end if

    ! Make a further reduction in ALPHA if necessary to preserve feasibility.
    alphm = alpha
    ad = -ONE
    ad(trueloc(resnew > 0)) = matprod(d, amat(:, trueloc(resnew > 0)))
    frac = alpha
    frac(trueloc(ad > 0)) = resnew(trueloc(ad > 0)) / ad(trueloc(ad > 0))
    frac(trueloc(is_nan(frac))) = alpha
    jsav = 0_IK
    if (any(frac < alpha)) then
        jsav = int(minloc(frac, dim=1), IK)
        alpha = frac(jsav)
    end if
    !----------------------------------------------------------------------------------------------!
    ! Alternatively, JSAV and ALPHA can be calculated as below.
    !JSAV = INT(MINLOC([ALPHA, FRAC], DIM=1), IK) - 1_IK
    !ALPHA = MINVAL([ALPHA, FRAC])  ! This line cannot be exchanged with the last.
    ! We prefer our implementation as the code is more explicit; in addition, it is more flexible:
    ! we can change the condition ANY(FRAC < ALPHA) to ANY(FRAC < (1 - EPS) * ALPHA) or
    ! ANY(FRAC < (1 + EPS) * ALPHA), depending on whether we believe a false positive or a false
    ! negative of JSAV > 0 is more harmful.
    !----------------------------------------------------------------------------------------------!

    alpha = min(max(alpha, alpbd), alphm)
    if (icount == nact) then
        alpha = min(alpha, ONE)
    end if

    ! Update S, G, RESNEW, RESACT and REDUCT.
    sold = s
    s = s + alpha * d
    ss = sum(s**2)
    if (.not. is_finite(ss)) then
        s = sold
        exit
    end if
    g = g + alpha * dw
    if (.not. is_finite(sum(abs(g)))) then
        exit
    end if
    restmp = resnew - alpha * ad  ! Only RESTMP(TRUELOC(RESNEW > 0)) is needed.
    resnew(trueloc(resnew > 0)) = max(TINYCV, restmp(trueloc(resnew > 0)))
    !!MATLAB: mask = (resnew > 0); resnew(mask) = max(TINYCV, resnew(mask) - alpha * ad(mask));
    if (icount == nact) then
        !resact(1:nact) = (ONE - bstep) * resact(1:nact)
        resact(1:nact) = (ONE - alpha*bstep) * resact(1:nact)
    end if
    reduct = reduct - alpha * (dg + HALF * alpha * dgd)
    if (reduct <= 0 .or. is_nan(reduct)) then
        s = sold
        exit
    end if

    ! Test for termination. Branch to a new loop if there is a new active constraint and if the
    ! distance from S to the trust region boundary is at least 0.2*DELTA.
    if (alpha >= alpht .or. -alphm * (dg + HALF * alphm * dgd) <= ctest * reduct) then
        exit
    end if
    if (jsav > 0) then
        !if (ss <= 0.64_RP * delsq) then
        get_act = .true.
        cycle
        !else
        !    exit
        !end if
    end if
    if (icount == n) then
        exit
    end if

    ! Calculate the next search direction, which is conjugate to the previous one if ICOUNT /= NACT.
    ! N.B.: NACT < 0 is impossible unless GETACT is buggy; NACT = 0 can happen, particularly if
    ! there is no constraint. In theory, the code for the second case below covers the first as well.
    if (nact <= 0) then
        gw = g
    else
        gw = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
        !!MATLAB: gw = qfac(:, nact+1:n) * (g' * qfac(:, nact+1:n))';
    end if

    if (icount == nact) then
        beta = ZERO
    else
        beta = inprod(gw, dw) / dgd
    end if
    d = -gw + beta * d
    alpbd = ZERO
    get_act = .false.
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(s) == n .and. all(is_finite(s)), 'SIZE(S) == N, S is finite', srname)
    call assert(norm(s) <= TWO * delta, '|S| <= 2*DELTA', srname)
    ! Due to rounding, it may happen that |S| > DELTA, but |S| > 2*DELTA is highly improbable.
    call assert(nact >= 0 .and. nact <= min(m, n), '0 <= NACT <= MIN(M, N)', srname)
    call assert(size(qfac, 1) == n .and. size(qfac, 2) == n, 'SIZE(QFAC) == [N, N]', srname)
    call assert(isorth(qfac, tol), 'QFAC is orthogonal', srname)
    call assert(size(rfac, 1) == n .and. size(rfac, 2) == n, 'SIZE(RFAC) == [N, N]', srname)
    call assert(istriu(rfac), 'RFAC is upper triangular', srname)
end if

end subroutine trstep


end module trustregion_mod

!--------------------------------------------------------------------------------------------------!
! Zaikun 20220417:
! For GW, the schemes below work evidently worse than the one above in a test on 20220417. Why?
!-----------------------------------------------------------------------!
! VERSION 1:
!!gw = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
!-----------------------------------------------------------------------!
! VERSION 2:
!!if (2 * nact < n) then
!!    gw = g - matprod(qfac(:, 1:nact), matprod(g, qfac(:, 1:nact)))
!!else
!!    gw = matprod(qfac(:, nact + 1:n), matprod(g, qfac(:, nact + 1:n)))
!!end if
!-----------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
