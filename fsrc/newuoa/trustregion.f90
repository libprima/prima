module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region iterations.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Wednesday, January 05, 2022 PM11:16:05
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trsapp, trrad

contains


subroutine trsapp(delta, gq, hq, pq, tol, x, xpt, crvmin, s, info)
!--------------------------------------------------------------------------------------------------!
! TRSAPP finds an approximate solution to the N-dimensional trust region subproblem
!
! min <X+S, GQ> + 0.5*<X+S, HESSIAN*(X+S)> s.t. ||S|| <= DELTA
!
! Note that the HESSIAN here is the sum of an explicit part HQ and an implicit part (PQ, XPT):
!
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*XPT(:, K)*XPT(:, K)' .
!
! The calculation of S begins with the truncated conjugate gradient method. If the boundary of the
! trust region is reached, then further changes to S may be made, each one being in the 2D space
! spanned by the current S and the corresponding gradient of Q. Thus S should provide a substantial
! reduction to Q within the trust region. See Section 5 of the NEWUOA paper.
!
! At return, S will be the approximate solution. CRVMIN will be set to the least curvature of
! HESSIAN along the conjugate directions that occur, except that it is set to ZERO if S goes all the
! way to the trust-region boundary. INFO is an exit flag:
! INFO = 0: an approximate solution satisfying one of the
! following conditions is found:
! 1. ||G+HS||/||GBEG|| <= TOL,
! 2. ||S|| = DELTA and <S, -(G+HS)> >= (1 - TOL)*||S||*||G+HS||,
! where TOL is a tolerance that is set to 1e-2 in NEWUOA.
! INFO = 1: the iteration is reducing Q only slightly;
! INFO = 2: the maximal number of iterations is attained;
! INFO = -1: too much rounding error to continue

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, ZERO, TENTH, DEBUGGING
use, non_intrinsic :: circle_mod, only : circle_min
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : inprod, issymmetric, norm, project, hess_mul

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: gq(:)   ! GQ(N)
real(RP), intent(in) :: hq(:, :)    ! HQ(N, N)
real(RP), intent(in) :: pq(:)   ! PQ(NPT)
real(RP), intent(in) :: tol
real(RP), intent(in) :: x(:)    ! X(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)

! Outputs
integer(IK), intent(out) :: info
real(RP), intent(out) :: crvmin
real(RP), intent(out) :: s(:)   ! S(N)

! Local variables
character(len=*), parameter :: srname = 'TRSAPP'
integer(IK) :: iter
integer(IK) :: itermax
integer(IK) :: n
integer(IK) :: npt
logical :: twod_search
real(RP) :: alpha
real(RP) :: angle
real(RP) :: args(4)
real(RP) :: bstep
real(RP) :: cth
real(RP) :: d(size(x))
real(RP) :: dd
real(RP) :: delsq
real(RP) :: dg
real(RP) :: dhd
real(RP) :: dhs
real(RP) :: ds
real(RP) :: g(size(x))
real(RP) :: gg
real(RP) :: gg0
real(RP) :: ggsav
real(RP) :: hd(size(x))
real(RP) :: hs(size(x))
real(RP) :: hx(size(x))
real(RP) :: hypt
real(RP) :: qadd
real(RP) :: qred
real(RP) :: reduc
real(RP) :: sg
real(RP) :: shs
real(RP) :: sold(size(x))
real(RP) :: ss
real(RP) :: sth

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(size(gq) == n, 'SIZE(GQ) = N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is an NxN symmetric matrix', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) = NPT', srname)
    call assert(size(x) == n .and. all(is_finite(x)), 'SIZE(X) == N, X is finite', srname)
    call assert(size(s) == n, 'SIZE(S) == N', srname)
end if

!====================!
! Calculation starts !
!====================!

s = ZERO
crvmin = ZERO
qred = ZERO
info = 2  ! Default exit flag is 2, i.e., itermax is attained

! Prepare for the first line search.
hx = hess_mul(hq, pq, xpt, x)
!--------------------------------------------------------------------------------------------------!
! N.B.: During the iterations, G is NOT updated, and it equals always GQ+HX, which is the gradient
! of the trust-region model at the trust-region center X. However, GG is updated: GG = |G + HS|^2,
! which is the norm square of the gradient at the current iterate.
g = gq + hx
gg = inprod(g, g)
!--------------------------------------------------------------------------------------------------!
gg0 = gg
d = -g
dd = gg
ds = ZERO
ss = ZERO
hs = ZERO
delsq = delta * delta
itermax = n

twod_search = .false.

! The truncated-CG iterations.
!
! The iteration will be terminated in 4 possible cases:
! 1. the maximal number of iterations is attained;
! 2. QADD <= TOL*QRED or ||G|| <= TOL*||GBEG||, where QADD is the reduction of Q due to the latest
! CG step, QRED is the reduction of Q since the beginning until the latest CG step, G is the
! current gradient, and GBEG is the initial gradient; see (5.13) of the NEWUOA paper;
! 3. DS <= 0
! 4. ||S|| = DELTA, i.e., CG path cuts the trust region boundary.
!
! In the 4th case, twod_search will be set to true, meaning that S will be improved by a sequence of
! two-dimensional search, the two-dimensional subspace at each iteration being span(S, -G).
do iter = 1, itermax
    ! Exit if GG is small. This must be done first; otherwise, DD can be 0 and BSTEP is not well
    ! defined. The inequality below must be non-strict so that GG = 0 will trigger the exit.
    if (gg <= (tol**2) * gg0) then
        info = 0_IK
        exit
    end if

    ! Set BSTEP to the step length such that |S + BSTEP*D| = DELTA.
    if (iter == 1) then
        bstep = delta / sqrt(dd)
    else
        hypt = sqrt(ds**2 + dd * (delsq - ss))
        ! Powell's code does not distinguish the following two cases, which have no difference in
        ! precise arithmetic. The following scheme stabilizes the calculation. Copied from LINCOA.
        if (ds <= 0) then
            bstep = (hypt - ds) / dd
        else
            bstep = (delsq - ss) / (ds + hypt)
        end if
    end if
    hd = hess_mul(hq, pq, xpt, d)
    dhd = inprod(d, hd)

    ! Set the step-length ALPHA and update CRVMIN.
    if (dhd <= 0) then
        alpha = bstep
    else
        alpha = min(bstep, gg / dhd)
        if (iter == 1) then
            crvmin = dhd / dd
        else
            crvmin = min(crvmin, dhd / dd)
        end if
    end if
    ! QADD is the reduction of Q due to the new CG step.
    qadd = alpha * (gg - HALF * alpha * dhd)
    ! QRED is the reduction of Q up to now.
    qred = qred + qadd
    ! QADD and QRED will be used in the 2D minimization if any.

    ! Update S, HS, and GG.
    sold = s
    s = s + alpha * d
    ss = inprod(s, s)
    hs = hs + alpha * hd
    ggsav = gg  ! Gradient norm square before this iteration
    gg = inprod(g + hs, g + hs)  ! Current gradient norm square
    ! We may record g+hs for later usage:
    ! gnew = g + hs
    ! Note that we should NOT set g = g + hs, because g contains the gradient of Q at x.

    ! Check whether to exit. This should be done after updating HS and GG, which will be used for
    ! the 2D minimization if any.
    ! Exit in case of Inf/NaN in S. This should come the first! Otherwise, we may return an S that
    ! contains NaN and fulfills other exit conditions.
    if (.not. is_finite(sum(abs(s)))) then
        s = sold
        info = -1_IK
        exit
    end if

    ! Exit if CG path cuts the boundary. It is the only possibility that TWOD_SEARCH is true.
    if (alpha >= bstep .or. ss >= delsq) then
        crvmin = ZERO
        !twod_search = (n >= 2)  ! TWOD_SEARCH should be FALSE if N = 1.
        twod_search = (n >= 2 .and. gg > (tol**2) * gg0)  ! TWOD_SEARCH should be FALSE if N = 1.
        exit
    end if

    ! Exit due to small QADD.
    if (qadd <= tol * qred) then
        info = 1_IK
        exit
    end if

    ! Prepare for the next CG iteration.
    d = (gg / ggsav) * d - g - hs  ! CG direction
    dd = inprod(d, d)
    ds = inprod(d, s)
    if (ds <= 0) then
        ! DS is positive in theory.
        info = -1_IK
        exit
    end if
end do

if (ss <= 0 .or. is_nan(ss)) then
    ! This may occur for ill-conditioned problems due to rounding.
    info = -1_IK
    twod_search = .false.
end if

if (twod_search) then
    ! At least 1 iteration of 2D minimization
    itermax = max(int(1, kind(itermax)), itermax - iter)
else
    itermax = 0_IK
end if

! The 2D minimization
! N.B.: During the iterations, G is NOT updated, and it equals always GQ+HX, which is the gradient
! of the trust-region model at the trust-region center X. However, GG is updated: GG = |G + HS|^2,
! which is the norm square of the gradient at the current iterate.
do iter = 1, itermax
    ! Exit if GG is small. The inequality must be non-strict so that GG = 0 can trigger the exit.
    if (gg <= (tol**2) * gg0) then
        info = 0_IK
        exit
    end if
    sg = inprod(s, g)
    shs = inprod(s, hs)

    !----------------------------------------------------------------------------------------------!
    !!----------------------------------------------!
    !!sgk = sg + shs
    !!if (sgk / sqrt(gg * delsq) <= tol - ONE) then
    !!    info = 0_IK
    !!    exit
    !!end if
    !!t = sqrt(delsq * gg - sgk**2)
    !!d = (delsq / t) * (g + hs) - (sgk / t) * s
    !!----------------------------------------------!
    ! The above is Powell's code. In precise arithmetic, INPROD(D, S) = 0 and |D| = |S|. However,
    ! when DELSQ*GG - SGK**2 is tiny, the error in D can be large and hence damage these
    ! inequalities significantly. This did happen in tests, especially when using the single
    ! precision. We calculate D as below. It did improve the performance of NEWUOA in our test.
    !----------------------------------------------------------------------------------------------!

    ! Begin the 2D minimization by calculating D and HD and some scalar products.
    !!sunit = s / norm(s)  ! Normalization seems to stabilize the projection below.
    !!d = (g + hs) - inprod(g + hs, sunit) * sunit
    d = (g + hs) - project(g + hs, s)
    ! N.B.:
    ! 1. The condition |D|<=SQRT(TOL*GG) below is equivalent to |INPROD(G+HS,S)|<=SQRT((1-TOL)*GG*SS)
    ! in theory. As given above, Powell's code triggers an exit if INPROD(G+HS,S)=SGK<=(TOL-1)*GG*SS.
    ! Since |SQRT(1-TOL) - (1-TOL)| <= TOL/2, our condition is close to |SGK| <= (TOL-1)*GG*SS.
    ! When |SGK| is tiny, S and G+HS are nearly parallel and hence the 2D search cannot continue.
    ! Note that SGK is unlikely positive if everything goes well.
    ! 2. SQRT(TOL)*SQRT(GG) is less likely to encounter underflow than SQRT(TOL*GG).
    ! 3. The condition below should be non-strict so that |D| = 0 can trigger the exit.
    if (norm(d) <= sqrt(tol) * sqrt(gg)) then
        info = 0_IK
        exit
    end if
    d = (d / norm(d)) * norm(s)
    ! In precise arithmetic, INPROD(D, S) = 0 and |D| = |S| = DELTA.
    if (abs(inprod(d, s)) >= TENTH * norm(d) * norm(s) .or. norm(d) >= TWO * delta) then
        info = -1_IK
        exit
    end if

    hd = hess_mul(hq, pq, xpt, d)

    ! Seek the value of the angle that minimizes Q.
    ! First, calculate the coefficients of the objective function on the circle.
    dg = inprod(d, g)
    dhd = inprod(hd, d)
    dhs = inprod(hd, s)
    args = [sg, HALF * (shs - dhd), dg, dhs]
    ! The 50 in the line below was chosen by Powell. It works the best in tests, magically. Larger
    ! (e.g., 60, 100) or smaller (e.g., 20, 40) values will worsen the performance of NEWUOA. Why??
    angle = circle_min(circle_fun_trsapp, args, 50_IK)

    ! Calculate the new S.
    cth = cos(angle)
    sth = sin(angle)
    sold = s
    s = cth * s + sth * d

    ! Exit in case of Inf/NaN in S.
    if (.not. is_finite(sum(abs(s)))) then
        s = sold
        info = -1_IK
        exit
    end if

    ! Test for convergence.
    reduc = circle_fun_trsapp(ZERO, args) - circle_fun_trsapp(angle, args)
    qred = qred + reduc
    if (reduc / qred <= tol) then
        info = 1_IK
        exit
    end if

    ! Calculate HS.
    hs = cth * hs + sth * hd
    gg = inprod(g + hs, g + hs)
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(s) == n .and. all(is_finite(s)), 'SIZE(S) == N, S is finite', srname)
    call assert(norm(s) <= TWO * delta, '|S| <= 2*DELTA', srname)
    ! Due to rounding, it may happen that |S| > DELTA, but |S| > 2*DELTA is highly improbable.
end if

end subroutine trsapp


function circle_fun_trsapp(theta, args) result(f)
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none
! Inputs
real(RP), intent(in) :: theta
real(RP), intent(in) :: args(:)

! Outputs
real(RP) :: f

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_FUN_TRSAPP'
real(RP) :: cth
real(RP) :: sth

! Preconditions
if (DEBUGGING) then
    call assert(size(args) == 4, 'SIZE(ARGS) == 4', srname)
end if

!====================!
! Calculation starts !
!====================!

cth = cos(theta)
sth = sin(theta)
f = (args(1) + args(2) * cth) * cth + (args(3) + args(4) * cth) * sth

!====================!
!  Calculation ends  !
!====================!
end function circle_fun_trsapp


function trrad(delta0, dnorm, eta1, eta2, gamma1, gamma2, ratio) result(delta)
!--------------------------------------------------------------------------------------------------!
! This function updates the trust region radius according to RATIO and DNORM.
!--------------------------------------------------------------------------------------------------!

! Generic module
use, non_intrinsic :: consts_mod, only : RP, HALF, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Input
real(RP), intent(in) :: delta0   ! Current trust-region radius
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
    call assert(delta0 >= dnorm .and. dnorm > 0, 'DELTA0 >= DNORM > 0', srname)
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
    delta = max(HALF * delta0, dnorm)
else
    delta = max(HALF * delta0, gamma2 * dnorm)
end if

! For noisy problems, the following may work better.
!!if (ratio <= eta1) then
!!    delta = gamma1 * dnorm
!!elseif (ratio <= eta2) then  ! Ensure DELTA >= DELTA0
!!    delta = delta0
!!else  ! Ensure DELTA > DELTA0 with a constant factor
!!    delta = max(delta0 * (1.0_RP + gamma2) / 2.0_RP, gamma2 * dnorm)
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
