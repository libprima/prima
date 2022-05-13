module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines concerning the trust-region calculations of UOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Friday, May 13, 2022 AM02:15:22
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trstep


contains


subroutine trstep(delta, g, h, tol, d, crvmin)
!--------------------------------------------------------------------------------------------------!
! This subroutine solves the trust-region subproblem
!
!     minimize <G, D> + 0.5 * <D, H*D> subject to ||D|| <= DELTA.
!
! The algorithm first tridiagonalizes H and then applies the More-Sorensen method in
! More, J. J., and Danny C. S., "Computing a trust region step", SIAM J. Sci. Stat. Comput. 4:
! 553-572, 1983.
! See Section 2 of the UOBYQA paper and
! Powell, M. J. D., "Trust region calculations revisited", Numerical Analysis 1997: Proceedings of
! the 17th Dundee Biennial Numerical Analysis Conference, 1997, 193--211,
! Powell, M. J. D., "The use of band matrices for second derivative approximations in trust region
! algorithms", Advances in Nonlinear Programming: Proceedings of the 96 International Conference on
! Nonlinear Programming, 1998, 3--28.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : issymmetric, inprod, hessenberg, eigmin, trueloc

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: g(:)  ! G(N)
real(RP), intent(in) :: h(:, :)  ! H(N, N)
real(RP), intent(in) :: tol

! In-outputs
real(RP), intent(out) :: d(:)  ! D(N)
real(RP), intent(out) :: crvmin

! Local variables
character(len=*), parameter :: srname = 'TRSTEP'
integer(IK) :: n
real(RP) :: gg(size(g))
real(RP) :: hh(size(g), size(g))
real(RP) :: piv(size(g))
real(RP) :: td(size(g))
real(RP) :: tn(size(g) - 1)
!real(RP) :: w(size(g))
real(RP) :: z(size(g))
real(RP) :: dold(size(g)) !!!
real(RP) :: dnewton(size(g))  ! Newton-Raphson step; only calculated when N = 1.
real(RP) :: delsq, dhd, dnorm, dsq, dtg, dtz, gam, gnorm,     &
&        gsq, hnorm, par, parl, parlest, paru,         &
&        paruest, phi, phil, phiu, &
&        slope, partmp, &
&        tnz, tempa, tempb, wsq, wwsq, zsq
integer(IK) :: iter, k, ksav, maxiter
logical :: posdef, negcrv, d_initialized

!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q. Only the upper triangular and diagonal
!       parts need be set. The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.
!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
!     CRVMIN will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step. Then CRVMIN will be positive, but otherwise it
!     will be set to ZERO.
!
!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| <= DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties
!
!             ||D|| <= DELTA  and  ACTRED <= (1-TOL)*MAXRED.

! Sizes.
n = int(size(g), kind(n))

if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(size(h, 1) == n .and. issymmetric(h), 'H is n-by-n and symmetric', srname)
    call assert(size(d) == n, 'SIZE(D) == N', srname)
end if

d = ZERO
crvmin = ZERO

gsq = sum(g**2)
gnorm = sqrt(gsq)

if (.not. any(abs(h) > 0)) then
    if (gnorm > 0) then
        d = -delta * (g / gnorm)
    end if
    return
end if

! Zaikun 20220301, 20220305:
! Powell's original code requires that N >= 2.  When N = 1, the code does not work (sometimes even
! encounter memory errors). This is indeed why the original UOBYQA code constantly terminates with
! "a trust region step has failed to reduce the quadratic model" when applied to univariate problems.
if (n == 1) then
    d = sign(delta, -g)  !!MATLAB: D_OUT = -DELTA * SIGN(G)
    if (h(1, 1) > 0) then
        dnewton = -g / h(1, 1)
        if (abs(dnewton(1)) <= delta) then
            d = dnewton
            crvmin = h(1, 1)
        end if
    end if
    return
end if

delsq = delta * delta

! Apply Householder transformations to get a tridiagonal matrix similar to H (i.e., the Hessenberg
! form of H), and put the elements of the Householder vectors in the lower triangular part of HH.
! Further, TD and TN will contain the diagonal and other nonzero elements of the tridiagonal matrix.
! In the comments hereafter, H indeed means this tridiagonal matrix.
hh = h
call hessenberg(hh, td, tn)  !!MATLAB: [P, h] = hess(h); td = diag(h); tn = diag(h, 1)

! Form GG by applying the similarity transformation to G.
gg = g
do k = 1, n - 1_IK
    gg(k + 1:n) = gg(k + 1:n) - inprod(gg(k + 1:n), hh(k + 1:n, k)) * hh(k + 1:n, k)
end do
!!MATLAB: gg = (g'*P)';  % gg = P'*g;

!--------------------------------------------------------------------------------------------------!
! Zaikun 20220303: Exit if GG, HH, TD, or TN is not finite. Otherwise, the behavior of this
! subroutine is not predictable. For example, if HNORM = GNORM = Inf, it is observed that the
! initial value of PARL defined below will change when we add code that should not affect PARL
! (e.g., print it, or add TD = 0, TN = 0, PIV = 0 at the beginning of this subroutine).
! This is probably because the behavior of MAX is undefined if it receives NaN (if GNORM and HNORM
! are both Inf, then GNORM/DELTA - HNORM = NaN).
!--------------------------------------------------------------------------------------------------!
if (.not. is_finite(sum(abs(gg)) + sum(abs(hh)) + sum(abs(td)) + sum(abs(tn)))) then
    return
end if

! Begin the trust region calculation with a tridiagonal matrix by calculating the L_1-norm of the
! Hessenberg form of H, which is an upper bound for the spectral norm of H.
hnorm = maxval(abs([ZERO, tn]) + abs(td) + abs([tn, ZERO]))

! Set the initial values of PAR and its bounds.
! N.B.: PAR is the parameter LAMBDA in More-Sorensen 1983 and Powell 1997, as well as the THETA in
! Section 2 of the UOBYQA paper. The algorithm looks for the optimal PAR characterized in Lemmas
! 2.1--2.3 of More-Sorensen 1983.
parl = maxval([ZERO, -minval(td), gnorm / delta - hnorm])  ! Lower bound for the optimal PAR
parlest = parl  ! Estimation for PARL
par = parl
paru = ZERO  ! Upper bound for the optimal PAR
paruest = ZERO  ! Estimation for PARU
posdef = .false.
iter = 0
maxiter = min(1000_IK, 100_IK * int(n, IK))  ! What is the theoretical bound of iterations?

140 continue

iter = iter + 1_IK

! Zaikun 26-06-2019: The original code can encounter infinite cycling, which did happen when testing
! the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS. Indeed, in all these cases, Inf and NaN
! appear in D due to extremely large values in the Hessian matrix (up to 10^219).
if (.not. is_finite(sum(abs(d)))) then
    d = dold
    goto 370
else
    dold = d
end if
if (iter > maxiter) then
    goto 370
end if

! Calculate the pivots of the Cholesky factorization of (H + PAR*I), i.e., the square of the diagonal
! of L in LL^T, or the diagonal of D in LDL^T). Note that the LDL factorization of H + PAR*I is
! L*diag(PIV)*L^T, where diag(PIV) is the diagonal matrix with PIV being the diagonal, and L is
! a lower triangular matrix with all the diagonal entries being 1, the subdiagonal being the vector
! TN/PIV(1:N-1), and all the other entries being 0.

piv = ZERO  ! PIV must be initialized, so that we know that any NaN in PIV is due to the loop below.
piv(1) = td(1) + par
! Powell implemented the loop by a GOTO, and K = N when the loop exits. It may not be true here.
do k = 1, n - 1_IK
    if (piv(k) > 0) then
        piv(k + 1) = td(k + 1) + par - tn(k)**2 / piv(k)
    elseif (abs(piv(k)) + abs(tn(k)) <= 0) then  ! PIV(K) == 0 == TN(K)
        piv(k + 1) = td(k + 1) + par
    else  ! PIV(K) < 0 .OR. (PIV(K) == 0 .AND. TN(K) /= 0)
        exit
    end if
end do

! Zaikun 20220509
if (any(is_nan(piv))) then
    goto 370  ! Better action to take???
end if

! NEGCRV is TRUE iff H + PAR*I has at least one negative eigenvalue (CRV means curvature).
negcrv = any(piv < 0 .or. (piv <= 0 .and. abs([tn, 0.0_RP]) > 0))

if (negcrv) then
    ! Set K to the first index corresponding to a negative curvature.
    k = minval(trueloc(piv < 0 .or. (piv <= 0 .and. abs([tn, 0.0_RP]) > 0)))
else
    ! Set K to the last index corresponding to a zero curvature; K = 0 if no such curvature exits.
    k = maxval([0_IK, trueloc(abs(piv) + abs([tn, 0.0_RP]) <= 0)])
end if

! At this point, K == 0 iff H + PAR*I is positive definite.

! Handle the case where H + PAR*I is positive semidefinite.
if (.not. negcrv) then
    ! Handle the case where the gradient at the trust region center is zero.
    if (gsq <= 0) then
        paru = par
        paruest = par
        if (par == ZERO) then  ! A rare case: the trust region center is optimal.
            goto 370
        end if
        if (k == 0) then  ! H + PAR*I is positive definite.
            goto 190
        end if
    end if

    ! Handle the case where the gradient at the trust region center is nonzero and H + PAR*I is
    ! positive definite.
    if (gsq > 0 .and. k == 0) then
        goto 230
    end if
end if


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! We arrive here only if 1 <= K <= N, when H + PAR*I has at least one nonpositive eigenvalue.
call assert(k >= 1 .and. k <= n, '1 <= K <= N', srname)
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


! Handle the case where H + PAR*I has at least one nonpositive eigenvalue.
! Set D to a direction of nonpositive curvature of the tridiagonal matrix, and thus revise PARLEST.
d = ZERO  ! Zaikun 20220512: Powell's code does not include this initialization. Is it correct???
d(k) = ONE  ! Zaikun 20220512: D(K+1:N) = ?

!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
! The following Boolean variable serves to check that D is always initialized in Powell's code
! whenever D is used.
d_initialized = (k == n)  ! Zaikun 20220512, TO BE REMOVED
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!

dhd = piv(k)

! In Fortran, the following two IFs CANNOT be merged into IF(K < N .AND. ABS(TN(K)) > ABS(PIV(K))).
! This is because Fortran may not perform a short-circuit evaluation of this logic expression, and
! hence TN(K) may be accessed even if K >= N, leading to an out-of-boundary index since SIZE(TN) is
! only N-1. This is not a problem in C, MATLAB, Python, Julia, or R, where short circuit is ensured.
if (k < n) then
    if (abs(tn(k)) > abs(piv(k))) then
        ! PIV(K+1) was named as "TEMP" in Powell's code. Is PIV(K+1) consistent with the meaning of PIV?
        piv(k + 1) = td(k + 1) + par
        if (piv(k + 1) <= abs(piv(k))) then
            d(k + 1) = sign(ONE, -tn(k))  !!MATLAB: d(k + 1) = -sing(tn(k))
            dhd = piv(k) + piv(k + 1) - TWO * abs(tn(k))
        else
            d(k + 1) = -tn(k) / piv(k + 1)
            dhd = piv(k) + tn(k) * d(k + 1)
        end if
    end if
end if

ksav = k

do k = ksav - 1_IK, 1, -1
    ! It may happen that TN(K) == 0 == PIV(K). Without checking TN(K), we will get D(K) = NaN.
    if (tn(k) /= ZERO) then
        d(k) = -tn(k) * d(k + 1) / piv(k)
    else
        d(1:k) = ZERO
        exit
    end if
end do

dsq = sum(d**2)
parl = par
parlest = par - dhd / dsq

190 continue

! Terminate with D set to a multiple of the current D if the following test suggests so.
if (gsq <= 0) then
    partmp = paruest * (ONE - tol)
else
    partmp = paruest
end if
if (paruest > 0 .and. parlest >= partmp) then

!----------------------------------------------------------------!
!----------------------------------------------------------------!
    call assert(d_initialized, 'D is initialized', srname)  ! Zaikun 20220512, TO BE REMOVED
!----------------------------------------------------------------!
!----------------------------------------------------------------!

    dtg = inprod(d, gg)
    d = -sign(delta / sqrt(dsq), dtg) * d  !!MATLAB: d = -sign(dtg) * (delta / sqrt(dsq)) * d
    goto 370
end if

220 continue

! Pick the value of PAR for the next iteration.
if (paru == ZERO) then
    par = TWO * parlest + gnorm / delta
else
    par = HALF * (parl + paru)
    par = max(par, parlest)
end if
if (paruest > 0) par = min(par, paruest)
goto 140

230 continue

! Calculate D = -(H + PAR*I)^{-1}*G for the current PAR in the positive definite case. The two loops
! below find D using the LDL factorization of the (tridiagonalized) H + PAR*I = L*diag(PIV)*L^T.
d(1) = -gg(1) / piv(1)
do k = 1, n - 1_IK  ! The loop sets D = -PIV^{-1}L^{-1}*GG
    d(k + 1) = -(gg(k + 1) + tn(k) * d(k)) / piv(k + 1)
end do
wsq = inprod(piv, d**2)  ! GG^T*(H+PAR*I)^{-1}*GG. Needed in the convergence test.
do k = n - 1_IK, 1, -1  ! The loop sets D = L^{-T}*D = -L^{-T}*PIV^{-1}*L^{-1}*GG = -(H+PAR*I)^{-1}*GG.
    d(k) = d(k) - tn(k) * d(k + 1) / piv(k)
end do

!----------------------------------------------------------------!
!----------------------------------------------------------------!
d_initialized = .true.  ! Zaikun 20220512, TO BE REMOVED
!----------------------------------------------------------------!
!----------------------------------------------------------------!


dsq = sum(d**2)

! Return if the Newton-Raphson step is feasible, setting CRVMIN to the least eigenvalue of Hessian.
if (par == ZERO .and. dsq <= delsq) then
    crvmin = eigmin(td, tn, 1.0E-2_RP)
    !!MATLAB:
    !!% It is critical for the efficiency to use `spdiags` to construct `tridh` in the sparse form.
    !!tridh = spdiags([[tn; 0], td, [0; tn]], -1:1, n, n);
    !!crvmin = eigs(tridh, 1, 'smallestreal');
    goto 370
end if

! Make the usual test for acceptability of a full trust region step.
dnorm = sqrt(dsq)
phi = ONE / dnorm - ONE / delta
if (tol * (ONE + par * dsq / wsq) - dsq * phi * phi >= 0) then
    d = (delta / dnorm) * d
    goto 370
end if
if (iter >= 2 .and. par <= parl) goto 370
if (paru > 0 .and. par >= paru) goto 370

! Complete the iteration when PHI is negative.
if (phi < 0) then
    parlest = par
    if (posdef) then
        if (phi <= phil) goto 370
        slope = (phi - phil) / (par - parl)
        parlest = par - phi / slope
    end if
    if (paru > 0) then
        slope = (phiu - phi) / (paru - par)
    else
        slope = ONE / gnorm
    end if
    partmp = par - phi / slope
    if (paruest > 0) then
        paruest = min(partmp, paruest)
    else
        paruest = partmp
    end if
    posdef = .true.
    parl = par
    phil = phi
    goto 220
end if

! If required, calculate Z for the alternative test for convergence.
! For information on Z, see the discussions below (16) in Section 2 of the UOBYQA paper (2002 version
! in Math. Program.; in the DAMTP 2000/NA14 report, it is below (2.8) in Section 2). The two loops
! below find Z using the LDL factorization of the (tridiagonalized) H + PAR*I.
if (.not. posdef) then
    z(1) = ONE / piv(1)
    do k = 1, n - 1_IK
        tnz = tn(k) * z(k)
        z(k + 1) = -(sign(ONE, tnz) + tnz) / piv(k + 1)
    end do
    wwsq = inprod(piv, z**2)  ! Needed in the convergence test.
    do k = n - 1_IK, 1, -1
        z(k) = z(k) - tn(k) * z(k + 1) / piv(k)
    end do

    zsq = sum(z**2)
    dtz = inprod(d, z)

    ! Apply the alternative test for convergence.
    tempa = abs(delsq - dsq)
    tempb = sqrt(dtz * dtz + tempa * zsq)
    gam = tempa / (sign(tempb, dtz) + dtz)
    if (tol * (wsq + par * delsq) - gam * gam * wwsq >= 0) then
        d = d + gam * z
        goto 370
    end if
    parlest = max(parlest, par - wwsq / zsq)
end if

! Complete the iteration when PHI is positive.
slope = ONE / gnorm
if (paru > 0) then
    if (phi >= phiu) goto 370
    slope = (phiu - phi) / (paru - par)
end if
parlest = max(parlest, par - phi / slope)
paruest = par
if (posdef) then
    slope = (phi - phil) / (par - parl)
    paruest = par - phi / slope
end if
paru = par
phiu = phi
goto 220

370 continue

! Apply the inverse Householder transformations to D.
do k = n - 1_IK, 1, -1
    d(k + 1:n) = d(k + 1:n) - inprod(d(k + 1:n), hh(k + 1:n, k)) * hh(k + 1:n, k)
end do
!!MATLAB: d = P*d;

end subroutine trstep


end module trustregion_mod
