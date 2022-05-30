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
! Last Modified: Tuesday, May 31, 2022 AM02:05:54
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
!
! The major calculations of the More-Sorensen method lie in the Cholesky factorization or LDL
! factorization of H + PAR*I with the iteratively selected values of PAR (in More-Sorensen 1983, the
! parameter is named LAMBDA; in Powell's UOBYQA paper, it is THETA). Powell's method in this code
! simplifies the calculations by first tridiagonalizing H with an orthogonal transformation. If a
! matrix T is tridiagonal, its LDL factorization, if exits, can be obtained easily:
!
!       T = L*diag(PIV)*L^T,
!
! where diag(PIV) is the diagonal matrix with the diagonal entries being PIV(1:N), i.e., "the pivots
! of the Cholesky factorization" in Powell's comments on his code, and L is the lower triangular
! matrix with all the diagonal entries being 1, the subdiagonal being the subdiagonal of T divided
! by PIV(1:N-1), and all the other entries being 0. PIV can be obtained by a simple recursion.
!
! For more information, see Section 2 of the UOBYQA paper and
! Powell, M. J. D., "Trust region calculations revisited", Numerical Analysis 1997: Proceedings of
! the 17th Dundee Biennial Numerical Analysis Conference, 1997, 193--211,
! Powell, M. J. D., "The use of band matrices for second derivative approximations in trust region
! algorithms", Advances in Nonlinear Programming: Proceedings of the 96 International Conference on
! Nonlinear Programming, 1998, 3--28.
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
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
integer(IK) :: i, iter, k, maxiter
logical :: posdef, negcrv
logical :: d_initialized  ! TO BE REMOVED.

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

if (is_nan(gsq)) then
    return
elseif (.not. any(abs(h) > 0)) then
    if (gnorm > 0) then
        d = -(delta / gnorm) * g
    end if
    return
end if

! Zaikun 20220301, 20220305:
! Powell's original code requires that N >= 2.  When N = 1, the code does not work (sometimes even
! encounters memory errors). This is indeed why the original UOBYQA code constantly terminates with
! "a trust region step has failed to reduce the quadratic model" when applied to univariate problems.
if (n == 1) then
    d = sign(delta, -g)  !!MATLAB: d = -delta * sign(g)
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
paru = ZERO  ! Upper bound for the optimal PAR ??? The initial value is less than PARL. Why?
paruest = ZERO  ! Estimation for PARU
posdef = .false.
iter = 0
maxiter = min(1000_IK, 100_IK * int(n, IK))  ! What is the theoretical bound of iterations?

do while (.true.)
    iter = iter + 1_IK

! Zaikun 26-06-2019: The original code can encounter infinite cycling, which did happen when testing
! the CUTEst problems GAUSS1LS, GAUSS2LS, and GAUSS3LS. Indeed, in all these cases, Inf and NaN
! appear in D due to extremely large values in the Hessian matrix (up to 10^219).
    if (.not. is_finite(sum(abs(d)))) then
        d = dold
        exit
    else
        dold = d
    end if
    if (iter > maxiter) then
        exit
    end if

! Calculate the pivots of the Cholesky factorization of (H + PAR*I), which correspond to the squares
! of the diagonal entries of L in the Cholesky factorization LL^T, or the diagonal matrix in the LDL
! factorization. After getting PIV, we can get the LDL factorization of H + PAR*I easily: it is
! L*diag(PIV)*L^T, where diag(PIV) is the diagonal matrix with PIV being the diagonal, and L is the
! lower triangular matrix with all the diagonal entries being 1, the subdiagonal being the vector
! TN/PIV(1:N-1) (entrywise division), and all the other entries being 0.
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
        exit  ! Better action to take???
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

! Handle the case where H + PAR*I is positive semidefinite.
    if (.not. negcrv) then
        ! Handle the case where the gradient at the trust region center is zero.
        if (gsq <= 0) then
            paru = par
            paruest = par
            if (par == ZERO) then  ! A rare case: the trust region center is optimal.
                exit
            end if
        end if
    end if


! At this point, K == 0 iff H + PAR*I is positive definite.

    if (k >= 1) then
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! We arrive here only if 1 <= K <= N, when H + PAR*I has at least one nonpositive eigenvalue.
        call assert(k >= 1 .and. k <= n, '1 <= K <= N', srname)
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


! Handle the case where H + PAR*I has at least one nonpositive eigenvalue.
! Set D to a direction of nonpositive curvature of the tridiagonal matrix, and thus revise PARLEST.

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! Zaikun 20220512: Powell's code does not include the following initialization. Consequently,
! D(KSAV+1:N) or D(KSAV+2:N) will not be initialized but inherit values from the previous iteration.
! Is this intended?
        d = ZERO
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

        d(k) = ONE  ! Zaikun 20220512: D(K+1:N) = ?

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! The code until "Terminate with D set to a multipe of the current D ..." sets only D(1:KSAV) or
! D(1:KSAV+1), with the KSAV defined later. D_INITIALIZED indicates whether D(1:N) is fully
! initialized in this process (TRUE) or not (FALSE). See the comments above
! CALL WASSERT(D_INITIALIZED, 'D IS INITIALIZED', SRNAME) for details.
        d_initialized = (k == n)  ! Zaikun 20220512, TO BE REMOVED
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

        dhd = piv(k)

! In Fortran, the following two IFs CANNOT be merged into IF(K < N .AND. ABS(TN(K)) > ABS(PIV(K))).
! This is because Fortran may not perform a short-circuit evaluation of this logic expression, and
! hence TN(K) may be accessed even if K >= N, leading to an out-of-boundary index since SIZE(TN) is
! only N-1. This is not a problem in C, MATLAB, Python, Julia, or R, where short circuit is ensured.
        if (k < n) then
            if (abs(tn(k)) > abs(piv(k))) then

                !---------------------------------------------------------------------!
                !---------------------------------------------------------------------!
                d_initialized = (k == n - 1)  ! Zaikun 20220512, TO BE REMOVED
                !---------------------------------------------------------------------!
                !---------------------------------------------------------------------!

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

        do i = k - 1_IK, 1, -1
            ! It may happen that TN(I) == 0 == PIV(I). Without checking TN(I), we will get D(I) = NaN.
            ! Once we encounter a zero TN(I), D(I) is set to zero, and D(1:I-1) will consequently be
            ! zero as well, because D(J) is a multiple of D(J+1) for each J.
            if (tn(i) /= ZERO) then
                d(i) = -tn(i) * d(i + 1) / piv(i)
            else
                d(1:i) = ZERO
                exit
            end if
        end do

        dsq = sum(d**2)
        parl = par
        parlest = par - dhd / dsq
    end if

    if (gsq <= 0 .or. k >= 1) then

        ! Terminate with D set to a multiple of the current D if the following test suggests so.
        if (gsq <= 0) then
            partmp = paruest * (ONE - tol)
        else
            partmp = paruest
        end if
        if (paruest > 0 .and. parlest >= partmp) then

            !----------------------------------------------------------------------------------------------!
            !----------------------------------------------------------------------------------------------!
            ! Zaikun 20220512, TO BE REMOVED
            ! The definition of D below requires that D is initialized. In Powell's code, it may happen that
            ! only D(1:KSAV) or D(1:KSAV+1) is initialized during the current iteration, but the other
            ! entries are inherited from the previous iteration OR from the initial value before the
            ! iterations start, which is 0. If such inheriting happens, D_INITIALIZED will be FALSE. In
            ! tests on 20220514, both cases did occur. Interestingly, in both cases, the inherited values
            ! were all zero or close to zero (1E-16), and hence not very different from the initial value
            ! zero that we set above. Is this intended?
            call wassert(d_initialized, 'D is initialized', srname)
            !----------------------------------------------------------------------------------------------!
            !----------------------------------------------------------------------------------------------!

            dtg = inprod(d, gg)
            if (dtg > 0) then
                d = -(delta / sqrt(dsq)) * d
            else  ! This ELSE covers the unlikely yet possible case where DTG is zero or even NaN.
                d = (delta / sqrt(dsq)) * d
            end if
            ! N.B.: As per Powell's code, the lines above would be D = -SIGN(DELTA/SQRT(DSQ), DTG) * D.
            ! However, our version here seems more reasonable in case DTG == 0, which is unlikely but did
            ! happen numerically. Note that SIGN(A, 0) = |A| /= -SIGN(A, 0).
            exit
        end if

! Pick the value of PAR for the next iteration.
        if (paru == ZERO) then
            par = TWO * parlest + gnorm / delta
        else
            par = HALF * (parl + paru)
            par = max(par, parlest)
        end if
        if (paruest > 0) par = min(par, paruest)
    else

! Handle the case where the gradient at the trust region center is nonzero and H + PAR*I is
! positive definite.
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
            exit
        end if

! Make the usual test for acceptability of a full trust region step.
        dnorm = sqrt(dsq)
        phi = ONE / dnorm - ONE / delta
        if (tol * (ONE + par * dsq / wsq) - dsq * phi * phi >= 0) then
            d = (delta / dnorm) * d
            exit
        end if
        if (iter >= 2 .and. par <= parl) exit
        if (paru > 0 .and. par >= paru) exit

! Complete the iteration when PHI is negative.
        if (phi < 0) then
            parlest = par
            if (posdef) then
                if (phi <= phil) exit
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
! Pick the value of PAR for the next iteration.
            if (paru == ZERO) then
                par = TWO * parlest + gnorm / delta
            else
                par = HALF * (parl + paru)
                par = max(par, parlest)
            end if
            if (paruest > 0) par = min(par, paruest)
            cycle
        end if

! If required, calculate Z for the alternative test for convergence.
! For information on Z, see the discussions below (16) in Section 2 of the UOBYQA paper (the 2002
! version in Math. Program.; in the DAMTP 2000/NA14 report, it is below (2.8) in Section 2). The two
! loops below find Z using the LDL factorization of the (tridiagonalized) H + PAR*I.
        if (.not. posdef) then
            z(1) = ONE / piv(1)
            do k = 1, n - 1_IK
                tnz = tn(k) * z(k)
                if (tnz > 0) then
                    z(k + 1) = -(ONE + tnz) / piv(k + 1)
                else
                    z(k + 1) = (ONE - tnz) / piv(k + 1)
                end if
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
            if (abs(dtz) > 0) then
                gam = tempa / (sign(tempb, dtz) + dtz)  !!MATLAB: gam = tempa / (sign(dtz)*tempb + dtz)
            else  ! This ELSE covers the unlikely yet possible case where DTZ is zero or even NaN.
                gam = sqrt(tempa / zsq)
            end if
            if (tol * (wsq + par * delsq) - gam * gam * wwsq >= 0) then
                d = d + gam * z
                exit
            end if
            parlest = max(parlest, par - wwsq / zsq)
        end if

! Complete the iteration when PHI is positive.
        slope = ONE / gnorm
        if (paru > 0) then
            if (phi >= phiu) exit
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
! Pick the value of PAR for the next iteration.
        if (paru == ZERO) then
            par = TWO * parlest + gnorm / delta
        else
            par = HALF * (parl + paru)
            par = max(par, parlest)
        end if
        if (paruest > 0) par = min(par, paruest)
    end if
end do

! Apply the inverse Householder transformations to D.
do k = n - 1_IK, 1, -1
    d(k + 1:n) = d(k + 1:n) - inprod(d(k + 1:n), hh(k + 1:n, k)) * hh(k + 1:n, k)
end do
!!MATLAB: d = P*d;

end subroutine trstep


end module trustregion_mod
