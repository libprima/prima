module geometry_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the geometry-improving of the interpolation set XPT.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Monday, September 26, 2022 PM09:48:13
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: setdrop_tr, geostep


contains


function setdrop_tr(idz, kopt, tr_success, bmat, d, delta, rho, xpt, zmat) result(knew)
!--------------------------------------------------------------------------------------------------!
! This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
! REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
! XPT(:, KNEW) is replaced by XNEW = XOPT + D, where D is the trust-region step.
! N.B.:
! 1. If TR_SUCCESS = TRUE, then KNEW > 0 so that XNEW is included into XPT. Otherwise, it is a bug.
! 2. If TR_SUCCESS = FALSE, then KNEW /= KOPT so that XPT(:, KOPT) stays. Otherwise, it is a bug.
! 3. It is tempting to take the function value into consideration when defining KNEW, for example,
! set KNEW so that FVAL(KNEW) = MAX(FVAL) as long as F(XNEW) < MAX(FVAL), unless there is a better
! choice. However, this is not a good idea, because the definition of KNEW should benefit the
! quality of the model that interpolates f at XPT. A set of points with low function values is not
! necessarily a good interpolation set. In contrast, a good interpolation set needs to include
! points with relatively high function values; otherwise, the interpolant will unlikely reflect the
! landscape of the function sufficiently.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: HDIAG(NPT), DENABS(NPT), SCORE(NPT) VLAG(N+NPT), XDIST(NPT)
! Size of local arrays: REAL(RP)*(4*NPT+N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TENTH, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : issymmetric
use, non_intrinsic :: powalg_mod, only : calden
implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: kopt
logical, intent(in) :: tr_success
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)
real(RP), intent(in) :: delta
real(RP), intent(in) :: rho
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
integer(IK) :: knew

! Local variables
character(len=*), parameter :: srname = 'SETDROP_TR'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: denabs(size(xpt, 2))
real(RP) :: distsq(size(xpt, 2))
real(RP) :: score(size(xpt, 2))
real(RP) :: weight(size(xpt, 2))

! Sizes
n = int(size(xpt, 1), kind(npt))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(delta >= rho .and. rho > 0, 'DELTA >= RHO > 0', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! Calculate the distance squares between the interpolation points and the "optimal point". When
! identifying the optimal point, as suggested in (7.5) of the NEWUOA paper, it is reasonable to
! take into account the new trust-region trial point XPT(:, KOPT) + D, which will become the optimal
! point in the next interpolation if TR_SUCCESS is TRUE. Strangely, considering this new point does
! not always lead to a better performance of NEWUOA. Here, we choose not to check TR_SUCCESS, as
! the performance of NEWUOA is better in this way.
! HOWEVER, THIS MAY WELL CHANGE IF THE OTHER PARTS OF NEWUOA ARE IMPLEMENTED DIFFERENTLY.
!if (tr_success) then
!    distsq = sum((xpt - spread(xpt(:, kopt) + d, dim=2, ncopies=npt))**2, dim=1)
!    !!MATLAB: distsq = sum((xpt - (xpt(:, kopt) + d)).^2)  % d should be a column!! Implicit expansion
!else
!    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
!    !!MATLAB: distsq = sum((xpt - xpt(:, kopt)).^2)  % Implicit expansion
!end if
distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)

weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**3
!--------------------------------------------------------------------------------------------------!
! Other possible definitions of WEIGHT.
!weight = max(ONE, distsq / rho**2)**3  ! This works almost the same as the above.
!weight = max(ONE, distsq / delta**2)**3  ! Code from BOBYQA. It does not work as well as the above.
!weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**4  ! It does not work as well as the above.
!weight = max(ONE, distsq / max(TENTH * delta, rho)**2)**2  ! This works poorly.
!--------------------------------------------------------------------------------------------------!

denabs = abs(calden(kopt, bmat, d, xpt, zmat, idz))
score = weight * denabs

! If the new F is not better than FVAL(KOPT), we set SCORE(KOPT) = -1 to avoid KNEW = KOPT.
if (.not. tr_success) then
    score(kopt) = -ONE
end if

! Changing the following IF to `IF (ANY(SCORE > 0)) THEN` does not lead to a better performance.
if (any(score > 1) .or. (tr_success .and. any(score > 0))) then
    ! See (7.5) of the NEWUOA paper for the definition of KNEW in this case.
    ! SCORE(K) is NaN implies DENABS(K) is NaN, but we want DENABS to be big. So we exclude such K.
    knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
    !!MATLAB: [~, knew] = max(score, [], 'omitnan');
elseif (tr_success) then
    ! Powell's code does not include the following instructions. With Powell's code, if DENABS
    ! consists of only NaN, then KNEW can be 0 even when TR_SUCCESS is TRUE.
    knew = int(maxloc(distsq, dim=1), IK)
else
    knew = 0_IK  ! We arrive here when TR_SUCCESS = FALSE and no entry of SCORE exceeds one.
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(knew /= kopt .or. tr_success, 'KNEW /= KOPT unless TR_SUCCESS = TRUE', srname)
    call assert(knew >= 1 .or. .not. tr_success, 'KNEW >= 1 unless TR_SUCCESS = FALSE', srname)
    ! KNEW >= 1 when TR_SUCCESS = TRUE unless NaN occurs in DISTSQ, which should not happen if the
    ! starting point does not contain NaN and the trust-region/geometry steps never contain NaN.
end if

end function setdrop_tr


function geostep(idz, knew, kopt, bmat, delbar, xpt, zmat) result(d)
!--------------------------------------------------------------------------------------------------!
! This subroutine finds a step D that intends to improve the geometry of the interpolation set when
! XPT(:, KNEW) is changed to XOPT + D, where XOPT = XPT(:, KOPT).
!
! XPT contains the current interpolation points.
! BMAT provides the last N ROWs of H.
! ZMAT and IDZ give a factorization of the first NPT by NPT sub-matrix of H.
! KNEW is the index of the interpolation point to be dropped.
! DELBAR is the trust region bound for the geometry step
! D will be set to the step from X to the new point.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: D(N), DDEN, PQLAG(NPT), VLAG(N+NPT), XOPT(N)
! Size of local arrays: REAL(RP)*(2*NPT+4*N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : issymmetric, norm
use, non_intrinsic :: powalg_mod, only : omega_col, calvlag, calbeta
implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: delbar
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP) :: d(size(xpt, 1))     ! D(N)

! Local variables
character(len=*), parameter :: srname = 'GEOSTEP'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: alpha
real(RP) :: beta
real(RP) :: dden(size(xpt, 1))
real(RP) :: denom
real(RP) :: denrat
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))
real(RP) :: xopt(size(xpt, 1))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt, 'KNEW /= KOPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(delbar > 0, 'DELBAR > 0', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

d = biglag(idz, knew, bmat, delbar, xopt, xpt, zmat)

! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC.
pqlag = omega_col(idz, zmat, knew)
alpha = pqlag(knew)  ! ALPHA is the KNEW-th diagonal entry of H, i.e., that of Omega.

! Calculate VLAG and BETA for D. Indeed, only VLAG(KNEW) is needed.
vlag = calvlag(kopt, bmat, d, xpt, zmat, idz)
beta = calbeta(kopt, bmat, d, xpt, zmat, idz)
denom = alpha * beta + vlag(knew)**2

! If the cancellation in DENOM is unacceptable, then BIGDEN calculates an alternative model step D.
! As in (6.17) of the NEWUOA paper, DENRAT = |ALPHA*BETA + TAU^2| / TAU^2 with TAU = VLAG(KNEW).
! Powell's code does not check whether VLAG(KNEW)**2 > 0, which holds in theory. VLAG(KNEW) can
! become 0 or NaN numerically, which did happen in tests, indicating a failure of BIGLAG, because
! BIGLAG should maximize |VLAG(KNEW)|. Upon this failure, it is reasonable to call BIGDEN. For the
! same reason, we check whether BETA is NaN. Why not check ALPHA? Because BIGDEN cannot improve ALPHA.
! Powell's code takes DDEN once it is calculated. We take it only if it renders a bigger denominator.
denrat = -ONE
if (vlag(knew)**2 > 0 .and. .not. is_nan(beta)) then
    denrat = abs(ONE + alpha * beta / vlag(knew)**2)
end if
! If DENRAT is NaN at this point, then ALPHA is NaN, and there is no need to call BIGDEN.
if (denrat <= 0.8_RP) then
    dden = bigden(idz, knew, kopt, bmat, d, xpt, zmat)
    vlag = calvlag(kopt, bmat, dden, xpt, zmat, idz)
    beta = calbeta(kopt, bmat, dden, xpt, zmat, idz)
    if (abs(alpha * beta + vlag(knew)**2) >= abs(denom) .or. is_nan(denom)) then
        d = dden
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(norm(d) <= TWO * delbar, '|D| <= 2*DELBAR', srname)
    ! Due to rounding, it may happen that |D| > DELBAR, but |D| > 2*DELBAR is highly improbable.
end if

end function geostep


function biglag(idz, knew, bmat, delbar, x, xpt, zmat) result(d)
!--------------------------------------------------------------------------------------------------!
! This subroutine calculates a D by approximately solving
!
! max |LFUNC(X + D)|, subject to ||D|| <= DELBAR,
!
! where LFUNC is the KNEW-th Lagrange function. See Section 6 of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: D(N), CF(5), DOLD(N), GC(N), GD(N), PQLAG(NPT), S(N), W(N)
! Size of local arrays: REAL(RP)*(5+6*N+NPT)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, QUART, TENTH, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, issymmetric, norm, project
use, non_intrinsic :: powalg_mod, only : omega_col, hess_mul
use, non_intrinsic :: univar_mod, only : circle_maxabs

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: delbar
real(RP), intent(in) :: x(:)        ! X(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP) :: d(size(xpt, 1))       ! D(N)

! Local variables
character(len=*), parameter :: srname = 'BIGLAG'
integer(IK) :: iter
integer(IK) :: maxiter
integer(IK) :: n
integer(IK) :: npt
real(RP) :: angle
real(RP) :: cf(5)
real(RP) :: cth
real(RP) :: dd
real(RP) :: dhd
real(RP) :: dold(size(x))
real(RP) :: gc(size(x))
real(RP) :: gd(size(x))
real(RP) :: gg
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: s(size(x))
real(RP) :: scaling
real(RP) :: sp
real(RP) :: ss
real(RP) :: sth
real(RP) :: t
real(RP) :: tau  ! LFUNC(X)
real(RP) :: tol
real(RP) :: w(size(x))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(delbar > 0, 'DELBAR > 0', srname)
    call assert(size(x) == n .and. all(is_finite(x)), 'SIZE(X) == N, X is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC.
pqlag = omega_col(idz, zmat, knew)

! Set the unscaled initial D. Form the gradient of LFUNC at X, and multiply D by the Hessian of LFUNC.
d = xpt(:, knew) - x
dd = inprod(d, d)
gd = hess_mul(d, xpt, pqlag)  ! GD = MATPROD(XPT, PQLAG * MATPROD(D, XPT))

!--------------------------------------------------------------------------------------------------!
gc = bmat(:, knew) + hess_mul(x, xpt, pqlag) ! GC = BMAT(:,KNEW) + MATPROD(XPT,PQLAG*MATPROD(X,XPT))
! The following is mathematically equivalent to the last but seems to work numerically better (why?)
!gc = Ax_plus_y(xpt, pqlag * matprod(x, xpt), bmat(:, knew))
!--------------------------------------------------------------------------------------------------!

! Scale D and GD, with a sign change if needed. Set S to another vector in the initial 2-D subspace.
gg = inprod(gc, gc)
sp = inprod(d, gc)
dhd = inprod(d, gd)
scaling = delbar / sqrt(dd)
if (sp * dhd < ZERO) then
    scaling = -scaling
end if
t = ZERO
if (sp**2 > 0.99_RP * dd * gg) then
    t = ONE
end if
tau = scaling * (abs(sp) + HALF * scaling * abs(dhd))
if (gg * delbar**2 < 1.0E-2_RP * tau**2) then
    t = ONE
end if
if (is_finite(sum(abs(scaling * d)))) then
    d = scaling * d
    gd = scaling * gd
    s = gc + t * gd
    maxiter = n
else
    maxiter = 0_IK  ! Return immediately to avoid producing a D containing NaN/Inf.
end if

! Begin the iteration by overwriting S with a vector that has the required length and direction,
! except that termination occurs if the given D and S are nearly parallel.
! TOL is the tolerance for telling whether S and D are nearly parallel. In Powell's code, the
! tolerance is 1.0D-4. We adapt it to the following value in case single precision is in use.
tol = min(1.0E-1_RP, max(EPS**QUART, 1.0E-4_RP))
do iter = 1, maxiter

    !----------------------------------------------------------------------------------------------!
    !!--------------------------------------------------!
    !!ds = inprod(d, s)
    !!ss = inprod(s, s)
    !!if (dd * ss - ds**2 <= 1.0E-8_RP * dd * ss) then
    !!    exit
    !!end if
    !!denom = sqrt(dd * ss - ds**2)
    !!s = (dd * s - ds * d) / denom
    !!--------------------------------------------------!
    ! The above is Powell's code. In precise arithmetic, INPROD(S, D) = 0 and |S| = |D|. However,
    ! when DD*SS - DS**2 is tiny, the error in S can be large and hence damage these inequalities
    ! significantly. This did happen in tests, especially when using the single precision. We
    ! calculate S as below. It did improve the performance of NEWUOA in our test.
    !----------------------------------------------------------------------------------------------!

    ss = inprod(s, s)
    s = s - project(s, d)  ! PROJECT(X, V) returns the projection of X to SPAN(V): X'*(V/|V|)*(V/|V|)
    ! N.B.:
    ! 1. The condition |S|<=TOL*SQRT(SS) below is equivalent to DS^2>=(1-TOL^2)*DD*SS in theory. As
    ! shown above, Powell's code triggers an exit if DS^2>=(1-1.0E-8)*DD*SS. So our condition is the
    ! same except that we take the machine epsilon into account in case single precision is in use.
    ! 2. The condition below should be non-strict so that |S| = 0 can trigger the exit.
    if (norm(s) <= tol * sqrt(ss)) then
        exit
    end if
    s = (norm(d) / norm(s)) * s

    ! In precise arithmetic, INPROD(S, D) = 0 and |S| = |D| = DELBAR.
    if (abs(inprod(d, s)) >= TENTH * norm(d) * norm(s) .or. norm(s) >= TWO * delbar) then
        exit
    end if

    w = hess_mul(s, xpt, pqlag)  ! W = MATPROD(XPT, PQLAG * MATPROD(S, XPT))

    ! Seek the value of the angle that maximizes |TAU|.
    ! First, calculate the coefficients of the objective function on the circle.
    cf(1) = HALF * inprod(s, w)
    cf(2) = inprod(d, gc)
    cf(3) = inprod(s, gc)
    cf(4) = HALF * inprod(d, gd) - cf(1)
    cf(5) = inprod(s, gd)
    ! The 50 in the line below was chosen by Powell. It works the best in tests, MAGICALLY. Larger
    ! (e.g., 60, 100) or smaller (e.g., 20, 40) values will worsen the performance of NEWUOA. Why??
    angle = circle_maxabs(circle_fun_biglag, cf, 50_IK)

    ! Calculate the new D and GD.
    cth = cos(angle)
    sth = sin(angle)
    dold = d
    d = cth * d + sth * s

    ! Exit in case of Inf/NaN in D.
    if (.not. is_finite(sum(abs(d)))) then
        d = dold
        exit
    end if

    ! Test for convergence.
    if (abs(circle_fun_biglag(angle, cf)) <= 1.1_RP * abs(circle_fun_biglag(ZERO, cf))) then
        exit
    end if

    ! Calculate GD and S.
    gd = cth * gd + sth * w
    s = gc + gd
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(norm(d) <= TWO * delbar, '|D| <= 2*DELBAR', srname)
    ! Due to rounding, it may happen that |D| > DELBAR, but |D| > 2*DELBAR is highly improbable.
end if

end function biglag


function bigden(idz, knew, kopt, bmat, d0, xpt, zmat) result(d)
!--------------------------------------------------------------------------------------------------!
! BIGDEN calculates a D by approximately solving
!
! max |SIGMA(XOPT + D)|, subject to ||D|| <= DELBAR,
!
! where SIGMA is the denominator sigma in the updating formula (4.11)--(4.12) for H, which is the
! inverse of the coefficient matrix for the interpolation system (see (3.12)). Indeed, each column
! of H corresponds to a Lagrange basis function of the interpolation problem.  See Section 6 of the
! NEWUOA paper.
! N.B.:
! In Powell's code, BIGDEN calculates also the VLAG and BETA for the selected D. Here, to reduce the
! coupling of code, we return only D but compute VLAG and BETA outside by calling VLAGBETA. It makes
! no difference mathematically, but the computed VLAG/BETA will change slightly due to rounding.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: D(N), DEN(9), DENEX(9), DOLD(N), DSTEMP(NPT), PQLAG(NPT), PAR(5), PROD(NPT+N, 5), &
!     & S(N), SSTEMP(NPT), V(NPT), VLAG(NPT+N), W(NPT+N, 5), X(N), XNEW(N), XPTEMP(N, NPT)
! Size of local arrays: REAL(RP)*(23+12*N+11*NPT+N*NPT) (TO BE REDUCED by removing XPTEMP)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, QUART, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod, issymmetric, norm, project
use, non_intrinsic :: powalg_mod, only : omega_col, omega_mul
use, non_intrinsic :: univar_mod, only : circle_maxabs

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: d0(:)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP) :: d(size(xpt, 1))     ! D(N)

! Local variable
character(len=*), parameter :: srname = 'BIGDEN'
integer(IK) :: iter
integer(IK) :: j
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: nw
real(RP) :: alpha
real(RP) :: angle
real(RP) :: dd
real(RP) :: den(9)
real(RP) :: denex(9)
real(RP) :: denmax
real(RP) :: densav
real(RP) :: dold(size(xpt, 1))
real(RP) :: ds
real(RP) :: dstemp(size(xpt, 2))
real(RP) :: dtest
real(RP) :: dxn
real(RP) :: pqlag(size(xpt, 2))
real(RP) :: par(5)
real(RP) :: prod(size(xpt, 1) + size(xpt, 2), 5)
real(RP) :: s(size(xpt, 1))
real(RP) :: ss
real(RP) :: sstemp(size(xpt, 2))
real(RP) :: tau
real(RP) :: tempa
real(RP) :: tempb
real(RP) :: tempc
real(RP) :: tol
real(RP) :: v(size(xpt, 2))
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))
real(RP) :: w(size(xpt, 1) + size(xpt, 2), 5)
real(RP) :: x(size(xpt, 1))
real(RP) :: xd
real(RP) :: xnew(size(xpt, 1))
real(RP) :: xnsq
real(RP) :: xptemp(size(xpt, 1), size(xpt, 2))
real(RP) :: xs
real(RP) :: xsq

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(knew /= kopt, 'KNEW /= KOPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(d0) == n .and. all(is_finite(d0)), 'SIZE(D0) == N, D0 is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
end if

!====================!
! Calculation starts !
!====================!

x = xpt(:, kopt) ! For simplicity, we use X to denote XOPT.

! PQLAG contains the leading NPT elements of the KNEW-th column of H, and it provides the second
! derivative parameters of LFUNC.
pqlag = omega_col(idz, zmat, knew)
alpha = pqlag(knew)  ! ALPHA is the KNEW-th diagonal entry of H, i.e., that of Omega.

! The initial search direction D is taken from the last call of BIGLAG, and the initial S is set
! below, usually to the direction from X to X_KNEW, but a different direction to an interpolation
! point may be chosen, in order to prevent S from being nearly parallel to D.
d = d0
dd = inprod(d, d)
s = xpt(:, knew) - x
ds = inprod(d, s)
ss = inprod(s, s)
xsq = inprod(x, x)

if (.not. (ds**2 <= 0.99_RP * dd * ss)) then
    ! `.NOT. (A <= B)` differs from `A > B`.  The former holds iff A > B or {A, B} contains NaN.
    dtest = ds**2 / ss
    xptemp = xpt - spread(x, dim=2, ncopies=npt)
    !!MATLAB: xptemp = xpt - x  % x should be a column!! Implicit expansion
    !----------------------------------------------------------------!
    !---------!dstemp = matprod(d, xpt) - inprod(x, d) !-------------!
    dstemp = matprod(d, xptemp)
    !----------------------------------------------------------------!
    sstemp = sum((xptemp)**2, dim=1)

    dstemp(kopt) = TWO * ds + ONE
    sstemp(kopt) = ss
    k = int(minloc(dstemp**2 / sstemp, dim=1), kind(k))
    ! K can be 0 due to NaN. In that case, set K = KNEW. Otherwise, memory errors will occur.
    if (k == 0) then
        k = knew
    end if
    if ((.not. (dstemp(k)**2 / sstemp(k) >= dtest)) .and. k /= kopt) then
        ! `.NOT. (A >= B)` differs from `A < B`.  The former holds iff A < B or {A, B} contains NaN.
        ! Although unlikely, if NaN occurs, it may happen that K = KOPT.
        s = xpt(:, k) - x
    end if
end if

densav = ZERO

! Begin the iteration by overwriting S with a vector that has the required length and direction.

! TOL is the tolerance for telling whether S and D are nearly parallel. In Powell's code, the
! tolerance is 1.0D-4. We adapt it to the following value in case single precision is in use.
tol = min(1.0E-1_RP, max(EPS**QUART, 1.0E-4_RP))
do iter = 1, n

    !----------------------------------------------------------------------------------------------!
    !!--------------------------------------------!
    !!ds = inprod(d, s)
    !!ss = inprod(s, s)
    !!ssden = dd * ss - ds**2
    !!if (ssden < 1.0E-8_RP * dd * ss) then
    !    exit
    !end if
    !!s = (ONE / sqrt(ssden)) * (dd * s - ds * d)
    !!--------------------------------------------!
    ! The above is Powell's code. In precise arithmetic, INPROD(S, D) = 0 and |S| = |D|. However,
    ! when DD*SS - DS**2 is tiny, the error in S can be large and hence damage these inequalities
    ! significantly. This did happen in tests, especially when using the single precision. We
    ! calculate S as below. It did improve the performance of NEWUOA in our test.
    !----------------------------------------------------------------------------------------------!

    ss = inprod(s, s)
    s = s - project(s, d)  ! PROJECT(X, V) returns the projection of X to SPAN(V): X'*(V/|V|)*(V/|V|)
    ! N.B.:
    ! 1. The condition |S|<=TOL*SQRT(SS) below is equivalent to DS^2>=(1-TOL^2)*DD*SS in theory. As
    ! shown above, Powell's code triggers an exit if DS^2>=(1-1.0E-8)*DD*SS. So our condition is the
    ! same except that we take the machine epsilon into account in case single precision is in use.
    ! 2. The condition below should be non-strict so that |S| = 0 can trigger the exit.
    if (norm(s) <= tol * sqrt(ss)) then
        exit
    end if
    s = (s / norm(s)) * norm(d)
    ! In precise arithmetic, INPROD(S, D) = 0 and |S| = |D| = DELBAR = |D0|.
    if (abs(inprod(d, s)) >= TENTH * norm(d) * norm(s) .or. norm(s) >= TWO * norm(d0)) then
        exit
    end if

    ! Set the coefficients of the first two terms of BETA.
    xd = inprod(x, d)
    xs = inprod(x, s)
    dd = inprod(d, d)
    tempa = HALF * xd * xd
    tempb = HALF * xs * xs
    den(1) = dd * (xsq + HALF * dd) + tempa + tempb
    den(2) = TWO * xd * dd
    den(3) = TWO * xs * dd
    den(4) = tempa - tempb
    den(5) = xd * xs
    den(6:9) = ZERO

    ! Put the coefficients of WCHECK in W.
    do k = 1, npt
        tempa = inprod(xpt(:, k), d)
        tempb = inprod(xpt(:, k), s)
        tempc = inprod(xpt(:, k), x)
        w(k, 1) = QUART * (tempa**2 + tempb**2)
        w(k, 2) = tempa * tempc
        w(k, 3) = tempb * tempc
        w(k, 4) = QUART * (tempa**2 - tempb**2)
        w(k, 5) = HALF * tempa * tempb
    end do
    w(npt + 1:npt + n, 1:5) = ZERO
    w(npt + 1:npt + n, 2) = d
    w(npt + 1:npt + n, 3) = s

    ! Put the coefficients of THETA*WCHECK in PROD.
    do j = 1, 5
        prod(1:npt, j) = omega_mul(idz, zmat, w(1:npt, j))
        nw = npt
        if (j == 2 .or. j == 3) then
            prod(1:npt, j) = prod(1:npt, j) + matprod(w(npt + 1:npt + n, j), bmat(:, 1:npt))
            nw = npt + n
        end if
        prod(npt + 1:npt + n, j) = matprod(bmat(:, 1:nw), w(1:nw, j))
    end do

    ! Include in DEN the part of BETA that depends on THETA.
    do k = 1, npt + n
        par(1:5) = HALF * prod(k, 1:5) * w(k, 1:5)
        den(1) = den(1) - par(1) - sum(par(1:5))
        tempa = prod(k, 1) * w(k, 2) + prod(k, 2) * w(k, 1)
        tempb = prod(k, 2) * w(k, 4) + prod(k, 4) * w(k, 2)
        tempc = prod(k, 3) * w(k, 5) + prod(k, 5) * w(k, 3)
        den(2) = den(2) - tempa - HALF * (tempb + tempc)
        den(6) = den(6) - HALF * (tempb - tempc)
        tempa = prod(k, 1) * w(k, 3) + prod(k, 3) * w(k, 1)
        tempb = prod(k, 2) * w(k, 5) + prod(k, 5) * w(k, 2)
        tempc = prod(k, 3) * w(k, 4) + prod(k, 4) * w(k, 3)
        den(3) = den(3) - tempa - HALF * (tempb - tempc)
        den(7) = den(7) - HALF * (tempb + tempc)
        tempa = prod(k, 1) * w(k, 4) + prod(k, 4) * w(k, 1)
        den(4) = den(4) - tempa - par(2) + par(3)
        tempa = prod(k, 1) * w(k, 5) + prod(k, 5) * w(k, 1)
        tempb = prod(k, 2) * w(k, 3) + prod(k, 3) * w(k, 2)
        den(5) = den(5) - tempa - HALF * tempb
        den(8) = den(8) - par(4) + par(5)
        tempa = prod(k, 4) * w(k, 5) + prod(k, 5) * w(k, 4)
        den(9) = den(9) - HALF * tempa
    end do

    par(1:5) = HALF * prod(knew, 1:5)**2
    denex(1) = alpha * den(1) + par(1) + sum(par(1:5))
    tempa = TWO * prod(knew, 1) * prod(knew, 2)
    tempb = prod(knew, 2) * prod(knew, 4)
    tempc = prod(knew, 3) * prod(knew, 5)
    denex(2) = alpha * den(2) + tempa + tempb + tempc
    denex(6) = alpha * den(6) + tempb - tempc
    tempa = TWO * prod(knew, 1) * prod(knew, 3)
    tempb = prod(knew, 2) * prod(knew, 5)
    tempc = prod(knew, 3) * prod(knew, 4)
    denex(3) = alpha * den(3) + tempa + tempb - tempc
    denex(7) = alpha * den(7) + tempb + tempc
    tempa = TWO * prod(knew, 1) * prod(knew, 4)
    denex(4) = alpha * den(4) + tempa + par(2) - par(3)
    tempa = TWO * prod(knew, 1) * prod(knew, 5)
    denex(5) = alpha * den(5) + tempa + prod(knew, 2) * prod(knew, 3)
    denex(8) = alpha * den(8) + par(4) - par(5)
    denex(9) = alpha * den(9) + prod(knew, 4) * prod(knew, 5)

    ! Seek the value of the angle that maximizes the |DENOM|.
    angle = circle_maxabs(circle_fun_bigden, denex, 50_IK)

    ! Calculate the new D.
    dold = d
    d = cos(angle) * d + sin(angle) * s

    ! Exit in case of Inf/NaN in D.
    if (.not. is_finite(sum(abs(d)))) then
        d = dold
        exit
    end if

    ! Test for convergence.
    if (iter > 1) then
        densav = max(densav, circle_fun_bigden(ZERO, denex))
    end if
    denmax = circle_fun_bigden(angle, denex)
    if (abs(denmax) <= 1.1_RP * abs(densav)) then
        exit
    end if
    densav = denmax

    ! Set S to HALF the gradient of the denominator with respect to D. First, calculate the new VLAG.
    par = [ONE, cos(angle), sin(angle), cos(2.0_RP * angle), sin(2.0_RP * angle)]
    vlag = matprod(prod, par)
    tau = vlag(knew)
    xnew = x + d
    dxn = inprod(d, xnew)
    xnsq = inprod(xnew, xnew)
    v = (tau * pqlag - alpha * vlag(1:npt)) * matprod(xnew, xpt)
    s = tau * bmat(:, knew) + alpha * (dxn * x + xnsq * d - vlag(npt + 1:npt + n))
    s = s + matprod(xpt, v)
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(norm(d) <= TWO * norm(d0), '|D| <= 2*DELBAR', srname)
    ! Due to rounding, it may happen that |D| > DELBAR, but |D| > 2*DELBAR is highly improbable.
end if

end function bigden


function circle_fun_biglag(theta, args) result(f)
!--------------------------------------------------------------------------------------------------!
! This function defines the objective function of the 2D search on a circle in BIGLAG.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): NONE
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none
! Inputs
real(RP), intent(in) :: theta
real(RP), intent(in) :: args(:)

! Outputs
real(RP) :: f

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_FUN_BIGLAG'
real(RP) :: cth
real(RP) :: sth

! Preconditions
if (DEBUGGING) then
    call assert(size(args) == 5, 'SIZE(ARGS) == 5', srname)
end if

!====================!
! Calculation starts !
!====================!

cth = cos(theta)
sth = sin(theta)
f = args(1) + (args(2) + args(4) * cth) * cth + (args(3) + args(5) * cth) * sth

!====================!
!  Calculation ends  !
!====================!
end function circle_fun_biglag


function circle_fun_bigden(theta, args) result(f)
!--------------------------------------------------------------------------------------------------!
! This function defines the objective function of the 2D search on a circle in BIGDEN.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack): NONE
! REAL(RP) :: PAR(9)
! Size of local arrays: REAL(RP)*9
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, ONE, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : inprod
implicit none
! Inputs
real(RP), intent(in) :: theta
real(RP), intent(in) :: args(:)

! Outputs
real(RP) :: f

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_FUN_BIGDEN'
real(RP) :: par(size(args))

! Preconditions
if (DEBUGGING) then
    call assert(size(args) == 9, 'SIZE(ARGS) == 9', srname)
end if

!====================!
! Calculation starts !
!====================!

par(1) = ONE
par(2:8:2) = cos(theta*[1.0_RP, 2.0_RP, 3.0_RP, 4.0_RP])
par(3:9:2) = sin(theta*[1.0_RP, 2.0_RP, 3.0_RP, 4.0_RP])
f = inprod(args, par)

!====================!
!  Calculation ends  !
!====================!
end function circle_fun_bigden


end module geometry_mod
