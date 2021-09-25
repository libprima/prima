! GEOMETRY_MOD is a module containing subroutines that are concerned with geometry-improving of the
! interpolation set XPT.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Last Modified: Friday, July 23, 2021 AM11:52:31

module geometry_mod

implicit none
private
public :: setdrop_tr, geostep


contains


function setdrop_tr(idz, kopt, beta, delta, ratio, rho, vlag, xopt, xpt, zmat) result(knew)
! SETDROP sets KNEW to the index of the interpolation point that will be deleted AFTER A TRUST
! REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
! XPT(:, KNEW) is replaced by XNEW = XOPT + D, where D is the trust-region step.  Note that the
! information of XNEW is included in VLAG and BETA, which are calculated according to D.
!
! N.B.: At the entry of this function is invoked, XOPT may differ from XPT(:, KOPT), because XOPT is
! updated but KOPT is not. See NEWUOB for details.

! Generic modules
use consts_mod, only : RP, IK, ONE, ZERO, TENTH, SRNLEN, DEBUGGING
use debug_mod, only : errstop, verisize

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: beta
real(RP), intent(in) :: delta
real(RP), intent(in) :: ratio
real(RP), intent(in) :: rho
real(RP), intent(in) :: vlag(:)  ! VLAG(NPT)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Output
integer(IK) :: knew

! Local variables
integer(IK) :: n
integer(IK) :: npt
real(RP) :: hdiag(size(zmat, 1))
real(RP) :: rhosq
real(RP) :: sigma(size(xpt, 2))
real(RP) :: xdsq(size(xpt, 2))
character(len=SRNLEN), parameter :: srname = 'SETDROP'


! Get and verify the sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

if (DEBUGGING) then
    if (n == 0 .or. npt < n + 2) then
        call errstop(srname, 'SIZE(XPT) is invalid')
    end if
    call verisize(vlag, npt)
    call verisize(xopt, n)
    call verisize(zmat, npt, int(npt - n - 1, kind(n)))
end if

rhosq = max(TENTH * delta, rho)**2
hdiag = -sum(zmat(:, 1:idz - 1)**2, dim=2) + sum(zmat(:, idz:npt - n - 1)**2, dim=2)
xdsq = sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1)
sigma = abs(beta * hdiag + vlag(1:npt)**2)
sigma = sigma * max(xdsq / rhosq, ONE)**3
if (ratio <= ZERO) then
    ! If the new F is not better than the current FVAL(KOPT), we set SIGMA(KOPT) = -1 to prevent
    ! KNEW from being KOPT.
    sigma(kopt) = -ONE
end if
if (maxval(sigma) > ONE .or. ratio > ZERO) then
    ! KNEW > 0 unless MAXVAL(SIGMA) <= 1 and RATIO <= ZERO. If RATIO > ZERO (the new F is less than
    ! the current FVAL(KOPT)), then we always set KNEW > 0, ensuring XNEW to be included in XPT.
    knew = int(maxloc(sigma, dim=1), kind(knew))
else
    knew = 0
end if
! It is tempting to take the function value into consideration when defining KNEW, for example,
! set KNEW so that FVAL(KNEW) = MAX(FVAL) as long as F(XNEW) < MAX(FVAL), unless there is a better
! choice. However, this is not a good idea, because the definition of KNEW should benefit the
! quality of the model that interpolates f at XPT. A set of points with low function values is not
! necessarily a good interplolation set. In contrast, a good interpolation set needs to include
! points with relatively high function values; otherwise, the interpolant will unlikely reflect the
! landscape of the function sufficiently.
end function setdrop_tr


function geostep(idz, knew, kopt, bmat, delbar, xpt, zmat) result(d)
! This subroutine finds a step D such that the geometry of the interplolation set is improved when
! XPT(:, KNEW) is changed to XOPT + D, where XOPT = XPT(:, KOPT)

! Generic modules
use consts_mod, only : RP, IK, ONE, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
use lina_mod, only : inprod

! Solver-specific module
use vlagbeta_mod, only : vlagbeta

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: delbar
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP) :: d(size(xpt, 1))     ! D(N)

! Local variables
integer(IK) :: n
integer(IK) :: npt
real(RP) :: alpha
real(RP) :: beta
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))
real(RP) :: xopt(size(xpt, 1))
real(RP) :: zknew(size(zmat, 2))
character(len=SRNLEN), parameter :: srname = 'GEOSTEP'


! Get and verify the sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

if (DEBUGGING) then
    if (n == 0 .or. npt < n + 2) then
        call errstop(srname, 'SIZE(XPT) is invalid')
    end if
    call verisize(bmat, n, npt + n)
    call verisize(zmat, npt, int(npt - n - 1, kind(n)))
    call verisize(vlag, npt + n)
end if

xopt = xpt(:, kopt)  ! Read XOPT.

d = biglag(idz, knew, delbar, bmat, xopt, xpt, zmat)

! ALPHA is the KNEW-th diagonal entry of H
zknew = zmat(knew, :)
zknew(1:idz - 1) = -zknew(1:idz - 1)
alpha = inprod(zmat(knew, :), zknew)

! Calculate VLAG and BETA for D.
call vlagbeta(idz, kopt, bmat, d, xpt, zmat, beta, vlag)

! If the cancellation in DENOM is unacceptable, then BIGDEN calculates an alternative model step D.
! VLAG and BETA for this D are calculated within BIGDEN.
if (abs(ONE + alpha * beta / vlag(knew)**2) <= 0.8_RP) then
    d = bigden(idz, knew, kopt, bmat, d, xpt, zmat)
end if

end function geostep


function biglag(idz, knew, delbar, bmat, x, xpt, zmat) result(d)
! BIGLAG calculates a D by approximately solving
!
! max |LFUNC(X + D)|, subject to ||D|| <= DELBAR,
!
! where LFUNC is the KNEW-th Lagrange function. See Setion 6 of the NEWUOA paper.

! Generic modules
use consts_mod, only : RP, IK, ONE, TWO, HALF, PI, ZERO, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
use lina_mod, only : Ax_plus_y, inprod, matprod

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
real(RP), intent(in) :: delbar
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: x(:)        ! X(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Output
real(RP) :: d(size(xpt, 1))       ! D(N)

! Local variables
integer(IK) :: i
integer(IK) :: isave
integer(IK) :: iterc
integer(IK) :: iu
integer(IK) :: n
integer(IK) :: npt
real(RP) :: angle
real(RP) :: cf(5)
real(RP) :: cth
real(RP) :: dd
real(RP) :: denom
real(RP) :: dhd
real(RP) :: gc(size(x))
real(RP) :: gd(size(x))
real(RP) :: gg
real(RP) :: hcol(size(xpt, 2))
real(RP) :: s(size(x))
real(RP) :: scaling
real(RP) :: sp
real(RP) :: ss
real(RP) :: step
real(RP) :: sth
real(RP) :: t
real(RP) :: tau
real(RP) :: taua
real(RP) :: taub
real(RP) :: taubeg
real(RP) :: taumax
real(RP) :: tauold
real(RP) :: unitang
real(RP) :: w(size(x))
real(RP) :: zknew(size(zmat, 2))
character(len=SRNLEN), parameter :: srname = 'BIGLAG'


! N is the number of variables.
! NPT is the number of interpolation equations.
! XPT contains the current interpolation points.
! BMAT provides the last N ROWs of H.
! ZMAT and IDZ give a factorization of the first NPT by NPT sub-matrix of H.
! KNEW is the index of the interpolation point to be dropped.
! DELBAR is the trust region bound for BIGLAG.
! D will be set to the step from X to the new point.
! HCOL, GC, GD, S and W will be used for working space.


! Get and verify the sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

if (DEBUGGING) then
    if (n == 0 .or. npt < n + 2) then
        call errstop(srname, 'SIZE(XPT) is invalid')
    end if
    call verisize(x, n)
    call verisize(bmat, n, npt + n)
    call verisize(zmat, npt, int(npt - n - 1, kind(n)))
end if

! Set HCOL to the leading NPT elements of the KNEW-th column of H.
zknew = zmat(knew, :)
zknew(1:idz - 1) = -zknew(1:idz - 1)
hcol = matprod(zmat, zknew)

! Set the unscaled initial direction D. Form the gradient of LFUNC at X, and multiply D by the
! Hessian of LFUNC.
d = xpt(:, knew) - x
dd = inprod(d, d)

gd = matprod(xpt, hcol * matprod(d, xpt))

!----------------------------------------------------------------!
gc = bmat(:, knew) + matprod(xpt, hcol * matprod(x, xpt)) !-!
!gc = Ax_plus_y(xpt, hcol * matprod(x, xpt), bmat(:, knew))
!----------------------------------------------------------------!

! Scale D and GD, with a sign change if required. Set S to another vector in the initial two
! dimensional subspace.
gg = inprod(gc, gc)
sp = inprod(d, gc)
dhd = inprod(d, gd)
scaling = delbar / sqrt(dd)
if (sp * dhd < ZERO) then
    scaling = -scaling
end if
t = ZERO
if (sp * sp > 0.99_RP * dd * gg) then
    t = ONE
end if
tau = scaling * (abs(sp) + HALF * scaling * abs(dhd))
if (gg * (delbar * delbar) < 1.0E-2_RP * tau * tau) then
    t = ONE
end if
d = scaling * d
gd = scaling * gd
s = gc + t * gd

! Begin the iteration by overwriting S with a vector that has the required length and direction,
! except that termination occurs if the given D and S are nearly parallel.
do iterc = 1, n
    dd = inprod(d, d)
    sp = inprod(d, s)
    ss = inprod(s, s)
    if (dd * ss - sp * sp <= 1.0E-8_RP * dd * ss) then
        exit
    end if
    denom = sqrt(dd * ss - sp * sp)
    s = (dd * s - sp * d) / denom

    w = matprod(xpt, hcol * matprod(s, xpt))

    ! Calculate the coefficients of the objective function on the circle, beginning with the
    ! multiplication of S by the second derivative matrix.
    cf(1) = inprod(s, w)
    cf(2) = inprod(d, gc)
    cf(3) = inprod(s, gc)
    cf(4) = inprod(d, gd)
    cf(5) = inprod(s, gd)
    cf(1) = HALF * cf(1)
    cf(4) = HALF * cf(4) - cf(1)

    ! Seek the value of the angle that maximizes |TAU|.
    taubeg = cf(1) + cf(2) + cf(4)
    taumax = taubeg
    tauold = taubeg
    isave = 0
    iu = 49
    unitang = (TWO * PI) / real(iu + 1, RP)

    do i = 1, iu
        angle = real(i, RP) * unitang
        cth = cos(angle)
        sth = sin(angle)
        tau = cf(1) + (cf(2) + cf(4) * cth) * cth + (cf(3) + cf(5) * cth) * sth
        if (abs(tau) > abs(taumax)) then
            taumax = tau
            isave = i
            taua = tauold
        else if (i == isave + 1) then
            taub = tau
        end if
        tauold = tau
    end do

    if (isave == 0) then
        taua = tau
    end if
    if (isave == iu) then
        taub = taubeg
    end if
    if (abs(taua - taub) > ZERO) then
        taua = taua - taumax
        taub = taub - taumax
        step = HALF * (taua - taub) / (taua + taub)
    else
        step = ZERO
    end if
    angle = unitang * (real(isave, RP) + step)

    ! Calculate the new D and GD. Then test for convergence.
    cth = cos(angle)
    sth = sin(angle)
    tau = cf(1) + (cf(2) + cf(4) * cth) * cth + (cf(3) + cf(5) * cth) * sth
    d = cth * d + sth * s
    gd = cth * gd + sth * w
    s = gc + gd
    if (abs(tau) <= 1.1_RP * abs(taubeg)) then
        exit
    end if
end do

end function biglag


function bigden(idz, knew, kopt, bmat, d0, xpt, zmat) result(d)
! BIGDEN calculates a D by approximately solving
!
! max |SIGMA(XOPT + D)|, subject to ||D|| <= DELBAR,
!
! where SIGMA is the denominator sigma in the updating formula (4.11)--(4.12) for H, which is the
! inverse of the coefficient matrix for the interplolation system (see (3.12)). Indeed, each column
! of H corresponds to a Lagrange basis function of the interpolation problem.  See Section 6 of the
! NEWUOA paper.
! N.B.:
! In Powell's code, BIGDEN calculates also the VLAG and BETA for the selected D. Here, to reduce the
! coupling of code, we return only D but computes VLAG and BETA outside by calling VLAGBETA. This
! does not change the mathematics, but the computed VLAG (BETA) will be numerically different due to
! rounding errors.

! Generic modules
use consts_mod, only : RP, IK, ONE, TWO, HALF, QUART, PI, ZERO, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
use lina_mod, only : Ax_plus_y, inprod, matprod

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: knew
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT+N)
real(RP), intent(in) :: d0(:)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Output
real(RP) :: d(size(xpt, 1))     ! D(N)

! Local variable
integer(IK) :: i
integer(IK) :: isave
integer(IK) :: iterc
integer(IK) :: iu
integer(IK) :: j
integer(IK) :: jc
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: nw
real(RP) :: alpha
real(RP) :: angle
real(RP) :: dd
real(RP) :: den(9)
real(RP) :: dena
real(RP) :: denb
real(RP) :: denex(9)
real(RP) :: denmax
real(RP) :: denold
real(RP) :: denom
real(RP) :: denomold
real(RP) :: densav
real(RP) :: ds
real(RP) :: dstemp(size(xpt, 2))
real(RP) :: dtest
real(RP) :: dxn
real(RP) :: hcol(size(xpt, 2))
real(RP) :: par(9)
real(RP) :: prod(size(xpt, 2) + size(xpt, 1), 5)
real(RP) :: s(size(xpt, 1))
real(RP) :: ss
real(RP) :: ssden
real(RP) :: sstemp(size(xpt, 2))
real(RP) :: step
real(RP) :: tau
real(RP) :: tempa
real(RP) :: tempb
real(RP) :: tempc
real(RP) :: unitang
real(RP) :: v(size(xpt, 2))
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))    ! VLAG(NPT + N)
real(RP) :: w(size(xpt, 2) + size(xpt, 1), 5)
real(RP) :: wz(size(zmat, 2))
real(RP) :: x(size(xpt, 1))
real(RP) :: xd
real(RP) :: xnew(size(xpt, 1))
real(RP) :: xnsq
real(RP) :: xptemp(size(xpt, 1), size(xpt, 2))
real(RP) :: xs
real(RP) :: xsq
real(RP) :: zknew(size(zmat, 2))
character(len=SRNLEN), parameter :: srname = 'BIGDEN'

! N is the number of variables.
! NPT is the number of interpolation equations.
! X is the best interpolation point so far.
! XPT contains the current interpolation points.
! BMAT provides the last N ROWs of H.
! ZMAT and IDZ give a factorization of the first NPT by NPT sub-matrix of H.
! NDIM is the second dimension of BMAT and has the value NPT + N.
! KOPT is the index of the optimal interpolation point.
! KNEW is the index of the interpolation point to be dropped.
! D will be set to the step from X to the new point, and on entry it should be the D that was
! calculated by the last call of BIGLAG. The length of the initial D provides a trust region bound
! on the final D.
!
! D is calculated in a way that should provide a denominator with a large modulus in the updating
! formula when the KNEW-th interpolation point is shifted to the new position X + D.


! Get and verify the sizes.
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

if (DEBUGGING) then
    if (n == 0 .or. npt < n + 2) then
        call errstop(srname, 'SIZE(XPT) is invalid')
    end if
    call verisize(x, n)
    call verisize(bmat, n, npt + n)
    call verisize(d0, n)
    call verisize(zmat, npt, int(npt - n - 1, kind(n)))
end if

x = xpt(:, kopt) ! For simplicity, use x to denote XOPT.

! Store the first NPT elements of the KNEW-th column of H in HCOL.
zknew = zmat(knew, :)
zknew(1:idz - 1) = -zknew(1:idz - 1)
hcol = matprod(zmat, zknew)
alpha = hcol(knew)

! The initial search direction D is taken from the last call of BIGLAG, and the initial S is set
! below, usually to the direction from X to X_KNEW, but a different direction to an interpolation
! point may be chosen, in order to prevent S from being nearly parallel to D.
d = d0
dd = inprod(d, d)
s = xpt(:, knew) - x
ds = inprod(d, s)
ss = inprod(s, s)
xsq = inprod(x, x)

if (.not. (ds * ds <= 0.99_RP * dd * ss)) then
    ! `.NOT. (A <= B)` differs from `A > B`.  The former holds iff A > B or {A, B} contains NaN.
    dtest = ds * ds / ss
    xptemp = xpt - spread(x, dim=2, ncopies=npt)
    !----------------------------------------------------------------!
    !---------!dstemp = matprod(d, xpt) - inprod(x, d) !-------------!
    dstemp = matprod(d, xptemp)
    !----------------------------------------------------------------!
    sstemp = sum((xptemp)**2, dim=1)

    dstemp(kopt) = TWO * ds + ONE
    sstemp(kopt) = ss
    k = int(minloc(dstemp * dstemp / sstemp, dim=1), kind(k))
    if ((.not. (dstemp(k) * dstemp(k) / sstemp(k) >= dtest)) .and. k /= kopt) then
        ! `.NOT. (A >= B)` differs from `A < B`.  The former holds iff A < B or {A, B} contains NaN.
        ! Although unlikely, if NaN occurs, it may happen that K = KOPT.
        s = xpt(:, k) - x
        ds = dstemp(k)
        ss = sstemp(k)
    end if
end if

ssden = dd * ss - ds * ds
densav = ZERO

! Begin the iteration by overwriting S with a vector that has the required length and direction.
do iterc = 1, n
    s = (ONE / sqrt(ssden)) * (dd * s - ds * d)
    xd = inprod(x, d)
    xs = inprod(x, s)

    ! Set the coefficients of the first two terms of BETA.
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
        w(k, 1) = QUART * (tempa * tempa + tempb * tempb)
        w(k, 2) = tempa * tempc
        w(k, 3) = tempb * tempc
        w(k, 4) = QUART * (tempa * tempa - tempb * tempb)
        w(k, 5) = HALF * tempa * tempb
    end do
    w(npt + 1:npt + n, 1:5) = ZERO
    w(npt + 1:npt + n, 2) = d
    w(npt + 1:npt + n, 3) = s

    ! Put the coefficents of THETA*WCHECK in PROD.
    do jc = 1, 5
        wz = matprod(w(1:npt, jc), zmat)
        wz(1:idz - 1) = -wz(1:idz - 1)
        prod(1:npt, jc) = matprod(zmat, wz)

        nw = npt
        if (jc == 2 .or. jc == 3) then
            prod(1:npt, jc) = prod(1:npt, jc) + matprod(w(npt + 1:npt + n, jc), bmat(:, 1:npt))
            nw = npt + n
        end if
        prod(npt + 1:npt + n, jc) = matprod(bmat(:, 1:nw), w(1:nw, jc))
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
    denom = denex(1) + denex(2) + denex(4) + denex(6) + denex(8)
    denold = denom
    denmax = denom
    isave = 0
    iu = 49
    unitang = (TWO * PI) / real(iu + 1, RP)
    par(1) = ONE
    do i = 1, iu
        angle = real(i, RP) * unitang
        par(2) = cos(angle)
        par(3) = sin(angle)
        do j = 4, 8, 2
            par(j) = par(2) * par(j - 2) - par(3) * par(j - 1)
            par(j + 1) = par(2) * par(j - 1) + par(3) * par(j - 2)
        end do
        denomold = denom
        denom = inprod(denex(1:9), par(1:9))
        if (abs(denom) > abs(denmax)) then
            denmax = denom
            isave = i
            dena = denomold
        else if (i == isave + 1) then
            denb = denom
        end if
    end do
    if (isave == 0) then
        dena = denom
    end if
    if (isave == iu) then
        denb = denold
    end if
    if (abs(dena - denb) > 0) then
        dena = dena - denmax
        denb = denb - denmax
        step = HALF * (dena - denb) / (dena + denb)
    else
        step = ZERO
    end if
    angle = unitang * (real(isave, RP) + step)

    ! Calculate the new parameters of the denominator, the new VLAG vector and the new D. Then test
    ! for convergence.
    par(2) = cos(angle)
    par(3) = sin(angle)
    do j = 4, 8, 2
        par(j) = par(2) * par(j - 2) - par(3) * par(j - 1)
        par(j + 1) = par(2) * par(j - 1) + par(3) * par(j - 2)
    end do

    !beta = inprod(den(1:9), par(1:9))  ! Not needed since we do not return BETA.

    denmax = inprod(denex(1:9), par(1:9))

    vlag = matprod(prod(:, 1:5), par(1:5))

    tau = vlag(knew)

    d = par(2) * d + par(3) * s
    dd = inprod(d, d)
    xnew = x + d
    dxn = inprod(d, xnew)
    xnsq = inprod(xnew, xnew)

    if (iterc > 1) then
        densav = max(densav, denold)
    end if
    if (abs(denmax) <= 1.1_RP * abs(densav)) then
        exit
    end if
    densav = denmax

    ! Set S to HALF the gradient of the denominator with respect to D.
    s = tau * bmat(:, knew) + alpha * (dxn * x + xnsq * d - vlag(npt + 1:npt + n))
    v = matprod(xnew, xpt)
    v = (tau * hcol - alpha * vlag(1:npt)) * v
    !------------------------------------------------------!
    s = s + matprod(xpt, v) !-------------------!
    !s = Ax_plus_y(xpt, v, s)
    !------------------------------------------------------!

    ss = inprod(s, s)
    ds = inprod(d, s)
    ssden = dd * ss - ds * ds
    if (ssden < 1.0E-8_RP * dd * ss) then
        exit
    end if
end do

!vlag(kopt) = vlag(kopt) + ONE  ! Note needed since we do not return VLAG.

end function bigden


end module geometry_mod
