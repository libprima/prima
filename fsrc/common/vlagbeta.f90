module vlagbeta_mod
!--------------------------------------------------------------------------------------------------!
! This module contains a subroutine that calculates VLAG and BETA for a given step D. Both VLAG and
! BETA are critical for the updating procedure of H, which is detailed formula (4.11) of the NEWUOA
! paper. See (4.12) for the definition of BETA, and VLAG is indeed Hw; (4.25)--(4.26) formulate the
! actual calculating scheme of VLAG and BETA.
!
! N.B.:
! 1. In languages like MATLAB/Python/Julia/R, CALVLAG and CALBETA should be implemented into one
! single function, as they share most of the calculation. We separate them in Fortran because
! Fortran functions can only have one output.
! 2. Given any D and t in {1, ..., NPT}, VLAG(t) = e_t^T Hw = L_t(XBASE + XOPT + D), where L_t is
! the t-th Lagrange basis function corresponding to the interpolation set that defines H. See
! (6.3)--(6.4) of the NEWUOA paper for details. As a consequence, SUM(VLAG(1:NPT)) = 1 in theory.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: July 2020
!
! Last Modified: Saturday, April 09, 2022 PM05:24:50
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: calvlag, calbeta


contains


function calvlag(kopt, bmat, d, xpt, zmat, idz) result(vlag)
!--------------------------------------------------------------------------------------------------!
! This function calculates VLAG = Hw for a given step D. See (4.25) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: VLAG(NPT+N), WCHECK(NPT), XOPT(N)
! Size of local arrays: REAL(RP)*(2*NPT+2*N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod, issymmetric
use, non_intrinsic :: powalg_mod, only : omega_mul

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
integer(IK), intent(in), optional :: idz  ! Absent in BOBYQA, being equivalent to IDZ = 1

! Outputs
real(RP) :: vlag(size(bmat, 2))    ! VLAG(NPT + N)

! Local variables
character(len=*), parameter :: srname = 'CALVLAG'
integer(IK) :: idz_loc
integer(IK) :: n
integer(IK) :: npt
real(RP) :: wcheck(size(zmat, 1))
real(RP) :: xopt(size(xpt, 1))
real(RP) :: tol  ! For debugging only

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Read IDZ, which is not present in BOBYQA, being equivalent to IDZ = 1.
idz_loc = 1_IK
if (present(idz)) then
    idz_loc = idz
end if

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz_loc >= 1 .and. idz_loc <= size(zmat, 2) + 1, '1 <= ID <= SIZE(ZMAT, 2) + 1', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

! WCHECK contains the first NPT entries of (w-v) for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, and also hat{w} in eq(6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
wcheck = matprod(d, xpt)
wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))
! Assume that the |D| ~ DELTA, |XPT| ~ |XOPT|, and DELTA < |XOPT|. Then WCHECK is in the order of
! DELTA*|XOPT|^3, which is can quickly become tiny.

!vlag(1:npt) = matprod(d, bmat(:, 1:npt))
!vlag(1:npt) = vlag(1:npt) + omega_mul(idz_loc, zmat, wcheck)

vlag(1:npt) = omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt))

!------------------------------------------------------!
vlag(npt + 1:npt + n) = matprod(bmat, [wcheck, d])
!vlag(npt + 1:npt + n) = matprod(bmat(:, 1:npt), wcheck) + matprod(bmat(:, npt + 1:npt + n), d)
!------------------------------------------------------!

vlag(kopt) = vlag(kopt) + ONE

!====================!
!  Calculation ends  !
!====================!

write (16, *) abs(sum(vlag(1:npt)) - 1.0_RP)

! Postconditions
if (DEBUGGING) then
    call assert(size(vlag) == npt + n, 'SIZE(VLAG) == NPT + N', srname)
    tol = max(1.0E-10_RP, min(1.0E-2_RP, 1.0E6_RP * EPS * real(npt + n, RP)))
    call assert(abs(sum(vlag(1:npt)) - ONE) <= tol, 'SUM(VLAG(1:NPT)) == 1', srname)
end if

end function calvlag


function calbeta(kopt, bmat, d, xpt, zmat, idz) result(beta)
!--------------------------------------------------------------------------------------------------!
! This function calculates BETA for a given step D. See (4.12) and (4.26) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
! List of local arrays (including function-output arrays; likely to be stored on the stack):
! REAL(RP) :: BW(N), BD(N), WCHECK(NPT), XOPT(N)
! Size of local arrays: REAL(RP)*(NPT+3*N)
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod, issymmetric
use, non_intrinsic :: powalg_mod, only : omega_inprod, omega_mul

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
integer(IK), intent(in), optional :: idz  ! Absent in BOBYQA, being equivalent to IDZ = 1

! Outputs
real(RP) :: beta

! Local variables
character(len=*), parameter :: srname = 'CALBETA'
integer(IK) :: idz_loc
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bw(size(bmat, 1))
real(RP) :: bd(size(bmat, 1))
real(RP) :: bsum
real(RP) :: dsq
real(RP) :: dxopt
real(RP) :: wcheck(size(zmat, 1))
real(RP) :: xopt(size(xpt, 1))
real(RP) :: xoptsq
real(RP) :: x(size(xpt, 1))
real(RP) :: wv(size(xpt, 1) + size(xpt, 2)), Hwv(size(xpt, 1) + size(xpt, 2)), wvHwv, beta1

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Read IDZ, which is absent from BOBYQA, being equivalent to IDZ = 1.
idz_loc = 1_IK
if (present(idz)) then
    idz_loc = idz
end if

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2', srname)
    call assert(idz_loc >= 1 .and. idz_loc <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

dxopt = inprod(d, xopt)
dsq = inprod(d, d)
xoptsq = inprod(xopt, xopt)

! WCHECK contains the first NPT entries of (w-v) for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, and also hat{w} in eq(6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
wcheck = matprod(d, xpt)
wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

!bw = matprod(bmat(:, 1:npt), wcheck)
!bd = matprod(bmat(:, npt + 1:npt + n), d)
!bsum = sum(bd * d + bw * d + bw * d)  ! VERSION 1

!bsum = inprod(bd + TWO * bw, d)  ! VERSION 2
!bsum = inprod(bd + bw + bw, d)  ! VERSION 5
!bsum = inprod(bw + bw + bd, d)  ! VERSION 5

!bw = matprod(bmat, [TWO * wcheck, d]); bsum = inprod(bw, d)  ! VERSION 3

!beta = dxopt**2 + dsq * (xoptsq + TWO * dxopt + HALF * dsq) - omega_inprod(idz_loc, zmat, wcheck, wcheck) - bsum  ! VERSIONa
!beta = dxopt**2 + dsq * (xoptsq + dxopt + dxopt + HALF * dsq) - omega_inprod(idz_loc, zmat, wcheck, wcheck) - bsum  ! VERSIONb


!beta = inprod(wcheck, omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt))) + inprod(d, matprod(bmat, [wcheck, d]))
!beta = inprod(wcheck, omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt))) + inprod(d, matprod(bmat, [wcheck, d]))


!----------------------------------------------------------------------------------------!
wv = [wcheck, d]
Hwv = [omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt)), matprod(bmat, wv)]
wvHwv = inprod(wv, Hwv)
!----------------------------------------------------------------------------------------!

!Hwv = [omega_mul(idz_loc, zmat, wcheck) + matprod(d, bmat(:, 1:npt)), &
!    & matprod(bmat(:, 1:npt), wcheck) + matprod(bmat(:, npt + 1:npt + n), d)]
!wvHwv = inprod(wcheck, Hwv(1:npt)) + inprod(d, Hwv(npt + 1:npt + n))

!x = xopt + d
!beta = HALF * (inprod(x, x)**2 + inprod(xopt, xopt)**2) - inprod(x, xopt)**2 - wvHwv

beta = dxopt**2 + dsq * (xoptsq + dxopt + dxopt + HALF * dsq) - wvHwv

!write (*, *) abs(beta1 - beta) / max(1.0_RP, abs(beta))

!====================!
!  Calculation ends  !
!====================!

end function calbeta

end module vlagbeta_mod
