module vlagbeta_mod
!--------------------------------------------------------------------------------------------------!
! This module contains a subroutine that calculates VLAG and BETA for a given step D. Both VLAG and
! BETA are critical for the updating procedure of H, which is detailed formula (4.11) of the NEWUOA
! paper. See (4.12) for the definition of BETA, and VLAG is indeed Hw.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Wednesday, November 03, 2021 PM09:54:25
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: calvlag, calbeta


contains


function calvlag(idz, kopt, bmat, d, xpt, zmat) result(vlag)
!--------------------------------------------------------------------------------------------------!
! This function calculates VLAG = Hw for a given step D. See (4.11) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : matprod

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP) :: vlag(size(bmat, 2))    ! VLAG(NPT + N)

! Local variables
character(len=*), parameter :: srname = 'CALVLAG'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: wcheck(size(zmat, 1))
real(RP) :: wz(size(zmat, 2))
real(RP) :: xopt(size(xpt, 1))

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
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

vlag(1:npt) = matprod(d, bmat(:, 1:npt))
wz = matprod(wcheck, zmat)
wz(1:idz - 1) = -wz(1:idz - 1)
vlag(1:npt) = vlag(1:npt) + matprod(zmat, wz)
vlag(kopt) = vlag(kopt) + ONE  ! The calculation of VLAG(1:NPT) finishes.
vlag(npt + 1:npt + n) = matprod(bmat, [wcheck, d]) ! The calculation of VLAG finishes.

!====================!
!  Calculation ends  !
!====================!

end function calvlag


function calbeta(idz, kopt, bmat, d, xpt, zmat) result(beta)
!--------------------------------------------------------------------------------------------------!
! This function calculates BETA for a given step D. See (4.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : inprod, matprod

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP) :: beta

! Local variables
character(len=*), parameter :: srname = 'CALBETA'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bw(size(bmat, 1))
real(RP) :: bd(size(bmat, 1))
real(RP) :: bsum
real(RP) :: dsq
real(RP) :: dx
real(RP) :: wcheck(size(zmat, 1))
real(RP) :: wz(size(zmat, 2))
real(RP) :: wzsav(size(wz))
real(RP) :: xopt(size(xpt, 1))
real(RP) :: xoptsq

! Sizes
n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N + 2', srname)
    call assert(idz >= 1 .and. idz <= npt - n, '1 <= IDZ <= NPT-N', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == N, D is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, &
        & 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

dx = inprod(d, xopt)
dsq = inprod(d, d)
xoptsq = inprod(xopt, xopt)

! WCHECK contains the first NPT entries of (w-v) for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, and also hat{w} in eq(6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
wcheck = matprod(d, xpt)
wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

wz = matprod(wcheck, zmat)
wzsav = wz
wz(1:idz - 1) = -wz(1:idz - 1)

bw = matprod(bmat(:, 1:npt), wcheck)
bd = matprod(bmat(:, npt + 1:npt + n), d)
bsum = sum(bd * d + bw * d + bw * d)

beta = dx**2 + dsq * (xoptsq + TWO * dx + HALF * dsq) - inprod(wzsav, wz) - bsum

!====================!
!  Calculation ends  !
!====================!

end function calbeta

end module vlagbeta_mod
