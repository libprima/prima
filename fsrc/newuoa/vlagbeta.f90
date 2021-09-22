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
! Last Modified: Wednesday, September 22, 2021 AM11:35:46
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: vlagbeta


contains


subroutine vlagbeta(idz, kopt, bmat, d, xpt, zmat, beta, vlag)
!--------------------------------------------------------------------------------------------------!
! VLAGBETA calculates VLAG = Hw and BETA for a given step D. See (4.11), (4.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: linalg_mod, only : Ax_plus_y, xA_plus_y, xpy_dot_z, inprod, matprod

implicit none

! Inputs
integer(IK), intent(in) :: idz
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(in) :: d(:)    ! D(N)
real(RP), intent(in) :: xpt(:, :)   ! XPT(N, NPT)
real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP), intent(out) :: beta
real(RP), intent(out) :: vlag(:)    ! VLAG(NPT + N)

! Local variables
character(len=*), parameter :: srname = 'VLAGBETA'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: bw(size(bmat, 1))
real(RP) :: bwvd
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
    call assert(size(vlag) == n + npt, 'SIZE(VLAG) == N + NPT', srname)
end if

!====================!
! Calculation starts !
!====================!

xopt = xpt(:, kopt)  ! Read XOPT.

! This is one of the two places where WCHECK is calculated, the other being in BIGDEN but removed
! WCHECK contains the first NPT entries of (w-v) for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, and also hat{w} in eq(6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
wcheck = matprod(d, xpt)
wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

vlag(1:npt) = matprod(d, bmat(:, 1:npt))

wz = matprod(wcheck, zmat)
wzsav = wz
wz(1:idz - 1) = -wz(1:idz - 1)
beta = -inprod(wzsav, wz)
!----------------------------------------------------------------------!
!-----!vlag(1 : npt) = vlag(1 : npt) + matprod(zmat, wz) !-------------!
vlag(1:npt) = Ax_plus_y(zmat, wz, vlag(1:npt))
!----------------------------------------------------------------------!

vlag(kopt) = vlag(kopt) + ONE  ! The calculation of VLAG(1:NPT) finishes.

bw = matprod(bmat(:, 1:npt), wcheck)
!----------------------------------------------------------------------!
!vlag(npt + 1 : npt + n) = bw + matprod(d, bmat(:, npt + 1 : npt + n)) !
vlag(npt + 1:npt + n) = xA_plus_y(d, bmat(:, npt + 1:npt + n), bw)
! The calculation of VLAG finishes.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!-----!bwvd = inprod(bw + vlag(npt + 1 : npt + n), d) !----------------!
bwvd = xpy_dot_z(bw, vlag(npt + 1:npt + n), d)
!----------------------------------------------------------------------!

dx = inprod(d, xopt)
dsq = inprod(d, d)
xoptsq = inprod(xopt, xopt)

! The final value of BETA is calculated as follows.
beta = dx**2 + dsq * (xoptsq + TWO * dx + HALF * dsq) + beta - bwvd

!====================!
!  Calculation ends  !
!====================!

end subroutine vlagbeta


end module vlagbeta_mod
