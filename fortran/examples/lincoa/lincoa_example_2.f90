!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: July 2020
!
! Last Modified: Friday, March 15, 2024 PM02:06:02
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALFUN -------------------------------!
module tetrahedron_mod

implicit none
private
public :: RP, IK, calfun, setup
integer, parameter :: RP = kind(0.0D0)
integer, parameter :: IK = kind(0)
! N.B.: We assume that PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0
! (default kind). Revise RP and IK if this is not the case.

contains

subroutine calfun(x, f)
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
integer :: j
integer, parameter :: np = 50
real(RP) :: theta(np), xp(np), yp(np), zp(np), xs, ys, zs, ss
real(RP) :: v12, v13, v14, v23, v24, v34, del1, del2, del3, del4, tmp

theta = 4.0_RP * atan(1.0_RP)*[(real(j - 1, RP) / real(np - 1, RP), j=1, np)]
xp = cos(theta) * cos(2.0_RP * theta)
yp = sin(theta) * cos(2.0_RP * theta)
zp = sin(2.0_RP * theta)
xp = xp - sum(xp) / real(np, RP)
yp = yp - sum(yp) / real(np, RP)
zp = zp - sum(zp) / real(np, RP)
xs = minval([0.0_RP, xp])
ys = minval([0.0_RP, yp])
zs = minval([0.0_RP, zp])
ss = maxval([0.0_RP, xp + yp + zp])
f = (ss - xs - ys - zs)**3 / 6.0_RP

v12 = x(1) * x(5) - x(4) * x(2)
v13 = x(1) * x(8) - x(7) * x(2)
v14 = x(1) * x(11) - x(10) * x(2)
v23 = x(4) * x(8) - x(7) * x(5)
v24 = x(4) * x(11) - x(10) * x(5)
v34 = x(7) * x(11) - x(10) * x(8)
del1 = v23 * x(12) - v24 * x(9) + v34 * x(6)
if (del1 <= 0) return
del2 = -v34 * x(3) - v13 * x(12) + v14 * x(9)
if (del2 <= 0) return
del3 = -v14 * x(6) + v24 * x(3) + v12 * x(12)
if (del3 <= 0) return
del4 = -v12 * x(9) + v13 * x(6) - v23 * x(3)
if (del4 <= 0) return
tmp = (del1 + del2 + del3 + del4)**3 / (del1 * del2 * del3 * del4)
f = min(tmp / 6.0_RP, f)

end subroutine calfun


subroutine setup(x0, Aineq, bineq)
implicit none

real(RP), intent(out) :: x0(:), Aineq(:, :), bineq(:)

integer, parameter :: np = 50
integer :: i, j
real(RP) :: theta(np), xp(np), yp(np), zp(np), xs, ys, zs, ss

theta = 4.0_RP * atan(1.0_RP)*[(real(j - 1, RP) / real(np - 1, RP), j=1, np)]
xp = cos(theta) * cos(2.0_RP * theta)
yp = sin(theta) * cos(2.0_RP * theta)
zp = sin(2.0_RP * theta)
xp = xp - sum(xp) / real(np, RP)
yp = yp - sum(yp) / real(np, RP)
zp = zp - sum(zp) / real(np, RP)
xs = minval([0.0_RP, xp])
ys = minval([0.0_RP, yp])
zs = minval([0.0_RP, zp])
ss = maxval([0.0_RP, xp + yp + zp])

x0(2:8) = 0.0_RP
x0(1) = 1.0_RP / xs
x0(5) = 1.0_RP / ys
x0(9) = 1.0_RP / zs
x0(10:12) = 1.0_RP / ss

Aineq = 0.0_RP
do j = 1, np
    do i = 1, 4
        Aineq(4 * j + i - 4, 3 * i - 2) = xp(j)
        Aineq(4 * j + i - 4, 3 * i - 1) = yp(j)
        Aineq(4 * j + i - 4, 3 * i) = zp(j)
    end do
end do

bineq = 1.0_RP

end subroutine setup

end module tetrahedron_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program lincoa_exmp

! The following line makes the solver available.
use lincoa_mod, only : lincoa

! The following line specifies which module provides CALFUN.
use tetrahedron_mod, only : RP, IK, calfun, setup

implicit none

integer, parameter :: n = 12, np = 50
integer :: nf, info
real(RP) :: f, cstrv, x(n), x0(n), Aineq(4 * np, n), bineq(4 * np)

! Set up X0 (starting point), Aineq, and bineq.
call setup(x0, Aineq, bineq)

! The following lines illustrates how to call the solver.
x = x0
call lincoa(calfun, x, f, cstrv, Aineq, bineq)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = x0
call lincoa(calfun, x, f, cstrv, Aineq, bineq, rhobeg=1.0_RP, iprint=1_IK, nf=nf, info=info)

end program lincoa_exmp
