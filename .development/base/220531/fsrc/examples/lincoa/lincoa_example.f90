!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code.
!
! Started: July 2020
!
! Last Modified: Monday, March 07, 2022 AM12:55:12
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALFUN -------------------------------!
module tetrahedron_mod

implicit none
private
public :: calfun, setup

contains

subroutine calfun(x, f)
implicit none

! Inputs
real(kind(0.0D0)), intent(in) :: x(:)

! Outputs
real(kind(0.0D0)), intent(out) :: f

! Local variables
integer :: j
integer, parameter :: np = 50
real(kind(0.0D0)) :: theta(np), xp(np), yp(np), zp(np), xs, ys, zs, ss
real(kind(0.0D0)) :: v12, v13, v14, v23, v24, v34, del1, del2, del3, del4, tmp

theta = 4.0D0 * atan(1.0D0)*[(real(j - 1, kind(0.0D0)) / real(np - 1, kind(0.0D0)), j=1, np)]
xp = cos(theta) * cos(2.0D0 * theta)
yp = sin(theta) * cos(2.0D0 * theta)
zp = sin(2.0D0 * theta)
xp = xp - sum(xp) / real(np, kind(0.0D0))
yp = yp - sum(yp) / real(np, kind(0.0D0))
zp = zp - sum(zp) / real(np, kind(0.0D0))
xs = minval([0.0D0, xp])
ys = minval([0.0D0, yp])
zs = minval([0.0D0, zp])
ss = maxval([0.0D0, xp + yp + zp])
f = (ss - xs - ys - zs)**3 / 6.0D0

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
f = min(tmp / 6.0D0, f)

end subroutine calfun


subroutine setup(x0, A, b)
implicit none

real(kind(0.0D0)), intent(out) :: x0(:), A(:, :), b(:)

integer, parameter :: np = 50
integer :: i, j
real(kind(0.0D0)) :: theta(np), xp(np), yp(np), zp(np), xs, ys, zs, ss

theta = 4.0D0 * atan(1.0D0)*[(real(j - 1, kind(0.0D0)) / real(np - 1, kind(0.0D0)), j=1, np)]
xp = cos(theta) * cos(2.0D0 * theta)
yp = sin(theta) * cos(2.0D0 * theta)
zp = sin(2.0D0 * theta)
xp = xp - sum(xp) / real(np, kind(0.0D0))
yp = yp - sum(yp) / real(np, kind(0.0D0))
zp = zp - sum(zp) / real(np, kind(0.0D0))
xs = minval([0.0D0, xp])
ys = minval([0.0D0, yp])
zs = minval([0.0D0, zp])
ss = maxval([0.0D0, xp + yp + zp])

x0(2:8) = 0.0D0
x0(1) = 1.0D0 / xs
x0(5) = 1.0D0 / ys
x0(9) = 1.0D0 / zs
x0(10:12) = 1.0D0 / ss

A = 0.0D0
do j = 1, np
    do i = 1, 4
        A(3 * i - 2, 4 * j + i - 4) = xp(j)
        A(3 * i - 1, 4 * j + i - 4) = yp(j)
        A(3 * i, 4 * j + i - 4) = zp(j)
    end do
end do

b = 1.0D0

end subroutine setup

end module tetrahedron_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program lincoa_exmp

!--------------------------------------------------------------------------------------------------!
! The following line makes the solver available.
use lincoa_mod, only : lincoa
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! The following line specifies which module provides CALFUN. If CALFUN is given by an external
! subroutine instead of a module, remove this line and uncomment the "external calfun" line below.
use tetrahedron_mod, only : calfun
!--------------------------------------------------------------------------------------------------!

use tetrahedron_mod, only : setup

implicit none

integer, parameter :: n = 12
integer, parameter :: np = 50
integer :: nf
real(kind(0.0D0)) :: x(n), f
real(kind(0.0D0)) :: x0(n), A(n, 4 * np), b(4 * np)

!--------------------------------------------------------------------------------------------------!
! If CALFUN is an external subroutine, then remove the line of  "use calfun_mod, only : calfun", and
! uncomment the following line.
!!external calfun
!--------------------------------------------------------------------------------------------------!

! Set up X0 (starting point), A, and b.
call setup(x0, A, b)

! The following lines illustrates how to call the solver.
x = x0
call lincoa(calfun, x, f, A=A, b=b, nf=nf)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = x0
call lincoa(calfun, x, f, A=A, b=b, rhobeg=0.5D0, iprint=1, nf=nf)
write (*, *) 'nf = ', nf
write (*, *) 'f = ', f
write (*, *) 'x = ', x

end program lincoa_exmp
