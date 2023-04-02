!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: July 2020
!
! Last Modified: Tuesday, February 28, 2023 PM08:46:50
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALFUN -------------------------------!
module calfun_mod

implicit none
private
public :: calfun

contains

subroutine calfun(x, f)
implicit none

! Inputs
real(kind(0.0D0)), intent(in) :: x(:)

! Outputs
real(kind(0.0D0)), intent(out) :: f

! Local variables
integer :: i
integer :: j
integer :: n
real(kind(0.0D0)) :: temp

n = int(size(x), kind(n))

f = 0.0D0
do i = 4, n, 2
    do j = 2, i - 2, 2
        temp = (x(i - 1) - x(j - 1))**2 + (x(i) - x(j))**2
        temp = max(temp, 1.0D-6)
        f = f + 1.0D0 / sqrt(temp)
    end do
end do

end subroutine calfun

end module calfun_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program bobyqa_exmp

! The following line makes the solver available.
use bobyqa_mod, only : bobyqa

! The following line specifies which module provides CALFUN.
use calfun_mod, only : calfun

implicit none

integer :: i
integer, parameter :: n = 20
integer :: nf
real(kind(0.0D0)) :: x(n)
real(kind(0.0D0)) :: f
real(kind(0.0D0)) :: x0(n), lb(n), ub(n), angle

! Define the starting point.
do i = 1, n / 2
    angle = real(i, kind(0.0D0)) * 8.0D0 * atan(1.0D0) / real(n / 2, kind(0.0D0))
    x0(2 * i - 1) = cos(angle)
    x0(2 * i) = sin(angle)
end do
lb = -1.0D0
ub = 1.0D0

! The following lines illustrates how to call the solver.
x = x0
call bobyqa(calfun, x, f, lb, ub)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = x0
call bobyqa(calfun, x, f, lb, ub, rhobeg=1.0D-1, iprint=2, nf=nf)


end program bobyqa_exmp
