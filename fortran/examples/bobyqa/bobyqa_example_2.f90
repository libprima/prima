!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: July 2020
!
! Last Modified: Tuesday, May 30, 2023 PM06:26:48
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALFUN -------------------------------!
module calfun_mod

implicit none
private
public :: RP, calfun
integer, parameter :: RP = kind(0.0D0)

contains

subroutine calfun(x, f)
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
integer :: i
integer :: j
integer :: n
real(RP) :: temp

n = int(size(x), kind(n))

f = 0.0_RP
do i = 4, n, 2
    do j = 2, i - 2, 2
        temp = (x(i - 1) - x(j - 1))**2 + (x(i) - x(j))**2
        temp = max(temp, 1.0D-6)
        f = f + 1.0_RP / sqrt(temp)
    end do
end do

end subroutine calfun

end module calfun_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program bobyqa_exmp

! The following line makes the solver available.
use bobyqa_mod, only : bobyqa

! The following line specifies which module provides CALFUN.
use calfun_mod, only : RP, calfun

implicit none

integer, parameter :: n = 20
integer :: i, nf, info
real(RP) :: f, x(n), x0(n), lb(n), ub(n), angle

! Define the starting point.
do i = 1, n / 2
    angle = real(i, RP) * 8.0_RP * atan(1.0_RP) / real(n / 2, RP)
    x0(2 * i - 1) = cos(angle)
    x0(2 * i) = sin(angle)
end do
lb = -1.0_RP
ub = 1.0_RP

! The following lines illustrates how to call the solver.
x = x0
call bobyqa(calfun, x, f, lb, ub)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = x0
call bobyqa(calfun, x, f, lb, ub, rhobeg=1.0D-1, iprint=1, nf=nf, info=info)


end program bobyqa_exmp
