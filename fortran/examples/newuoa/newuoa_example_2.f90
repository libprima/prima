!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: July 2020
!
! Last Modified: Tuesday, March 05, 2024 AM09:42:53
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALFUN -------------------------------!
module calfun_mod

implicit none
private
public :: RP, IK, calfun
integer, parameter :: RP = kind(0.0D0)
integer, parameter :: IK = kind(0)
! N.B.: We assume that PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0
! (default kind). Revise RP and IK if this is not the case.

contains

subroutine calfun(x, f)
! The Chebyquad test problem (Fletcher, 1965)
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

integer :: i, n
real(RP) :: y(size(x) + 1, size(x) + 1), tmp

n = size(x)

y(1:n, 1) = 1.0_RP
y(1:n, 2) = 2.0_RP * x - 1.0_RP
do i = 2, n
    y(1:n, i + 1) = 2.0_RP * y(1:n, 2) * y(1:n, i) - y(1:n, i - 1)
end do

f = 0.0_RP
do i = 1, n + 1
    tmp = sum(y(1:n, i)) / real(n, RP)
    if (modulo(i, 2) /= 0) then
        tmp = tmp + 1.0_RP / real(i * i - 2 * i, RP)
    end if
    f = f + tmp * tmp
end do
end subroutine calfun

end module calfun_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program newuoa_exmp

! The following line makes the solver available.
use newuoa_mod, only : newuoa

! The following line specifies which module provides CALFUN.
use calfun_mod, only : RP, IK, calfun

implicit none

integer, parameter :: n = 6
integer :: i, nf, info
real(RP) :: f, x(n)

! The following lines illustrates how to call the solver to solve the Chebyquad problem.
x = [(real(i, RP) / real(n + 1, RP), i=1, n)]  ! Define the starting point.
call newuoa(calfun, x, f)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = [(real(i, RP) / real(n + 1, RP), i=1, n)]  ! Define the starting point.
call newuoa(calfun, x, f, rhobeg=0.2_RP * x(1), iprint=1_IK, nf=nf, info=info)

end program newuoa_exmp
