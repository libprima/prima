!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code.
!
! Started: July 2020
!
! Last Modified: Saturday, January 01, 2022 PM12:41:42
!--------------------------------------------------------------------------------------------------!

!!!!!! THE MODULE THAT IMPLEMENTS CALFUN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module calfun_mod

implicit none
private
public :: calfun

contains

subroutine calfun(x, f)
! The Chebyquad test problem (Fletcher, 1965)
implicit none

real(kind(0.0D0)), intent(in) :: x(:)
real(kind(0.0D0)), intent(out) :: f

integer :: i, n
real(kind(0.0D0)) :: y(size(x) + 1, size(x) + 1), tmp

n = size(x)

y(1:n, 1) = 1.0D0
y(1:n, 2) = 2.0D0 * x - 1.0D0
do i = 2, n
    y(1:n, i + 1) = 2.0D0 * y(1:n, 2) * y(1:n, i) - y(1:n, i - 1)
end do

f = 0.0D0
do i = 1, n + 1
    tmp = sum(y(1:n, i)) / real(n, kind(0.0D0))
    if (mod(i, 2) /= 0) then
        tmp = tmp + 1.0D0 / real(i * i - 2 * i, kind(0.0D0))
    end if
    f = f + tmp * tmp
end do
end subroutine calfun

end module calfun_mod


!!!!!!! THE MAIN PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program newuoa_exmp

! The following line makes the solver available.
!--------------------------------------------------------------------------------------------------!
use newuoa_mod, only : newuoa
!--------------------------------------------------------------------------------------------------!

! The following line specifies which module provides CALFUN. If CALFUN is given by an external
! subroutine instead of a module, remove this line and uncomment the "external calfun" line below.
!--------------------------------------------------------------------------------------------------!
use calfun_mod, only : calfun
!--------------------------------------------------------------------------------------------------!

implicit none

integer :: i
integer, parameter :: n = 6
real(kind(0.0D0)) :: x(n)
real(kind(0.0D0)) :: f

! If CALFUN is an external subroutine, then remove the line of  "use calfun_mod, only : calfun", and
! uncomment the following line.
!--------------------------------------------------------------------------------------------------!
!external calfun
!--------------------------------------------------------------------------------------------------!

! The following lines illustrates how to call the solver to solve the Chebyquad problem.
!----------------------------------------------------------------------------------------------!
x = [(real(i, kind(0.0D0)) / real(n + 1, kind(0.0D0)), i=1, n)]
call newuoa(calfun, x, f)  ! This call will not print anything.
!----------------------------------------------------------------------------------------------!
! In addition to the compulsory argument, and M, the following illustration specifies also
! RHOBEG and IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN,
! etc.) will take their default values coded in the solver.
!----------------------------------------------------------------------------------------------!
call newuoa(calfun, x, f, rhobeg=0.2D0 * x(1), iprint=1)
!----------------------------------------------------------------------------------------------!

end program newuoa_exmp
