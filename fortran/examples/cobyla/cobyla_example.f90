!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: July 2020
!
! Last Modified: Friday, August 04, 2023 PM06:01:56
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALCFC -------------------------------!
module calcfc_mod

implicit none
private
public :: RP, calcfc_chebyquad, calcfc_hexagon
integer, parameter :: RP = kind(0.0D0)

contains

subroutine calcfc_chebyquad(x, f, constr)
! The Chebyquad test problem (Fletcher, 1965)
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)

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
constr = 0.0_RP  ! Without this line, compilers may complain that CONSTR is not set.
end subroutine calcfc_chebyquad

subroutine calcfc_hexagon(x, f, constr)
! Test problem 10 (Hexagon area) in Powell's original package.
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_HEXAGON'
real(kind(1.0_RP)), intent(in) :: x(:)
real(kind(1.0_RP)), intent(out) :: constr(:)
real(kind(1.0_RP)), intent(out) :: f

call assert(size(x) == 9 .and. size(constr) == 14, 'SIZE(X) == 9, SIZE(CONSTR) == 14', srname)

f = -0.5_RP * (x(1) * x(4) - x(2) * x(3) + x(3) * x(9) - x(5) * x(9) + x(5) * x(8) - x(6) * x(7))
constr(1) = -1.0_RP + x(3)**2 + x(4)**2
constr(2) = -1.0_RP + x(9)**2
constr(3) = -1.0_RP + x(5)**2 + x(6)**2
constr(4) = -1.0_RP + x(1)**2 + (x(2) + x(9))**2
constr(5) = -1.0_RP + (x(1) + x(5))**2 + (x(2) + x(6))**2
constr(6) = -1.0_RP + (x(1) + x(7))**2 + (x(2) + x(8))**2
constr(7) = -1.0_RP + (x(3) + x(5))**2 + (x(4) + x(6))**2
constr(8) = -1.0_RP + (x(3) + x(7))**2 + (x(4) + x(8))**2
constr(9) = -1.0_RP + x(7)**2 + (x(8) + x(9))**2
constr(10) = -x(1) * x(4) + x(2) * x(3)
constr(11) = -x(3) * x(9)
constr(12) = -x(5) * x(9)
constr(13) = -x(5) * x(8) + x(6) * x(7)
constr(14) = -x(9)
end subroutine calcfc_hexagon

end module calcfc_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program cobyla_exmp

! The following line makes the solver available.
use cobyla_mod, only : cobyla

! The following line specifies which module provides CALCFC.
use calcfc_mod, only : RP, calcfc_chebyquad, calcfc_hexagon

implicit none

integer, parameter :: n_chebyquad = 6
real(RP) :: x_chebyquad(n_chebyquad)
integer, parameter :: n_hexagon = 9
real(RP) :: x_hexagon(n_hexagon)
real(RP) :: f, cstrv
real(RP), allocatable :: constr(:)
integer :: m, i, nf, info

! The following lines illustrates how to call the solver to solve the Chebyquad problem.
x_chebyquad = [(real(i, RP) / real(n_chebyquad + 1, RP), i=1, n_chebyquad)] ! Starting point
m = 0  ! Dimension of constraints. M must the specified correctly, or the program will crash!
call cobyla(calcfc_chebyquad, m, x_chebyquad, f, cstrv)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also CONSTR, RHOBEG,
! and IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x_chebyquad = [(real(i, RP) / real(n_chebyquad + 1, RP), i=1, n_chebyquad)] ! Starting point
call cobyla(calcfc_chebyquad, m, x_chebyquad, f, cstrv, rhobeg=0.2_RP * x_chebyquad(1), iprint=1, nf=nf, info=info)

! The following lines illustrates how to call the solver to solve the Hexagon problem.
x_hexagon = 2.0_RP  ! Starting point.
m = 14  ! Dimension of constraints. M must the specified correctly, or the program will crash!
call cobyla(calcfc_hexagon, m, x_hexagon, f, cstrv)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also CONSTR, RHOBEG,
! and IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver. Note that CONSTR is an output, which will be set to
! the value of CONSTR(X_HEXAGON) when the solver returns.
x_hexagon = 2.0_RP  ! Starting point.
call cobyla(calcfc_hexagon, m, x_hexagon, f, cstrv, nlconstr=constr, rhobeg=1.0_RP, rhoend=1.0D-3, iprint=1, nf=nf, info=info)

deallocate (constr) ! Deallocate the array CONSTR, which is allocated by the solver. Otherwise, memory leaks.

end program cobyla_exmp
