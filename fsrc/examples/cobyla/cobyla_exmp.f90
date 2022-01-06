!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code.
!
! Started: July 2020
!
! Last Modified: Thursday, January 06, 2022 PM04:03:59
!--------------------------------------------------------------------------------------------------!


!-------------------------------- THE MODULE THAT IMPLEMENTS CALCFC -------------------------------!
module calcfc_mod

implicit none
private
public :: calcfc_chebyquad, calcfc_hexagon

contains

subroutine calcfc_chebyquad(x, f, constr)
! The Chebyquad test problem (Fletcher, 1965)
implicit none

real(kind(0.0D0)), intent(in) :: x(:)
real(kind(0.0D0)), intent(out) :: f
real(kind(0.0D0)), intent(out) :: constr(:)

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
constr = 0.0D0  ! Without this line, compilers may complain that CONSTR is not set.
end subroutine calcfc_chebyquad

subroutine calcfc_hexagon(x, f, constr)
! Test problem 10 (Hexagon area) in Powell's original package.
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_HEXAGON'
real(kind(1.0D0)), intent(in) :: x(:)
real(kind(1.0D0)), intent(out) :: constr(:)
real(kind(1.0D0)), intent(out) :: f

call assert(size(x) == 9 .and. size(constr) == 14, 'SIZE(X) == 9, SIZE(CONSTR) == 14', srname)

f = -0.5D0 * (x(1) * x(4) - x(2) * x(3) + x(3) * x(9) - x(5) * x(9) + x(5) * x(8) - x(6) * x(7))
constr(1) = 1.0D0 - x(3)**2 - x(4)**2
constr(2) = 1.0D0 - x(9)**2
constr(3) = 1.0D0 - x(5)**2 - x(6)**2
constr(4) = 1.0D0 - x(1)**2 - (x(2) - x(9))**2
constr(5) = 1.0D0 - (x(1) - x(5))**2 - (x(2) - x(6))**2
constr(6) = 1.0D0 - (x(1) - x(7))**2 - (x(2) - x(8))**2
constr(7) = 1.0D0 - (x(3) - x(5))**2 - (x(4) - x(6))**2
constr(8) = 1.0D0 - (x(3) - x(7))**2 - (x(4) - x(8))**2
constr(9) = 1.0D0 - x(7)**2 - (x(8) - x(9))**2
constr(10) = x(1) * x(4) - x(2) * x(3)
constr(11) = x(3) * x(9)
constr(12) = -x(5) * x(9)
constr(13) = x(5) * x(8) - x(6) * x(7)
constr(14) = x(9)
end subroutine calcfc_hexagon

end module calcfc_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program cobyla_exmp

! The following line makes the solver available.
use cobyla_mod, only : cobyla

! The following line specifies which module provides CALCFC. If CALCFC is given by an external
! subroutine instead of a module, remove this line and uncomment the "external" line below.
use calcfc_mod, only : calcfc_chebyquad, calcfc_hexagon

implicit none

real(kind(0.0D0)) :: f
integer, parameter :: n_chebyquad = 6
real(kind(0.0D0)) :: x_chebyquad(n_chebyquad)
integer, parameter :: n_hexagon = 9
real(kind(0.0D0)) :: x_hexagon(n_hexagon)
real(kind(0.0D0)), allocatable :: constr(:)
integer :: m, i

! If CALCFC is an external subroutine, then remove the line of  "use calcfc_mod", and uncomment the
! following line.
!--------------------------------------------------------------------------------------------------!
!!external calcfc_chebyquad, calcfc_hexagon
!--------------------------------------------------------------------------------------------------!

! The following lines illustrates how to call the solver to solve the Chebyquad problem.
x_chebyquad = [(real(i, kind(0.0D0)) / real(n_chebyquad + 1, kind(0.0D0)), i=1, n_chebyquad)]
m = 0  ! Dimension of constraints defined by CALCFC_CHEBYQUAD (there is none).
call cobyla(calcfc_chebyquad, x_chebyquad, f, m)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also CONSTR, RHOBEG,
! and IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
call cobyla(calcfc_chebyquad, x_chebyquad, f, m, rhobeg=0.2D0 * x_chebyquad(1), iprint=1)

! The following lines illustrates how to call the solver to solve the Hexagon problem.
x_hexagon = 1.0D0
m = 14  ! Dimension of constraints defined by CALCFC_HEXAGON.
call cobyla(calcfc_hexagon, x_hexagon, f, m)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also CONSTR, RHOBEG,
! and IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver. Note that CONSTR is an output, which will be set to
! the value of CONSTR(X_HEXAGON) when the solver returns.
call cobyla(calcfc_hexagon, x_hexagon, f, m, constr=constr, rhobeg=1.0D0, iprint=1)

end program cobyla_exmp
