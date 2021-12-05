module fun_mod
implicit none
private
public :: FUN
abstract interface
    subroutine FUN(x, f)
    implicit none
    real, intent(in) :: x(:)
    real, intent(out) :: f
    end subroutine FUN
end interface
end module fun_mod


module test_mod
use fun_mod, only : FUN
implicit none
private
public :: PROBLEM, quadratic

type PROBLEM
    integer :: n
    real, allocatable :: x0(:)
    procedure(FUN), nopass, pointer :: objective => null()
end type PROBLEM

contains

subroutine quadratic(x, f)
implicit none
real, intent(in) :: x(:)
real, intent(out) :: f
f = sum(x**2)
end subroutine quadratic

end module test_mod


program main

use test_mod, only : PROBLEM, quadratic
implicit none

type(PROBLEM) :: p
procedure(FUN) :: quadratic
real :: y

integer :: i

p % n = 4
allocate (p % x0(p % n))
p % x0 = [(i, i=1, p % n)]
p % objective => quadratic
call p % objective(p % x0, y)

write (*, *) y

end program main
