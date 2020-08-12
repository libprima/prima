!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module proc_mod 
! PROC_MOD provides a subroutine PROC(X, Y, F) that invokes a subroutine
! F at (X, Y) and then increment Y by 1.

implicit none
private 
public :: proc


contains


subroutine proc(x, y, f)

implicit none

real, intent(in) :: x(:)
real, intent(out) :: y

interface
    subroutine f(x, y)
    real, intent(in) :: x(:)
    real, intent(out) :: y
    end subroutine f
end interface

call f(x, y)
y = y + 1.0

end subroutine proc 


end module proc_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fun_mod
! FUN_MOD provides a subroutine FUN(X, Y) that sets Y = SUM(X**2).

implicit none
private 
public :: fun 


contains 


subroutine fun(x, y)

implicit none
real, intent(in) :: x(:)
real, intent(out)  :: y 

y = sum(x**2)

end subroutine fun


end module fun_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main 
! MAIN tests PROC(X, Y, FUN) with PROC from PROC_MOD and FUN from FUN_MOD.

use proc_mod, only : proc
use fun_mod, only : fun
implicit none

real :: x(3), y

call random_number(x)
print *, x

call proc(x, y, fun)
print *, y

print *, sum(x**2) + 1.0

end program main 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
