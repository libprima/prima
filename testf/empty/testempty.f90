!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! testempty.f90

!!!!!! A module that offends AF90/AF95 !!!!!!!!!!!!!!!!!!!!!!!!
module empty_mod

implicit none
private
public :: foo

contains

subroutine foo(n)
implicit none
integer, intent(in) :: n
integer :: a(0)
integer :: b(n - 1)
integer, allocatable :: c(:)
allocate (c(n - 1))
call bar(c)  ! Update: AF90/AF95 is happy with this line.
call bar(a)  ! AF90/AF95 is happy with this line.
call bar(b)  ! AF90/AF95 is angry with this line.
end subroutine foo

subroutine bar(x)
implicit none
integer, intent(out) :: x(:)
x = 1  ! BAR(B) annoys AF90/AF95 regardless of this line.
end subroutine bar

end module empty_mod
!!!!!! Module ends !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!! Main program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program testempty

use empty_mod, only : foo

implicit none

call foo(2)  ! AF90/AF95 is happy with this line.
call foo(1)  ! AF90/AF95 is angry with this line.
write (*, *) 'Succeed!'  ! Declare victory when arriving here.

end program testempty
!!!!!! Main program ends !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
