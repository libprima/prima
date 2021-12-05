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
public :: PROBLEM, init


type PROBLEM
    integer :: n
    real, allocatable :: x0(:)
    procedure(FUN), nopass, pointer :: objective => null()
!contains
!    procedure, pass :: init  ! This is a type-bounding procedure
end type PROBLEM


contains


subroutine init(prob, probname, n)
implicit none
!class(PROBLEM), intent(inout) :: prob
type(PROBLEM) :: prob
character(len=*), intent(in) :: probname
integer, intent(in) :: n
integer :: i
select case (probname)
case ('quad')
    prob % objective => quad
! There can be other cases ...
case default
    prob % objective => null()
end select
prob % n = n
allocate (prob % x0(n))  ! Needed by Absoft Fortran 22.0
prob % x0 = [(i, i=1, n)]
end subroutine init

subroutine quad(x, f)
implicit none
real, intent(in) :: x(:)
real, intent(out) :: f
f = sum(x**2)
end subroutine quad

end module test_mod


program main

use test_mod!, only : PROBLEM, init
implicit none

type(PROBLEM) :: p
real :: y

!call p % init('quad', 4)
call init(p, 'quad', 4)
call p % objective(p % x0, y)

write (*, *) y

end program main
