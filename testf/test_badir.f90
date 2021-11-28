module test_mod

contains

subroutine test
implicit none
integer :: k
real :: a(1) = 0.0
call bar(foo([(a(k), k=1, 1)]))
end subroutine test

pure function foo(x) result(y)
implicit none
real, intent(in) :: x(:)
real :: y
y = x(1)
end function foo

subroutine bar(y)
real, intent(in) :: y
print *, y
end subroutine bar

end module test_mod

program testprog
use test_mod
call test
end program testprog
