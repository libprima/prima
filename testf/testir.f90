module test_mod
contains

subroutine assert(condition)
logical, intent(in) :: condition

if (condition) then
    print *, 'Condition true.'
else
    print *, "Condition false."
end if
end subroutine assert

function lt(x, y) result(l)
real, intent(in) :: x(:), y(:)
logical :: l
l = all(x < y)
end function lt

end module test_mod

program testcrash
use, non_intrinsic :: test_mod, only : assert, lt
implicit none
integer :: i, n
real :: x(2), y(2)

n = size(x); 
x = 1.0
y = 2.0

if (.true.) then
    call assert(all([(.true., i=1, 1)]))
    call assert(all([(lt([x(i), x(i + 1)], [y(i), y(i + 1)]), i=1, n - 1)]))
    call assert(maxval([(x(i) * y(i), i=1, 2)]) >= 1)
end if

end program testcrash
