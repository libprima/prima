module string_mod
implicit none
private
public :: int2str

contains
function int2str(i) result(res)
character, allocatable :: res(:)
integer, intent(in) :: i
character(range(i) + 2) :: tmp
write (tmp, '(i0)') i
write (*, *) tmp, len(tmp), trim(tmp), len(trim(tmp))
allocate (res(len(trim(tmp))))
res = trim(tmp)
end function
end module string_mod

program test
use, non_intrinsic :: string_mod, only : int2str
print *, int2str(-1)
end program test
