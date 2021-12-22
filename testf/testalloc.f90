!        This is file : testalloc
! Author= zaikunzhang
! Started at: 21.12.2021
! Last Modified: Tuesday, December 21, 2021 PM03:21:42

module alloc_mod
implicit none
private
public :: alloc

contains

subroutine alloc(x, n)
integer, intent(out), allocatable :: x(:)
integer, intent(in) :: n
allocate (x(n))
end subroutine alloc

end module alloc_mod

program testalloc
use, non_intrinsic :: alloc_mod, only : alloc
implicit none
integer, allocatable :: x(:)
integer :: n

n = 0

call alloc(x, n)

print *, x, size(x)

end program testalloc
