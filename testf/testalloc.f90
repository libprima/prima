!        This is file : testalloc
! Author= zaikunzhang
! Started at: 21.12.2021
! Last Modified: Friday, January 28, 2022 PM12:18:51

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

print *, 'Test ALLOCATE(X(0))'
if (allocated(x)) deallocate (x); allocate (x(0))
print *, 'Finished with SIZE(X) = ', size(x)

print *, 'Test ALLOCATE(X(-1))'
if (allocated(x)) deallocate (x); allocate (x(-1))
print *, 'Finished with SIZE(X) = ', size(x)

end program testalloc
