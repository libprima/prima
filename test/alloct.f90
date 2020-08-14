program alloct
real, allocatable :: a(:)
real :: b(2)
a = (/ b(1:0), b(2:1), b(3:2)/) 
print *, allocated(a), size(a)
if (allocated(a)) deallocate(a)
a = b(1:0)
print *, allocated(a), size(a)
end program alloct
