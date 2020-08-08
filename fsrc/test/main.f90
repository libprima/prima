program main 

use proc_mod, only : proc
use fun_mod, only : calfun
implicit none

real :: x(2), y

call random_number(x)
print *, x

call proc(x, y, calfun)
print *, y

print *, sum(x**2) + 1.0

end program main 
