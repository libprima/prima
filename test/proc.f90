module proc_mod 

implicit none
private 
public :: proc


contains


subroutine proc(x, y, f)

use intf_mod, only : func
use subproc_mod, only : subproc
implicit none

real, intent(in) :: x(:)
real, intent(out) :: y
procedure (func) :: f

call subproc(x, y, f)
y = y + 1.0

end subroutine proc 


end module proc_mod
