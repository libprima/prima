module subproc_mod

implicit none
private
public :: subproc


contains


subroutine subproc(x, y, f)

use intf_mod, only : func
implicit none

real, intent(in) :: x(:)
real, intent(out) :: y
procedure(func) :: f


call f(x, y)

end subroutine subproc


end module subproc_mod
