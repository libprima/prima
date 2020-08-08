module proc_mod 

implicit none
private 
public :: proc


contains


subroutine proc(x, y, f)

implicit none

real, intent(in) :: x(:)
real, intent(out) :: y

interface
    subroutine f(x, y)
    real, intent(in) :: x(:)
    real, intent(out) :: y
    end subroutine f
end interface

call f(x, y)
y = y + 1.0

end subroutine proc 


end module proc_mod
