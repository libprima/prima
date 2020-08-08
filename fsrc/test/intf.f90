module intf_mod 

implicit none
private 
public :: func

interface
    subroutine func(x, y)
    real, intent(in) :: x(:)
    real, intent(out) :: y
    end subroutine func
end interface

end module intf_mod
