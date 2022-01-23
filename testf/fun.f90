module fun_mod

implicit none
private
public :: calfun


contains


subroutine calfun(x, y)

implicit none
real, intent(in) :: x(:)
real, intent(out) :: y

y = sum(x**2)

end subroutine calfun


end module fun_mod
