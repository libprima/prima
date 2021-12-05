module test
implicit none
type problem
    integer :: n
    real, allocatable :: x0(:)
    procedure(FUN), nopass, pointer :: objective => null()
contains
    procedure, pass :: init !this is type bound procedure
end type problem

abstract interface
    subroutine FUN(x, f)
    implicit none
    real, intent(in) :: x(:)
    real, intent(out) :: f
    end subroutine FUN
end interface
contains
subroutine init(this, func, num)
class(problem), intent(inout) :: this
procedure(FUN) :: func
integer, intent(in) :: num
this % objective => func
this % n = num
!....
end subroutine init

end module test
